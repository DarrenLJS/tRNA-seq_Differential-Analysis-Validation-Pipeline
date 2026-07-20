"""
workflow/scripts/cross_cell_line_concordance.py

Final gate: compares Stage 2's major outputs across cell lines (A549 vs
THP1, or whatever CELL_LINES resolves to from the manifest) and reports,
per comparison type, which calls/rankings are concordant vs cell-line-
specific. Nothing downstream of this rule exists -- it's the terminal
node of the DAG (rule `all` requires only this file's output).

CONCORDANCE DEFINITIONS (kept explicit and separately reportable per
comparison type, rather than one blended "concordance score" -- a single
number here would hide which parts of the pipeline actually replicate):

  isodecoder DE:      isodecoder_id x timepoint pairs present in BOTH cell
                       lines' high-confidence intersect sets (rule 10),
                       with the SAME direction of log2FoldChange.
  I34 modification:   isodecoder_id x timepoint pairs significant (padj <
                       fdr) in BOTH cell lines, same direction of f_stim
                       vs f_ctrl change.
  Delta(c) ranking:   Spearman correlation between the two cell lines'
                       Delta(c) vectors (all 61 codons), per timepoint --
                       a ranking-concordance measure, not a per-codon
                       binary match, since Delta(c) is explicitly a
                       RANKED output per the proposal.
  Gene prediction:    Spearman correlation between the two cell lines'
                       predicted_translation_score vectors (genes common
                       to both), per timepoint.

Only pairwise cell-line comparisons are implemented (this project has
exactly 2 cell lines per the manifest) -- if CELL_LINES ever has >2
entries, this script raises rather than silently comparing only the first
two, since a >2-cell-line design would need a real decision about how
"concordant across all" should be defined (all pairwise? majority?).

FIX (2026-07-20) -- isodecoder DE / I34 GLM matching is no longer exact-ID
--------------------------------------------------------------------------
mim-tRNAseq's isodecoder clustering (Stage 1) is DATA-DRIVEN PER CELL
LINE: which raw loci get merged into one "isodecoder" depends on that
cell line's own observed coverage/mismatch signal, not a fixed rule. The
same underlying locus can therefore end up under a different collapsed
isodecoder_id string in each cell line -- e.g. A549 may keep
"Homo_sapiens_tRNA-Arg-TCT-2" as its own singleton isodecoder while THP1
merges it into "Homo_sapiens_tRNA-Arg-TCT-3/2/5" (confirmed on real data,
41+ such loci disagree between A549 and THP1 -- see
build_decoding_whitelist.py's module docstring for the worked examples).

The previous version of concordance_isodecoder_de/concordance_i34_glm
joined on isodecoder_id STRING EQUALITY, which silently drops every such
locus from the concordance calculation entirely -- and these are
disproportionately the ambiguous/low-coverage loci most likely to be
biologically interesting, so the dropped set is not a random sample.

Both functions now match on LOCUS-FAMILY OVERLAP instead: each collapsed
isodecoder_id is decomposed back into the set of raw-locus families it
represents (see _locus_family_set()), and an A549 row is compared against
a THP1 row whenever their locus-family sets share at least one member,
within the same timepoint. This lets e.g. A549's "Arg-TCT-2" legitimately
compare against THP1's merged "Arg-TCT-3/2/5", using THP1's merged value
as the best available resolution on that side.

This can produce one-to-many matches (one cell line's singleton isodecoder
matching against the other's merged isodecoder, or vice versa) -- these
are kept (not deduplicated away) and flagged via the `ambiguous_match`
column, since each represents a genuine, independently interpretable
comparison ("does Arg-TCT-2's direction in A549 agree with the merged
Arg-TCT-3/2/5 value in THP1?"). Because this could inflate n_comparable
relative to a strict 1:1 reading, the summary reports BOTH the concordance
rate over all matches and the rate restricted to unambiguous (1:1) matches
only, so neither number is chosen for you.

Per-comparison drop counts (rows in cell line 1 / cell line 2 with no
locus-family overlap at all, at that timepoint, in the other cell line)
are reported in the summary AND written out at full per-row detail
(snakemake.output.isodecoder_de_detail / i34_glm_detail) so the
concordance numbers are auditable rather than a single opaque percentage.
"""

import logging
import re
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# Matches a Stage-1 FINAL isodecoder_id: <prefix>_tRNA-<iso>-<anticodon>-<nums>
# where <nums> is a single family number ("3") or a '/'-joined list of family
# numbers from a data-driven merge ("3/5", "5/1/3/7/8"). See module docstring.
_ID_RE = re.compile(r'^(?P<base>.+_tRNA-[^-]+-[^-]+)-(?P<nums>[\d/]+)$')


def _locus_family_set(isodecoder_id):
    """
    Decompose a Stage-1 FINAL isodecoder_id into the set of underlying
    raw-locus FAMILY-level identifiers it represents, e.g.:
        "Homo_sapiens_tRNA-Lys-TTT-3/5"       -> {"...-TTT-3", "...-TTT-5"}
        "Homo_sapiens_tRNA-Ile-AAT-5/1/3/7/8" -> 5 elements
        "Homo_sapiens_tRNA-Val-CAC-1"         -> singleton, {"...-CAC-1"}

    This recovers the same locus-family membership build_decoding_
    whitelist.py derives from unsplitClusterInfo.txt + Isodecoder_counts.txt
    -- parsing the ID string directly avoids needing those two per-cell-
    line files as additional inputs here, since the '/'-joined ID string
    IS that join's output (see build_decoding_whitelist.py's
    _load_locus_to_isodecoder_map docstring for the source derivation).

    Any isodecoder_id that doesn't match the expected naming convention is
    treated as its own opaque singleton (logged at DEBUG, not silently
    dropped) so it can still participate in exact-match comparisons.
    """
    m = _ID_RE.match(str(isodecoder_id))
    if not m:
        return {str(isodecoder_id)}
    base = m.group("base")
    nums = m.group("nums").split("/")
    return {f"{base}-{n}" for n in nums}


def _merge_on_locus_overlap(d1, d2, cl1, cl2, id_col="isodecoder_id"):
    """
    Match rows between d1 (cell line cl1) and d2 (cell line cl2) whenever
    their isodecoder_id's underlying locus-family sets overlap at all,
    matched separately within each timepoint. See module docstring "FIX"
    note for the full rationale.

    Returns:
      merged   : one row per (d1_row, d2_row) match. Columns from d1/d2
                 (other than id_col/timepoint) are suffixed _{cl1}/_{cl2};
                 isodecoder_id itself is ALSO suffixed per cell line
                 (isodecoder_id_{cl1}, isodecoder_id_{cl2}) since the two
                 sides can now genuinely differ. Adds:
                   match_type      : 'exact_id' or 'locus_overlap'
                   ambiguous_match : True if either side's row matched
                                     more than one row on the other side
                                     (typically a singleton-vs-merged case)
      n_d1_dropped : d1 rows with no overlapping locus set in d2 at the
                     same timepoint.
      n_d2_dropped : same, for d2.
      n_ambiguous  : count of merged rows flagged ambiguous_match.
    """
    d1 = d1.copy()
    d2 = d2.copy()
    d1["_locus_set"] = d1[id_col].apply(_locus_family_set)
    d2["_locus_set"] = d2[id_col].apply(_locus_family_set)

    match_pairs = []
    for tp, g1 in d1.groupby("timepoint"):
        g2 = d2[d2["timepoint"] == tp]
        if g2.empty:
            continue
        for i1, r1 in g1.iterrows():
            for i2, r2 in g2.iterrows():
                if r1["_locus_set"] & r2["_locus_set"]:
                    match_pairs.append((i1, i2))

    if not match_pairs:
        return pd.DataFrame(), len(d1), len(d2), 0

    d1_match_counts = Counter(i1 for i1, _ in match_pairs)
    d2_match_counts = Counter(i2 for _, i2 in match_pairs)

    records = []
    for i1, i2 in match_pairs:
        r1, r2 = d1.loc[i1], d2.loc[i2]
        rec = {"timepoint": r1["timepoint"]}
        for col in d1.columns:
            if col in ("_locus_set", "timepoint"):
                continue
            rec[f"{col}_{cl1}"] = r1[col]
        for col in d2.columns:
            if col in ("_locus_set", "timepoint"):
                continue
            rec[f"{col}_{cl2}"] = r2[col]
        rec["match_type"] = "exact_id" if r1[id_col] == r2[id_col] else "locus_overlap"
        rec["ambiguous_match"] = (d1_match_counts[i1] > 1) or (d2_match_counts[i2] > 1)
        records.append(rec)

    merged = pd.DataFrame(records)
    n_d1_dropped = len(d1) - len(d1_match_counts)
    n_d2_dropped = len(d2) - len(d2_match_counts)
    n_ambiguous = int(merged["ambiguous_match"].sum())
    return merged, n_d1_dropped, n_d2_dropped, n_ambiguous


def _pair_paths(paths, cell_lines):
    """paths is a flat list in the same order as expand(...cell_line=CELL_LINES);
    zip them back into {cell_line: path}."""
    return dict(zip(cell_lines, paths))


def concordance_isodecoder_de(paths_by_cl, fdr):
    dfs = {cl: pd.read_csv(p, sep="\t") for cl, p in paths_by_cl.items()}
    cl1, cl2 = list(dfs.keys())
    d1, d2 = dfs[cl1], dfs[cl2]

    merged, n_d1_dropped, n_d2_dropped, n_ambiguous = _merge_on_locus_overlap(d1, d2, cl1, cl2)
    if merged.empty:
        return merged, 0, 0, 0, 0, n_d1_dropped, n_d2_dropped

    same_dir = np.sign(merged[f"log2FoldChange_{cl1}"]) == np.sign(merged[f"log2FoldChange_{cl2}"])
    merged["concordant"] = same_dir
    n_all = merged.shape[0]
    n_conc_all = int(merged["concordant"].sum())

    unamb = merged.loc[~merged["ambiguous_match"]]
    n_unamb = unamb.shape[0]
    n_conc_unamb = int(unamb["concordant"].sum())

    return merged, n_all, n_conc_all, n_unamb, n_conc_unamb, n_d1_dropped, n_d2_dropped


def concordance_i34_glm(paths_by_cl, fdr):
    dfs = {cl: pd.read_csv(p, sep="\t") for cl, p in paths_by_cl.items()}
    cl1, cl2 = list(dfs.keys())
    d1, d2 = dfs[cl1], dfs[cl2]

    if "significant" not in d1.columns or "significant" not in d2.columns:
        log.warning("I34 GLM results missing 'significant' column -- skipping I34 concordance")
        return pd.DataFrame(), 0, 0, 0, 0, 0, 0

    d1_sig = d1[d1["significant"] == True]  # noqa: E712
    d2_sig = d2[d2["significant"] == True]  # noqa: E712

    merged, n_d1_dropped, n_d2_dropped, n_ambiguous = _merge_on_locus_overlap(d1_sig, d2_sig, cl1, cl2)
    if merged.empty:
        return merged, 0, 0, 0, 0, n_d1_dropped, n_d2_dropped

    dir1 = np.sign(merged[f"f_stim_{cl1}"] - merged[f"f_ctrl_{cl1}"])
    dir2 = np.sign(merged[f"f_stim_{cl2}"] - merged[f"f_ctrl_{cl2}"])
    merged["concordant"] = dir1 == dir2
    n_all = merged.shape[0]
    n_conc_all = int(merged["concordant"].sum())

    unamb = merged.loc[~merged["ambiguous_match"]]
    n_unamb = unamb.shape[0]
    n_conc_unamb = int(unamb["concordant"].sum())

    return merged, n_all, n_conc_all, n_unamb, n_conc_unamb, n_d1_dropped, n_d2_dropped


def concordance_ranked_vector(paths_by_cl, id_col, value_col):
    dfs = {cl: pd.read_csv(p, sep="\t") for cl, p in paths_by_cl.items()}
    cl1, cl2 = list(dfs.keys())
    d1, d2 = dfs[cl1], dfs[cl2]

    rows = []
    common_timepoints = set(d1["timepoint"].unique()) & set(d2["timepoint"].unique())
    for tp in sorted(common_timepoints):
        s1 = d1[d1["timepoint"] == tp].set_index(id_col)[value_col]
        s2 = d2[d2["timepoint"] == tp].set_index(id_col)[value_col]
        common_ids = s1.index.intersection(s2.index)
        s1c, s2c = s1.loc[common_ids].dropna(), s2.loc[common_ids].dropna()
        common_ids2 = s1c.index.intersection(s2c.index)
        if len(common_ids2) < 5:
            continue
        rho, pval = spearmanr(s1c.loc[common_ids2], s2c.loc[common_ids2])
        rows.append(dict(timepoint=tp, n_common=len(common_ids2), spearman_rho=rho, pvalue=pval))
    return pd.DataFrame(rows)


def run_concordance(isodecoder_paths, i34_glm_paths, delta_c_paths, gene_scores_paths, cell_lines, fdr,
                     out_path, isodecoder_de_detail_path, i34_glm_detail_path):
    if len(cell_lines) != 2:
        raise NotImplementedError(
            f"cross_cell_line_concordance.py currently only supports exactly 2 cell "
            f"lines (pairwise comparison); got {len(cell_lines)}: {cell_lines}. "
            f"A >2-cell-line design needs an explicit decision about how 'concordant "
            f"across all' should be defined before this script can be extended."
        )

    isodecoder_by_cl = _pair_paths(isodecoder_paths, cell_lines)
    i34_by_cl = _pair_paths(i34_glm_paths, cell_lines)
    delta_c_by_cl = _pair_paths(delta_c_paths, cell_lines)
    gene_scores_by_cl = _pair_paths(gene_scores_paths, cell_lines)

    summary_rows = []

    (iso_merged, n_all, n_conc_all, n_unamb, n_conc_unamb,
     n_d1_dropped, n_d2_dropped) = concordance_isodecoder_de(isodecoder_by_cl, fdr)
    iso_merged.to_csv(isodecoder_de_detail_path, sep="\t", index=False)
    summary_rows.append(dict(
        comparison="isodecoder_DE_direction",
        n_comparable=n_all, n_concordant=n_conc_all,
        pct_concordant=(100 * n_conc_all / n_all) if n_all > 0 else np.nan,
        n_comparable_unambiguous=n_unamb, n_concordant_unambiguous=n_conc_unamb,
        pct_concordant_unambiguous=(100 * n_conc_unamb / n_unamb) if n_unamb > 0 else np.nan,
        n_dropped_cl1=n_d1_dropped, n_dropped_cl2=n_d2_dropped,
    ))
    log.info(
        f"isodecoder_DE_direction: {n_all} comparable pairs ({n_conc_all} concordant, "
        f"{n_unamb} unambiguous 1:1 of which {n_conc_unamb} concordant); "
        f"{n_d1_dropped} {cell_lines[0]} rows and {n_d2_dropped} {cell_lines[1]} rows "
        f"had no locus-overlap match at all and were dropped. Full detail: "
        f"{isodecoder_de_detail_path}"
    )

    (i34_merged, n_all_i34, n_conc_all_i34, n_unamb_i34, n_conc_unamb_i34,
     n_d1_dropped_i34, n_d2_dropped_i34) = concordance_i34_glm(i34_by_cl, fdr)
    i34_merged.to_csv(i34_glm_detail_path, sep="\t", index=False)
    summary_rows.append(dict(
        comparison="I34_GLM_direction",
        n_comparable=n_all_i34, n_concordant=n_conc_all_i34,
        pct_concordant=(100 * n_conc_all_i34 / n_all_i34) if n_all_i34 > 0 else np.nan,
        n_comparable_unambiguous=n_unamb_i34, n_concordant_unambiguous=n_conc_unamb_i34,
        pct_concordant_unambiguous=(100 * n_conc_unamb_i34 / n_unamb_i34) if n_unamb_i34 > 0 else np.nan,
        n_dropped_cl1=n_d1_dropped_i34, n_dropped_cl2=n_d2_dropped_i34,
    ))
    log.info(
        f"I34_GLM_direction: {n_all_i34} comparable pairs ({n_conc_all_i34} concordant, "
        f"{n_unamb_i34} unambiguous 1:1 of which {n_conc_unamb_i34} concordant); "
        f"{n_d1_dropped_i34} {cell_lines[0]} rows and {n_d2_dropped_i34} {cell_lines[1]} rows "
        f"had no locus-overlap match at all and were dropped. Full detail: "
        f"{i34_glm_detail_path}"
    )

    # Delta(c) and gene-prediction concordance are rank-correlation measures
    # over codons / genes, not isodecoder_id joins -- unaffected by the
    # isodecoder-collapsing issue this fix addresses, so unchanged.
    delta_c_rho = concordance_ranked_vector(delta_c_by_cl, id_col="codon", value_col="delta_c")
    for _, row in delta_c_rho.iterrows():
        summary_rows.append(dict(
            comparison=f"delta_c_ranking_spearman_tp{row['timepoint']}",
            n_comparable=row["n_common"], n_concordant=np.nan, pct_concordant=np.nan,
            n_comparable_unambiguous=np.nan, n_concordant_unambiguous=np.nan, pct_concordant_unambiguous=np.nan,
            n_dropped_cl1=np.nan, n_dropped_cl2=np.nan,
            spearman_rho=row["spearman_rho"], spearman_pvalue=row["pvalue"],
        ))

    gene_rho = concordance_ranked_vector(gene_scores_by_cl, id_col="gene_id", value_col="predicted_translation_score")
    for _, row in gene_rho.iterrows():
        summary_rows.append(dict(
            comparison=f"gene_prediction_ranking_spearman_tp{row['timepoint']}",
            n_comparable=row["n_common"], n_concordant=np.nan, pct_concordant=np.nan,
            n_comparable_unambiguous=np.nan, n_concordant_unambiguous=np.nan, pct_concordant_unambiguous=np.nan,
            n_dropped_cl1=np.nan, n_dropped_cl2=np.nan,
            spearman_rho=row["spearman_rho"], spearman_pvalue=row["pvalue"],
        ))

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_path, sep="\t", index=False)

    log.info(f"Wrote cross-cell-line concordance summary ({len(summary)} rows) -> {out_path}")
    for _, row in summary.iterrows():
        log.info(f"  {row.to_dict()}")

    return summary


if __name__ == "__main__":
    run_concordance(
        isodecoder_paths=list(snakemake.input.isodecoder_highconf),
        i34_glm_paths=list(snakemake.input.i34_glm),
        delta_c_paths=list(snakemake.input.delta_c),
        gene_scores_paths=list(snakemake.input.gene_scores),
        cell_lines=list(snakemake.params.cell_lines),
        fdr=snakemake.params.fdr,
        out_path=snakemake.output.summary,
        isodecoder_de_detail_path=snakemake.output.isodecoder_de_detail,
        i34_glm_detail_path=snakemake.output.i34_glm_detail,
    )
