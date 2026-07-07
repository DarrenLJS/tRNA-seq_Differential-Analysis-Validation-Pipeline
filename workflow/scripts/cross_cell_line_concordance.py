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
"""

import logging
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def _pair_paths(paths, cell_lines):
    """paths is a flat list in the same order as expand(...cell_line=CELL_LINES);
    zip them back into {cell_line: path}."""
    return dict(zip(cell_lines, paths))


def concordance_isodecoder_de(paths_by_cl, fdr):
    dfs = {cl: pd.read_csv(p, sep="\t") for cl, p in paths_by_cl.items()}
    cl1, cl2 = list(dfs.keys())
    d1, d2 = dfs[cl1], dfs[cl2]

    merged = d1.merge(d2, on=["isodecoder_id", "timepoint"], suffixes=(f"_{cl1}", f"_{cl2}"))
    if merged.empty:
        return pd.DataFrame(), 0, 0

    same_dir = np.sign(merged[f"log2FoldChange_{cl1}"]) == np.sign(merged[f"log2FoldChange_{cl2}"])
    merged["concordant"] = same_dir
    n_both_highconf = merged.shape[0]
    n_concordant = merged["concordant"].sum()
    return merged, n_both_highconf, n_concordant


def concordance_i34_glm(paths_by_cl, fdr):
    dfs = {cl: pd.read_csv(p, sep="\t") for cl, p in paths_by_cl.items()}
    cl1, cl2 = list(dfs.keys())
    d1, d2 = dfs[cl1], dfs[cl2]

    if "significant" not in d1.columns or "significant" not in d2.columns:
        log.warning("I34 GLM results missing 'significant' column -- skipping I34 concordance")
        return pd.DataFrame(), 0, 0

    d1_sig = d1[d1["significant"] == True]
    d2_sig = d2[d2["significant"] == True]
    merged = d1_sig.merge(d2_sig, on=["isodecoder_id", "timepoint"], suffixes=(f"_{cl1}", f"_{cl2}"))
    if merged.empty:
        return pd.DataFrame(), 0, 0

    dir1 = np.sign(merged[f"f_stim_{cl1}"] - merged[f"f_ctrl_{cl1}"])
    dir2 = np.sign(merged[f"f_stim_{cl2}"] - merged[f"f_ctrl_{cl2}"])
    merged["concordant"] = dir1 == dir2
    return merged, merged.shape[0], merged["concordant"].sum()


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


def run_concordance(isodecoder_paths, i34_glm_paths, delta_c_paths, gene_scores_paths, cell_lines, fdr, out_path):
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

    iso_merged, n_both, n_conc = concordance_isodecoder_de(isodecoder_by_cl, fdr)
    summary_rows.append(dict(
        comparison="isodecoder_DE_direction",
        n_comparable=n_both, n_concordant=n_conc,
        pct_concordant=(100 * n_conc / n_both) if n_both > 0 else np.nan,
    ))

    i34_merged, n_both_i34, n_conc_i34 = concordance_i34_glm(i34_by_cl, fdr)
    summary_rows.append(dict(
        comparison="I34_GLM_direction",
        n_comparable=n_both_i34, n_concordant=n_conc_i34,
        pct_concordant=(100 * n_conc_i34 / n_both_i34) if n_both_i34 > 0 else np.nan,
    ))

    delta_c_rho = concordance_ranked_vector(delta_c_by_cl, id_col="codon", value_col="delta_c")
    for _, row in delta_c_rho.iterrows():
        summary_rows.append(dict(
            comparison=f"delta_c_ranking_spearman_tp{row['timepoint']}",
            n_comparable=row["n_common"], n_concordant=np.nan,
            pct_concordant=np.nan, spearman_rho=row["spearman_rho"], spearman_pvalue=row["pvalue"],
        ))

    gene_rho = concordance_ranked_vector(gene_scores_by_cl, id_col="gene_id", value_col="predicted_translation_score")
    for _, row in gene_rho.iterrows():
        summary_rows.append(dict(
            comparison=f"gene_prediction_ranking_spearman_tp{row['timepoint']}",
            n_comparable=row["n_common"], n_concordant=np.nan,
            pct_concordant=np.nan, spearman_rho=row["spearman_rho"], spearman_pvalue=row["pvalue"],
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
    )
