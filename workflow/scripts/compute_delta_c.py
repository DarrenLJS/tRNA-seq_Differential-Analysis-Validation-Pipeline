"""
workflow/scripts/compute_delta_c.py

Computes the per-codon score Delta(c) for all 61 sense codons, per
timepoint, from:
  - decoding_whitelist.tsv   (rule 09 -- structural reachability + term_type)
  - isodecoder DESeq2 results (rule 10 -- FC(i), per timepoint)
  - I34 GLM results           (rule 11 -- f_stim/f_ctrl, per timepoint)
  - Q34 GLM results           (rule 11 -- f_stim/f_ctrl, per timepoint;
                                exploratory / low-confidence per proposal)
  - kappa                     (Q34 confidence dial, swept in rule 16)

FORMULA (see project derivation -- this is the fully-expanded, corrected
version, not the original 3-term proposal formula):

  Delta(c) = sum_i contribution(i, c)

  where, per whitelist term_type:
    canonical    -> log2( FC(i) )
    both_I       -> log2( FC(i) )                               [I34, U-ending]
    mod_only_I   -> log2( FC(i) * f_stim_I(i) / f_ctrl_I(i) )    [I34, C/A-ending]
    both_Q_C     -> log2( FC(i) * (f_stim_Q(i)/f_ctrl_Q(i))**kappa )  [Q34, C-ending]
    both_Q_U     -> log2( FC(i) )                               [Q34, U-ending]

MISSING DATA HANDLING
----------------------
An isodecoder can be whitelist-reachable for codon c but have no DESeq2
result (filtered for low counts) or no GLM result (insufficient position-34
coverage). Rather than silently treating a missing FC(i) as FC=1 (i.e. "no
change", which would be a real assumption, not a neutral default), missing
terms are DROPPED from the sum for that (codon, timepoint) and counted --
Delta(c)'s output row reports how many of the whitelist's contributing
isodecoders were actually usable, so a codon with 1/6 isodecoders covered
is visibly less trustworthy than one with 6/6, rather than both looking
like an equally solid number.

f_ctrl(i) == 0 or f_stim(i) == 0 edge case: the log ratio is undefined.
Rather than crash or silently substitute a pseudocount that would bias
small values, these are dropped with a logged note per isodecoder --
pseudocount handling should happen upstream in the GLM fit (rule 11), not
be invented here.
"""

import logging
from itertools import product

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

STOP_CODONS_DNA = {"TAA", "TAG", "TGA"}
SENSE_CODONS = sorted(
    a + b + c
    for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if a + b + c not in STOP_CODONS_DNA
)


def _load_fc_table(path):
    """
    Expects long format: isodecoder_id, timepoint, log2FoldChange
    (DESeq2's native output column name -- kept as-is rather than
    renamed, so this function is a thin, traceable wrapper, not a
    silent transform). FC(i) = 2**log2FoldChange is computed here.
    """
    df = pd.read_csv(path, sep="\t")
    required = {"isodecoder_id", "timepoint", "log2FoldChange"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"FC table {path} missing required columns: {missing}. Found: {list(df.columns)}")
    df = df.copy()
    df["FC"] = 2 ** df["log2FoldChange"]
    return df


def _load_glm_table(path):
    """
    Expects long format: isodecoder_id, timepoint, f_stim, f_ctrl.
    Rows with f_stim<=0 or f_ctrl<=0 are flagged (not dropped here --
    dropped per-use in compute_delta_c so the reason is attributable to
    the specific (codon, timepoint) it would have contributed to).
    """
    df = pd.read_csv(path, sep="\t")
    required = {"isodecoder_id", "timepoint", "f_stim", "f_ctrl"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GLM table {path} missing required columns: {missing}. Found: {list(df.columns)}")
    return df


def compute_delta_c(whitelist_path, fc_path, i34_glm_path, q34_glm_path, kappa, out_path, log_path=None):
    whitelist = pd.read_csv(whitelist_path, sep="\t")
    fc = _load_fc_table(fc_path)
    i34_glm = _load_glm_table(i34_glm_path)
    q34_glm = _load_glm_table(q34_glm_path)

    timepoints = sorted(fc["timepoint"].unique())
    log.info(f"Computing Delta(c) for {len(timepoints)} timepoint(s), kappa={kappa}")

    fc_idx      = fc.set_index(["isodecoder_id", "timepoint"])["FC"]
    i34_idx     = i34_glm.set_index(["isodecoder_id", "timepoint"])[["f_stim", "f_ctrl"]]
    q34_idx     = q34_glm.set_index(["isodecoder_id", "timepoint"])[["f_stim", "f_ctrl"]]

    results = []
    skip_log = []

    for tp in timepoints:
        for codon in SENSE_CODONS:
            rows = whitelist[whitelist["codon"] == codon]
            contributions = []
            n_reachable = len(rows)
            n_used = 0

            for _, wr in rows.iterrows():
                iso_id = wr["isodecoder_id"]
                term_type = wr["term_type"]

                key = (iso_id, tp)
                if key not in fc_idx.index:
                    skip_log.append(f"[{tp}] {codon}: {iso_id} ({term_type}) skipped -- no DESeq2 FC result (filtered/low count)")
                    continue
                FC_i = fc_idx.loc[key]
                if not np.isfinite(FC_i) or FC_i <= 0:
                    skip_log.append(f"[{tp}] {codon}: {iso_id} ({term_type}) skipped -- FC non-finite or <=0 ({FC_i})")
                    continue

                if term_type in ("canonical", "both_I", "both_Q_U"):
                    contributions.append(np.log2(FC_i))
                    n_used += 1

                elif term_type == "mod_only_I":
                    if key not in i34_idx.index:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (mod_only_I) skipped -- no I34 GLM result")
                        continue
                    f_stim, f_ctrl = i34_idx.loc[key]
                    if not (np.isfinite(f_stim) and np.isfinite(f_ctrl)) or f_stim <= 0 or f_ctrl <= 0:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (mod_only_I) skipped -- f_stim/f_ctrl invalid ({f_stim}, {f_ctrl})")
                        continue
                    contributions.append(np.log2(FC_i * (f_stim / f_ctrl)))
                    n_used += 1

                elif term_type == "both_Q_C":
                    if key not in q34_idx.index:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (both_Q_C) skipped -- no Q34 GLM result")
                        continue
                    f_stim, f_ctrl = q34_idx.loc[key]
                    if not (np.isfinite(f_stim) and np.isfinite(f_ctrl)) or f_stim <= 0 or f_ctrl <= 0:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (both_Q_C) skipped -- f_stim/f_ctrl invalid ({f_stim}, {f_ctrl})")
                        continue
                    if kappa == 0:
                        contributions.append(np.log2(FC_i))  # (ratio)**0 == 1, collapses to plain FC term
                    else:
                        contributions.append(np.log2(FC_i * (f_stim / f_ctrl) ** kappa))
                    n_used += 1

                else:
                    raise ValueError(f"Unrecognised term_type '{term_type}' in whitelist row for {iso_id}/{codon}")

            delta_c = float(np.sum(contributions)) if contributions else np.nan
            results.append(dict(
                timepoint=tp, codon=codon, delta_c=delta_c,
                n_reachable_isodecoders=n_reachable, n_used_isodecoders=n_used,
                kappa=kappa,
            ))

    out = pd.DataFrame(results).sort_values(["timepoint", "delta_c"], ascending=[True, False])
    out.to_csv(out_path, sep="\t", index=False)

    n_undefined = out["delta_c"].isna().sum()
    n_partial = ((out["n_used_isodecoders"] < out["n_reachable_isodecoders"]) & (out["n_used_isodecoders"] > 0)).sum()
    log.info(f"Wrote Delta(c): {len(out)} (timepoint, codon) rows -> {out_path}")
    log.info(f"Undefined Delta(c) (zero usable isodecoders): {n_undefined}")
    log.info(f"Partially-covered Delta(c) (some but not all whitelist isodecoders usable): {n_partial}")

    if log_path:
        with open(log_path, "w") as fh:
            fh.write(f"Delta(c) computation summary (kappa={kappa})\n")
            fh.write(f"Timepoints: {timepoints}\n")
            fh.write(f"Undefined Delta(c) rows: {n_undefined}\n")
            fh.write(f"Partially-covered Delta(c) rows: {n_partial}\n\n")
            fh.write("Per-isodecoder skip reasons:\n")
            fh.write("\n".join(skip_log[:2000]))  # cap to avoid unbounded log files
            if len(skip_log) > 2000:
                fh.write(f"\n... ({len(skip_log) - 2000} more skip entries truncated)\n")

    return out


if __name__ == "__main__":
    compute_delta_c(
        whitelist_path=snakemake.input.whitelist,
        fc_path=snakemake.input.isodecoder_fc,
        i34_glm_path=snakemake.input.i34_glm,
        q34_glm_path=snakemake.input.q34_glm,
        kappa=float(snakemake.wildcards.kappa),
        out_path=snakemake.output.delta_c,
        log_path=snakemake.log[0] if len(snakemake.log) else None,
    )
