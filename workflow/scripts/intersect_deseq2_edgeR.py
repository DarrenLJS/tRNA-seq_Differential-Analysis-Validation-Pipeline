"""
workflow/scripts/intersect_deseq2_edgeR.py

Intersects DESeq2 and edgeR significant calls (both FDR < threshold, same
direction of log-fold-change) into a "high-confidence" isodecoder set per
timepoint, and applies the replicate-Pearson-r QC gate carried over from
Stage 1 -- timepoints whose within-condition replicate correlation falls
below diff_abundance.min_replicate_r are flagged `qc_flag=exploratory`
rather than dropped, so the information isn't lost, just clearly marked.
"""

import logging
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def compute_replicate_r(counts_path, coldata_path):
    """
    Per timepoint x condition group, mean pairwise Pearson r across
    replicates (log2(count+1) space, matching Stage 1's QC convention).
    Returns dict[(timepoint, condition)] = mean_r.
    """
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    coldata = pd.read_csv(coldata_path, sep="\t", index_col=0)
    coldata = coldata.loc[counts.columns]

    log_counts = np.log2(counts + 1)
    out = {}
    for (tp, cond), sub in coldata.groupby(["timepoint", "condition"]):
        samples = sub.index.tolist()
        if len(samples) < 2:
            out[(tp, cond)] = np.nan
            continue
        rs = []
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                r, _ = pearsonr(log_counts[samples[i]], log_counts[samples[j]])
                rs.append(r)
        out[(tp, cond)] = float(np.mean(rs)) if rs else np.nan
    return out


def intersect_calls(deseq2_path, edgeR_path, counts_path, coldata_path, fdr, min_replicate_r, out_path):
    deseq2 = pd.read_csv(deseq2_path, sep="\t")
    edgeR  = pd.read_csv(edgeR_path, sep="\t")

    # Only pairwise (per-timepoint) calls are intersected -- the LRT row
    # (timepoint == "all_LRT") has no edgeR counterpart by design (see
    # deseq2_isodecoder.R docstring) and is passed through separately.
    deseq2_pw = deseq2[deseq2["test"] == "pairwise_stim_vs_ctrl"].copy()
    deseq2_lrt = deseq2[deseq2["test"] == "LRT_condition_x_timepoint"].copy()

    merged = deseq2_pw.merge(
        edgeR, on=["isodecoder_id", "timepoint"], suffixes=("_deseq2", "_edgeR"), how="inner"
    )

    deseq2_sig = merged["padj"] < fdr
    edgeR_sig  = merged["FDR"] < fdr
    same_dir   = np.sign(merged["log2FoldChange"]) == np.sign(merged["logFC"])

    merged["highconf"] = deseq2_sig & edgeR_sig & same_dir

    replicate_r = compute_replicate_r(counts_path, coldata_path)
    def qc_flag(row):
        r_stim = replicate_r.get((row["timepoint"], "stimulated"), np.nan)
        r_ctrl = replicate_r.get((row["timepoint"], "control"), np.nan)
        worst = np.nanmin([r_stim, r_ctrl]) if not (np.isnan(r_stim) and np.isnan(r_ctrl)) else np.nan
        if np.isnan(worst):
            return "unknown_replicate_r"
        return "pass" if worst >= min_replicate_r else "exploratory_low_replicate_r"

    merged["qc_flag"] = merged.apply(qc_flag, axis=1)

    highconf_only = merged[merged["highconf"]].copy()
    highconf_only.to_csv(out_path, sep="\t", index=False)

    n_total = merged.shape[0]
    n_highconf = highconf_only.shape[0]
    n_exploratory = (highconf_only["qc_flag"] == "exploratory_low_replicate_r").sum()
    log.info(f"Intersected {n_total} isodecoder x timepoint pairs")
    log.info(f"High-confidence (DESeq2 + edgeR agree, same direction): {n_highconf}")
    log.info(f"...of which flagged exploratory (replicate r below {min_replicate_r}): {n_exploratory}")

    return highconf_only


if __name__ == "__main__":
    intersect_calls(
        deseq2_path=snakemake.input.deseq2_results,
        edgeR_path=snakemake.input.edgeR_results,
        counts_path=snakemake.input.counts,
        coldata_path=snakemake.input.coldata,
        fdr=snakemake.params.fdr,
        min_replicate_r=snakemake.params.min_replicate_r,
        out_path=snakemake.output.highconf,
    )
