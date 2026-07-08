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
An isodecoder can be whitelist-reachable for codon c but have no DESeq2 row
in the high-confidence intersect (rule 10) -- either because it failed the
DESeq2/edgeR significance-agreement filter, or was low-count-filtered
upstream. This IS treated as FC(i)=1 (log2FoldChange=0, i.e. "no detected
change") rather than dropped from the sum. This is a deliberate modelling
assumption, not a neutral default -- "not significant" is being read as
"assume no change" rather than "insufficient power to know." It is tracked
separately from observed contributions: each Delta(c) row reports
n_observed_isodecoders (real FC(i) used) vs n_imputed_isodecoders (FC=1
assumed) vs n_reachable_isodecoders (whitelist total), so a codon resting
entirely on imputed zeros is distinguishable from one backed by real
fold-changes, rather than both silently looking like equally solid numbers.
(Only FC(i) is imputed this way. I34/Q34 GLM f_stim/f_ctrl gaps are NOT
imputed -- see below -- those reflect insufficient position-34 read
coverage, a different failure mode from "not significant.")

f_ctrl(i) == 0 or f_stim(i) == 0 edge case, and missing GLM rows generally:
the log ratio is undefined, or the modification-rate measurement simply
doesn't exist for that isodecoder/timepoint. Rather than crash or silently
substitute a pseudocount that would bias small values, these terms are
DROPPED (not imputed) with a logged note per isodecoder -- pseudocount
handling should happen upstream in the GLM fit (rule 11), not be invented
here. Exception: for both_Q_C at kappa=0, the term mathematically reduces
to log2(FC(i)) regardless of f_stim/f_ctrl (ratio**0 == 1), so it does NOT
require a GLM row at all in that case -- see both_Q_C branch below.
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
    impute_log = []

    for tp in timepoints:
        for codon in SENSE_CODONS:
            rows = whitelist[whitelist["codon"] == codon]
            contributions = []
            n_reachable = len(rows)
            n_observed = 0
            n_imputed = 0

            for _, wr in rows.iterrows():
                iso_id = wr["isodecoder_id"]
                term_type = wr["term_type"]

                key = (iso_id, tp)
                fc_imputed = key not in fc_idx.index
                if fc_imputed:
                    # Not in the high-confidence DESeq2/edgeR intersect (either
                    # failed the significance-agreement filter, or low-count
                    # filtered upstream). Imputed as FC=1 ("no detected
                    # change") rather than dropped -- see MISSING DATA
                    # HANDLING docstring for the reasoning and caveats.
                    impute_log.append(f"[{tp}] {codon}: {iso_id} ({term_type}) FC imputed=1.0 -- not in high-confidence DESeq2/edgeR intersect")
                    FC_i = 1.0
                else:
                    FC_i = fc_idx.loc[key]
                    if not np.isfinite(FC_i) or FC_i <= 0:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} ({term_type}) skipped -- FC non-finite or <=0 ({FC_i})")
                        continue

                if term_type in ("canonical", "both_I", "both_Q_U"):
                    contributions.append(np.log2(FC_i))
                    n_imputed += 1 if fc_imputed else 0
                    n_observed += 0 if fc_imputed else 1

                elif term_type == "mod_only_I":
                    if key not in i34_idx.index:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (mod_only_I) skipped -- no I34 GLM result")
                        continue
                    f_stim, f_ctrl = i34_idx.loc[key]
                    if not (np.isfinite(f_stim) and np.isfinite(f_ctrl)) or f_stim <= 0 or f_ctrl <= 0:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (mod_only_I) skipped -- f_stim/f_ctrl invalid ({f_stim}, {f_ctrl})")
                        continue
                    contributions.append(np.log2(FC_i * (f_stim / f_ctrl)))
                    n_imputed += 1 if fc_imputed else 0
                    n_observed += 0 if fc_imputed else 1

                elif term_type == "both_Q_C":
                    # kappa=0 short-circuit FIRST: the term mathematically
                    # reduces to log2(FC(i)) when kappa=0 ((ratio)**0 == 1)
                    # and does not depend on f_stim/f_ctrl at all -- so it
                    # must not require a Q34 GLM row to be present. Gating
                    # on q34_idx membership before this check was the bug:
                    # it silently dropped kappa=0 contributions purely for
                    # lacking GLM coverage they never needed.
                    if kappa == 0:
                        contributions.append(np.log2(FC_i))
                        n_imputed += 1 if fc_imputed else 0
                        n_observed += 0 if fc_imputed else 1
                        continue
                    if key not in q34_idx.index:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (both_Q_C) skipped -- no Q34 GLM result")
                        continue
                    f_stim, f_ctrl = q34_idx.loc[key]
                    if not (np.isfinite(f_stim) and np.isfinite(f_ctrl)) or f_stim <= 0 or f_ctrl <= 0:
                        skip_log.append(f"[{tp}] {codon}: {iso_id} (both_Q_C) skipped -- f_stim/f_ctrl invalid ({f_stim}, {f_ctrl})")
                        continue
                    contributions.append(np.log2(FC_i * (f_stim / f_ctrl) ** kappa))
                    n_imputed += 1 if fc_imputed else 0
                    n_observed += 0 if fc_imputed else 1

                else:
                    raise ValueError(f"Unrecognised term_type '{term_type}' in whitelist row for {iso_id}/{codon}")

            n_used = n_observed + n_imputed
            delta_c = float(np.sum(contributions)) if contributions else np.nan
            results.append(dict(
                timepoint=tp, codon=codon, delta_c=delta_c,
                n_reachable_isodecoders=n_reachable,
                n_observed_isodecoders=n_observed,
                n_imputed_isodecoders=n_imputed,
                n_used_isodecoders=n_used,
                kappa=kappa,
            ))

    out = pd.DataFrame(results).sort_values(["timepoint", "delta_c"], ascending=[True, False])
    out.to_csv(out_path, sep="\t", index=False)

    n_undefined = out["delta_c"].isna().sum()
    n_partial = ((out["n_used_isodecoders"] < out["n_reachable_isodecoders"]) & (out["n_used_isodecoders"] > 0)).sum()
    n_fully_imputed = ((out["n_observed_isodecoders"] == 0) & (out["n_imputed_isodecoders"] > 0)).sum()
    log.info(f"Wrote Delta(c): {len(out)} (timepoint, codon) rows -> {out_path}")
    log.info(f"Undefined Delta(c) (zero usable isodecoders): {n_undefined}")
    log.info(f"Partially-covered Delta(c) (some but not all whitelist isodecoders usable): {n_partial}")
    log.info(f"Fully-imputed Delta(c) (no real FC observed, resting entirely on FC=1 imputation): {n_fully_imputed}")
    log.info(f"FC(i) values imputed as 1.0 (not in high-confidence DESeq2/edgeR intersect): {len(impute_log)}")

    if log_path:
        with open(log_path, "w") as fh:
            fh.write(f"Delta(c) computation summary (kappa={kappa})\n")
            fh.write(f"Timepoints: {timepoints}\n")
            fh.write(f"Undefined Delta(c) rows: {n_undefined}\n")
            fh.write(f"Partially-covered Delta(c) rows: {n_partial}\n")
            fh.write(f"Fully-imputed Delta(c) rows (no real FC observed at all): {n_fully_imputed}\n")
            fh.write(f"Total FC(i) imputations (FC=1, not in high-confidence intersect): {len(impute_log)}\n\n")
            fh.write("Per-isodecoder skip reasons (GLM missing/invalid -- NOT imputed):\n")
            fh.write("\n".join(skip_log[:2000]))  # cap to avoid unbounded log files
            if len(skip_log) > 2000:
                fh.write(f"\n... ({len(skip_log) - 2000} more skip entries truncated)\n")
            fh.write("\n\nPer-isodecoder FC imputations (FC=1 assumed -- not in high-confidence intersect):\n")
            fh.write("\n".join(impute_log[:2000]))
            if len(impute_log) > 2000:
                fh.write(f"\n... ({len(impute_log) - 2000} more imputation entries truncated)\n")

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
