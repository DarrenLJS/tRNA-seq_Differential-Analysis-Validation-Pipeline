"""
workflow/scripts/compute_delta_c_v2.py

Second, count-based per-codon score -- companion to compute_delta_c.py's
rate/GLM-based Delta(c) (rule 14), NOT a replacement for it. Both are
computed and both feed rule 16/17.

Unlike compute_delta_c.py, there is no whitelist-weighted term-type sum to
do here: rule build_codon_count_matrix has already produced real,
DESeq2/edgeR-tested codon-level counts (canonical + I34/Q34-variant
anticodon contributions already summed in at the count stage, per
build_codon_count_matrix.py), so this script's only job is reshaping
codon_highconf_intersect_kappa{kappa}.tsv into the SAME output schema
compute_delta_c.py uses (timepoint, codon, delta_c, ..., kappa), so rule
16/17 can consume delta_c_v1 and delta_c_v2 uniformly without
special-casing which one they're looking at.

delta_c(v2) for a given (codon, timepoint) is simply that row's DESeq2
log2FoldChange from the codon-level fit -- "the fold change between
control and treated in the anticodons [that decode this codon]", per
supervisor spec, now estimated directly from real per-codon counts and
replicate variance (DESeq2 + edgeR agreement, same FDR/replicate-r gating
as rule 10) rather than combined post-hoc from separate FC(i) and f(i)
estimates the way v1 does.

MISSING DATA: a codon absent from the high-confidence intersect (failed
DESeq2/edgeR significance-agreement, or filtered out upstream) gets
delta_c = NaN here, NOT imputed to 0/FC=1 the way v1 imputes missing
isodecoder FCs -- v1's imputation is a deliberate whitelist-completeness
choice (every whitelist-reachable isodecoder must contribute something);
v2 has no such completeness obligation since it is not summing
per-isodecoder terms, so an untested codon is left honestly undefined
rather than assigned a value it doesn't have evidence for.
"""

import logging
import os
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


def compute_delta_c_v2(codon_highconf_path, kappa, out_path, log_path=None):
    highconf = pd.read_csv(codon_highconf_path, sep="\t")

    required = {"isodecoder_id", "timepoint", "log2FoldChange"}
    missing = required - set(highconf.columns)
    if missing:
        raise ValueError(
            f"{codon_highconf_path} missing required columns: {missing}. "
            f"Found: {list(highconf.columns)}"
        )

    # NOTE: 'isodecoder_id' is the column name deseq2_isodecoder.R always
    # writes regardless of level (see rule deseq2_codon docstring for why
    # -- reused script, hardcoded rowname column) -- at codon level it
    # actually holds the codon string. Renamed here at the boundary so
    # nothing downstream has to remember this quirk.
    highconf = highconf.rename(columns={"isodecoder_id": "codon"})

    # Only pairwise per-timepoint rows carry a real codon-level contrast;
    # the LRT row's "timepoint" is the literal string "all_LRT" and isn't
    # part of the per-timepoint SENSE_CODONS grid below -- pass it through
    # separately rather than forcing it into the same completeness table.
    pairwise = highconf[highconf["timepoint"] != "all_LRT"].copy()
    lrt = highconf[highconf["timepoint"] == "all_LRT"].copy()

    pairwise["timepoint"] = pairwise["timepoint"].astype(int)
    timepoints = sorted(pairwise["timepoint"].unique())

    fc_idx = pairwise.set_index(["codon", "timepoint"])["log2FoldChange"]

    results = []
    for tp in timepoints:
        for codon in SENSE_CODONS:
            key = (codon, tp)
            delta_c = float(fc_idx.loc[key]) if key in fc_idx.index else np.nan
            results.append(dict(
                timepoint=tp, codon=codon, delta_c=delta_c,
                tested=key in fc_idx.index,
                kappa=kappa,
            ))

    out = pd.DataFrame(results).sort_values(["timepoint", "delta_c"], ascending=[True, False])
    out.to_csv(out_path, sep="\t", index=False)

    n_tested = out["tested"].sum()
    n_undefined = out["delta_c"].isna().sum()
    log.info(f"Wrote Delta(c) v2: {len(out)} (timepoint, codon) rows -> {out_path}")
    log.info(f"Codons with a real high-confidence codon-level FC: {n_tested}/{len(out)}")
    log.info(f"Undefined (not in high-confidence intersect): {n_undefined}")
    if not lrt.empty:
        log.info(f"LRT rows present but not included in the per-timepoint grid: {len(lrt)}")

    if log_path:
        with open(log_path, "a") as fh:
            fh.write(f"kappa={kappa}\n")
            fh.write(f"Tested: {n_tested}/{len(out)}, undefined: {n_undefined}\n")

    return out


if __name__ == "__main__":
    log_path = snakemake.log[0] if len(snakemake.log) else None
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()

    compute_delta_c_v2(
        codon_highconf_path=snakemake.input.codon_highconf,
        kappa=float(snakemake.wildcards.kappa),
        out_path=snakemake.output.delta_c_v2,
        log_path=log_path,
    )
