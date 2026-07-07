"""
workflow/scripts/gene_translation_prediction.py

Per-gene predicted translational efficiency score = dot product of the
gene's codon usage frequency vector (rule 09) with Delta(c) (rule 14),
per timepoint. A gene enriched for codons with high positive Delta(c)
(decoding capacity increasing under stimulation) gets a high score;
enriched for negative-Delta(c) codons gets a low score.

Codons with undefined Delta(c) (no usable contributing isodecoder -- see
compute_delta_c.py) are treated as 0 contribution for that gene's dot
product, NOT dropped from the gene's codon usage normalisation -- i.e. a
gene's score is computed over whatever codons Delta(c) actually covers,
weighted by that gene's real usage of them, rather than silently
renormalising the codon usage vector to only the covered codons (which
would inflate the contribution of well-covered codons in a
gene-composition-dependent way). This is a modelling choice worth a
methods footnote, not a hidden default.
"""

import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def predict_gene_scores(delta_c_path, codon_usage_path, out_path):
    delta_c = pd.read_csv(delta_c_path, sep="\t")
    codon_usage = pd.read_csv(codon_usage_path, sep="\t", index_col=0)

    timepoints = sorted(delta_c["timepoint"].unique())
    results = []

    for tp in timepoints:
        dc_tp = delta_c[delta_c["timepoint"] == tp].set_index("codon")["delta_c"]
        # Reindex to the full codon_usage column set, filling undefined
        # Delta(c) codons with 0 (explicit "no contribution", per docstring).
        dc_vec = dc_tp.reindex(codon_usage.columns).fillna(0.0)

        n_covered = dc_tp.notna().sum()
        n_total = len(codon_usage.columns)
        if n_covered < n_total:
            log.warning(f"Timepoint {tp}: only {n_covered}/{n_total} codons have a defined Delta(c); "
                        f"remaining codons contribute 0 to every gene's score this timepoint.")

        scores = codon_usage.values @ dc_vec.values
        tp_result = pd.DataFrame({
            "gene_id": codon_usage.index,
            "timepoint": tp,
            "predicted_translation_score": scores,
        })
        results.append(tp_result)

    final = pd.concat(results, ignore_index=True)
    final["rank_within_timepoint"] = final.groupby("timepoint")["predicted_translation_score"] \
        .rank(ascending=False, method="min")
    final = final.sort_values(["timepoint", "predicted_translation_score"], ascending=[True, False])
    final.to_csv(out_path, sep="\t", index=False)
    log.info(f"Wrote {len(final)} gene x timepoint prediction rows -> {out_path}")
    return final


if __name__ == "__main__":
    predict_gene_scores(
        delta_c_path=snakemake.input.delta_c,
        codon_usage_path=snakemake.input.codon_usage,
        out_path=snakemake.output.scores,
    )
