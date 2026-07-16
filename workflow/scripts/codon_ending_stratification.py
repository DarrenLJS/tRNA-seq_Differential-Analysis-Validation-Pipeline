"""
workflow/scripts/codon_ending_stratification.py

Directly tests the project's central hypothesis (dissertation Section 1.6):
that poly(I:C)-induced tRNA pool remodelling selectively favours decoding
of G/C-ending codons and disfavours decoding of A/U-ending codons -- using
ONLY Delta(c) (rule 14), already computed. This does NOT depend on the
Watson et al. polysome data (rule 16) or on rules 15/17 -- it is an
internal-consistency test of the direction of Delta(c) itself, not a
validation of Delta(c) against an external ground truth. Keep this
distinction explicit in any downstream write-up: this rule shows Delta(c)
points the direction the hypothesis predicts; rule 16 is the separate,
still-pending test of whether that direction actually predicts translation.

WHY STRATIFY BY CODON, NOT BY ISODECODER
------------------------------------------
"Codon ending" is a property of a (isodecoder, codon) pair in the decoding
whitelist, not of an isodecoder alone -- a single I34-eligible isodecoder
decodes both a U-ending codon (both_I term) and C/A-ending codons
(mod_only_I term) of the same split box, so labelling the ISODECODER as
"G/C-decoding" or "A/U-decoding" is not well-defined in general. Delta(c)
(rule 14 output) is already aggregated to exactly one row per (timepoint,
codon), so it is the natural, unambiguous unit for this test.

CODON-ENDING DEFINITION
------------------------
DNA-letter convention (matching the whitelist/compute_delta_c.py
convention, T not U): a codon's wobble (3rd) position determines its
"ending". G/C-ending = codon[-1] in {G, C}. A/U-ending = codon[-1] in
{A, T}. This is computed directly from the `codon` column already present
in delta_c_kappa{kappa}.tsv -- no new upstream input is required.

STATISTICAL TEST
-----------------
Two independent groups of codons (G/C-ending vs A/U-ending) -- NOT paired
observations -- so a Mann-Whitney U test (the two-sample rank-sum test;
equivalent to an unpaired Wilcoxon rank-sum test, distinct from the PAIRED
Wilcoxon signed-rank test) is the correct non-parametric choice, rather
than a paired Wilcoxon test. Reported per timepoint (kappa is fixed by
the wildcard, matching rule 14's own per-kappa outputs):

  - two-sided p-value: is there ANY directional difference between the
    two groups' Delta(c) distributions?
  - one-sided p-value (alternative: G/C-ending > A/U-ending): the exact
    directional prediction in Section 1.6 -- reported alongside, not
    instead of, the two-sided result, so the test isn't silently
    re-framed as one-sided only when it happens to support the hypothesis.
  - rank-biserial correlation as a non-parametric effect size (range
    -1 to 1; positive = G/C-ending ranks higher, consistent with the
    hypothesis direction).

MISSING DATA
------------
Codons with undefined Delta(c) (no whitelist-reachable isodecoder at all
for that codon/timepoint -- see compute_delta_c.py) are excluded from
both groups before testing, not imputed. The count excluded from each
group is reported so a codon-ending imbalance in *coverage* (as opposed
to in Delta(c) itself) is visible rather than silently folded into the
test.
"""

import logging

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

GC_ENDING = {"G", "C"}
AU_ENDING = {"A", "T"}


def _rank_biserial(u_stat, n1, n2):
    """
    Rank-biserial correlation for a Mann-Whitney U test: r = 1 - 2U/(n1*n2).
    Positive r means group 1 (G/C-ending, passed first to mannwhitneyu)
    tends to rank higher than group 2 (A/U-ending) -- i.e. consistent with
    the hypothesis direction.
    """
    if n1 == 0 or n2 == 0:
        return np.nan
    return 1.0 - (2.0 * u_stat) / (n1 * n2)


def stratify_by_codon_ending(delta_c_path, alpha, out_path, log_path=None):
    delta_c = pd.read_csv(delta_c_path, sep="\t")
    required = {"timepoint", "codon", "delta_c", "kappa"}
    missing = required - set(delta_c.columns)
    if missing:
        raise ValueError(f"{delta_c_path} missing required columns: {missing}. Found: {list(delta_c.columns)}")

    kappa = delta_c["kappa"].iloc[0] if len(delta_c) else float("nan")
    ending = delta_c["codon"].str[-1]
    is_gc = ending.isin(GC_ENDING)
    is_au = ending.isin(AU_ENDING)
    n_other = int((~is_gc & ~is_au).sum())
    if n_other:
        log.warning(f"{n_other} codon row(s) have a wobble base outside {{A,C,G,T}} -- excluded from both groups.")

    results = []
    for tp, sub in delta_c.groupby("timepoint"):
        gc_all = sub[is_gc.loc[sub.index]]
        au_all = sub[is_au.loc[sub.index]]

        gc_defined = gc_all["delta_c"].dropna()
        au_defined = au_all["delta_c"].dropna()

        n_excluded_gc = len(gc_all) - len(gc_defined)
        n_excluded_au = len(au_all) - len(au_defined)

        row = dict(
            timepoint=tp,
            kappa=kappa,
            n_gc_ending=len(gc_defined),
            n_au_ending=len(au_defined),
            n_excluded_undefined_gc=n_excluded_gc,
            n_excluded_undefined_au=n_excluded_au,
            median_delta_c_gc=float(gc_defined.median()) if len(gc_defined) else np.nan,
            median_delta_c_au=float(au_defined.median()) if len(au_defined) else np.nan,
        )

        if len(gc_defined) < 2 or len(au_defined) < 2:
            log.warning(
                f"Timepoint '{tp}': too few defined Delta(c) codons in one or both groups "
                f"(G/C n={len(gc_defined)}, A/U n={len(au_defined)}) -- skipping test."
            )
            row.update(
                median_diff=np.nan, u_statistic=np.nan,
                pvalue_two_sided=np.nan, pvalue_one_sided_gc_greater=np.nan,
                rank_biserial=np.nan, significant_two_sided=False, significant_one_sided=False,
            )
            results.append(row)
            continue

        row["median_diff"] = row["median_delta_c_gc"] - row["median_delta_c_au"]

        u_two, p_two = mannwhitneyu(gc_defined, au_defined, alternative="two-sided")
        u_one, p_one = mannwhitneyu(gc_defined, au_defined, alternative="greater")

        row["u_statistic"] = float(u_two)
        row["pvalue_two_sided"] = float(p_two)
        row["pvalue_one_sided_gc_greater"] = float(p_one)
        row["rank_biserial"] = _rank_biserial(u_two, len(gc_defined), len(au_defined))
        row["significant_two_sided"] = bool(p_two < alpha)
        row["significant_one_sided"] = bool(p_one < alpha)

        results.append(row)

    out = pd.DataFrame(results).sort_values("timepoint")
    out.to_csv(out_path, sep="\t", index=False)

    n_sig_one_sided = int(out["significant_one_sided"].sum())
    log.info(f"Wrote {len(out)} timepoint row(s) -> {out_path}")
    log.info(
        f"Timepoints where G/C-ending Delta(c) is significantly greater than A/U-ending "
        f"(one-sided, alpha={alpha}): {n_sig_one_sided}/{len(out)}"
    )

    if log_path:
        with open(log_path, "w") as fh:
            fh.write(f"Codon-ending stratification of Delta(c) (kappa={kappa}, alpha={alpha})\n")
            fh.write(
                "This tests whether Delta(c) itself already points in the direction the "
                "central hypothesis predicts (G/C-ending codons favoured, A/U-ending "
                "disfavoured). It does NOT validate Delta(c) against external polysome "
                "data (rule 16) -- that is a separate, still-pending test.\n\n"
            )
            fh.write(out.to_string(index=False))
            fh.write("\n")

    return out


if __name__ == "__main__":
    stratify_by_codon_ending(
        delta_c_path=snakemake.input.delta_c,
        alpha=snakemake.params.alpha,
        out_path=snakemake.output.summary,
        log_path=snakemake.log[0] if len(snakemake.log) else None,
    )
