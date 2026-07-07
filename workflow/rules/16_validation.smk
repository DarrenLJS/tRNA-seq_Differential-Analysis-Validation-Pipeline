# =============================================================================
# workflow/rules/16_validation.smk
# Validation against Watson et al. (2020) polysome profiling data --
# proposal Section 4, criteria 1-4. Fisher's exact (ISG/housekeeping
# representation among Watson et al. UP/DOWN sets) + Spearman correlation
# (predicted score vs observed polysome log2FC), per timepoint. Also runs
# the kappa sweep across rule 14/15's outputs to empirically justify (or
# not) a nonzero Q34 weight, rather than assuming one.
# =============================================================================

rule validate_fisher_spearman:
    """
    Per-cell-line validation at kappa_default (config wobble_glm.kappa_default,
    i.e. the value actually used for the "headline" result). The full
    kappa_sweep comparison is a separate rule (validate_kappa_sweep) so
    the headline result and the sweep diagnostics are clearly distinguished
    outputs, not mixed into one file.
    """
    input:
        gene_scores = f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_kappa" + str(config["wobble_glm"]["kappa_default"]) + ".tsv",
        watson_fc   = f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",
        gene_sets   = f"{STAGE2_ROOT}/references/isg_housekeeping_gene_sets.tsv",
    output:
        summary = f"{STAGE2_ROOT}/validation/{{cell_line}}/validation_summary.tsv",
    params:
        fisher_alpha   = config["validation"]["fisher_alpha"],
        spearman_alpha = config["validation"]["spearman_alpha"],
    log:
        f"{STAGE2_ROOT}/logs/16_validation/{{cell_line}}_validation.log",
    resources:
        sge_extra = sge_extra("validation"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/validate_fisher_spearman.R"


rule validate_kappa_sweep:
    """
    Sweeps kappa (config wobble_glm.kappa_sweep) across all cell lines and
    timepoints, recomputing the Spearman correlation (predicted vs Watson
    et al. observed) at each value. This is the empirical justification
    step for whatever kappa ends up used as the headline value -- if
    kappa=0 (i.e. ignoring Q34 entirely) turns out to correlate best, that
    is itself a reportable, legitimate finding (Q34 signal not adding
    predictive value at current sequencing depth), not a failure of the
    sweep.
    """
    input:
        gene_scores = expand(
            f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),
        watson_fc = f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",
    output:
        sweep_summary = f"{STAGE2_ROOT}/validation/kappa_sweep_summary.tsv",
    params:
        kappa_values = config["wobble_glm"]["kappa_sweep"],
        cell_lines   = CELL_LINES,
        stage2_root  = STAGE2_ROOT,
    log:
        f"{STAGE2_ROOT}/logs/16_validation/kappa_sweep.log",
    resources:
        sge_extra = sge_extra("validation"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/kappa_sweep_summary.R"
