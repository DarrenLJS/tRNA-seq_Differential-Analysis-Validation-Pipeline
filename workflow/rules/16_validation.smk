# =============================================================================
# workflow/rules/16_validation.smk
# Validation against Watson et al. (2020) polysome profiling data --
# proposal Section 4, criteria 1-4. Fisher's exact (ISG/housekeeping
# representation among Watson et al. UP/DOWN sets) + Spearman correlation
# (predicted score vs observed polysome log2FC), per timepoint. Also runs
# the kappa sweep across rule 14/15's outputs to empirically justify (or
# not) a nonzero Q34 weight, rather than assuming one.
#
# EXTENDED: also validates delta_c_v2 (count/FC-based companion score,
# rule 14/15's v2 outputs), as a fully parallel set of rules -- v2 is
# compared against the same Watson et al. benchmark, not folded into the
# v1 numbers. validate_fisher_spearman.R is fully generic (just consumes
# whatever gene_scores file it's pointed at) and is reused UNCHANGED for
# validate_fisher_spearman_v2. kappa_sweep_summary.R needed an actual
# edit (it built its own v1-only file paths internally, unlike the other
# validation scripts which take pre-resolved Snakemake inputs) -- it now
# takes a required score_version param ("v1"/"v2") so both sweep rules
# below call the SAME edited script rather than diverging into two copies.
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


rule validate_fisher_spearman_v2:
    """
    v2 (count/FC-based delta_c_v2) sibling of validate_fisher_spearman
    above, at the same kappa_default. Same script, unchanged -- it only
    ever reads whatever gene_scores path it's given, with no v1-specific
    assumptions baked in.
    """
    input:
        gene_scores = f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_v2_kappa" + str(config["wobble_glm"]["kappa_default"]) + ".tsv",
        watson_fc   = f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",
        gene_sets   = f"{STAGE2_ROOT}/references/isg_housekeeping_gene_sets.tsv",
    output:
        summary = f"{STAGE2_ROOT}/validation/{{cell_line}}/validation_summary_v2.tsv",
    params:
        fisher_alpha   = config["validation"]["fisher_alpha"],
        spearman_alpha = config["validation"]["spearman_alpha"],
    log:
        f"{STAGE2_ROOT}/logs/16_validation/{{cell_line}}_validation_v2.log",
    resources:
        sge_extra = sge_extra("validation"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/validate_fisher_spearman.R"   # reused unchanged, see module docstring


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
        kappa_values  = config["wobble_glm"]["kappa_sweep"],
        cell_lines    = CELL_LINES,
        stage2_root   = STAGE2_ROOT,
        score_version = "v1",
    log:
        f"{STAGE2_ROOT}/logs/16_validation/kappa_sweep.log",
    resources:
        sge_extra = sge_extra("validation"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/kappa_sweep_summary.R"


rule validate_kappa_sweep_v2:
    """
    v2 sibling of validate_kappa_sweep above -- same script
    (kappa_sweep_summary.R), score_version="v2" switches which
    gene_translation_scores filename pattern it reads (see script
    docstring/edit note). `gene_scores` here is listed only so Snakemake
    schedules the v2 gene-prediction rule first; the script itself still
    resolves its own file paths internally per cell_line/kappa (unchanged
    from v1's design), it doesn't consume this input list directly.
    """
    input:
        gene_scores = expand(
            f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_v2_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),
        watson_fc = f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",
    output:
        sweep_summary = f"{STAGE2_ROOT}/validation/kappa_sweep_summary_v2.tsv",
    params:
        kappa_values  = config["wobble_glm"]["kappa_sweep"],
        cell_lines    = CELL_LINES,
        stage2_root   = STAGE2_ROOT,
        score_version = "v2",
    log:
        f"{STAGE2_ROOT}/logs/16_validation/kappa_sweep_v2.log",
    resources:
        sge_extra = sge_extra("validation"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/kappa_sweep_summary.R"
