# =============================================================================
# workflow/rules/17_cross_cell_line_concordance.smk
# Final gate: nothing gets called a "robust biological signal" unless it
# reproduces in both A549 and THP1. Compares rule 10 (isodecoder DE), rule
# 11 (I34 GLM), rule 14 (Delta(c) ranking), and rule 15 (gene prediction
# ranking) across cell lines.
# =============================================================================

rule cross_cell_line_concordance:
    input:
        isodecoder_highconf = expand(
            f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_highconf_intersect.tsv",
            cell_line=CELL_LINES,
        ),
        i34_glm = expand(
            f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/I34_glm_results.tsv",
            cell_line=CELL_LINES,
        ),
        delta_c = expand(
            f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_kappa" + str(config["wobble_glm"]["kappa_default"]) + ".tsv",
            cell_line=CELL_LINES,
        ),
        gene_scores = expand(
            f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_kappa" + str(config["wobble_glm"]["kappa_default"]) + ".tsv",
            cell_line=CELL_LINES,
        ),
    output:
        summary = f"{STAGE2_ROOT}/concordance/cross_cell_line_concordance_summary.tsv",
    params:
        cell_lines = CELL_LINES,
        fdr = config["diff_abundance"]["fdr_threshold"],
    log:
        f"{STAGE2_ROOT}/logs/17_concordance/cross_cell_line_concordance.log",
    resources:
        sge_extra = sge_extra("cross_cell_line_concordance"),
    conda:
        os.path.join(STAGE1_ENV_DIR, "environment.yaml")
    script:
        "../scripts/cross_cell_line_concordance.py"
