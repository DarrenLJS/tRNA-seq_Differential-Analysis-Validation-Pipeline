# =============================================================================
# workflow/rules/17_cross_cell_line_concordance.smk
# Final gate: nothing gets called a "robust biological signal" unless it
# reproduces in both A549 and THP1. Compares rule 10 (isodecoder DE), rule
# 11 (I34 GLM), rule 14 (Delta(c) ranking), and rule 15 (gene prediction
# ranking) across cell lines.
#
# FIX (2026-07-20): isodecoder DE / I34 GLM matching now joins on
# locus-family overlap rather than exact isodecoder_id string equality,
# since mim-tRNAseq's per-cell-line data-driven clustering means the same
# underlying locus can collapse into differently-named isodecoders in each
# cell line (see cross_cell_line_concordance.py module docstring for the
# full rationale and worked example). Two new detail outputs
# (isodecoder_de_detail, i34_glm_detail) report every matched pair at full
# resolution -- including match_type (exact_id vs locus_overlap) and
# ambiguous_match flags -- so the concordance percentages in `summary` are
# auditable rather than a single opaque number.
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
        isodecoder_de_detail = f"{STAGE2_ROOT}/concordance/isodecoder_DE_matched_detail.tsv",
        i34_glm_detail       = f"{STAGE2_ROOT}/concordance/I34_glm_matched_detail.tsv",
    params:
        cell_lines = CELL_LINES,
        fdr = config["diff_abundance"]["fdr_threshold"],
    log:
        f"{STAGE2_ROOT}/logs/17_concordance/cross_cell_line_concordance.log",
    resources:
        sge_extra = sge_extra("cross_cell_line_concordance"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/cross_cell_line_concordance.py"
