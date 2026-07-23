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
#
# EXTENDED: also runs concordance against delta_c_v2 / v2 gene scores.
# cross_cell_line_concordance.py is fully generic (delta_c_paths and
# gene_scores_paths are plain file lists resolved by Snakemake, columns
# codon/delta_c and gene_id/predicted_translation_score respectively --
# identical schema for v1 and v2) and is reused UNCHANGED. The
# isodecoder_de_detail/i34_glm_detail outputs are NOT score-version-
# specific (they only depend on rule 10/11, not on Delta(c) at all), so
# rule cross_cell_line_concordance_v2 necessarily recomputes byte-
# identical copies of those two files under _v2-suffixed names -- flagged
# here as intentional, low-cost duplication (this rule is cheap, see
# config resources.cross_cell_line_concordance) rather than something to
# "optimize away" by splitting the script, which isn't worth the
# complexity for two small detail tables.
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


rule cross_cell_line_concordance_v2:
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
            f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_v2_kappa" + str(config["wobble_glm"]["kappa_default"]) + ".tsv",
            cell_line=CELL_LINES,
        ),
        gene_scores = expand(
            f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_v2_kappa" + str(config["wobble_glm"]["kappa_default"]) + ".tsv",
            cell_line=CELL_LINES,
        ),
    output:
        summary = f"{STAGE2_ROOT}/concordance/cross_cell_line_concordance_summary_v2.tsv",
        isodecoder_de_detail = f"{STAGE2_ROOT}/concordance/isodecoder_DE_matched_detail_v2.tsv",
        i34_glm_detail       = f"{STAGE2_ROOT}/concordance/I34_glm_matched_detail_v2.tsv",
    params:
        cell_lines = CELL_LINES,
        fdr = config["diff_abundance"]["fdr_threshold"],
    log:
        f"{STAGE2_ROOT}/logs/17_concordance/cross_cell_line_concordance_v2.log",
    resources:
        sge_extra = sge_extra("cross_cell_line_concordance"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/cross_cell_line_concordance.py"   # reused unchanged, see module docstring
