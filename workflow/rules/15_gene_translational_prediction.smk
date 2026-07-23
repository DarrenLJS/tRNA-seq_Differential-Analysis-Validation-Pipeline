# =============================================================================
# workflow/rules/15_gene_translational_prediction.smk
# Gene-level predicted translational efficiency score: dot product of
# Delta(c) with each gene's codon usage vector (rule 09). Computed across
# the same kappa_sweep as rule 14 so rule 16 can evaluate which kappa
# value actually improves agreement with Watson et al.'s observed data.
#
# EXTENDED: also runs against delta_c_v2 (rule 14's count/FC-based
# companion score). gene_translation_prediction.py is fully generic (only
# needs a delta_c-shaped input: columns timepoint, codon, delta_c) and is
# reused UNCHANGED -- delta_c_v2_kappa{kappa}.tsv already has exactly
# those columns (see compute_delta_c_v2.py), so no script edit was needed
# here, only a new rule pointing at the v2 input/output paths.
# =============================================================================

rule gene_translation_prediction:
    input:
        delta_c     = f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_kappa{{kappa}}.tsv",
        codon_usage = f"{STAGE2_ROOT}/references/gene_codon_usage.tsv",
    output:
        scores = f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_kappa{{kappa}}.tsv",
    log:
        f"{STAGE2_ROOT}/logs/15_gene_prediction/{{cell_line}}_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("gene_translation_prediction"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/gene_translation_prediction.py"


rule gene_translation_prediction_v2:
    input:
        delta_c     = f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_v2_kappa{{kappa}}.tsv",
        codon_usage = f"{STAGE2_ROOT}/references/gene_codon_usage.tsv",
    output:
        scores = f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_v2_kappa{{kappa}}.tsv",
    log:
        f"{STAGE2_ROOT}/logs/15_gene_prediction/{{cell_line}}_v2_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("gene_translation_prediction"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/gene_translation_prediction.py"   # reused unchanged, see module docstring
