# =============================================================================
# workflow/rules/15_gene_translational_prediction.smk
# Gene-level predicted translational efficiency score: dot product of
# Delta(c) with each gene's codon usage vector (rule 09). Computed across
# the same kappa_sweep as rule 14 so rule 16 can evaluate which kappa
# value actually improves agreement with Watson et al.'s observed data.
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
