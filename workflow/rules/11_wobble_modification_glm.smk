# =============================================================================
# workflow/rules/11_wobble_modification_glm.smk
# Wobble position-34 modification estimation -- proposal Section 3.6,
# split I34 (confirmatory) vs Q34 (exploratory), per the assay-sensitivity
# discussion: standard-RT misincorporation is a near-complete reporter of
# inosine (~100% ceiling) but only ~1% above baseline for queuosine even at
# full modification, so the two are NOT treated with the same statistical
# confidence, despite sharing the same underlying pos34_coverage/mismatch
# matrices from Stage 1.
#
# Isodecoder partitioning uses the SAME explicit whitelist (isotype +
# structural N34/N35 check) as rule 09's decoding_whitelist.tsv, not a
# separate ad hoc filter -- both consume config wobble_glm.i34_isotypes /
# q34_isotypes so the GLM's isodecoder set and Delta(c)'s whitelist buckets
# can never silently drift apart.
# =============================================================================

rule wobble_glm_i34:
    """
    Binomial GLM (successes=mismatch count, trials=coverage) at position
    34, restricted to I34-eligible isodecoders (A34, whitelisted isotype),
    stimulated vs control, per timepoint. Full confirmatory treatment --
    CI-based, FDR-corrected.
    """
    input:
        coverage = f"{SCRATCH}/deseq2_input/{{cell_line}}/pos34_coverage_matrix.tsv",
        mismatch = f"{SCRATCH}/deseq2_input/{{cell_line}}/pos34_mismatch_matrix.tsv",
        coldata  = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
        whitelist = f"{STAGE2_ROOT}/references/{{cell_line}}/decoding_whitelist.tsv",
    output:
        results = f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/I34_glm_results.tsv",
    params:
        fdr = config["wobble_glm"]["fdr_threshold"],
        bucket = "I34",
        confidence_tier = "confirmatory",
    log:
        f"{STAGE2_ROOT}/logs/11_wobble_glm/{{cell_line}}_I34.log",
    resources:
        sge_extra = sge_extra("wobble_glm"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/wobble_glm.R"


rule wobble_glm_q34:
    """
    Same binomial GLM machinery, restricted to Q34-eligible isodecoders
    (G34, N35=T, whitelisted isotype). Output is explicitly labeled
    "exploratory" -- reported with wide-uncertainty framing, NOT treated
    as a CI-based confirmatory result, per the assay-sensitivity discussion
    (standard-RT misincorporation for Q34 is close to background noise).
    The prob/f values here still feed the kappa-weighted Delta(c) term
    (rule 14), but kappa defaults to 0 (config wobble_glm.kappa_default)
    until the rule-16 sweep empirically justifies otherwise.
    """
    input:
        coverage = f"{SCRATCH}/deseq2_input/{{cell_line}}/pos34_coverage_matrix.tsv",
        mismatch = f"{SCRATCH}/deseq2_input/{{cell_line}}/pos34_mismatch_matrix.tsv",
        coldata  = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
        whitelist = f"{STAGE2_ROOT}/references/{{cell_line}}/decoding_whitelist.tsv",
    output:
        results = f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/Q34_glm_results_exploratory.tsv",
    params:
        fdr = config["wobble_glm"]["fdr_threshold"],
        bucket = "Q34",
        confidence_tier = "exploratory",
    log:
        f"{STAGE2_ROOT}/logs/11_wobble_glm/{{cell_line}}_Q34.log",
    resources:
        sge_extra = sge_extra("wobble_glm"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/wobble_glm.R"


rule wobble_glm_tgt_expression_check:
    """
    Companion check for Q34: does TGT (queuine tRNA-ribosyltransferase,
    QTRT1/QTRT2) expression itself change with poly(I:C) stimulation.
    Since direct Q34 sequencing signal is unreliable (see rule
    wobble_glm_q34 docstring), this orthogonal readout is the more
    defensible way to argue Q34 regulation is changing at all -- it does
    NOT feed into Delta(c) directly, it is reported alongside the Q34 GLM
    results in rule 16's validation summary as corroborating (or
    non-corroborating) evidence.

    FIX: requires QTRT1/QTRT2 to be present in whatever gene-level
    expression quantification is available -- Stage 1's pipeline is
    tRNA-focused and does NOT produce a genome-wide mRNA expression
    matrix. This rule currently assumes such a matrix exists somewhere
    (e.g. if total RNA-seq was also generated for these samples) -- if it
    does not, this rule cannot run as designed and TGT expression would
    need a separate qPCR/RNA-seq data source not currently in either
    stage's scope. Flagging this now rather than silently building a rule
    with no real input.
    """
    input:
        # PLACEHOLDER -- no Stage 1 output currently provides genome-wide
        # gene expression. This input will not resolve until a real source
        # is identified; left explicit rather than fabricated so Snakemake
        # fails loudly here instead of silently downstream.
        mrna_expression = f"{SCRATCH}/PLACEHOLDER_genome_wide_expression_matrix.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/TGT_expression_check.tsv",
    log:
        f"{STAGE2_ROOT}/logs/11_wobble_glm/{{cell_line}}_TGT_expression.log",
    resources:
        sge_extra = sge_extra("wobble_glm"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/tgt_expression_check.R"
