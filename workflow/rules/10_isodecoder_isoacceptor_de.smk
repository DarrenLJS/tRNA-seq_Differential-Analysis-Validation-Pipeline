# =============================================================================
# workflow/rules/10_isodecoder_isoacceptor_de.smk
# Isodecoder/isoacceptor differential abundance -- proposal Section 3.5.
# DESeq2 pairwise (stimulated vs control per timepoint) + LRT across
# timepoints, per cell line, at both isodecoder and isoacceptor level.
# edgeR run in parallel as a sensitivity check; intersection = high-
# confidence set. Runs directly on Stage 1's already-built count matrices
# -- no new count generation happens here.
# =============================================================================

rule deseq2_isodecoder:
    input:
        counts  = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results  = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_DESeq2_results.tsv",
        rds      = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_dds.rds",
    params:
        lrt_model = config["diff_abundance"]["lrt_model"],
        fdr       = config["diff_abundance"]["fdr_threshold"],
        level     = "isodecoder",
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_deseq2_isodecoder.log",
    resources:
        sge_extra = sge_extra("deseq2_isodecoder"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/deseq2_isodecoder.R"


rule deseq2_isoacceptor:
    input:
        counts  = f"{SCRATCH}/deseq2_input/{{cell_line}}/isoacceptor_counts_matrix.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results  = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isoacceptor_DESeq2_results.tsv",
        rds      = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isoacceptor_dds.rds",
    params:
        lrt_model = config["diff_abundance"]["lrt_model"],
        fdr       = config["diff_abundance"]["fdr_threshold"],
        level     = "isoacceptor",
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_deseq2_isoacceptor.log",
    resources:
        sge_extra = sge_extra("deseq2_isodecoder"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/deseq2_isodecoder.R"   # same script, level param switches counts/output naming


rule edgeR_isodecoder_sensitivity:
    """
    edgeR run on the same isodecoder count matrix as an independent
    sensitivity check (proposal 3.5). Not a replacement for DESeq2 -- the
    intersection of the two calls is what feeds the "high-confidence" set
    used downstream in rule 14 (isodecoder FC values still come from
    DESeq2, per Stage 1's convention; edgeR here only gates which calls
    are trusted, it doesn't supply the FC(i) values itself).
    """
    input:
        counts  = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_edgeR_results.tsv",
    params:
        fdr = config["diff_abundance"]["fdr_threshold"],
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_edgeR_isodecoder.log",
    resources:
        sge_extra = sge_extra("edgeR_sensitivity_check"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/edgeR_sensitivity_check.R"


rule intersect_deseq2_edgeR:
    """
    Intersection of DESeq2 and edgeR significant calls (both FDR <
    diff_abundance.fdr_threshold, same direction of effect) = the
    high-confidence isodecoder set. Also applies the replicate-Pearson-r
    QC gate (diff_abundance.min_replicate_r) carried over from Stage 1 --
    timepoints failing it are flagged exploratory in the output, not
    silently dropped.
    """
    input:
        deseq2_results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_DESeq2_results.tsv",
        edgeR_results   = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_edgeR_results.tsv",
        counts          = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
        coldata         = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        highconf = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_highconf_intersect.tsv",
    params:
        fdr           = config["diff_abundance"]["fdr_threshold"],
        min_replicate_r = config["diff_abundance"]["min_replicate_r"],
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_intersect.log",
    resources:
        sge_extra = sge_extra("edgeR_sensitivity_check"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/intersect_deseq2_edgeR.py"
