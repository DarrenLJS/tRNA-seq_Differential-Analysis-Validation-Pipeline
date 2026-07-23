# =============================================================================
# workflow/rules/10_isodecoder_isoacceptor_de.smk
# Isodecoder/isoacceptor differential abundance -- proposal Section 3.5.
# DESeq2 pairwise (stimulated vs control per timepoint) + LRT across
# timepoints, per cell line, at both isodecoder and isoacceptor level.
# edgeR run in parallel as a sensitivity check; intersection = high-
# confidence set. Runs directly on Stage 1's already-built count matrices
# -- no new count generation happens here.
#
# EXTENDED (count-based modification score branch): this file also builds
# and tests the anticodon- and codon-level count matrices for the second,
# count/FC-based per-codon score (delta_c_v2, rule 14) -- a companion to
# the existing rate/GLM-based Delta(c), not a replacement for it. This is
# the first place in file 10 that depends on rule 09's
# decoding_whitelist.tsv (the original 4 rules above only ever consume
# Stage 1 outputs) -- monotonic with the existing 09 -> 10/11 -> 14
# convention, just a dependency edge that didn't exist here before.
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

    FIX (feedback item 8b): `highconf` alone is pre-filtered to only the
    passing rows (the `highconf` boolean is all-True by construction in
    that file), so no downstream plot reading only `highconf` can ever show
    a failing isodecoder -- they were dropped before the file was written,
    not filtered in the plotting script. `highconf_all` is a second output
    carrying the FULL DESeq2 x edgeR merge (every isodecoder x timepoint
    pair with results in both tools), with the `highconf` boolean column
    intact but NOT filtered on, so "who passed" and "who didn't" both live
    in one well-defined, versioned file.
    """
    input:
        deseq2_results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_DESeq2_results.tsv",
        edgeR_results   = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_edgeR_results.tsv",
        counts          = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
        coldata         = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        highconf     = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_highconf_intersect.tsv",
        highconf_all = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_highconf_intersect_all.tsv",
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


# =============================================================================
# Count-based modification score branch (delta_c_v2, see rule 14).
# Chain: isodecoder_mismatch_table -> anticodon_count_matrix (kappa) ->
#        DESeq2/edgeR/intersect @ anticodon level -> codon_count_matrix ->
#        DESeq2/edgeR/intersect @ codon level -> (rule 14) compute_delta_c_v2
# =============================================================================

rule build_isodecoder_mismatch_table:
    """
    Collated isodecoder mismatch table: one row per (isodecoder, sample),
    with total_count_isd (Stage 1 isodecoder_counts_matrix.tsv),
    mismatch_A/T/G/C (per-nucleotide-type pos-34 misincorporation
    proportions, read from mim-tRNAseq's RAW mismatch/ directory -- NOT
    Stage 1's pos34_mismatch_matrix.tsv, which sums across type and would
    lose exactly the per-letter breakdown this needs), and actual_count_isd
    (canonical/reference-matching read-count estimate). See
    build_isodecoder_mismatch_table.py docstring for the full formula and
    the approximation it makes explicit (position-34 misread rate
    extrapolated across the whole isodecoder's read population).

    No Stage 1 code touched or rerun -- reads the same raw mismatch/
    directory Stage 1's own build_mismatch_matrices rule already consumes,
    as a second, independent parallel consumer.
    """
    input:
        mismatch_dir      = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/mismatch",
        isodecoder_counts = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
    output:
        mismatch_table = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/isodecoder_mismatch_table.tsv",
    params:
        cell_line = "{cell_line}",
        manifest  = config["manifest"],
    log:
        f"{STAGE2_ROOT}/logs/10_count_mod_score/{{cell_line}}_isodecoder_mismatch_table.log",
    resources:
        sge_extra = sge_extra("build_isodecoder_mismatch_table"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/build_isodecoder_mismatch_table.py"


rule build_anticodon_count_matrix:
    """
    Preprocessed anticodon count table (canonical + I34/Q34-modified
    variant rows, e.g. AGC / IGC / QGC), plus the wide count matrix for
    DESeq2/edgeR below. Reuses decoding_whitelist.tsv (rule 09) for
    isotype/eligibility/position34_base/codon-reach -- no wobble logic
    reimplemented here.

    I34 signal: full weight, not kappa-weighted (confirmatory, mirrors
    rule 11's treatment). Q34 signal: kappa-weighted, swept over the same
    wobble_glm.kappa_sweep values as the existing rate-based branch;
    kappa=0 means no Q34-variant rows at all. See
    build_anticodon_count_matrix.py docstring for the full I34/Q34 signal
    letter assumptions (config count_mod_score) and why they're flagged
    as an assumption rather than an established fact for Q34 specifically.
    """
    input:
        mismatch_table = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/isodecoder_mismatch_table.tsv",
        whitelist      = f"{STAGE2_ROOT}/references/{{cell_line}}/decoding_whitelist.tsv",
    output:
        anticodon_table  = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_count_table_kappa{{kappa}}.tsv",
        anticodon_matrix = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_counts_matrix_kappa{{kappa}}.tsv",
    params:
        i34_ref_base     = config["count_mod_score"]["i34_ref_base"],
        i34_signal_base  = config["count_mod_score"]["i34_signal_base"],
        q34_ref_base     = config["count_mod_score"]["q34_ref_base"],
        q34_signal_bases = config["count_mod_score"]["q34_signal_bases"],
    log:
        f"{STAGE2_ROOT}/logs/10_count_mod_score/{{cell_line}}_anticodon_matrix_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("build_anticodon_count_matrix"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/build_anticodon_count_matrix.py"


rule deseq2_anticodon:
    input:
        counts  = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_counts_matrix_kappa{{kappa}}.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_DESeq2_results_kappa{{kappa}}.tsv",
        rds     = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_dds_kappa{{kappa}}.rds",
    params:
        lrt_model = config["diff_abundance"]["lrt_model"],
        fdr       = config["diff_abundance"]["fdr_threshold"],
        level     = "anticodon",
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_deseq2_anticodon_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("deseq2_isodecoder"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/deseq2_isodecoder.R"   # reused unchanged; level param only affects logging/labels


rule edgeR_anticodon_sensitivity:
    input:
        counts  = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_counts_matrix_kappa{{kappa}}.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_edgeR_results_kappa{{kappa}}.tsv",
    params:
        fdr = config["diff_abundance"]["fdr_threshold"],
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_edgeR_anticodon_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("edgeR_sensitivity_check"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/edgeR_sensitivity_check.R"


rule intersect_anticodon_deseq2_edgeR:
    input:
        deseq2_results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_DESeq2_results_kappa{{kappa}}.tsv",
        edgeR_results   = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_edgeR_results_kappa{{kappa}}.tsv",
        counts          = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_counts_matrix_kappa{{kappa}}.tsv",
        coldata         = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        highconf     = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_highconf_intersect_kappa{{kappa}}.tsv",
        highconf_all = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_highconf_intersect_all_kappa{{kappa}}.tsv",
    params:
        fdr             = config["diff_abundance"]["fdr_threshold"],
        min_replicate_r = config["diff_abundance"]["min_replicate_r"],
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_intersect_anticodon_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("edgeR_sensitivity_check"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/intersect_deseq2_edgeR.py"


rule build_codon_count_matrix:
    """
    Expands the anticodon-level long table into a codon-level wide count
    matrix (many-to-many: one anticodon's count is duplicated into every
    codon it reaches; several anticodons reaching the same codon sum).
    No whitelist re-parsing -- read_codon reach is already carried per row
    by build_anticodon_count_matrix.py. See build_codon_count_matrix.py
    docstring.
    """
    input:
        anticodon_table = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_count_table_kappa{{kappa}}.tsv",
    output:
        codon_matrix = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/codon_counts_matrix_kappa{{kappa}}.tsv",
    log:
        f"{STAGE2_ROOT}/logs/10_count_mod_score/{{cell_line}}_codon_matrix_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("build_anticodon_count_matrix"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/build_codon_count_matrix.py"


rule deseq2_codon:
    input:
        counts  = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/codon_counts_matrix_kappa{{kappa}}.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_DESeq2_results_kappa{{kappa}}.tsv",
        rds     = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_dds_kappa{{kappa}}.rds",
    params:
        lrt_model = config["diff_abundance"]["lrt_model"],
        fdr       = config["diff_abundance"]["fdr_threshold"],
        level     = "codon",
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_deseq2_codon_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("deseq2_isodecoder"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/deseq2_isodecoder.R"   # reused unchanged, see rule deseq2_anticodon


rule edgeR_codon_sensitivity:
    input:
        counts  = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/codon_counts_matrix_kappa{{kappa}}.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_edgeR_results_kappa{{kappa}}.tsv",
    params:
        fdr = config["diff_abundance"]["fdr_threshold"],
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_edgeR_codon_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("edgeR_sensitivity_check"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/edgeR_sensitivity_check.R"


rule intersect_codon_deseq2_edgeR:
    input:
        deseq2_results = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_DESeq2_results_kappa{{kappa}}.tsv",
        edgeR_results   = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_edgeR_results_kappa{{kappa}}.tsv",
        counts          = f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/codon_counts_matrix_kappa{{kappa}}.tsv",
        coldata         = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        highconf     = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_highconf_intersect_kappa{{kappa}}.tsv",
        highconf_all = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_highconf_intersect_all_kappa{{kappa}}.tsv",
    params:
        fdr             = config["diff_abundance"]["fdr_threshold"],
        min_replicate_r = config["diff_abundance"]["min_replicate_r"],
    log:
        f"{STAGE2_ROOT}/logs/10_diff_abundance/{{cell_line}}_intersect_codon_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("edgeR_sensitivity_check"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/intersect_deseq2_edgeR.py"
