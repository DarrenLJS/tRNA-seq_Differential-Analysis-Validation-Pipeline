# =============================================================================
# workflow/rules/13_trf_differential_abundance.smk
# tRF differential abundance -- proposal Section 3.8. Parses TRAX's
# per-cell-line count tables into 5'-tRF / 3'-tRF / i-tRF / tiRNA classes,
# DESeq2 per class, stimulated vs control, per timepoint. tiRNA class
# flagged specially (stress-induced translation-initiation inhibitor).
#
# This rule (along with 09's Watson/GEO fetch) is one of the two places
# explicitly flagged as designed against documentation rather than a real
# file in hand -- see Snakefile docstring. Expect a correction pass once
# run against real TRAX output, the same way Stage 1's 07_trax.smk needed
# several FIX iterations.
# =============================================================================

rule parse_trax_trf_classes:
    """
    Parse TRAX's per-cell-line output into 5'/3'/i-tRF/tiRNA class count
    matrices. TRAX's exact column/file naming for tRF-class labelling
    varies by version -- see parse_trax_tRF_classes.py docstring for the
    defensive multi-pattern parsing this uses instead of a single hardcoded
    filename.
    """
    input:
        # FIX (2026-07-XX): was pointed at -normalizedreadcounts.txt (TRAX's
        # DESeq2/CPM-normalized, fractional output) -- DESeq2 requires raw
        # integer counts (it does its own normalization internally), so
        # feeding it already-normalized floats failed every class/timepoint
        # with "some values in assay are not integers". -readcounts.txt is
        # TRAX's raw (unnormalized) counterpart, written unconditionally
        # alongside the normalized file by the same processsamples.py run --
        # already present on disk from the existing rule 07 run, no rerun
        # needed. Confirmed integer-valued on real data (awk sweep, every
        # row, both cell lines: A549 34540/34540, THP1 40340/40340) -- see
        # parse_trax_tRF_classes.py for the strict integer check applied
        # on read-in.
        trax_readcounts = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-readcounts.txt",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        class_matrices_dir = directory(f"{STAGE2_ROOT}/trf_diff_abundance/{{cell_line}}/class_matrices"),
    params:
        trf_classes = config["trf_diff_abundance"]["trf_classes"],
        min_unique_cov = config["trf_diff_abundance"]["min_unique_cov"],
    log:
        f"{STAGE2_ROOT}/logs/13_trf_diff_abundance/{{cell_line}}_parse_trax.log",
    resources:
        sge_extra = sge_extra("trf_diff_abundance"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/parse_trax_tRF_classes.py"


rule deseq2_trf:
    """
    DESeq2 differential abundance per tRF class, stimulated vs control,
    per timepoint. One combined output across all classes (class is a
    column, not a separate file per class) so downstream consumers don't
    need to know the class list in advance.
    """
    input:
        class_matrices_dir = f"{STAGE2_ROOT}/trf_diff_abundance/{{cell_line}}/class_matrices",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/trf_diff_abundance/{{cell_line}}/trf_DESeq2_results.tsv",
    params:
        fdr = config["trf_diff_abundance"]["fdr_threshold"],
        trf_classes = config["trf_diff_abundance"]["trf_classes"],
    log:
        f"{STAGE2_ROOT}/logs/13_trf_diff_abundance/{{cell_line}}_deseq2_trf.log",
    resources:
        sge_extra = sge_extra("deseq2_isodecoder"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/deseq2_trf.R"
