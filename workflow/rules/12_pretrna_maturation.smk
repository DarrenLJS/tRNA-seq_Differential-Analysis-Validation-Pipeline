# =============================================================================
# workflow/rules/12_pretrna_maturation.smk
# Pre-tRNA:mature tRNA ratio -- proposal Section 3.7. The one place Stage 2
# needs genuinely new quantification: Stage 1's isodecoder counts are
# 97%-identity clusters, but this analysis needs PER-LOCUS resolution.
# featureCounts against Pass-1 functional BAMs (mature) and Pass-2 pre-tRNA
# BAMs (pre-tRNA), per locus, then a linear model. Reported as
# supplementary/exploratory per the proposal's own framing.
# =============================================================================

rule pretrna_locus_counts:
    """
    Per-locus read counts, mature (Pass-1 functional BAM) vs pre-tRNA
    (Pass-2 BAM), for every sample of a cell line.

    Mature-locus annotation is derived directly from functional_bams'
    own BAM header (one whole-length SAF feature per @SQ contig), NOT
    from the GtRNAdb genomic BED (config references.gtrndb_bed) -- that
    BED describes hg38 genomic coordinates, but functional.bam is aligned
    to a per-sequence mature-tRNA reference (mim-tRNAseq-style: one
    contig per tRNA locus), so the two coordinate spaces never overlap.
    Confirmed via samtools view -H against a real functional.bam: using
    gtrndb_bed silently produced 0.0% assigned alignments for every
    locus, every sample (featureCounts exits cleanly even when nothing
    is assigned). See pretrna_locus_counts.py docstring for detail.

    Pre-tRNA loci use pretRNA_fasta_spliced, treated as its own set of
    "loci" via a SAF derived from the same FASTA headers (each FASTA
    entry = one whole-sequence feature).

    pretRNA.bam is genuinely paired-end (confirmed via samtools flagstat)
    while functional.bam is single-end -- featureCounts is run with
    -p --countReadPairs for the pre-tRNA pass only.
    """
    input:
        functional_bams = expand(
            f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam",
            sample=SAMPLES,
        ),
        pretrna_bams = expand(
            f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA.bam",
            sample=SAMPLES,
        ),
        pretrna_fasta = config["references"]["pretRNA_fasta_spliced"],
    output:
        locus_counts = f"{STAGE2_ROOT}/pretrna_ratio/{{cell_line}}/locus_counts.tsv",
    params:
        cell_line = "{cell_line}",
        samples   = lambda wc: samples_for(wc.cell_line),
        scratch   = SCRATCH,
        min_coverage = config["pretrna_ratio"]["min_locus_coverage"],
    log:
        f"{STAGE2_ROOT}/logs/12_pretrna_maturation/{{cell_line}}_locus_counts.log",
    resources:
        sge_extra = sge_extra("pretrna_locus_counts"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/pretrna_locus_counts.py"


rule pretrna_ratio_lm:
    """
    Linear model (pre-tRNA:mature ratio ~ condition, per locus, per
    timepoint) on the per-locus counts from pretrna_locus_counts.
    Supplementary/exploratory output per proposal 3.7's own framing --
    NOT gated into the high-confidence isodecoder set rule 14 consumes.
    """
    input:
        locus_counts = f"{STAGE2_ROOT}/pretrna_ratio/{{cell_line}}/locus_counts.tsv",
        coldata = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    output:
        results = f"{STAGE2_ROOT}/pretrna_ratio/{{cell_line}}/pretrna_mature_ratio_lm.tsv",
    log:
        f"{STAGE2_ROOT}/logs/12_pretrna_maturation/{{cell_line}}_ratio_lm.log",
    resources:
        sge_extra = sge_extra("pretrna_ratio_lm"),
    conda:
        "../../envs/r_stats.yaml"
    script:
        "../scripts/pretrna_ratio_lm.R"
