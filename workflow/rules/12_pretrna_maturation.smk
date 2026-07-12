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
        # FIX (2026-07-XX): was expand(..., sample=SAMPLES) -- the GLOBAL
        # 30-sample list (manifest order: all 15 A549 rows, then all 15
        # THP1 rows), not this rule's own {cell_line} wildcard. Confirmed
        # on real data: for cell_line=THP1, featureCounts' own output
        # header showed A549 BAM paths in the first 15 sample-column
        # positions. rename_sample_cols() in pretrna_locus_counts.py then
        # zip()s those 30 columns against params.samples =
        # samples_for(wildcards.cell_line) (only 15, correctly THP1) --
        # zip() truncates to the shorter list, so THP1's locus_counts.tsv
        # was silently getting A549's mature/pre-tRNA counts relabeled as
        # THP1 sample IDs, while THP1's real BAM columns (positions 15-29,
        # never renamed) were dropped entirely by the later
        # melt(value_vars=samples) call. A549 happened to run correctly
        # only because it's first in the manifest, not because the logic
        # was right. Scoping these to samples_for(wildcards.cell_line) --
        # the same call already used for params.samples -- guarantees the
        # BAM list and the sample-ID list passed to rename_sample_cols are
        # always the same length, same order, same cell line.
        functional_bams = lambda wildcards: expand(
            f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam",
            sample=samples_for(wildcards.cell_line),
        ),
        pretrna_bams = lambda wildcards: expand(
            f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA.bam",
            sample=samples_for(wildcards.cell_line),
        ),
        pretrna_fasta = config["references"]["pretRNA_fasta_spliced"],
        # FIX (2026-07-XX): needed to redirect pre-tRNA loci belonging to
        # an isodecoder family that mim-tRNAseq's --cluster-id 0.97
        # collapsing absorbed into another family's cluster (and which
        # therefore has no contig of its own in the mature BAM) to the
        # Parent family that does. Same file rule 09's
        # build_decoding_whitelist already reads for the same underlying
        # reason -- see pretrna_locus_counts.py's
        # _load_absorbed_family_redirect for the real-data-confirmed
        # column format.
        unsplit_cluster_info = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/_run/annotation/{{cell_line}}_tRNAseq_unsplitClusterInfo.txt",
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
