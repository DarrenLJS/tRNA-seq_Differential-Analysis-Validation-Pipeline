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
    (Pass-2 BAM), for every sample of a cell line. Uses the GtRNAdb BED
    (mature loci) and the pre-tRNA reference (pretRNA_fasta_spliced,
    treated as its own set of "loci" via a BED derived from the same
    FASTA headers) as featureCounts annotation.

    FIX-check before first run: pretRNA_fasta_spliced's headers need a
    BED-equivalent (whole-sequence-as-one-feature) annotation for
    featureCounts -- not currently guaranteed to exist as a BED anywhere
    in Stage 1's reference set. build_pretrna_locus_bed in
    pretrna_locus_counts.py derives one from the FASTA index (.fai) if a
    real BED isn't found; confirm coordinates make sense against a real
    pre-tRNA reference before trusting the counts.
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
        gtrndb_bed  = config["references"]["gtrndb_bed"],
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
