# =============================================================================
# workflow/rules/09_reference_prep_stage2.smk
# One-off reference builds specific to Stage 2. Nothing here re-derives or
# duplicates a Stage-1 reference -- these all consume Stage-1 outputs
# read-only (anticodon_map, GRCh38 genome/GTF) or fetch genuinely new
# external data (Watson et al. polysome data, ISG/housekeeping gene lists).
#
# Outputs:
#   decoding_whitelist.tsv           -- isodecoder x codon reachability +
#                                        term_type/gamma (rule 14 core input)
#   gene_codon_usage.tsv             -- per-gene 61-codon usage frequency
#   isg_housekeeping_gene_sets.tsv   -- curated gene ID sets for rule 16
#   watson_polysome_foldchange.tsv   -- per-gene translational log2FC, per
#                                        timepoint, from Watson et al. 2020
# =============================================================================

rule build_decoding_whitelist:
    """
    Build the explicit isodecoder -> reachable-codon whitelist (I34/Q34/
    default buckets) that Delta(c) is computed from. See
    workflow/scripts/build_decoding_whitelist.py docstring for the full
    derivation -- this is the single most load-bearing reference in Stage 2.

    FIXED -- now built per {cell_line}, not once globally, AND joined
    through the correct mim-tRNAseq output file. anticodon_map's raw
    per-locus IDs must be relabeled to Stage 1's FINAL, cell-line-specific
    collapsed isodecoder IDs (used by pos34_coverage_matrix.tsv,
    isodecoder_highconf_intersect.tsv, etc.) via
    {cell_line}_tRNAseq_unsplitClusterInfo.txt -- NOT
    {cell_line}_tRNAseqclusterInfo.txt, which reflects a different,
    earlier structural/alignment clustering pass unrelated to isodecoder
    identity (confirmed on real data: clusterInfo.txt grouped loci that
    Isodecoder_counts.txt does not actually merge). unsplitClusterInfo.txt
    is the genuine data-driven, per-cell-line coverage/mismatch-based
    collapsing (rows carry reasons like "insufficient coverage at
    mismatch X" / "potential mod at mismatch Y") -- confirmed A549 and
    THP1 disagree on 41+ loci, so a single shared whitelist cannot be
    correct for both. See build_decoding_whitelist.py module docstring
    "FIXED" / "CORRECTED" notes for the full derivation and worked
    real-data examples.
    """
    input:
        anticodon_map = config["references"]["anticodon_map"],
        unsplit_cluster_info = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/_run/annotation/{{cell_line}}_tRNAseq_unsplitClusterInfo.txt",
        isodecoder_counts = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isodecoder_counts.txt",
    output:
        whitelist = f"{STAGE2_ROOT}/references/{{cell_line}}/decoding_whitelist.tsv",
    params:
        i34_isotypes = config["wobble_glm"]["i34_isotypes"],
        q34_isotypes = config["wobble_glm"]["q34_isotypes"],
    log:
        f"{STAGE2_ROOT}/logs/09_reference_prep/{{cell_line}}_build_decoding_whitelist.log",
    resources:
        sge_extra = sge_extra("build_decoding_whitelist"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/build_decoding_whitelist.py"


rule build_codon_usage_table:
    """
    Per-gene codon usage frequency vector (61 sense codons) from GRCh38
    CDS sequences, used by rule 15 to dot-product against Delta(c).

    FIX-flag before first real run: needs a CDS FASTA, not just the genome
    + GTF -- calls gffread (or equivalent) to extract per-transcript CDS
    from config references.genome_fasta + references.ensembl_gtf, then
    picks one representative transcript per gene (longest CDS, matching
    the convention most codon-usage tools default to) before counting
    codons. Confirm gffread is available in envs/environment.yaml (it is
    not currently declared there -- add it if missing) before running.
    """
    input:
        genome_fasta = config["references"]["genome_fasta"],
        gtf          = config["references"]["ensembl_gtf"],
    output:
        codon_usage = f"{STAGE2_ROOT}/references/gene_codon_usage.tsv",
    params:
        outdir = f"{STAGE2_ROOT}/references",
    log:
        f"{STAGE2_ROOT}/logs/09_reference_prep/build_codon_usage_table.log",
    resources:
        sge_extra = sge_extra("build_codon_usage_table"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/build_codon_usage_table.py"


rule fetch_isg_housekeeping_lists:
    """
    Curated ISG (Interferome) and housekeeping (Eisenberg & Levanon 2013)
    gene ID lists for the rule-16 Fisher's exact enrichment test. Output
    is a single long-format TSV: gene_id, gene_set ("ISG" | "housekeeping").

    FIX-flag: Interferome requires either an account-gated bulk export or
    per-query scraping of their web interface -- this script currently
    targets Interferome's documented REST-ish query endpoint; confirm it
    still returns bulk results before trusting the output. If Interferome
    access is blocked/changed, fall back to a static, literature-cited ISG
    list (e.g. the ISG core list from Schoggins lab overexpression
    screens) -- flagged as a TODO in the script itself.
    """
    output:
        gene_sets = f"{STAGE2_ROOT}/references/isg_housekeeping_gene_sets.tsv",
    params:
        isg_source = config["stage2_references"]["isg_list_source"],
        hk_source  = config["stage2_references"]["housekeeping_list_source"],
    log:
        f"{STAGE2_ROOT}/logs/09_reference_prep/fetch_isg_housekeeping_lists.log",
    resources:
        sge_extra = sge_extra("reference_prep_stage2"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/fetch_isg_housekeeping_lists.py"


rule fetch_watson_polysome_data:
    """
    Fetch and parse Watson et al. (2020) polysome-profiling data -- per
    gene, per timepoint, translational log2 fold change (poly(I:C) vs
    control), used by rule 16 as the validation ground truth.

    Data source: GEO GSE130618 (superseries; polysome-seq + total RNA-seq
    subseries). Tries the GEO processed series matrix first; if that
    doesn't directly carry a gene-level translational-efficiency column
    (GEO superseries structure is not guaranteed to expose this cleanly),
    falls back to parsing the ~5.7MB supplementary ZIP hosted at NAR
    Online alongside the article, which most likely contains the authors'
    own DESeq2 output tables.

    NEEDS A SCOPING CHECK AGAINST REAL OUTPUT before the first full run --
    this is designed against GEO/NAR documentation, not a file in hand
    (see Snakefile docstring / project conversation). The NAR supplementary
    URL itself is not yet confirmed -- config stage2_references.watson_nar_supp_url
    is a placeholder; fill it in once located, or let the GEO path succeed
    on its own.
    """
    output:
        polysome_fc = f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",
    params:
        geo_accession = config["stage2_references"]["watson_geo_accession"],
        nar_supp_url  = config["stage2_references"]["watson_nar_supp_url"],
        outdir        = f"{STAGE2_ROOT}/references",
    log:
        f"{STAGE2_ROOT}/logs/09_reference_prep/fetch_watson_polysome_data.log",
    resources:
        sge_extra = sge_extra("reference_prep_stage2"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/fetch_watson_polysome_data.py"
