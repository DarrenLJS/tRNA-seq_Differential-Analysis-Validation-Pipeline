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
    Curated ISG (Interferome v2.0) and housekeeping (Eisenberg & Levanon
    2013) gene ID lists for the rule-16 Fisher's exact enrichment test.
    Output is a single long-format TSV: gene_id, gene_set ("ISG" |
    "housekeeping").

    Neither source has a documented bulk-download API -- both are staged
    manually as local export files (see fetch_isg_housekeeping_lists.py
    docstring for the exact download procedure for each). Config
    stage2_references.isg_list_path / housekeeping_list_path point at
    those staged files; this rule fails loudly with instructions if
    they're not there yet, rather than silently substituting a small
    placeholder list mislabelled as either source.
    """
    output:
        gene_sets = f"{STAGE2_ROOT}/references/isg_housekeeping_gene_sets.tsv",
    params:
        isg_list_path          = config["stage2_references"]["isg_list_path"],
        housekeeping_list_path = config["stage2_references"]["housekeeping_list_path"],
        isg_source             = config["stage2_references"]["isg_list_source"],
        housekeeping_source    = config["stage2_references"]["housekeeping_list_source"],
        isg_gene_col           = config["stage2_references"].get("isg_gene_col"),
        housekeeping_gene_col  = config["stage2_references"].get("housekeeping_gene_col"),
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
    Parse Watson, Bellora & Macias (2020, NAR) polysome-profiling data --
    per-gene translational log2 fold change (poly(I:C) vs control), used
    by rule 16 as the validation ground truth.

    CONFIRMED against the real supplementary file (2026-07-17): reads a
    LOCAL copy of "Supplementary Excel File 3.xlsx" (from
    gkz1060_supplemental_files.zip on the NAR article page), sheet
    "Poly siMock + p(IC) vs siMock" -- confirmed to be the polysome-
    fraction siMock+p(I:C) vs siMock comparison (Figure 4A / Suppl Fig
    S3B), not the siILF3 comparisons and not the total-RNA file (File 2).
    No download attempted -- GSE130618 provides raw SRA reads only, and
    NAR serves the ZIP behind a signed/expiring CDN URL, so this is a
    one-time manual staging step (see script docstring).

    Two real structural quirks of the source sheet, both handled in the
    script (see its docstring for detail): (1) it's two side-by-side
    down-/up-regulated blocks, not one contiguous table; (2) it only
    contains genes Watson et al. called significant (FDR<0.05), not a
    full per-gene log2FC table -- report this as a methods limitation.

    This is a SINGLE 4h poly(I:C) timepoint in the source data, not a
    timecourse -- output has no timepoint column. Rule 16's
    validate_fisher_spearman.R / kappa_sweep_summary.R validate each
    tRNA-seq timepoint's Delta(c) prediction against this one external
    benchmark (join on gene_id only), rather than looping over a
    timepoint dimension that doesn't exist in the source data.
    """
    output:
        polysome_fc = f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",
    params:
        nar_supp_path = config["stage2_references"]["watson_nar_supp_path"],
        sheet_name    = config["stage2_references"]["watson_nar_supp_sheet"],
        outdir        = f"{STAGE2_ROOT}/references",
    log:
        f"{STAGE2_ROOT}/logs/09_reference_prep/fetch_watson_polysome_data.log",
    resources:
        sge_extra = sge_extra("reference_prep_stage2"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/fetch_watson_polysome_data.py"
