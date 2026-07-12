"""
workflow/scripts/pretrna_locus_counts.py

Per-locus read counts, mature (Pass-1 functional BAM) vs pre-tRNA (Pass-2
BAM), for every sample of a cell line -- output feeds rule
pretrna_ratio_lm's linear model.

Uses featureCounts (subread) rather than a hand-rolled pysam counter,
since it's the same tool class Stage 1 already depends on conceptually
(samtools/bedtools family) and handles multi-mapping/strand options
correctly out of the box. Requires `featureCounts` on PATH -- part of the
`subread` bioconda package; NOT currently declared in Stage 1's
environment.yaml, add it if missing.

FIXED (was: FIX before first real run)
---------------------------------------
Originally this built the mature-locus SAF from `gtrndb_bed`
(hg38-tRNAs_nochr.bed), which is a GENOMIC hg38 coordinate BED. But
Pass-1's `functional.bam` is NOT aligned to the hg38 genome -- its @SQ
header shows one contig per individual mature tRNA sequence (e.g.
`Homo_sapiens_tRNA-Ala-TGC-1-1`, length ~72-76bp), i.e. mim-tRNAseq-style
mature-tRNA-space alignment. Genomic BED coordinates can never overlap
that BAM's contigs -- featureCounts silently "succeeded" while assigning
0.0% of every alignment, for every sample, every locus. Confirmed via
`samtools view -H` on a real functional.bam vs the SAF actually produced.

Fix: derive the mature-locus SAF directly from the functional BAM's own
header (`_derive_saf_from_bam_header`), the same whole-sequence-as-one-
-feature approach `_derive_saf_from_fasta_index` already used correctly
for the pre-tRNA side. This guarantees the annotation always matches
whatever space the BAM was actually aligned to, and removes the need for
`gtrndb_bed` as an input to this rule entirely (it was never the right
annotation source for functional.bam; if a true GtRNAdb genomic-BED-based
count is wanted elsewhere, that's a separate rule against a
genome-aligned BAM, not this one).

Also fixed: Pass-2 `pretRNA.bam` files are genuinely paired-end
(confirmed via `samtools flagstat` -- ~117M paired-in-sequencing reads,
58.7M read1/58.7M read2), while Pass-1 `functional.bam` files are
single-end (0 paired-in-sequencing). featureCounts was being run in
single-end mode for both, which made it abort on the paired BAMs
("Paired-end reads were detected in single-end read library"). Now passed
`-p --countReadPairs` for the pre-tRNA (paired) pass only.
"""

import logging
import os
import re
import subprocess
import shutil

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def _derive_saf_from_fasta_index(fasta_path, out_saf):
    """
    Build a minimal SAF (Simplified Annotation Format) treating each FASTA
    sequence as one whole-length feature, from samtools faidx's .fai index
    (built if not already present).
    """
    fai_path = fasta_path + ".fai"
    if not os.path.exists(fai_path):
        if shutil.which("samtools") is None:
            raise RuntimeError("samtools not found on PATH -- required to index the pre-tRNA FASTA.")
        subprocess.run(["samtools", "faidx", fasta_path], check=True)

    with open(fai_path) as fh, open(out_saf, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            seq_id, length = fields[0], int(fields[1])
            out.write(f"{seq_id}\t{seq_id}\t1\t{length}\t+\n")
    return out_saf


def _derive_saf_from_bam_header(bam_path, out_saf):
    """
    Build a minimal SAF treating each BAM reference contig (@SQ SN) as one
    whole-length feature. Used for mature tRNA counting because
    functional.bam is aligned to a mature-tRNA-sequence reference (one
    contig per locus, e.g. mim-tRNAseq-style), NOT the hg38 genome --
    genomic BED coordinates (gtrndb_bed) do not correspond to this BAM's
    contig space at all, which previously produced 0 assigned alignments
    for every locus/sample silently (featureCounts exits 0 even when it
    assigns nothing).
    """
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found on PATH -- required to read the BAM header.")
    header = subprocess.run(
        ["samtools", "view", "-H", bam_path], check=True, capture_output=True, text=True
    ).stdout

    with open(out_saf, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in header.splitlines():
            if not line.startswith("@SQ"):
                continue
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            seq_id, length = fields["SN"], int(fields["LN"])
            out.write(f"{seq_id}\t{seq_id}\t1\t{length}\t+\n")
    return out_saf


def _bed_to_saf(bed_path, out_saf):
    with open(bed_path) as fh, open(out_saf, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}"
            strand = fields[5] if len(fields) > 5 else "+"
            out.write(f"{name}\t{chrom}\t{start + 1}\t{end}\t{strand}\n")  # BED is 0-based, SAF is 1-based
    return out_saf


def run_featurecounts(bams, saf_path, out_path, threads=4, paired=False):
    if shutil.which("featureCounts") is None:
        raise RuntimeError(
            "featureCounts not found on PATH. Add `subread` to envs/environment.yaml "
            "(bioconda channel) -- required by pretrna_locus_counts.py."
        )
    cmd = [
        "featureCounts", "-F", "SAF", "-a", saf_path,
        "-o", out_path, "-T", str(threads),
        "-M", "--fraction",  # count multimappers fractionally -- pre-tRNA loci
                              # are highly homologous, dropping multimappers
                              # entirely would bias counts low
    ]
    if paired:
        # pretRNA.bam is genuinely paired-end (confirmed via samtools
        # flagstat) -- without this, featureCounts detects PE reads in
        # what it assumes is an SE library and aborts (exit 255).
        cmd += ["-p", "--countReadPairs"]
    cmd += bams
    log.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def _strip_copy_suffix(raw_locus):
    """
    Strip the trailing '-<copy_number>' segment from a locus/family string
    (e.g. "Homo_sapiens_tRNA-Val-CAC-11-1" -> family
    "Homo_sapiens_tRNA-Val-CAC-11"). Same convention as
    build_decoding_whitelist.py's _strip_copy_suffix -- kept as a separate
    copy here (not imported) since this module and build_decoding_whitelist
    each already declare their own logging/config setup and importing
    across sibling Snakemake `script:` entrypoints is fragile; if these
    drift, that's a sign they should be factored into a shared module.
    """
    return re.sub(r"-\d+$", "", str(raw_locus))


def _normalize_pretrna_locus_id(raw_header):
    """
    FIX (2026-07-XX): bedtools getfasta -name (rule 00b, Stage 1) appends
    '::chrom:start-end(strand)' to every header in this pipeline's bedtools
    version when -name is combined with -s (confirmed on real data: 100%
    of hg38_pretRNA_spliced.fa headers carry this suffix, not just
    disambiguated duplicates). GtRNAdb's own locus names also lack the
    'Homo_sapiens_' prefix that mim-tRNAseq's BAM contigs carry (e.g. BED
    name "tRNA-Val-CAC-11-1" vs BAM contig "Homo_sapiens_tRNA-Val-CAC-11-1").
    Strip the coordinate suffix and add the prefix so pre-tRNA locus IDs
    are directly comparable to mature BAM contig names before any
    family-level merge-cluster reconciliation happens.
    """
    bare = str(raw_header).split("::", 1)[0]
    return f"Homo_sapiens_{bare}"


def _load_family_to_bam_contig(bam_path):
    """
    Ground-truth family -> real BAM contig name map, built directly from
    the mature functional.bam header (NOT from Isodecoder_counts.txt's
    isodecoder/parent columns -- those are family-level strings with no
    copy number, e.g. "Homo_sapiens_tRNA-Val-CAC-11", which do NOT match
    the real BAM contig "Homo_sapiens_tRNA-Val-CAC-11-1"; confirmed on
    real A549 data). mim-tRNAseq keeps exactly one representative contig
    per isodecoder cluster (its own family + whichever copy number it
    happened to retain) -- stripping that contig's own copy suffix
    recovers its family, giving a direct, unambiguous family->contig map
    with no intermediate naming convention to get wrong.
    """
    header = subprocess.run(
        ["samtools", "view", "-H", bam_path], check=True, capture_output=True, text=True
    ).stdout
    mapping = {}
    collisions = []
    for line in header.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
        contig = fields["SN"]
        family = _strip_copy_suffix(contig)
        if family in mapping and mapping[family] != contig:
            collisions.append((family, mapping[family], contig))
        mapping[family] = contig
    if collisions:
        log.warning(f"{len(collisions)} family->contig collision(s) in BAM header "
                     f"(more than one contig strips to the same family) -- first few: {collisions[:5]}")
    return mapping


def _load_absorbed_family_redirect(unsplit_cluster_info_path):
    """
    Maps an ABSORBED family (one merged into another isodecoder's cluster
    during mim-tRNAseq's --cluster-id 0.97 collapsing, and therefore with
    NO contig of its own in the mature BAM) to the Parent family that
    holds the real, shared contig.

    unsplitClusterInfo.txt columns: Parent (raw Family-Copy locus, full
    "Homo_sapiens_tRNA-..." prefix), "Unsplit transcripts" (comma-sep,
    FAMILY-ONLY names, same prefix, no copy suffix -- confirmed against
    the real A549_tRNAseq_unsplitClusterInfo.txt, e.g. Parent=
    "Homo_sapiens_tRNA-Ile-AAT-5-4", Unsplit="Homo_sapiens_tRNA-Ile-AAT-7,
    Homo_sapiens_tRNA-Ile-AAT-3, Homo_sapiens_tRNA-Ile-AAT-8,
    Homo_sapiens_tRNA-Ile-AAT-1").
    """
    unsplit = pd.read_csv(unsplit_cluster_info_path, sep="\t")
    required = {"Parent", "Unsplit transcripts"}
    if not required.issubset(unsplit.columns):
        raise ValueError(f"{unsplit_cluster_info_path}: expected columns {required}, "
                          f"found {list(unsplit.columns)}.")

    redirect = {}
    for _, row in unsplit.iterrows():
        parent_family = _strip_copy_suffix(row["Parent"])
        for member in str(row["Unsplit transcripts"]).split(","):
            member_family = _strip_copy_suffix(member.strip())
            if member_family:
                redirect[member_family] = parent_family
    return redirect


def _resolve_pretrna_locus_ids(pretrna_long, family_to_contig, absorbed_redirect):
    """
    Reconciles raw pre-tRNA locus IDs (post _normalize_pretrna_locus_id)
    against the real mature BAM contig namespace, then sums pretrna_count
    across every raw locus that resolves to the same final contig
    (this is where multiple raw loci belonging to one mim-tRNAseq
    isodecoder cluster get correctly combined, instead of only the
    Parent's own raw count silently standing in for the whole cluster).

    Any locus whose family neither has its own contig nor a redirect
    entry is dropped (logged with a count and example names) rather than
    silently producing a 0/NaN join -- this is the explicit version of
    what the old direct-name join was doing accidentally and silently
    for every single locus.
    """
    df = pretrna_long.copy()
    df["family"] = df["locus_id"].apply(_strip_copy_suffix)

    def resolve(family):
        if family in family_to_contig:
            return family_to_contig[family]
        redirected = absorbed_redirect.get(family)
        if redirected is not None and redirected in family_to_contig:
            return family_to_contig[redirected]
        return None

    df["resolved_locus_id"] = df["family"].apply(resolve)

    n_total = df["locus_id"].nunique()
    unresolved = df[df["resolved_locus_id"].isna()]
    n_unresolved = unresolved["locus_id"].nunique()
    if n_unresolved:
        log.warning(
            f"{n_unresolved}/{n_total} raw pre-tRNA loci did not resolve to any mature "
            f"BAM contig (dropped) -- expected for GtRNAdb low-confidence/predicted loci "
            f"(tRX-*, Und-NNN-* families) that mim-tRNAseq's reference build excludes; "
            f"confirmed on real data these do not appear in functional.bam at all. "
            f"Examples: {sorted(unresolved['locus_id'].unique())[:10]}"
        )

    df = df[df["resolved_locus_id"].notna()].copy()
    df["locus_id"] = df["resolved_locus_id"]
    df = df.groupby(["locus_id", "sample_id"], as_index=False)["pretrna_count"].sum()
    return df


def build_locus_counts(functional_bams, pretrna_bams, pretrna_fasta, unsplit_cluster_info,
                        samples, min_coverage, out_path, workdir):
    os.makedirs(workdir, exist_ok=True)

    # Mature SAF comes from functional_bams' own header, not gtrndb_bed --
    # see module docstring "FIXED" note. All functional_bams share the
    # same reference (same Pass-1 alignment index), so the first one's
    # header is representative of all of them.
    mature_saf = _derive_saf_from_bam_header(functional_bams[0], os.path.join(workdir, "mature_loci.saf"))
    # NOTE: pretrna_saf still uses the LITERAL FASTA headers (including
    # the '::coords' suffix) -- this SAF is matched against pretrna_bams'
    # real @SQ names by featureCounts, and those BAM contigs were built by
    # bowtie2-build directly from this same FASTA's literal headers (rule
    # 00c, Stage 1), so it MUST stay unmodified here. Locus-ID
    # normalization and reconciliation against the mature namespace happens
    # afterwards, on the featureCounts *output* -- see
    # _resolve_pretrna_locus_ids below.
    pretrna_saf = _derive_saf_from_fasta_index(pretrna_fasta, os.path.join(workdir, "pretrna_loci.saf"))

    mature_counts_raw = os.path.join(workdir, "mature_featurecounts.txt")
    pretrna_counts_raw = os.path.join(workdir, "pretrna_featurecounts.txt")

    run_featurecounts(functional_bams, mature_saf, mature_counts_raw, paired=False)
    run_featurecounts(pretrna_bams, pretrna_saf, pretrna_counts_raw, paired=True)

    mature = pd.read_csv(mature_counts_raw, sep="\t", comment="#")
    pretrna = pd.read_csv(pretrna_counts_raw, sep="\t", comment="#")

    # featureCounts sample columns are the BAM paths -- rename to sample_id
    # via the order `samples` was passed in (matches expand() order in the
    # .smk rule, which iterates SAMPLES in manifest order for both BAM lists).
    def rename_sample_cols(df, bam_paths, samples):
        bam_cols = df.columns[6:]  # first 6 cols are GeneID/Chr/Start/End/Strand/Length
        mapping = dict(zip(bam_cols, samples))
        return df.rename(columns=mapping)

    mature = rename_sample_cols(mature, functional_bams, samples)
    pretrna = rename_sample_cols(pretrna, pretrna_bams, samples)

    mature_long = mature.melt(id_vars=["Geneid"], value_vars=samples,
                               var_name="sample_id", value_name="mature_count").rename(columns={"Geneid": "locus_id"})
    pretrna_long = pretrna.melt(id_vars=["Geneid"], value_vars=samples,
                                 var_name="sample_id", value_name="pretrna_count").rename(columns={"Geneid": "locus_id"})

    # FIX (2026-07-XX): the old direct name join was confirmed on real data
    # to silently produce pretrna_count=0 for essentially every row -- two
    # stacked causes: (1) bedtools getfasta -name appended '::chrom:start-
    # end(strand)' to every pre-tRNA header (not just disambiguated
    # duplicates) and lacked the 'Homo_sapiens_' prefix the mature BAM
    # contigs carry, and (2) mim-tRNAseq's --cluster-id 0.97 collapsing
    # means several raw pre-tRNA loci can share ONE mature contig (646 raw
    # pre-tRNA loci vs 329 mature BAM contigs on A549, confirmed), so even
    # a corrected direct name join would still miss every absorbed
    # (non-Parent) cluster member. _resolve_pretrna_locus_ids normalizes
    # names, redirects absorbed families to their Parent's real contig,
    # and sums pretrna_count across every raw locus that maps to the same
    # final contig, dropping (with a logged count) any locus that
    # genuinely doesn't exist in mim-tRNAseq's reference (GtRNAdb's
    # low-confidence tRX-*/Und-NNN-* loci -- confirmed absent from the
    # real functional.bam header entirely).
    pretrna_long["locus_id"] = pretrna_long["locus_id"].apply(_normalize_pretrna_locus_id)
    family_to_contig = _load_family_to_bam_contig(functional_bams[0])
    absorbed_redirect = _load_absorbed_family_redirect(unsplit_cluster_info)
    pretrna_long = _resolve_pretrna_locus_ids(pretrna_long, family_to_contig, absorbed_redirect)

    merged = mature_long.merge(pretrna_long, on=["locus_id", "sample_id"], how="outer")
    merged["mature_count"] = merged["mature_count"].fillna(0)
    merged["pretrna_count"] = merged["pretrna_count"].fillna(0)
    merged["total_coverage"] = merged["mature_count"] + merged["pretrna_count"]
    merged = merged[merged["total_coverage"] >= min_coverage].copy()
    merged["pretrna_mature_ratio"] = merged["pretrna_count"] / merged["mature_count"].replace(0, pd.NA)

    merged.to_csv(out_path, sep="\t", index=False)
    log.info(f"Wrote {len(merged)} locus x sample rows (min_coverage={min_coverage}) -> {out_path}")
    return merged


if __name__ == "__main__":
    build_locus_counts(
        functional_bams=list(snakemake.input.functional_bams),
        pretrna_bams=list(snakemake.input.pretrna_bams),
        pretrna_fasta=snakemake.input.pretrna_fasta,
        unsplit_cluster_info=snakemake.input.unsplit_cluster_info,
        samples=list(snakemake.params.samples),
        min_coverage=snakemake.params.min_coverage,
        out_path=snakemake.output.locus_counts,
        workdir=os.path.join(snakemake.params.scratch, "stage2", "pretrna_ratio", snakemake.params.cell_line, "_work"),
    )
