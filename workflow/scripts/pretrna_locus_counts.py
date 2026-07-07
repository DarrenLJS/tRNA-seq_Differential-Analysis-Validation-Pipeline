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

FIX before first real run
--------------------------
`pretRNA_fasta_spliced` (the Pass-2 reference) is a FASTA, and Stage 1's
config does not confirm a matching BED/SAF annotation exists for it.
`_derive_bed_from_fasta_index` below builds a whole-sequence-as-one-
feature SAF from the FASTA's .fai index as a fallback -- this treats each
pre-tRNA reference sequence as a single countable locus, which is
reasonable IF each FASTA entry already corresponds to one tRNA locus
(likely, given how Stage 1 built this reference), but confirm against the
real file rather than assuming.
"""

import logging
import os
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


def run_featurecounts(bams, saf_path, out_path, threads=4):
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
    ] + bams
    log.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def build_locus_counts(functional_bams, pretrna_bams, gtrndb_bed, pretrna_fasta,
                        samples, min_coverage, out_path, workdir):
    os.makedirs(workdir, exist_ok=True)

    mature_saf = _bed_to_saf(gtrndb_bed, os.path.join(workdir, "mature_loci.saf"))
    pretrna_saf = _derive_saf_from_fasta_index(pretrna_fasta, os.path.join(workdir, "pretrna_loci.saf"))

    mature_counts_raw = os.path.join(workdir, "mature_featurecounts.txt")
    pretrna_counts_raw = os.path.join(workdir, "pretrna_featurecounts.txt")

    run_featurecounts(functional_bams, mature_saf, mature_counts_raw)
    run_featurecounts(pretrna_bams, pretrna_saf, pretrna_counts_raw)

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

    # NOTE: mature loci (GtRNAdb BED) and pre-tRNA loci (spliced pre-tRNA
    # FASTA headers) may not share identical locus_id naming conventions --
    # FIX-check this join produces sensible overlap on a real file; if
    # naming differs, add an explicit locus_id mapping table rather than
    # relying on an outer join silently producing all-NaN pretrna_counts.
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
        gtrndb_bed=snakemake.input.gtrndb_bed,
        pretrna_fasta=snakemake.input.pretrna_fasta,
        samples=list(snakemake.params.samples),
        min_coverage=snakemake.params.min_coverage,
        out_path=snakemake.output.locus_counts,
        workdir=os.path.join(snakemake.params.scratch, "stage2", "pretrna_ratio", snakemake.params.cell_line, "_work"),
    )
