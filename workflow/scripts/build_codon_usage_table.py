"""
workflow/scripts/build_codon_usage_table.py

Builds a per-gene codon usage frequency table (61 sense codons) from
GRCh38 CDS sequences, for the rule-15 dot product against Delta(c).

APPROACH
--------
1. Extract per-transcript CDS sequences from the genome FASTA + GTF using
   gffread (must be on PATH -- add to envs/environment.yaml if missing;
   NOT currently declared there, see rule docstring).
2. Collapse to one representative transcript per gene: longest CDS,
   ties broken by transcript_id (deterministic), matching the convention
   most codon-usage tools default to when a "the gene's codon usage"
   single-vector answer is needed.
3. Count codons in frame (length must be a multiple of 3; sequences that
   aren't are logged and skipped rather than silently truncated).
4. Normalise counts to per-gene frequency (each gene's 61-codon vector
   sums to 1) so the rule-15 dot product with Delta(c) reflects codon
   COMPOSITION, not raw CDS length.

CAVEATS TO CHECK AT FIRST REAL RUN
-----------------------------------
- "Longest CDS per gene" is a reasonable default but not the only
  defensible choice (APPRIS principal isoform would be more rigorous if
  available in the GTF attributes) -- flagged here rather than silently
  assumed to be uncontroversial.
- Genes with multiple in-frame stop codons (readthrough, selenoprotein
  UGA-as-Sec) will have their internal stop codon excluded from the
  count automatically (STOP_CODONS_DNA filtered out below), which is the
  standard convention but worth a footnote if any selenoprotein genes are
  of specific interest downstream.
"""

import subprocess
import shutil
import logging
import gzip
import os
from collections import Counter

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

STOP_CODONS_DNA = {"TAA", "TAG", "TGA"}
SENSE_CODONS = sorted(
    a + b + c
    for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if a + b + c not in STOP_CODONS_DNA
)


def extract_cds_fasta(genome_fasta, gtf, out_fasta):
    if shutil.which("gffread") is None:
        raise RuntimeError(
            "gffread not found on PATH. Add `gffread` to envs/environment.yaml "
            "(bioconda channel) -- it is required by build_codon_usage_table.py "
            "and is not currently declared in Stage 1's shared environment file."
        )
    # gffread 0.12.9 does not reliably auto-decompress .gtf.gz input --
    # passing the gzipped file directly causes it to read raw compressed
    # bytes as text (silent corruption: "unexpected tab character" /
    # "invalid start coordinate" warnings, then a hard parse error).
    # Decompress to a plain-text temp GTF first and feed that to gffread.
    gtf_for_gffread = gtf
    tmp_gtf = None
    if gtf.endswith(".gz"):
        tmp_gtf = os.path.join(
            os.path.dirname(out_fasta) or ".",
            "_tmp_" + os.path.basename(gtf)[:-3],
        )
        log.info(f"Decompressing {gtf} -> {tmp_gtf} (gffread .gz support unreliable)")
        with gzip.open(gtf, "rb") as f_in, open(tmp_gtf, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gtf_for_gffread = tmp_gtf

    cmd = ["gffread", "-x", out_fasta, "-g", genome_fasta, gtf_for_gffread]
    log.info(f"Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    finally:
        if tmp_gtf is not None and os.path.exists(tmp_gtf):
            os.remove(tmp_gtf)


def parse_fasta(path):
    """Minimal FASTA parser -> dict[header] = sequence (str, uppercase)."""
    seqs = {}
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks).upper()
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks).upper()
    return seqs


def gene_id_from_transcript_header(header, gtf_tx2gene):
    """
    gffread's -x output headers are transcript IDs. Map back to gene_id
    via a transcript_id -> gene_id table parsed from the GTF (see
    `_build_tx2gene`). Transcripts absent from the map are skipped with a
    warning (should not happen for a well-formed GTF/CDS pair, but GTFs
    are not always perfectly clean).
    """
    return gtf_tx2gene.get(header)


def _build_tx2gene(gtf_path):
    """Parse a GTF for transcript_id -> gene_id, from any CDS/transcript line."""
    import gzip
    opener = gzip.open if gtf_path.endswith(".gz") else open
    tx2gene = {}
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "CDS"):
                continue
            attrs = fields[8]
            tx_id, gene_id = None, None
            for kv in attrs.strip().split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                if kv.startswith("transcript_id"):
                    tx_id = kv.split('"')[1] if '"' in kv else kv.split()[-1]
                elif kv.startswith("gene_id"):
                    gene_id = kv.split('"')[1] if '"' in kv else kv.split()[-1]
            if tx_id and gene_id:
                tx2gene[tx_id] = gene_id
    return tx2gene


def build_codon_usage(genome_fasta, gtf, out_path, tmp_fasta):
    extract_cds_fasta(genome_fasta, gtf, tmp_fasta)
    cds_seqs = parse_fasta(tmp_fasta)
    log.info(f"Extracted {len(cds_seqs)} CDS transcript sequences")

    tx2gene = _build_tx2gene(gtf)

    # Group by gene, keep longest CDS per gene
    gene_best = {}   # gene_id -> (length, seq)
    n_skipped_no_gene = 0
    n_skipped_not_multiple_of_3 = 0
    for tx_id, seq in cds_seqs.items():
        gene_id = tx2gene.get(tx_id)
        if gene_id is None:
            n_skipped_no_gene += 1
            continue
        if len(seq) % 3 != 0:
            n_skipped_not_multiple_of_3 += 1
            continue
        if gene_id not in gene_best or len(seq) > gene_best[gene_id][0]:
            gene_best[gene_id] = (len(seq), seq)

    log.info(f"Genes with usable CDS: {len(gene_best)}")
    if n_skipped_no_gene:
        log.warning(f"Skipped {n_skipped_no_gene} transcripts with no gene_id mapping")
    if n_skipped_not_multiple_of_3:
        log.warning(f"Skipped {n_skipped_not_multiple_of_3} transcripts with CDS length not a multiple of 3")

    rows = []
    for gene_id, (length, seq) in gene_best.items():
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        counts = Counter(c for c in codons if c in SENSE_CODONS)  # drop stop codon(s), any Ns
        total = sum(counts.values())
        if total == 0:
            continue
        row = {"gene_id": gene_id}
        row.update({c: counts.get(c, 0) / total for c in SENSE_CODONS})
        rows.append(row)

    usage = pd.DataFrame(rows).set_index("gene_id")
    usage.to_csv(out_path, sep="\t")
    log.info(f"Wrote codon usage table: {usage.shape[0]} genes x {usage.shape[1]} codons -> {out_path}")
    return usage


if __name__ == "__main__":
    build_codon_usage(
        genome_fasta=snakemake.input.genome_fasta,
        gtf=snakemake.input.gtf,
        out_path=snakemake.output.codon_usage,
        tmp_fasta=f"{snakemake.params.outdir}/_tmp_cds.fa",
    )
