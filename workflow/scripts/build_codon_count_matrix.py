"""
workflow/scripts/build_codon_count_matrix.py

Expands the anticodon-level long table (build_anticodon_count_matrix.py's
anticodon_count_table_kappa{kappa}.tsv) into a codon-level wide count
matrix (rows=codon, cols=sample), for DESeq2/edgeR at the codon level
(rule deseq2_codon etc.) -- the endpoint that feeds compute_delta_c_v2.

No whitelist re-parsing needed here: each row of the anticodon long table
already carries its own read_codon list (comma-separated), pre-derived
from decoding_whitelist.tsv by build_anticodon_count_matrix.py. This
script's only job is the many-to-many expansion:

  - One read_anticodon can decode several codons (wobble) -> its count
    is DUPLICATED into every codon it reaches, not split/divided between
    them. This is intentional, not double-counting in the "same read
    counted twice" sense -- a codon's count is asking "how much decoding
    capacity for me exists", and the same charged tRNA population
    genuinely represents capacity for every codon type it can pair with.
    Exactly the same convention the existing whitelist-based Delta(c)
    (rule 14 / compute_delta_c.py) already uses at the isodecoder level;
    this just applies it one level up, at the anticodon level.
  - Several read_anticodons (canonical + I34/Q34 variants, or distinct
    anticodons in the same box) can reach the same codon -> counts SUM.
"""

import logging
import os

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)


def build_codon_count_matrix(anticodon_table_path, out_matrix_path, log_path=None):
    long_tbl = pd.read_csv(anticodon_table_path, sep="\t", dtype={"read_codon": str})

    long_tbl = long_tbl[long_tbl["read_codon"].notna() & (long_tbl["read_codon"] != "")]

    exploded = long_tbl.assign(
        codon=long_tbl["read_codon"].str.split(",")
    ).explode("codon")
    exploded["codon"] = exploded["codon"].str.strip()

    codon_sample = (
        exploded.groupby(["codon", "sample"])["total_count_isa"]
        .sum()
        .reset_index()
    )

    matrix = codon_sample.pivot_table(
        index="codon", columns="sample", values="total_count_isa", fill_value=0,
    )

    # Drop all-zero rows -- matches Stage 1's collate_counts.py convention
    # of dropping all-zero features before handing a matrix to DESeq2.
    all_zero = (matrix == 0).all(axis=1)
    if all_zero.any():
        log.info(f"Dropping {all_zero.sum()} all-zero codon rows")
    matrix = matrix[~all_zero].round().astype(int)

    matrix.to_csv(out_matrix_path, sep="\t")
    log.info(
        f"Wrote codon count matrix: {matrix.shape[0]} codons x "
        f"{matrix.shape[1]} samples -> {out_matrix_path}"
    )

    if log_path:
        with open(log_path, "a") as fh:
            fh.write(f"Codon matrix: {matrix.shape[0]} x {matrix.shape[1]} -> {out_matrix_path}\n")
            fh.write(f"Codons covered: {sorted(matrix.index.tolist())}\n")

    return matrix


if __name__ == "__main__":
    log_path = snakemake.log[0] if len(snakemake.log) else None
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()

    build_codon_count_matrix(
        anticodon_table_path=snakemake.input.anticodon_table,
        out_matrix_path=snakemake.output.codon_matrix,
        log_path=log_path,
    )
