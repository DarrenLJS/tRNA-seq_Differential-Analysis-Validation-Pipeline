"""
workflow/scripts/fetch_watson_polysome_data.py

Parses Watson, Bellora & Macias (2020, NAR, doi:10.1093/nar/gkz1060)
polysome-profiling data into a per-gene translational log2FC table for
rule 16 validation.

CONFIRMED AGAINST THE REAL SUPPLEMENTARY FILE (2026-07-17) -- reads a
LOCAL xlsx file, no download attempted. The paper ships three
supplementary Excel files + a PDF (gkz1060_supplemental_files.zip from
the NAR article page). The one needed here is:

  "Supplementary Excel File 3.xlsx", sheet "Poly siMock + p(IC) vs siMock"

confirmed from that file's own "Suppl_table 3 Information" sheet to be
"Polysome associated RNA sequencing results" (i.e. the polysome-fraction
DE comparison, not total RNA) for exactly siMock+p(I:C) vs siMock --
matching Figure 4A / Supplementary Figure S3B, and NOT the siILF3
comparisons (File 2 is the total-RNA equivalent; ignored here per the
earlier scoping note that the polysome-fraction DE result IS the
translational-regulation readout, not a ratio to be computed).

TWO STRUCTURAL QUIRKS OF THIS SHEET, HANDLED BELOW
----------------------------------------------------
1. NOT ONE TABLE -- TWO SIDE-BY-SIDE BLOCKS. Row 1 is a summary header
   ("# N down-regulated" / "# M up-regulated"); row 2 has column headers
   repeated twice; columns A-D hold the down-regulated genes
   (Gene_id, AveExpr, l2FC, FDR), columns F-I hold the up-regulated genes
   (same 4 columns), column E is blank. The two blocks are different
   lengths (confirmed 1040 down / 1830 up rows, matching the row-1
   summary counts exactly) and must be parsed and concatenated
   separately, not read as one contiguous table. Confirmed zero gene_id
   overlap between the two blocks.
2. SIGNIFICANT GENES ONLY, NOT A FULL GENE UNIVERSE. Per the workbook's
   own "Suppl_table 3 Information" sheet, this table already is filtered
   to "Significantly differentially polysome associated genes (FDR<0.05)"
   -- it is NOT a full per-gene log2FC table with all tested genes. This
   is a real limitation worth a methods footnote: (a) the Spearman
   correlation in rule 16 is computed only over genes Watson et al.
   called significant, which will tend to inflate |rho| relative to a
   correlation against the full tested gene set (restriction to
   significant hits truncates the null-ish middle of the distribution);
   (b) the Fisher's exact test's background/denominator is these ~2870
   genes, not the whole transcriptome -- both tests should be reported
   with this caveat explicit, not silently treated as if this were an
   unfiltered DE table.

OUTPUT
------
watson_polysome_foldchange.tsv, columns: gene_id, log2FC, padj,
source ("NAR_supplementary_File3_polysome"). No timepoint column --
Watson et al.'s poly(I:C) stimulation is a single 4h timepoint (see
validate_fisher_spearman.R / kappa_sweep_summary.R for how rule 16
handles this against this pipeline's own 2h/4h/8h timepoints).
"""

import logging

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

SHEET_NAME = "Poly siMock + p(IC) vs siMock"
# Column layout within that sheet (0-indexed), confirmed against the real file:
#   down-regulated block: columns 0-3 = Gene_id, AveExpr, l2FC, FDR
#   blank separator:      column 4
#   up-regulated block:   columns 5-8 = Gene_id, AveExpr, l2FC, FDR
# Header row is index 1 (row 1 is the "# N down / # M up" summary line).
HEADER_ROW_IDX = 1
DOWN_COLS = [0, 1, 2, 3]
UP_COLS = [5, 6, 7, 8]
BLOCK_COLNAMES = ["gene_id", "AveExpr", "log2FC", "padj"]


def _parse_block(raw, cols, data_start_row):
    """Extract one side-by-side block (down- or up-regulated) from the
    raw (header=None) sheet, stopping at the first row where the block's
    gene_id column is empty."""
    id_col = cols[0]
    rows = []
    for row in raw[data_start_row:]:
        gene_id = row[id_col]
        # pandas gives NaN (float), not None, for empty cells when read via
        # pd.read_excel -- checking only `is None` silently fails to find
        # the block boundary and pulls in the padding rows below the
        # shorter block (confirmed: this bug made the down-regulated block
        # come out at 1830 rows instead of the correct 1040).
        if gene_id is None or (isinstance(gene_id, float) and pd.isna(gene_id)):
            break
        rows.append([row[c] for c in cols])
    df = pd.DataFrame(rows, columns=BLOCK_COLNAMES)
    return df


def parse_watson_supplementary(xlsx_path, sheet_name=SHEET_NAME):
    raw_df = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=None)
    raw = raw_df.values.tolist()

    down = _parse_block(raw, DOWN_COLS, HEADER_ROW_IDX + 1)
    up = _parse_block(raw, UP_COLS, HEADER_ROW_IDX + 1)

    log.info(f"Parsed down-regulated block: {len(down)} genes (all log2FC < 0 expected)")
    log.info(f"Parsed up-regulated block: {len(up)} genes (all log2FC > 0 expected)")

    n_down_wrong_sign = (down["log2FC"] >= 0).sum()
    n_up_wrong_sign = (up["log2FC"] <= 0).sum()
    if n_down_wrong_sign or n_up_wrong_sign:
        log.warning(
            f"Sign sanity check failed: {n_down_wrong_sign} 'down' rows with log2FC>=0, "
            f"{n_up_wrong_sign} 'up' rows with log2FC<=0 -- re-check DOWN_COLS/UP_COLS "
            "against the real sheet layout before trusting this output."
        )

    combined = pd.concat([down, up], ignore_index=True)
    overlap = set(down["gene_id"]) & set(up["gene_id"])
    if overlap:
        log.warning(
            f"{len(overlap)} gene_id(s) appear in BOTH the down- and up-regulated blocks "
            "-- unexpected, inspect the source sheet before trusting this output."
        )

    combined = combined[["gene_id", "log2FC", "padj"]].drop_duplicates(subset="gene_id")
    combined["source"] = "NAR_supplementary_File3_polysome"
    return combined


def fetch_watson_polysome_data(xlsx_path, out_path, sheet_name=SHEET_NAME):
    if not xlsx_path:
        raise RuntimeError(
            "config stage2_references.watson_nar_supp_path is empty. Stage "
            "'Supplementary Excel File 3.xlsx' from the NAR article's supplementary "
            "ZIP (gkz1060_supplemental_files.zip, from "
            "https://academic.oup.com/nar/article/48/1/116/5614571) locally and "
            "point stage2_references.watson_nar_supp_path at it."
        )
    result = parse_watson_supplementary(xlsx_path, sheet_name)

    result.to_csv(out_path, sep="\t", index=False)
    log.info(f"Wrote Watson et al. polysome data ({len(result)} rows) -> {out_path}")
    log.warning(
        "This table contains ONLY genes Watson et al. called significant "
        "(FDR<0.05) in the polysome-fraction siMock+p(I:C) vs siMock comparison -- "
        "not a full per-gene log2FC table. rule 16's Spearman correlation and "
        "Fisher's exact test are computed against this restricted gene set; report "
        "this as a methods limitation, not an unfiltered DE comparison."
    )
    return result


if __name__ == "__main__":
    fetch_watson_polysome_data(
        xlsx_path=snakemake.params.nar_supp_path,
        out_path=snakemake.output.polysome_fc,
        sheet_name=snakemake.params.sheet_name,
    )
