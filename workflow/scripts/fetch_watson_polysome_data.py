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
watson_polysome_foldchange.tsv, columns: gene_id, gene_symbol, log2FC,
padj, source ("NAR_supplementary_File3_polysome"). No timepoint column --
Watson et al.'s poly(I:C) stimulation is a single 4h timepoint (see
validate_fisher_spearman.R / kappa_sweep_summary.R for how rule 16
handles this against this pipeline's own 2h/4h/8h timepoints).

GENE ID FIX (2026-07-18)
-------------------------
The source sheet's "Gene_id" column is actually HGNC gene SYMBOLS
(e.g. AASDHPPT, ABCA8), not Ensembl gene IDs, despite the column name.
Confirmed via rule 16 returning n=0 overlapping genes at every timepoint
when joined directly against gene_translation_scores_kappa*.tsv (which
is keyed on Ensembl gene IDs, e.g. ENSG00000189060). This script now
builds a symbol -> Ensembl gene_id map from references.ensembl_gtf (the
same GTF build_codon_usage_table.py already parses) and maps the Watson
symbols onto it before writing output, so gene_id here means the same
thing it means everywhere else in this pipeline. The original symbol is
kept as gene_symbol for traceability. Symbols that don't map (not found
in the GTF, or ambiguous) are dropped and counted/logged rather than
silently included with a null gene_id.
"""

import logging
import gzip
import re

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


def build_symbol_to_ensembl_map(gtf_path):
    """Parse a gzipped Ensembl GTF's 'gene' feature lines into a
    gene_name (symbol) -> gene_id (Ensembl, version stripped) map.

    Ambiguous symbols (same gene_name mapping to >1 distinct gene_id)
    are dropped from the map entirely rather than resolved arbitrarily,
    and counted/logged -- silently picking one Ensembl ID for an
    ambiguous symbol would be a bigger correctness risk than dropping
    those few genes from the validation.
    """
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')

    symbol_to_ids = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            id_match = gene_id_re.search(attrs)
            name_match = gene_name_re.search(attrs)
            if not id_match or not name_match:
                continue
            gene_id = id_match.group(1).split(".")[0]  # strip version suffix
            gene_name = name_match.group(1)
            symbol_to_ids.setdefault(gene_name, set()).add(gene_id)

    ambiguous = {sym: ids for sym, ids in symbol_to_ids.items() if len(ids) > 1}
    if ambiguous:
        log.warning(
            f"{len(ambiguous)} gene symbol(s) in the GTF map to more than one "
            "Ensembl gene_id -- these are dropped from the symbol->Ensembl map "
            "rather than resolved arbitrarily."
        )

    symbol_to_ensembl = {
        sym: next(iter(ids)) for sym, ids in symbol_to_ids.items() if len(ids) == 1
    }
    log.info(
        f"Built symbol->Ensembl map from {gtf_path}: "
        f"{len(symbol_to_ensembl)} unambiguous symbols "
        f"({len(ambiguous)} ambiguous symbols dropped)"
    )
    return symbol_to_ensembl


def parse_watson_supplementary(xlsx_path, symbol_to_ensembl, sheet_name=SHEET_NAME):
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

    # The sheet's "gene_id" column is actually an HGNC symbol (see module
    # docstring) -- rename it and map onto real Ensembl gene IDs so this
    # output joins correctly against the rest of the pipeline.
    combined = combined.rename(columns={"gene_id": "gene_symbol"})
    combined["gene_id"] = combined["gene_symbol"].map(symbol_to_ensembl)

    n_unmapped = combined["gene_id"].isna().sum()
    if n_unmapped:
        log.warning(
            f"{n_unmapped} of {len(combined)} Watson et al. gene symbols did not "
            "map to an Ensembl gene_id via the GTF (not found, or ambiguous) -- "
            "dropped from output rather than kept with a null gene_id."
        )
    combined = combined.dropna(subset=["gene_id"])

    combined = combined[["gene_id", "gene_symbol", "log2FC", "padj"]]
    combined["source"] = "NAR_supplementary_File3_polysome"
    return combined


def fetch_watson_polysome_data(xlsx_path, gtf_path, out_path, sheet_name=SHEET_NAME):
    if not xlsx_path:
        raise RuntimeError(
            "config stage2_references.watson_nar_supp_path is empty. Stage "
            "'Supplementary Excel File 3.xlsx' from the NAR article's supplementary "
            "ZIP (gkz1060_supplemental_files.zip, from "
            "https://academic.oup.com/nar/article/48/1/116/5614571) locally and "
            "point stage2_references.watson_nar_supp_path at it."
        )
    symbol_to_ensembl = build_symbol_to_ensembl_map(gtf_path)
    result = parse_watson_supplementary(xlsx_path, symbol_to_ensembl, sheet_name)

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
        gtf_path=snakemake.input.gtf,
        out_path=snakemake.output.polysome_fc,
        sheet_name=snakemake.params.sheet_name,
    )
