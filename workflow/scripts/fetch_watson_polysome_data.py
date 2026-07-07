"""
workflow/scripts/fetch_watson_polysome_data.py

Fetches and parses Watson et al. (2020) polysome-profiling data into a
per-gene, per-timepoint translational log2FC table for rule 16 validation.

Source: GEO GSE130618 (superseries; polysome-seq + total RNA-seq
subseries for the antiviral/ILF3 translational program). A supplementary
data ZIP (~5.7MB) is also hosted at NAR Online alongside the article and
most likely contains the authors' own DESeq2-derived tables.

STRATEGY (in order)
--------------------
1. Try GEO's processed series matrix / supplementary files directly via
   the GEO FTP/HTTP layout for GSE130618. GEO superseries structure is
   NOT guaranteed to expose a single clean "translational efficiency"
   column -- it may instead provide separate polysome-fraction and
   total-RNA count matrices that require computing the poly:total ratio
   fold change ourselves (polysome_log2FC - total_log2FC, i.e. the
   standard translational-efficiency definition). Both cases are handled
   below; which one actually applies needs confirming against the real
   GEO download (FIX-flag).
2. If GEO's data isn't directly usable (missing gene-level FC, or the
   subseries structure doesn't parse as expected), fall back to
   downloading and parsing the NAR supplementary ZIP -- URL currently a
   placeholder in config (stage2_references.watson_nar_supp_url); locate
   and fill in before relying on this path.

OUTPUT
------
watson_polysome_foldchange.tsv, columns: gene_id, timepoint, log2FC,
padj (if available from the source), source ("GEO" | "NAR_supplementary").

NEEDS A REAL SCOPING PASS -- this script is designed against GEO/NAR
documentation, not a downloaded file. Expect a correction cycle here,
the same as Stage 1's TRAX rule needed against real TRAX output.
"""

import io
import gzip
import logging
import zipfile

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _geo_series_matrix_url(accession):
    # GEO's path convention: GSE130nnn/GSE130618/matrix/GSE130618_series_matrix.txt.gz
    stub = accession[:-3] + "nnn"
    return f"{GEO_FTP_BASE}/{stub}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def try_fetch_geo(accession):
    url = _geo_series_matrix_url(accession)
    log.info(f"Attempting GEO series matrix fetch: {url}")
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        raw = gzip.decompress(resp.content).decode("utf-8", errors="replace")

        # GEO series matrix files are metadata-heavy with a
        # "!series_matrix_table_begin" / "_end" delimited expression block.
        # This is unlikely to directly contain per-gene translational log2FC
        # (that's normally a DERIVED quantity, not raw expression) -- more
        # likely this path only gets us sample-level normalized counts,
        # requiring us to compute polysome-vs-total log2FC ourselves per
        # gene per timepoint from the sample metadata (!Sample_title /
        # !Sample_characteristics rows to identify polysome vs total,
        # timepoint, replicate). That mapping is NOT implemented here yet --
        # FIX once the real file's !Sample_* rows are inspected, since GEO
        # sample naming conventions vary by submitting lab and can't be
        # reliably guessed in advance.
        if "!series_matrix_table_begin" not in raw:
            raise ValueError("Series matrix file did not contain the expected data table markers")

        log.warning(
            "GEO series matrix fetched successfully, but this script does not yet "
            "implement the polysome-vs-total log2FC derivation from raw sample-level "
            "counts (requires inspecting real !Sample_* metadata rows to identify "
            "which samples are polysome vs total RNA, and which timepoint each "
            "belongs to). Falling back to the NAR supplementary path, which is more "
            "likely to contain the authors' own precomputed fold-change table."
        )
        return None
    except Exception as e:
        log.warning(f"GEO series matrix fetch/parse failed ({e})")
        return None


def try_fetch_nar_supplementary(nar_supp_url):
    if not nar_supp_url:
        log.warning(
            "config stage2_references.watson_nar_supp_url is empty -- cannot attempt "
            "the NAR supplementary fallback. Locate the exact supplementary ZIP URL "
            "from the article page (Nucleic Acids Research, Watson et al. 2020) and "
            "fill it into config_stage2.yaml."
        )
        return None
    log.info(f"Attempting NAR supplementary ZIP fetch: {nar_supp_url}")
    try:
        resp = requests.get(nar_supp_url, timeout=120)
        resp.raise_for_status()
        zf = zipfile.ZipFile(io.BytesIO(resp.content))
        log.info(f"Supplementary ZIP contents: {zf.namelist()}")

        # FIX: the exact filename(s) inside the ZIP containing per-gene
        # translational fold-change data are unknown until the real ZIP is
        # inspected. Try common patterns (xlsx/csv/tsv with "fold" or "DE"
        # or "translat" in the name) and log what's found rather than
        # guessing a single hardcoded filename.
        candidates = [
            n for n in zf.namelist()
            if any(kw in n.lower() for kw in ("fold", "de_", "deseq", "translat", "polysome"))
            and n.lower().endswith((".csv", ".tsv", ".xlsx", ".txt"))
        ]
        if not candidates:
            log.warning("No obviously-named fold-change table found in the ZIP -- manual inspection required.")
            return None

        log.info(f"Candidate table(s) in ZIP: {candidates}")
        # Use the first candidate; log clearly which one so a human can verify.
        target = candidates[0]
        with zf.open(target) as fh:
            if target.lower().endswith(".xlsx"):
                df = pd.read_excel(fh)
            else:
                sep = "," if target.lower().endswith(".csv") else "\t"
                df = pd.read_csv(fh, sep=sep)

        log.info(f"Parsed '{target}': columns = {list(df.columns)}")
        # FIX: column name mapping to (gene_id, timepoint, log2FC, padj) is
        # unknown until the real file is seen -- this needs a human pass to
        # confirm which columns correspond to which quantity before this
        # function returns a usable long-format table. Currently returns the
        # raw parsed frame with a warning rather than guessing column names.
        log.warning(
            "Returning RAW parsed supplementary table without column mapping to "
            "(gene_id, timepoint, log2FC, padj) -- inspect columns above and "
            "update the mapping logic in this function before trusting rule 16."
        )
        return df
    except Exception as e:
        log.warning(f"NAR supplementary fetch/parse failed ({e})")
        return None


def fetch_watson_polysome_data(geo_accession, nar_supp_url, out_path):
    result = try_fetch_geo(geo_accession)
    source = "GEO"
    if result is None:
        result = try_fetch_nar_supplementary(nar_supp_url)
        source = "NAR_supplementary"

    if result is None:
        raise RuntimeError(
            "Both the GEO and NAR supplementary fetch paths failed to produce a "
            "usable table. This rule needs a manual scoping pass against the real "
            "data source before it can run unattended -- see docstring and inline "
            "FIX comments in fetch_watson_polysome_data.py."
        )

    result["source"] = source
    result.to_csv(out_path, sep="\t", index=False)
    log.info(f"Wrote Watson et al. polysome data ({len(result)} rows, source={source}) -> {out_path}")
    return result


if __name__ == "__main__":
    fetch_watson_polysome_data(
        geo_accession=snakemake.params.geo_accession,
        nar_supp_url=snakemake.params.nar_supp_url,
        out_path=snakemake.output.polysome_fc,
    )
