"""
workflow/scripts/parse_trax_tRF_classes.py

Parses TRAX's per-cell-line normalized read count table into separate
5'-tRF / 3'-tRF / i-tRF count matrices for rule deseq2_trf.

CONFIRMED against the real file (A549-normalizedreadcounts.txt /
THP1-normalizedreadcounts.txt). Two things the original guessed-pattern
version got wrong, both now fixed:

1. FILE FORMAT: TRAX writes this file R-style -- the header row has one
   fewer field than the data rows (feature ID column is unlabeled row
   names). pandas' default read_csv auto-detects this and puts the
   feature ID into the DataFrame's INDEX, not into any column. The
   previous version searched df.columns for an ID-like column name,
   found none, and silently fell back to df.columns[0] -- which is
   actually the first SAMPLE's read counts, not feature IDs. That's why
   classification failed on every row: it was pattern-matching numbers,
   not strings.

2. CLASS VOCABULARY: the feature ID itself carries the class as a fixed
   suffix, e.g. "tRNA-Ala-AGC-1_fiveprime". The real suffix vocabulary is
   exactly five values -- wholecounts, fiveprime, threeprime, other,
   antisense -- not the "trf5"/"5-half"/"itrf"/"tirna"-style substrings
   originally guessed. Mapping to this pipeline's class names
   (config trf_diff_abundance.trf_classes):
     fiveprime   -> 5prime_tRF
     threeprime  -> 3prime_tRF
     other       -> i_tRF        (TRAX's "other" bucket is internal-fragment
                                   reads that don't fit the 5'/3'-end
                                   windows -- the closest real-data match
                                   to i-tRF)
     wholecounts -> excluded (full-length tRNA reads, not a fragment)
     antisense   -> excluded (antisense-mapping reads, QC/noise category)

   TRAX's default output has NO dedicated tiRNA (stress-induced tRNA
   half) category -- tiRNAs are ~30-40nt halves and tRFs are ~14-30nt
   fragments, a length distinction TRAX does not encode in this file at
   all (it would require joining against a read-length distribution,
   e.g. *-readlengths.txt, which is a separate design decision outside
   this fix's scope). If "tiRNA" is requested in
   trf_diff_abundance.trf_classes, this script writes an empty matrix for
   it (consistent with the pre-existing "no rows classified" fallback)
   and logs a clear, loud warning rather than silently guessing a length
   cutoff -- decide separately whether a length-based tiRNA split is
   actually needed before trusting any tiRNA-labelled output downstream.
"""

import logging

import os

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# TRAX's real, fixed suffix vocabulary (confirmed against
# A549-normalizedreadcounts.txt / THP1-normalizedreadcounts.txt), mapped
# to this pipeline's class names. Not all TRAX suffixes are tRF classes:
# wholecounts and antisense are deliberately excluded (see module
# docstring). tiRNA has no TRAX-native suffix at all -- it is NOT a key
# in this map, and is handled as an explicit empty-matrix case below if
# requested in config.
TRAX_SUFFIX_TO_CLASS = {
    "fiveprime":  "5prime_tRF",
    "threeprime": "3prime_tRF",
    "other":      "i_tRF",
}
EXCLUDED_TRAX_SUFFIXES = {"wholecounts", "antisense"}


def _classify_by_suffix(feature_id):
    """
    Deterministic classification: TRAX feature IDs are
    "{locus}_{suffix}" where suffix is exactly one of wholecounts /
    fiveprime / threeprime / other / antisense. Split on the last
    underscore and look up the suffix directly -- no regex guessing.
    """
    suffix = str(feature_id).rsplit("_", 1)[-1]
    if suffix in TRAX_SUFFIX_TO_CLASS:
        return TRAX_SUFFIX_TO_CLASS[suffix]
    return None  # covers EXCLUDED_TRAX_SUFFIXES and any unrecognised suffix


def parse_trf_classes(readcounts_path, trf_classes, min_unique_cov, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # TRAX writes this file R-style: the header row has one fewer field
    # than the data rows, so pandas auto-detects the feature ID as the
    # DataFrame's index rather than a column. Read explicitly with
    # index_col=0 to make that intentional rather than relying on
    # pandas' implicit fallback (which is what silently broke the
    # previous version -- it never noticed the ID wasn't in df.columns
    # at all).
    df = pd.read_csv(readcounts_path, sep="\t", index_col=0)
    log.info(f"Loaded TRAX readcounts: {df.shape[0]} rows x {df.shape[1]} samples. "
             f"Sample columns: {list(df.columns)[:5]}...")
    log.info(f"Feature ID index sample: {df.index[:5].tolist()}")

    class_series = df.index.to_series().apply(_classify_by_suffix)
    n_classified = class_series.notna().sum()
    n_excluded = df.shape[0] - n_classified

    if n_classified == 0:
        raise RuntimeError(
            "Could not classify any rows by TRAX suffix (_fiveprime/_threeprime/"
            "_other). This file's feature-ID naming convention does not match "
            f"the expected TRAX format -- sample IDs seen: {df.index[:10].tolist()}. "
            "Inspect the real file and update TRAX_SUFFIX_TO_CLASS in "
            "parse_trax_tRF_classes.py accordingly. Refusing to proceed with an "
            "unclassified or garbage classification."
        )

    log.info(f"Classified {n_classified}/{df.shape[0]} rows by suffix; "
             f"{n_excluded} rows excluded (wholecounts/antisense/unrecognised)")
    log.info(f"Class breakdown:\n{class_series.value_counts()}")

    df["_trf_class"] = class_series
    sample_cols = [c for c in df.columns if c != "_trf_class"]


    for cls in trf_classes:
        if cls == "tiRNA":
            log.warning(
                "'tiRNA' was requested in trf_diff_abundance.trf_classes, but TRAX's "
                "normalizedreadcounts.txt has NO tiRNA-specific category -- its "
                "suffixes are wholecounts/fiveprime/threeprime/other/antisense only. "
                "tiRNAs (~30-40nt halves) are not distinguished from tRFs (~14-30nt "
                "fragments) here; that split would need a read-length-based filter "
                "(e.g. against *-readlengths.txt), which this script does NOT "
                "attempt. Writing an empty tiRNA matrix -- decide separately whether "
                "a length-based tiRNA definition is actually needed before trusting "
                "any tiRNA-labelled output downstream."
            )
            mat = pd.DataFrame(columns=sample_cols)
            out_path = os.path.join(out_dir, f"{cls}_counts_matrix.tsv")
            mat.to_csv(out_path, sep="\t")
            log.info(f"Wrote {cls}: 0 features x {len(sample_cols)} samples (empty, by design) -> {out_path}")
            continue

        sub = df[df["_trf_class"] == cls]
        if sub.empty:
            log.warning(f"No rows classified as '{cls}' -- writing empty matrix (deseq2_trf.R must handle this gracefully)")
        mat = sub[sample_cols]
        # min_unique_cov filter -- drop features with < threshold summed
        # count across all samples (proxy for "uniquely mapping reads"
        # threshold from Stage 1's trax.min_unique_cov; TRAX's own
        # normalizedreadcounts table may already reflect unique-only
        # counts depending on version -- FIX-verify this isn't double-filtering).
        mat = mat[mat.sum(axis=1) >= min_unique_cov]
        out_path = os.path.join(out_dir, f"{cls}_counts_matrix.tsv")
        mat.to_csv(out_path, sep="\t")
        log.info(f"Wrote {cls}: {mat.shape[0]} features x {mat.shape[1]} samples -> {out_path}")


if __name__ == "__main__":
    parse_trf_classes(
        readcounts_path=snakemake.input.trax_readcounts,
        trf_classes=snakemake.params.trf_classes,
        min_unique_cov=snakemake.params.min_unique_cov,
        out_dir=snakemake.output.class_matrices_dir,
    )
