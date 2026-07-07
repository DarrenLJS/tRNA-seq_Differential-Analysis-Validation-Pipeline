"""
workflow/scripts/parse_trax_tRF_classes.py

Parses TRAX's per-cell-line normalized read count table into separate
5'-tRF / 3'-tRF / i-tRF / tiRNA count matrices for rule deseq2_trf.

DESIGNED AGAINST DOCUMENTATION, NOT A REAL FILE -- flagged explicitly in
the rule docstring. TRAX's read-count output typically includes a
fragment-type or feature-name column that encodes the tRF class (e.g. via
a suffix or prefix in the feature name -- conventions have varied across
TRAX versions between something like "trf5", "trf3", "trf-1", "itrf",
"tirna" substrings in a `Feature`/`type` column, or embedded in the
feature ID itself). Since the exact TRAX version's output schema was not
available to inspect while building this pipeline, this parser tries
several detection strategies in order and logs which one succeeded, so a
human can verify the classification is actually correct against real
output before trusting deseq2_trf's results.

If ALL detection strategies fail to find a usable class signal, this
script raises rather than silently producing an empty/wrong classification
-- a loud failure here is much cheaper to fix than a downstream DESeq2 run
on garbage class assignments.
"""

import logging
import os
import re

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# Class name -> regex patterns tried against column names / feature-ID
# strings, in order. FIX once real TRAX output is inspected -- these are
# best-guess patterns based on TRAX's documented naming conventions, not
# confirmed against a live file.
CLASS_PATTERNS = {
    "5prime_tRF": [r"\btrf.?5\b", r"\b5.?trf\b", r"5[-_]?half", r"5p.?tRF"],
    "3prime_tRF": [r"\btrf.?3\b", r"\b3.?trf\b", r"3[-_]?half", r"3p.?tRF"],
    "i_tRF":      [r"\bi.?trf\b", r"internal.?trf"],
    "tiRNA":      [r"\btirna\b", r"tRNA.?halves?", r"stress.?induced"],
}


def _detect_class_column(df):
    """Look for a column whose name itself suggests fragment type/class."""
    candidates = [c for c in df.columns if re.search(r"type|class|feature|fragment", c, re.IGNORECASE)]
    return candidates


def _classify_by_feature_id(feature_id):
    for cls, patterns in CLASS_PATTERNS.items():
        for pat in patterns:
            if re.search(pat, str(feature_id), re.IGNORECASE):
                return cls
    return None


def parse_trf_classes(readcounts_path, trf_classes, min_unique_cov, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    df = pd.read_csv(readcounts_path, sep="\t")
    log.info(f"Loaded TRAX readcounts: {df.shape[0]} rows x {df.shape[1]} cols. Columns: {list(df.columns)[:10]}...")

    id_col_candidates = [c for c in df.columns if re.search(r"^(feature|id|name|gene)", c, re.IGNORECASE)]
    if not id_col_candidates:
        id_col = df.columns[0]
        log.warning(f"No obvious ID column found by name; falling back to first column '{id_col}'")
    else:
        id_col = id_col_candidates[0]
    log.info(f"Using '{id_col}' as feature ID column")

    # ---- Strategy 1: dedicated class/type column ----
    class_cols = _detect_class_column(df)
    class_series = None
    strategy_used = None

    if class_cols:
        col = class_cols[0]
        log.info(f"Strategy 1: found candidate class column '{col}', attempting pattern match on its values")
        mapped = df[col].apply(_classify_by_feature_id)
        if mapped.notna().sum() > 0:
            class_series = mapped
            strategy_used = f"class_column:{col}"

    # ---- Strategy 2: classify from the feature ID string itself ----
    if class_series is None or class_series.notna().sum() == 0:
        log.info("Strategy 2: attempting classification directly from feature ID strings")
        mapped = df[id_col].apply(_classify_by_feature_id)
        if mapped.notna().sum() > 0:
            class_series = mapped
            strategy_used = f"feature_id_pattern:{id_col}"

    if class_series is None or class_series.notna().sum() == 0:
        raise RuntimeError(
            "Could not classify any rows into 5'-tRF/3'-tRF/i-tRF/tiRNA using either "
            "a dedicated class column or feature-ID pattern matching. This TRAX "
            "output's naming convention does not match any pattern in "
            "CLASS_PATTERNS -- inspect the real file's columns/feature IDs "
            f"(sample: {df[id_col].head(10).tolist()}) and update CLASS_PATTERNS "
            "in parse_trax_tRF_classes.py accordingly. Refusing to proceed with an "
            "unclassified or garbage classification."
        )

    log.info(f"Classification succeeded via: {strategy_used}")
    df["_trf_class"] = class_series
    n_classified = df["_trf_class"].notna().sum()
    n_unclassified = df["_trf_class"].isna().sum()
    log.info(f"Classified {n_classified}/{len(df)} rows; {n_unclassified} rows unclassified (dropped)")
    log.info(f"Class breakdown:\n{df['_trf_class'].value_counts()}")

    sample_cols = [c for c in df.columns if c not in (id_col, "_trf_class") and c not in class_cols]

    for cls in trf_classes:
        sub = df[df["_trf_class"] == cls]
        if sub.empty:
            log.warning(f"No rows classified as '{cls}' -- writing empty matrix (deseq2_trf.R must handle this gracefully)")
        mat = sub[[id_col] + sample_cols].set_index(id_col)
        # min_unique_cov filter -- drop features with < threshold summed
        # count across all samples (proxy for "uniquely mapping reads"
        # threshold from Stage 1's trax.min_unique_cov; TRAX's own
        # normalizedreadcounts table may already reflect unique-only
        # counts depending on version -- FIX-verify this isn't double-filtering).
        mat = mat[mat.sum(axis=1) >= min_unique_cov]
        out_path = os.path.join(out_dir, f"{cls}_counts_matrix.tsv")
        mat.to_csv(out_path, sep="\t")
        log.info(f"Wrote {cls}: {mat.shape[0]} features x {mat.shape[1]} samples -> {out_path}")

    return strategy_used


if __name__ == "__main__":
    parse_trf_classes(
        readcounts_path=snakemake.input.trax_readcounts,
        trf_classes=snakemake.params.trf_classes,
        min_unique_cov=snakemake.params.min_unique_cov,
        out_dir=snakemake.output.class_matrices_dir,
    )
