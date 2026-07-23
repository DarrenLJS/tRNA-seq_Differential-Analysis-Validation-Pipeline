"""
workflow/scripts/build_isodecoder_mismatch_table.py

Builds the collated isodecoder mismatch table (count-based modification
score branch, supervisor spec): one row per (isodecoder, sample), columns

    isodecoder, total_count_isd, mismatch_A, mismatch_T, mismatch_G,
    mismatch_C, actual_count_isd, sample

in that order, tab-separated.

WHY THIS READS THE RAW mismatch/ DIRECTORY, NOT Stage 1's
pos34_mismatch_matrix.tsv
------------------------------------------------------------------------
Stage 1's own collation script (collate_mismatch.py, rule
build_mismatch_matrices) deliberately SUMS the per-type ('type' column:
observed base A/C/G/T) mismatch proportions into one "total mismatch"
value per (isodecoder, sample), for the existing binomial GLM (rule 11)
which only needs "any mismatch vs none". That collapse throws away
exactly the per-letter breakdown this table needs, so this script goes
back to mim-tRNAseq's raw per-(isodecoder, pos, type, sample) mismatch
table (moved by Stage 1 rule `mimtrnaseq` to
{SCRATCH}/pass1_mimtrnaseq/{cell_line}/mismatch/, same directory Stage
1's own build_mismatch_matrices rule consumes) instead of the collated
matrix. No Stage 1 code is touched or rerun -- this is a second, parallel
consumer of an already-existing raw output.

mismatch_A/T/G/C SEMANTICS
------------------------------------------------------------------------
Each column is the PROPORTION (rate in [0, 1], not a raw count) of
position-34 reads for that isodecoder/sample observed as that base,
following the same "recover proportion from whatever mim-tRNAseq
actually stored" logic as Stage 1's collate_mismatch.py (mismatch data
may be stored as a raw count or as a pre-computed proportion depending on
mim-tRNAseq version; either way this script normalises to a proportion:
count / coverage). The base matching the isodecoder's own genomic
reference at position 34 will have a mismatch_X column that is ~0 by
construction (mim-tRNAseq does not usually emit a "mismatch to the
reference base" row at all; if it does, it is treated the same as any
other type here and will simply come out near zero).

actual_count_isd FORMULA (supervisor spec, applied literally)
------------------------------------------------------------------------
    actual_count_isd = total_count_isd
                        - total_count_isd * mismatch_A
                        - total_count_isd * mismatch_T
                        - total_count_isd * mismatch_G
                        - total_count_isd * mismatch_C

total_count_isd is the WHOLE-ISODECODER read count from Stage 1's
isodecoder_counts_matrix.tsv (not the position-34-specific coverage in
the mismatch file). This formula therefore extrapolates a locus-specific,
position-34 misread rate across the isodecoder's entire read population --
a deliberate approximation (per supervisor spec, not invented here), not
a claim that every read was actually observed at position 34. Flagging
this once, here, rather than re-deriving it silently downstream.
"""

import glob
import logging
import os
import re

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

MISMATCH_TYPES = ["A", "T", "G", "C"]


def _match_sample(label, cl_samples):
    """
    Map a mim-tRNAseq sample label (plain id or full BAM path) to a
    manifest sample_id. Identical logic to Stage 1's collate_mismatch.py
    match_sample() -- kept in sync deliberately so a sample that resolves
    there resolves the same way here.
    """
    label = str(label)
    if os.sep in label or "/" in label:
        label = os.path.basename(label)
        for suffix in (
            ".fq.gz.unpaired_uniq.bam",
            ".fastq.gz.unpaired_uniq.bam",
            ".fq.unpaired_uniq.bam",
            ".fastq.unpaired_uniq.bam",
            ".bam",
            ".fq.gz",
            ".fastq.gz",
            ".fq",
            ".fastq",
        ):
            if label.endswith(suffix):
                label = label[: -len(suffix)]
                break
        label = re.sub(r"_val_\d+$", "", label)
    for s in cl_samples:
        if label == s or label.startswith(s):
            return s
    return None


def _load_raw_mismatch(mismatch_dir):
    """
    Locate and parse mim-tRNAseq's raw mismatch table, same candidate-file
    search and flexible column detection as Stage 1's collate_mismatch.py,
    EXCEPT the 'type' column is deliberately kept (not summed away).
    """
    candidates = (
        glob.glob(os.path.join(mismatch_dir, "*mismatch*.csv")) +
        glob.glob(os.path.join(mismatch_dir, "*mismatch*.txt")) +
        glob.glob(os.path.join(mismatch_dir, "*.csv")) +
        glob.glob(os.path.join(mismatch_dir, "*.txt"))
    )
    seen = set()
    candidates = [p for p in candidates if not (p in seen or seen.add(p))]
    if not candidates:
        raise FileNotFoundError(f"No files found in mismatch dir: {mismatch_dir}")

    for fpath in candidates:
        try:
            df = pd.read_csv(fpath, sep=None, engine="python", dtype={"canon_pos": str})
            df.columns = [c.lower().strip() for c in df.columns]

            if "canon_pos" in df.columns and "pos" not in df.columns:
                df = df.rename(columns={"canon_pos": "pos"})
            elif "canon_pos" in df.columns:
                df["pos"] = df["canon_pos"].astype(str)

            # FIX: force 'pos' to string regardless of which source column it
            # came from. If the file's own 'pos' column (rather than
            # 'canon_pos') was used, pandas infers int64 from a
            # numeric-looking column before we ever get a chance to compare
            # against the string "34" below -- silently returning zero
            # position-34 rows instead of failing loudly. canon_pos is read
            # with dtype=str explicitly (see pd.read_csv call above) because
            # it can contain non-numeric variable-loop labels ("20a", "e8");
            # plain 'pos' has no such guarantee, so cast defensively here.
            if "pos" in df.columns:
                df["pos"] = df["pos"].astype(str).str.strip()

            if not {"pos", "cov", "type"}.issubset(df.columns):
                log.info(
                    f"Skipping {os.path.basename(fpath)}: missing pos/cov/type "
                    f"(found: {list(df.columns)})"
                )
                continue

            s_col = next(
                (c for c in df.columns
                 if c in ("sample", "id", "library", "name", "bam", "path", "file")),
                None,
            )
            if s_col is None:
                log.info(f"Skipping {os.path.basename(fpath)}: no sample column")
                continue

            i_col = next(
                (c for c in df.columns
                 if c in ("isodecoder", "cluster", "trna", "gene", "feature")
                 and c != s_col),
                None,
            )
            if i_col is None:
                reserved = {s_col, "pos", "cov", "type"}
                i_col = next(
                    (c for c in df.columns if c not in reserved and df[c].dtype == object),
                    None,
                )
            if i_col is None:
                log.info(f"Skipping {os.path.basename(fpath)}: no isodecoder column")
                continue

            m_col = next(
                (c for c in df.columns
                 if c in ("mm", "mismatch", "mismatch_count", "count",
                          "proportion", "rate", "mismatch_rate",
                          "mismatch_proportion", "fraction")),
                None,
            )
            if m_col is None:
                log.info(f"Skipping {os.path.basename(fpath)}: no mismatch column")
                continue

            log.info(
                f"Parsed {os.path.basename(fpath)}: isodecoder_col='{i_col}', "
                f"sample_col='{s_col}', mismatch_col='{m_col}', type_col='type'"
            )
            return df, i_col, s_col, m_col

        except Exception as e:
            log.warning(f"Could not parse {os.path.basename(fpath)}: {e}")
            continue

    raise RuntimeError(
        f"Could not parse any mismatch file in {mismatch_dir} with pos/cov/type "
        f"columns present. Check log for column names found."
    )


def build_isodecoder_mismatch_table(
    mismatch_dir, isodecoder_counts_path, manifest_path, cell_line, out_path, log_path=None,
):
    manifest = pd.read_csv(manifest_path, sep="\t", index_col="sample_id")
    cl_samples = manifest[manifest["cell_line"] == cell_line].index.tolist()
    log.info(f"Cell line {cell_line}: {len(cl_samples)} samples: {cl_samples}")

    df, iso_col, sample_col, mm_col = _load_raw_mismatch(mismatch_dir)

    # Normalise 'type' to uppercase single-letter base, restrict to A/T/G/C
    # (drop any 'N'/ambiguous type rows mim-tRNAseq may emit).
    df["type"] = df["type"].astype(str).str.upper().str.strip()
    df = df[df["type"].isin(MISMATCH_TYPES)].copy()

    pos34 = df[df["pos"] == "34"].copy()
    if pos34.empty:
        raise RuntimeError(
            f"No position-34 rows found. Positions present (first 20): "
            f"{list(df['pos'].unique()[:20])}"
        )

    pos34["sample_id"] = pos34[sample_col].apply(lambda x: _match_sample(x, cl_samples))
    unmatched = pos34["sample_id"].isna()
    if unmatched.any():
        log.warning(
            f"{unmatched.sum()} rows could not be matched to a sample_id -- "
            f"labels: {pos34.loc[unmatched, sample_col].unique().tolist()}"
        )
    pos34 = pos34[~unmatched].copy()

    # Recover proportion if the source column is a raw count rather than a
    # pre-computed rate -- same detection Stage 1 uses.
    is_proportion = mm_col in (
        "proportion", "rate", "mismatch_rate", "mismatch_proportion", "fraction",
    )
    if is_proportion:
        pos34["prop"] = pos34[mm_col].astype(float)
    else:
        pos34["prop"] = (pos34[mm_col].astype(float) / pos34["cov"].astype(float)).fillna(0.0)
    pos34["prop"] = pos34["prop"].clip(lower=0.0, upper=1.0).fillna(0.0)

    # Collapse any duplicate (isodecoder, sample_id, type) rows (e.g. two
    # raw mim-tRNAseq sample labels resolving to the same manifest
    # sample_id) by averaging the proportion -- same convention as Stage
    # 1's collate_mismatch.py uses when collapsing duplicates.
    collapsed = (
        pos34.groupby([iso_col, "sample_id", "type"])["prop"]
        .mean()
        .reset_index()
    )

    wide = collapsed.pivot_table(
        index=[iso_col, "sample_id"], columns="type", values="prop", fill_value=0.0,
    ).reset_index()
    for t in MISMATCH_TYPES:
        if t not in wide.columns:
            wide[t] = 0.0
    wide = wide.rename(columns={iso_col: "isodecoder", "sample_id": "sample"})
    wide = wide.rename(columns={t: f"mismatch_{t}" for t in MISMATCH_TYPES})

    # ---- join total_count_isd from Stage 1's isodecoder count matrix ----
    iso_counts = pd.read_csv(isodecoder_counts_path, sep="\t", index_col=0)
    iso_counts = iso_counts.reindex(columns=cl_samples)

    counts_long = (
        iso_counts.reset_index()
        .melt(id_vars=iso_counts.index.name or "index", var_name="sample", value_name="total_count_isd")
        .rename(columns={(iso_counts.index.name or "index"): "isodecoder"})
    )

    out = wide.merge(counts_long, on=["isodecoder", "sample"], how="inner")

    n_mismatch_only = len(wide) - len(out)
    if n_mismatch_only > 0:
        log.warning(
            f"{n_mismatch_only} (isodecoder, sample) mismatch rows had no matching "
            f"isodecoder in isodecoder_counts_matrix.tsv -- dropped (likely an "
            f"isodecoder that was mismatch-called but zero-count-filtered out of "
            f"Stage 1's count matrix)."
        )

    out["actual_count_isd"] = out["total_count_isd"] - out["total_count_isd"] * (
        out["mismatch_A"] + out["mismatch_T"] + out["mismatch_G"] + out["mismatch_C"]
    )
    out["actual_count_isd"] = out["actual_count_isd"].clip(lower=0.0).round().astype(int)

    out = out[
        ["isodecoder", "total_count_isd", "mismatch_A", "mismatch_T", "mismatch_G",
         "mismatch_C", "actual_count_isd", "sample"]
    ].sort_values(["isodecoder", "sample"]).reset_index(drop=True)

    out.to_csv(out_path, sep="\t", index=False)
    log.info(f"Wrote {len(out)} (isodecoder, sample) rows -> {out_path}")

    if log_path:
        with open(log_path, "a") as fh:
            fh.write(f"Wrote {len(out)} (isodecoder, sample) rows -> {out_path}\n")
            fh.write(f"Isodecoders: {out['isodecoder'].nunique()}, samples: {out['sample'].nunique()}\n")

    return out


if __name__ == "__main__":
    log_path = snakemake.log[0] if len(snakemake.log) else None
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()

    build_isodecoder_mismatch_table(
        mismatch_dir=snakemake.input.mismatch_dir,
        isodecoder_counts_path=snakemake.input.isodecoder_counts,
        manifest_path=snakemake.params.manifest,
        cell_line=snakemake.params.cell_line,
        out_path=snakemake.output.mismatch_table,
        log_path=log_path,
    )
