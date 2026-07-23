"""
workflow/scripts/build_anticodon_count_matrix.py

Builds the preprocessed anticodon count table (count-based modification
score branch, supervisor spec): one row per (read_anticodon, sample),
columns

    read_anticodon, read_codon, total_count_isa, sample

read_codon is a comma-separated list -- one read_anticodon can decode more
than one codon (wobble decoding), so this is genuinely a list column, not
a single value.

Also writes a wide count MATRIX (rows=read_anticodon, cols=sample) for
downstream DESeq2/edgeR (rule deseq2_anticodon etc.), analogous to Stage
1's isodecoder_counts_matrix.tsv.

INPUT
-----
isodecoder_mismatch_table.tsv  (this branch's own table #1, see
                                 build_isodecoder_mismatch_table.py)
decoding_whitelist.tsv         (rule 09 -- REUSED, not rebuilt: supplies
                                 isotype, anticodon, position34_base,
                                 bucket per isodecoder, and the
                                 isodecoder -> reachable-codon term_type
                                 map this script re-groups by anticodon)

MODIFIED-VARIANT ROW CONSTRUCTION
----------------------------------
For each anticodon, up to three rows are emitted:

  canonical row      read_anticodon = anticodon (unchanged)
  I34-variant row     read_anticodon = "I" + anticodon[1:]   (only if the
                       anticodon has >=1 I34-eligible isodecoder)
  Q34-variant row     read_anticodon = "Q" + anticodon[1:]   (only if the
                       anticodon has >=1 Q34-eligible isodecoder)

Per isodecoder, per sample:

  canonical_count   = actual_count_isd                       (table #1,
                        already nets out ALL FOUR mismatch types --
                        this is a stricter subtraction than "just the
                        modeled signal", see build_isodecoder_mismatch_table.py
                        docstring; kept as-is per supervisor spec)

  i34_mod_count     = round( total_count_isd * mismatch_<I34_SIGNAL_BASE> )
                        only for isodecoders with position34_base ==
                        I34_REF_BASE and bucket == "I34"; else 0.
                        NOT kappa-weighted -- I34 is the confirmatory
                        signal (mirrors rule 11's treatment: full weight,
                        no confidence dial), same as the rate-based branch.

  q34_mod_count     = round( kappa * total_count_isd *
                              sum(mismatch_<t> for t in Q34_SIGNAL_BASES) )
                        only for isodecoders with position34_base ==
                        Q34_REF_BASE and bucket == "Q34"; else 0.
                        kappa-weighted -- mirrors the existing rate-based
                        branch's Q34 confidence dial (rule 11 / 14),
                        swept over the same wobble_glm.kappa_sweep values.
                        kappa=0 -> Q34-variant counts are zero everywhere
                        (no Q34-variant row emitted at all), matching the
                        existing pipeline's default of ignoring Q34 signal
                        until empirically justified (rule 16 sweep).

I34_SIGNAL_BASE / Q34_SIGNAL_BASES ASSUMPTIONS (flagged, not silently
buried): I34 uses a single-letter signature (A34 read as G -- the
well-established inosine RT signature). Q34 has no comparable clean
single-letter literature signature (consistent with rule 11's framing
of Q34 detection as near-background), so the default here sums ALL
non-reference letters at a G34 position (A+T+C) as the Q34 proxy --
the same "any mismatch" logic rule 11's GLM already uses, just scoped to
Q34-eligible isodecoders specifically. Both are configurable via
config["count_mod_score"] rather than hardcoded, so this assumption can
be swapped without touching this script.

NOTE ON CONSERVATION: canonical_count is computed independently of
kappa/i34_mod_count/q34_mod_count (it is simply table #1's
actual_count_isd), so canonical + I34-variant + Q34-variant counts for a
given isodecoder/sample do NOT necessarily sum to total_count_isd -- some
reads (mismatches not attributed to either modeled signal, e.g. the
T/C "noise" at an A34 position) are excluded from all three rows by
design, and Q34-variant counts additionally scale down with kappa without
canonical scaling back up to compensate. This is intentional (each row
is only ever populated by the signal it claims to represent), not a
bookkeeping bug -- documented here so it isn't rediscovered as one.

read_codon RECONSTRUCTION FROM THE WHITELIST
----------------------------------------------
The whitelist is keyed per (isodecoder_id, codon) with a term_type label.
term_type depends on anticodon + bucket, not on isodecoder identity per
se, so codon-reach is re-derived here by GROUPING whitelist rows by
anticodon (not rebuilding wobble logic):

  canonical row's codons  <- term_type in {canonical, both_I, both_Q_C}
  I34-variant row's codons <- term_type == mod_only_I
  Q34-variant row's codons <- term_type == mod_only_Q

If isodecoders sharing the same anticodon disagree on their term_type/
codon set (should not happen -- term_type is a function of anticodon +
bucket -- but not structurally impossible if two isotypes ever shared an
anticodon spelling), this is logged as a warning and the union of codons
is used rather than silently picking one.
"""

import logging
import os

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

FLAT_TERM_TYPES = {"canonical", "both_I", "both_Q_C"}
I34_TERM_TYPE = "mod_only_I"
Q34_TERM_TYPE = "mod_only_Q"


def _codon_reach_by_anticodon(whitelist):
    """
    Returns three dict[anticodon] -> sorted list of codons: canonical,
    i34_variant, q34_variant. Logs (does not raise) if isodecoders sharing
    an anticodon disagree on term_type/codon membership for that group.
    """
    canonical, i34, q34 = {}, {}, {}

    for anticodon, grp in whitelist.groupby("anticodon"):
        flat_codons = sorted(grp.loc[grp["term_type"].isin(FLAT_TERM_TYPES), "codon"].unique())
        i34_codons = sorted(grp.loc[grp["term_type"] == I34_TERM_TYPE, "codon"].unique())
        q34_codons = sorted(grp.loc[grp["term_type"] == Q34_TERM_TYPE, "codon"].unique())

        # Consistency check: every isodecoder row contributing to this
        # anticodon's flat group should reach the same flat codon set as
        # every other isodecoder row for this anticodon (per isodecoder,
        # not per group) -- spot-check by isodecoder rather than assuming.
        per_iso_flat = (
            grp[grp["term_type"].isin(FLAT_TERM_TYPES)]
            .groupby("isodecoder_id")["codon"]
            .apply(lambda s: tuple(sorted(s.unique())))
        )
        if per_iso_flat.nunique() > 1:
            log.warning(
                f"Anticodon {anticodon}: isodecoders disagree on canonical codon "
                f"reach ({dict(per_iso_flat)}). Using the union: {flat_codons}."
            )

        canonical[anticodon] = flat_codons
        i34[anticodon] = i34_codons
        q34[anticodon] = q34_codons

    return canonical, i34, q34


def build_anticodon_count_matrix(
    mismatch_table_path, whitelist_path, kappa,
    i34_ref_base, i34_signal_base,
    q34_ref_base, q34_signal_bases,
    out_long_path, out_matrix_path, log_path=None,
):
    mm = pd.read_csv(mismatch_table_path, sep="\t")
    wl = pd.read_csv(whitelist_path, sep="\t")

    iso_info = wl[["isodecoder_id", "isotype", "anticodon", "position34_base", "bucket"]].drop_duplicates()

    merged = mm.merge(iso_info, left_on="isodecoder", right_on="isodecoder_id", how="left")
    n_unmatched = merged["anticodon"].isna().sum()
    if n_unmatched:
        missing_isodecoders = merged.loc[merged["anticodon"].isna(), "isodecoder"].unique().tolist()
        log.warning(
            f"{n_unmatched} (isodecoder, sample) rows have no entry in "
            f"decoding_whitelist.tsv -- excluded from the anticodon matrix "
            f"(cannot determine anticodon/codon reach without it). "
            f"Isodecoders affected: {missing_isodecoders[:10]}"
            + (" ... (truncated)" if len(missing_isodecoders) > 10 else "")
        )
    merged = merged.dropna(subset=["anticodon"]).copy()

    # ---- per-isodecoder/sample counts ----
    merged["canonical_count"] = merged["actual_count_isd"]

    is_i34 = (merged["bucket"] == "I34") & (merged["position34_base"] == i34_ref_base)
    merged["i34_mod_count"] = 0.0
    merged.loc[is_i34, "i34_mod_count"] = (
        merged.loc[is_i34, "total_count_isd"] * merged.loc[is_i34, f"mismatch_{i34_signal_base}"]
    )
    merged["i34_mod_count"] = merged["i34_mod_count"].round().clip(lower=0).astype(int)

    is_q34 = (merged["bucket"] == "Q34") & (merged["position34_base"] == q34_ref_base)
    q34_rate = merged[[f"mismatch_{t}" for t in q34_signal_bases]].sum(axis=1)
    merged["q34_mod_count"] = 0.0
    merged.loc[is_q34, "q34_mod_count"] = (
        kappa * merged.loc[is_q34, "total_count_isd"] * q34_rate.loc[is_q34]
    )
    merged["q34_mod_count"] = merged["q34_mod_count"].round().clip(lower=0).astype(int)

    # ---- aggregate to anticodon level ----
    canonical_by_ac = (
        merged.groupby(["anticodon", "sample"])["canonical_count"].sum().reset_index()
        .rename(columns={"anticodon": "read_anticodon", "canonical_count": "total_count_isa"})
    )
    i34_by_ac = (
        merged.groupby(["anticodon", "sample"])["i34_mod_count"].sum().reset_index()
    )
    i34_by_ac = i34_by_ac[i34_by_ac["i34_mod_count"] > 0].rename(
        columns={"i34_mod_count": "total_count_isa"}
    )
    i34_by_ac["read_anticodon"] = "I" + i34_by_ac["anticodon"].str[1:]

    q34_by_ac = (
        merged.groupby(["anticodon", "sample"])["q34_mod_count"].sum().reset_index()
    )
    q34_by_ac = q34_by_ac[q34_by_ac["q34_mod_count"] > 0].rename(
        columns={"q34_mod_count": "total_count_isa"}
    )
    q34_by_ac["read_anticodon"] = "Q" + q34_by_ac["anticodon"].str[1:]

    canonical_codons, i34_codons, q34_codons = _codon_reach_by_anticodon(wl)

    canonical_by_ac["read_codon"] = canonical_by_ac["read_anticodon"].map(
        lambda a: ",".join(canonical_codons.get(a, []))
    )
    i34_by_ac["read_codon"] = i34_by_ac["anticodon"].map(lambda a: ",".join(i34_codons.get(a, [])))
    q34_by_ac["read_codon"] = q34_by_ac["anticodon"].map(lambda a: ",".join(q34_codons.get(a, [])))

    long_out = pd.concat(
        [
            canonical_by_ac[["read_anticodon", "read_codon", "total_count_isa", "sample"]],
            i34_by_ac[["read_anticodon", "read_codon", "total_count_isa", "sample"]],
            q34_by_ac[["read_anticodon", "read_codon", "total_count_isa", "sample"]],
        ],
        ignore_index=True,
    )
    # Drop rows with no reachable codon at all (whitelist gave nothing to
    # decode) -- these can't contribute to any per-codon score and would
    # otherwise show up as an unexplained empty read_codon field.
    n_no_codon = (long_out["read_codon"] == "").sum()
    if n_no_codon:
        log.warning(
            f"{n_no_codon} anticodon/sample rows have zero reachable codons "
            f"(whitelist coverage gap) -- dropped."
        )
    long_out = long_out[long_out["read_codon"] != ""].copy()

    long_out = long_out.sort_values(["read_anticodon", "sample"]).reset_index(drop=True)
    long_out.to_csv(out_long_path, sep="\t", index=False)
    log.info(f"Wrote {len(long_out)} (read_anticodon, sample) rows -> {out_long_path}")

    # ---- wide count matrix for DESeq2/edgeR ----
    matrix = long_out.pivot_table(
        index="read_anticodon", columns="sample", values="total_count_isa",
        aggfunc="sum", fill_value=0,
    )
    matrix = matrix.round().astype(int)
    matrix.to_csv(out_matrix_path, sep="\t")
    log.info(
        f"Wrote anticodon count matrix: {matrix.shape[0]} read_anticodons x "
        f"{matrix.shape[1]} samples -> {out_matrix_path} (kappa={kappa})"
    )

    if log_path:
        with open(log_path, "a") as fh:
            fh.write(f"kappa={kappa}\n")
            fh.write(f"Long table: {len(long_out)} rows -> {out_long_path}\n")
            fh.write(f"Matrix: {matrix.shape[0]} x {matrix.shape[1]} -> {out_matrix_path}\n")
            fh.write(f"I34-variant anticodons: {i34_by_ac['read_anticodon'].nunique()}\n")
            fh.write(f"Q34-variant anticodons: {q34_by_ac['read_anticodon'].nunique()}\n")

    return long_out, matrix


if __name__ == "__main__":
    log_path = snakemake.log[0] if len(snakemake.log) else None
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()

    build_anticodon_count_matrix(
        mismatch_table_path=snakemake.input.mismatch_table,
        whitelist_path=snakemake.input.whitelist,
        kappa=float(snakemake.wildcards.kappa),
        i34_ref_base=snakemake.params.i34_ref_base,
        i34_signal_base=snakemake.params.i34_signal_base,
        q34_ref_base=snakemake.params.q34_ref_base,
        q34_signal_bases=snakemake.params.q34_signal_bases,
        out_long_path=snakemake.output.anticodon_table,
        out_matrix_path=snakemake.output.anticodon_matrix,
        log_path=log_path,
    )
