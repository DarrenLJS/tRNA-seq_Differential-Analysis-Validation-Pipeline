"""
workflow/scripts/build_decoding_whitelist.py

Builds the explicit isodecoder -> reachable-codon decoding whitelist that
Delta(c) (rule 14) is computed from. This is the single source of truth for
"which isodecoder, in which modification state, can decode which codon" --
every downstream Delta(c) number depends on this table being right, so the
logic here follows the case-by-case derivation worked through against the
GtRNAdb anticodon-count tables (Four Box / Six Box / Two Box / Two Box and
Other), not a generic anticodon-shape pattern match.

WHY NOT PATTERN-MATCH ON ANTICODON SHAPE
-----------------------------------------
"Any G34 tRNA in a two-codon Y-ending box" is NOT sufficient to identify Q34
substrates: Cys-GCA is G34, sits in a superficially Asp/Asn/His/Tyr-shaped
two-box context, but position 35 = C, not U -- it fails the G-U-N consensus
and is NOT a queuosine substrate in humans. Symmetrically, "any A34 tRNA" is
not sufficient to identify I34 substrates: Cys-ACA is A34 but Cys is not one
of the 8 canonical I34 amino-acid families. Ile is a further trap: Ile-AAT
(A34, I34-eligible) and Ile-TAT (U34, NOT I34-eligible, plain Watson-Crick
reader of AUA) must not be merged despite sharing an amino acid.

This script therefore requires BOTH an isotype whitelist (config
wobble_glm.i34_isotypes / q34_isotypes) AND a structural check (N34 base,
and for Q34 additionally N35 == T) before assigning an isodecoder to the
I34 or Q34 bucket. Any isotype in the config list whose anticodon fails the
structural check is logged as a warning and demoted to the default bucket
rather than silently trusted -- this is meant to catch anticodon_map data
errors, not just enforce the rule.

DECODING / WOBBLE RULES IMPLEMENTED
------------------------------------
Anticodon is read 5'->3' as N34 N35 N36 (confirmed against the GtRNAdb
table's own worked examples, e.g. Ala-AGC / Ile-AAT / Asp-GTC).
Codon (5'->3') = complement(N36) + complement(N35) + wobble(N34), i.e.
codon position 1 pairs N36, position 2 pairs N35 (always strict
Watson-Crick, position 2 is never a wobble position), position 3 (the
codon's own wobble position) pairs N34 per the rules below. DNA-letter
alphabet (T, not U) is used throughout to match the anticodon_map's own
convention and the GtRNAdb table's anticodon spelling.

  N34 = C            -> reads codon3 = G only (strict WC, no wobble at all)
  N34 = A, unmodified -> reads codon3 = T only (WC A:T)
  N34 = A, I34-edited -> reads codon3 in {T, C, A} (inosine wobble)
                          ONLY for isotypes in i34_isotypes
  N34 = G, unmodified -> reads codon3 in {C, T} (WC G:C + native G:T wobble)
                          -- this applies to EVERY G34 isodecoder, whether or
                          not it is Q34-eligible; unmodified wobble pairing
                          is a baseline structural fact, not something that
                          requires the Q34 whitelist.
  N34 = G, Q-modified -> SAME codon set {C, T} as unmodified -- queuosine
                          does not add a new reachable codon, it modulates
                          decoding EFFICIENCY on the C-ending codon only.
                          ONLY for isotypes in q34_isotypes AND N35 == T.
  N34 = T, unmodified -> reads codon3 in {A, G} (WC U:A + classical Crick
                          wobble U:G). No modification pathway is modeled
                          for U34 in this pipeline -- six-box / NAR-type
                          U34 chemistry (mcm5s2U etc.) is explicitly OUT OF
                          SCOPE per the proposal's stated boundaries, so
                          every U34 isodecoder falls into the default
                          bucket regardless of isotype.

TERM TYPES AND GAMMA (see rule 14 / compute_delta_c.py for how these are
consumed -- gamma is the exponent applied to the f_stim/f_ctrl ratio):

  canonical      gamma=0   plain log2[FC(i)]                  (default bucket)
  both_I         gamma=0   plain log2[FC(i)]                  (I34, T-ending codon:
                                                                 editing is strictly
                                                                 additive, so the
                                                                 U-ending codon is
                                                                 decoded by the WHOLE
                                                                 pool regardless of
                                                                 editing status)
  mod_only_I     gamma=1   log2[FC(i)*(f_stim/f_ctrl)]         (I34, C/A-ending codon:
                                                                 reachable ONLY via
                                                                 editing)
  both_Q_C       gamma=kappa log2[FC(i)*(f_stim/f_ctrl)^kappa] (Q34, C-ending codon:
                                                                 efficiency-modulated)
  both_Q_U       gamma=0   plain log2[FC(i)]                   (Q34, U-ending codon:
                                                                 no Q-dependence modeled)

Exactly one term type fires per (isodecoder, codon) pair that is structurally
reachable at all; unreachable pairs are simply absent from the whitelist
(this is the Sum #2 completeness property discussed in the proposal, not a
zero-valued row -- absence, not an explicit zero, keeps the whitelist small:
~400+ isodecoders x 61 codons would be a mostly-empty matrix otherwise).

STOP CODONS / SUPPRESSOR & SELENOCYSTEINE tRNAs
-------------------------------------------------
The 3 stop codons (TAA, TAG, TGA in DNA-letter convention) are excluded from
the 61-sense-codon output entirely. Any isodecoder whose only structural
"codon" would be a stop codon (SelCys, Supres suppressor tRNAs) is dropped
from the whitelist with a note in the QC log -- Delta(c) is defined over the
61 sense codons per the proposal and these are out of scope by definition,
not an oversight.

INPUT
-----
anticodon_map (Stage 1 reference, read-only): expected to provide, per
tRNA locus/isodecoder, at minimum an isotype/amino-acid column and an
anticodon column (3-letter, DNA alphabet, 5'->3', e.g. "AGC"). Column
names are auto-detected across a few common variants since the exact
header wasn't confirmed against a real file when this script was written
-- see `_detect_columns()`. FIX at first real run if the file's actual
header doesn't match any of the candidates tried.

OUTPUT
------
decoding_whitelist.tsv, one row per (isodecoder_id, codon) pair, columns:
  isodecoder_id, isotype, anticodon, position34_base, bucket,
  codon, term_type, gamma_expr, notes
"""

import sys
import logging
from itertools import product

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
STOP_CODONS_DNA = {"TAA", "TAG", "TGA"}

# Standard genetic code (DNA-letter codons -> 3-letter amino acid), used as
# a cross-check filter on every wobble-expanded codon candidate. This is
# NOT optional decoration -- without it, blanket wobble rules (e.g.
# "unmodified U34 reads {A,G}") will happily generate biologically
# impossible reach, such as an Ile-TAT tRNA "reading" ATG (Met). AUG is a
# split-box exception reserved for the dedicated Met-tRNA and is never
# read by any Ile isoacceptor -- caught by testing this script against the
# Ile three-way split before this table existed. Any wobble-rule-derived
# codon whose translated amino acid does not match the tRNA's own isotype
# is dropped, not just for U34 but uniformly across all buckets, since this
# is a general correctness property, not an Ile-specific patch.
STANDARD_GENETIC_CODE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "*",   "TAG": "*",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "*",   "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

# Isotype label normalisation -- anticodon_map may spell some isotypes
# differently from the standard 3-letter genetic-code table (e.g.
# "iMet/Met" for the initiator, "Sec"/"SelCys" for selenocysteine, which
# reads a UGA stop codon and is dropped by the stop-codon filter anyway).
# Anything not in this map is compared to the genetic code table as-is;
# if that also fails to match, the isotype is treated as "unmappable" and
# the genetic-code cross-check is SKIPPED with a logged note rather than
# incorrectly zeroing out a real isodecoder's reach.
ISOTYPE_NORMALISE = {
    "iMet/Met": "Met",
    "iMet": "Met",
    "fMet": "Met",
}
UNMAPPABLE_ISOTYPES = {"SelCys", "Sec", "Supres"}

# Candidate column-name variants to try when reading the Stage-1
# anticodon_map.tsv -- FIX to a single confirmed name once the real file
# header is checked (see build_anticodon_map.py in Stage 1 rule 00d).
ISODECODER_ID_CANDIDATES = ["isodecoder_id", "locus", "locus_id", "tRNA_id", "gene_id", "name"]
ISOTYPE_CANDIDATES       = ["isotype", "amino_acid", "aa", "AA"]
ANTICODON_CANDIDATES     = ["anticodon", "anticodon_seq", "AC"]


def _detect_columns(df):
    """Auto-detect isodecoder_id / isotype / anticodon columns by name."""
    def pick(candidates, label):
        for c in candidates:
            if c in df.columns:
                return c
        raise ValueError(
            f"Could not find a '{label}' column in anticodon_map among "
            f"candidates {candidates}. Actual columns present: {list(df.columns)}. "
            f"Update ISODECODER_ID_CANDIDATES/ISOTYPE_CANDIDATES/ANTICODON_CANDIDATES "
            f"in build_decoding_whitelist.py to match the real file."
        )
    id_col   = pick(ISODECODER_ID_CANDIDATES, "isodecoder_id")
    iso_col  = pick(ISOTYPE_CANDIDATES, "isotype")
    ac_col   = pick(ANTICODON_CANDIDATES, "anticodon")
    return id_col, iso_col, ac_col


def _codon_from_anticodon(anticodon, n34_override=None):
    """
    Given a 3-letter anticodon string (5'->3', N34 N35 N36), return the
    codon position-1/position-2 (fixed, strict WC) as a 2-letter prefix,
    and separately the set of possible codon-position-3 letters reachable
    given n34 (or n34_override, used to force the "modified" wobble base
    behaviour without mutating the anticodon string itself).
    """
    n34, n35, n36 = anticodon[0], anticodon[1], anticodon[2]
    codon_pos1 = COMPLEMENT[n36]
    codon_pos2 = COMPLEMENT[n35]
    prefix = codon_pos1 + codon_pos2

    base34 = n34_override if n34_override is not None else n34

    if base34 == "C":
        pos3_unmod = {"G"}
    elif base34 == "A":
        pos3_unmod = {"T"}
    elif base34 == "G":
        pos3_unmod = {"C", "T"}
    elif base34 == "T":
        pos3_unmod = {"A", "G"}
    else:
        raise ValueError(f"Unrecognised N34 base '{base34}' in anticodon '{anticodon}'")

    return prefix, pos3_unmod


def _codon_matches_isotype(codon, isotype):
    """
    True if `codon` translates (standard genetic code) to the same amino
    acid as `isotype`, OR if `isotype` can't be mapped to a standard
    3-letter code at all (SelCys/Supres edge cases) -- in which case the
    check is skipped rather than incorrectly zeroing out real reach.
    Returns (matches: bool, was_skipped: bool).
    """
    if isotype in UNMAPPABLE_ISOTYPES:
        return True, True
    norm = ISOTYPE_NORMALISE.get(isotype, isotype)
    expected_aa = STANDARD_GENETIC_CODE.get(codon)
    if expected_aa is None:
        return True, True  # stop codon or malformed -- handled elsewhere
    if norm not in set(STANDARD_GENETIC_CODE.values()):
        return True, True  # isotype label doesn't map cleanly -- skip, don't guess
    return (expected_aa == norm), False


def build_whitelist(anticodon_map_path, i34_isotypes, q34_isotypes, out_path, log_path=None):
    df = pd.read_csv(anticodon_map_path, sep="\t")
    id_col, iso_col, ac_col = _detect_columns(df)
    log.info(f"Using columns: isodecoder_id='{id_col}', isotype='{iso_col}', anticodon='{ac_col}'")

    i34_set = set(i34_isotypes)
    q34_set = set(q34_isotypes)

    rows = []
    dropped_stop = []
    demoted_warnings = []
    genetic_code_rejections = []
    genetic_code_skipped = []

    for _, r in df.iterrows():
        isodecoder_id = r[id_col]
        isotype       = str(r[iso_col]).strip()
        anticodon     = str(r[ac_col]).strip().upper()

        if len(anticodon) != 3 or any(b not in "ACGT" for b in anticodon):
            log.warning(f"Skipping {isodecoder_id}: malformed anticodon '{anticodon}'")
            continue

        n34 = anticodon[0]
        n35 = anticodon[1]

        is_i34_candidate = (n34 == "A") and (isotype in i34_set)
        is_q34_candidate = (n34 == "G") and (isotype in q34_set)

        # Defensive structural check -- demote rather than trust the config
        # list blindly, since a data error in anticodon_map should not
        # silently propagate into Delta(c).
        if is_q34_candidate and n35 != "T":
            demoted_warnings.append(
                f"{isodecoder_id} (isotype={isotype}, anticodon={anticodon}) is in "
                f"q34_isotypes but fails the G-U-N consensus (N35='{n35}' != 'T') -- "
                f"demoted to default bucket. Check anticodon_map for a data error, "
                f"or remove this isotype from config q34_isotypes if it's a genuine "
                f"non-Q34 paralog."
            )
            is_q34_candidate = False

        if is_i34_candidate and isotype not in i34_set:  # unreachable given above, kept for clarity
            is_i34_candidate = False

        bucket = "I34" if is_i34_candidate else ("Q34" if is_q34_candidate else "default")

        # ---- reachable codons + term assignment ----
        if bucket == "I34":
            # Unmodified A34 reaches only T-ending (prefix computed once;
            # I34's wobble range {T,C,A} is hardcoded below per the
            # inosine wobble rule rather than re-derived, since inosine is
            # not one of the four standard bases _codon_from_anticodon()
            # switches on).
            prefix, _pos3_unmod_unused = _codon_from_anticodon(anticodon)
            # T-ending: decoded by whole pool regardless of editing status -> both_I
            candidates = [("T", "both_I", "0",
                           "U-ending codon: decoded by entire pool (edited + unedited); editing is strictly additive for I34."),
                          ("C", "mod_only_I", "1",
                           "Reachable only by the I34-edited fraction of the pool (inosine wobble)."),
                          ("A", "mod_only_I", "1",
                           "Reachable only by the I34-edited fraction of the pool (inosine wobble).")]
            for pos3, term_type, gamma_expr, note in candidates:
                codon = prefix + pos3
                matches, skipped = _codon_matches_isotype(codon, isotype)
                if skipped:
                    genetic_code_skipped.append(f"{isodecoder_id} / {codon} (isotype '{isotype}' not mappable)")
                if not matches:
                    genetic_code_rejections.append(
                        f"{isodecoder_id} (isotype={isotype}, anticodon={anticodon}): wobble rule proposed "
                        f"codon {codon}, but genetic code says {codon}={STANDARD_GENETIC_CODE.get(codon)}, "
                        f"not {isotype} -- REJECTED (I34 bucket)."
                    )
                    continue
                rows.append(dict(
                    isodecoder_id=isodecoder_id, isotype=isotype, anticodon=anticodon,
                    position34_base=n34, bucket=bucket, codon=codon,
                    term_type=term_type, gamma_expr=gamma_expr, notes=note,
                ))

        elif bucket == "Q34":
            prefix, pos3_unmod = _codon_from_anticodon(anticodon)          # {'C','T'}
            candidates = [("C", "both_Q_C", "kappa",
                           "C-ending codon: reachable unmodified AND Q-modified; Q efficiency effect applied via kappa exponent on f_stim/f_ctrl."),
                          ("T", "both_Q_U", "0",
                           "U-ending codon: reachable regardless of Q status; no Q-dependence modeled (literature effect is C-ending-specific -- modeling simplification, see rule 09 docstring).")]
            for pos3, term_type, gamma_expr, note in candidates:
                codon = prefix + pos3
                matches, skipped = _codon_matches_isotype(codon, isotype)
                if skipped:
                    genetic_code_skipped.append(f"{isodecoder_id} / {codon} (isotype '{isotype}' not mappable)")
                if not matches:
                    genetic_code_rejections.append(
                        f"{isodecoder_id} (isotype={isotype}, anticodon={anticodon}): wobble rule proposed "
                        f"codon {codon}, but genetic code says {codon}={STANDARD_GENETIC_CODE.get(codon)}, "
                        f"not {isotype} -- REJECTED (Q34 bucket)."
                    )
                    continue
                rows.append(dict(
                    isodecoder_id=isodecoder_id, isotype=isotype, anticodon=anticodon,
                    position34_base=n34, bucket=bucket, codon=codon,
                    term_type=term_type, gamma_expr=gamma_expr, notes=note,
                ))

        else:
            if isotype in UNMAPPABLE_ISOTYPES:
                # These isotypes (SelCys/Sec, Supres) are non-standard-code
                # readers -- their real target is a recoded/suppressed stop
                # codon, already excluded by the stop-codon filter below.
                # The genetic-code cross-check is intentionally SKIPPED for
                # them (can't validate against a standard amino acid), which
                # means an ordinary wobble-rule candidate could otherwise
                # slip through unchecked (e.g. SelCys-TCA's T34 wobble rule
                # proposes both TGA [stop, correctly dropped] and TGG [Trp,
                # WRONG -- SelCys does not decode UGG in vivo]). Rather than
                # let an unvalidated candidate through, these isotypes are
                # excluded from the whitelist entirely -- caught by
                # inspecting real build output, not anticipated in advance.
                demoted_warnings.append(
                    f"{isodecoder_id} (isotype={isotype}) excluded entirely from whitelist "
                    f"-- non-standard-code isotype, genetic-code cross-check cannot validate "
                    f"its wobble-rule-derived reach, so no sense-codon assignment is trusted."
                )
                continue
            prefix, pos3_set = _codon_from_anticodon(anticodon)
            for pos3 in sorted(pos3_set):
                codon = prefix + pos3
                matches, skipped = _codon_matches_isotype(codon, isotype)
                if skipped:
                    genetic_code_skipped.append(f"{isodecoder_id} / {codon} (isotype '{isotype}' not mappable)")
                if not matches:
                    genetic_code_rejections.append(
                        f"{isodecoder_id} (isotype={isotype}, anticodon={anticodon}): wobble rule proposed "
                        f"codon {codon}, but genetic code says {codon}={STANDARD_GENETIC_CODE.get(codon)}, "
                        f"not {isotype} -- REJECTED (default bucket). This is exactly the Ile-TAT/AUG-type "
                        f"box-splitting trap the genetic-code cross-check exists to catch."
                    )
                    continue
                rows.append(dict(
                    isodecoder_id=isodecoder_id, isotype=isotype, anticodon=anticodon,
                    position34_base=n34, bucket=bucket, codon=codon,
                    term_type="canonical", gamma_expr="0",
                    notes="No modification pathway modeled for this isodecoder (out of I34/Q34 scope, or non-eligible paralog of a modified isotype).",
                ))

    wl = pd.DataFrame(rows)

    # Drop any row whose codon is a stop codon (SelCys / suppressor tRNAs,
    # or any structural artifact) -- Delta(c) is defined over the 61 sense
    # codons only.
    is_stop = wl["codon"].isin(STOP_CODONS_DNA)
    if is_stop.any():
        dropped_stop = wl.loc[is_stop, "isodecoder_id"].unique().tolist()
        wl = wl.loc[~is_stop].copy()

    wl = wl.sort_values(["codon", "bucket", "isodecoder_id"]).reset_index(drop=True)
    wl.to_csv(out_path, sep="\t", index=False)

    # ---- completeness check: every one of the 61 sense codons should have
    # at least one contributing isodecoder (Sum #2 completeness property) ----
    all_dna_codons = {"".join(p) for p in product("ACGT", repeat=3)}
    sense_codons = sorted(all_dna_codons - STOP_CODONS_DNA)
    covered = set(wl["codon"].unique())
    missing = sorted(set(sense_codons) - covered)

    n_i34 = wl.loc[wl.bucket == "I34", "isodecoder_id"].nunique()
    n_q34 = wl.loc[wl.bucket == "Q34", "isodecoder_id"].nunique()
    n_default = wl.loc[wl.bucket == "default", "isodecoder_id"].nunique()

    summary_lines = [
        f"Whitelist rows written: {len(wl)}",
        f"Isodecoders assigned: I34={n_i34}, Q34={n_q34}, default={n_default}",
        f"Sense codons covered: {len(covered)}/61",
    ]
    if missing:
        summary_lines.append(f"WARNING - codons with ZERO contributing isodecoders: {missing}")
    if dropped_stop:
        summary_lines.append(f"Dropped (stop-codon-only) isodecoders: {dropped_stop}")
    if demoted_warnings:
        summary_lines.append("Q34 structural-check demotions:")
        summary_lines.extend(f"  - {w}" for w in demoted_warnings)
    if genetic_code_rejections:
        summary_lines.append(f"Genetic-code cross-check REJECTIONS ({len(genetic_code_rejections)}) -- wobble-rule candidates dropped for crossing into a different amino acid's codon:")
        summary_lines.extend(f"  - {w}" for w in genetic_code_rejections)
    if genetic_code_skipped:
        summary_lines.append(f"Genetic-code cross-check SKIPPED ({len(genetic_code_skipped)}) -- isotype label not mappable to standard 3-letter code, filter not applied for these:")
        summary_lines.extend(f"  - {w}" for w in genetic_code_skipped)

    for line in summary_lines:
        log.info(line)

    if log_path:
        with open(log_path, "w") as fh:
            fh.write("\n".join(summary_lines) + "\n")

    if missing:
        log.warning(
            "Codon coverage is incomplete -- Delta(c) will be undefined (silently "
            "zero by omission) for the codons listed above. This most likely means "
            "the anticodon_map is missing loci for a codon family, or a bucket "
            "assignment bug excluded a valid isodecoder. Check before trusting "
            "downstream Delta(c) output."
        )

    return wl


if __name__ == "__main__":
    # Snakemake `script:` directive entry point.
    build_whitelist(
        anticodon_map_path=snakemake.input.anticodon_map,
        i34_isotypes=snakemake.params.i34_isotypes,
        q34_isotypes=snakemake.params.q34_isotypes,
        out_path=snakemake.output.whitelist,
        log_path=snakemake.log[0] if len(snakemake.log) else None,
    )
