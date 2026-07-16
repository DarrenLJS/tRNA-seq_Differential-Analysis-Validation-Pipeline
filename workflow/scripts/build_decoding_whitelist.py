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
  N34 = A, I34-edited -> reads codon3 in {T, C} (inosine wobble to U,C at
                          meaningful efficiency) ONLY for isotypes in
                          i34_isotypes. REVISED 2026-07-16: I34:A pairing
                          (the third leg of Crick's 1966 I34->{U,C,A} rule)
                          is dropped -- it is structurally real but occurs
                          at low efficiency in vivo (Curran 1995, PMC306738:
                          three independent in vivo assays in E. coli show
                          slow/unstable decoding at CGA, an I34:A-read
                          codon). This whitelist tracks functionally
                          meaningful decoding capacity, not every
                          structurally-possible pairing.
  N34 = G, unmodified -> reads codon3 = C ONLY (WC G:C). REVISED
                          2026-07-16: native G34:U wobble is dropped for
                          the same reason as I34:A above -- it is real but
                          low-efficiency in vivo (Fields et al. 2022,
                          Front. Mol. Biosci.: "Unmodified G34 tRNAs
                          efficiently decode C-ending codons through
                          traditional Watson-Crick interactions and wobble
                          to U-ending codons with reduced efficiency").
                          This applies to EVERY G34 isodecoder, whether or
                          not it is Q34-eligible -- it is a baseline
                          structural/efficiency fact, not something that
                          requires the Q34 whitelist.
  N34 = G, Q-modified -> ADDS codon3 = T as newly reachable at meaningful
                          efficiency (Meier et al. 1985, via Fields et al.
                          2022: Q34 tRNAs "read U- and C-ending codons
                          equally well"). ONLY for isotypes in q34_isotypes
                          AND N35 == T. This is a genuine reach-expansion
                          under this whitelist's efficiency standard (unlike
                          the pre-2026-07-16 model, which treated native
                          G34 as already reaching both C and U and Q34 as
                          purely efficiency-modulating).
  N34 = T, unmodified -> reads codon3 in {A, G} (WC U:A + classical Crick
                          wobble U:G). No modification pathway is modeled
                          for U34 in this pipeline -- six-box / NAR-type
                          U34 chemistry (mcm5s2U etc.) is explicitly OUT OF
                          SCOPE per the proposal's stated boundaries, so
                          every U34 isodecoder falls into the default
                          bucket regardless of isotype.

TERM TYPES AND GAMMA (see rule 14 / compute_delta_c.py for how these are
consumed -- gamma is the exponent applied to the f_stim/f_ctrl ratio).
REVISED 2026-07-16: Q34's two terms now mirror I34's structure exactly
(one flat/native term, one edit-gated term), per supervisor instruction;
only the edit-gated term's exponent differs (kappa, not full weight 1),
since Q34 detection is exploratory/low-sensitivity (rule 11 wobble_glm.R):

  canonical      gamma=0     plain log2[FC(i)]                  (default bucket)
  both_I         gamma=0     plain log2[FC(i)]                  (I34, T-ending codon:
                                                                   editing is strictly
                                                                   additive, so the
                                                                   U-ending codon is
                                                                   decoded by the WHOLE
                                                                   pool regardless of
                                                                   editing status)
  mod_only_I     gamma=1     log2[FC(i)*(f_stim/f_ctrl)]         (I34, C-ending codon:
                                                                   reachable at meaningful
                                                                   efficiency ONLY via
                                                                   editing; A-ending
                                                                   dropped, see above)
  both_Q_C       gamma=0     plain log2[FC(i)]                   (Q34, C-ending codon:
                                                                   native WC pairing for
                                                                   G34, decoded by the
                                                                   WHOLE pool regardless
                                                                   of editing status --
                                                                   mirrors both_I)
  mod_only_Q     gamma=kappa log2[FC(i)*(f_stim/f_ctrl)^kappa]   (Q34, U-ending codon:
                                                                   reachable at meaningful
                                                                   efficiency ONLY via
                                                                   Q-editing -- mirrors
                                                                   mod_only_I, but
                                                                   kappa-weighted rather
                                                                   than full weight 1)

Exactly one term type fires per (isodecoder, codon) pair that is structurally
reachable at all; unreachable pairs are simply absent from the whitelist
(this is the Sum #2 completeness property discussed in the proposal, not a
zero-valued row -- absence, not an explicit zero, keeps the whitelist small:
~400+ isodecoders x 61 codons would be a mostly-empty matrix otherwise).

MUTUAL EXCLUSIVITY CONSTRAINT (formal statement, finalised 2026-07-16)
------------------------------------------------------------------------
For a given isodecoder i, at most one of the five term-type indicators
below is 1 for a given codon c -- i.e. i contributes to Delta(c) through
exactly one mechanism, never two at once (no double-counting the same
decoding capacity):

  I(i_can, c) + I(i_mod^I34, c_C) + I(i_can, i_mod^I34, c_U)
             + I(i_can, i_mod^Q34, c_C) + I(i_mod^Q34, c_U)  <=  1

  I(i_can, c)                -- canonical: no fixed ending (fires on
                                 whichever codon ordinary WC pairing, or
                                 default-bucket wobble, reaches -- G-ending
                                 for C34, U-ending for non-I34 A34,
                                 C-ending for non-Q34/consensus-failing
                                 G34, A/G-ending for T34)
  I(i_mod^I34, c_C)           -- mod_only_I: edit-gated, fires on the
                                 C-ending codon only
  I(i_can, i_mod^I34, c_U)    -- both_I: flat/native, fires on the
                                 U-ending codon only
  I(i_can, i_mod^Q34, c_C)    -- both_Q_C: flat/native, fires on the
                                 C-ending codon only
  I(i_mod^Q34, c_U)           -- mod_only_Q: edit-gated, fires on the
                                 U-ending codon only

Reading the c_C / c_U subscripts: I34's edit-gated term is pinned to c_C
and its flat term to c_U; Q34's edit-gated term is pinned to c_U and its
flat term to c_C -- i.e. the two mirror each other with the C/U roles
swapped, which is the direct notational consequence of I34 gating on the
C-ending codon while Q34 (post-2026-07-16 revision) gates on the U-ending
codon. Pinning the subscript to a specific codon per term (rather than
using a shared unsubscripted c for the four modification-dependent terms)
means this expression can be summed across both codons of an isoacceptor's
box in one pass without ambiguity -- c_C and c_U are two distinct terms in
that sum, not two hidden instantiations of the same symbol. I(i_can, c)
is left unsubscripted deliberately: canonical isn't tied to one ending at
all (see above), so it should not carry a C/U subscript it doesn't have.

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
anticodon_map (Stage 1 reference, read-only): CONFIRMED against the real
file (Stage 1 rule 00d, build_anticodon_map.py) to have exactly two
columns -- `locus` and `anticodon` -- and no isotype/amino-acid column.
Column names for these two are still auto-detected across a few candidate
spellings for robustness (see `_detect_columns()`), but isotype is NOT
read from a column at all: it is derived from the `locus` string itself
via `_isotype_from_locus()`, mirroring the exact two parsing branches
Stage 1's build_anticodon_map.py used to build this same locus/anticodon
pairing in the first place:

  New format (GtRNAdb 2.0): "tRNA-{AminoAcid}-{Anticodon}-{Family}-{Copy}"
    e.g. "tRNA-Tyr-GTA-1-1" -> isotype = "Tyr"

  Old format (fallback): "...-{AminoAcid}{Anticodon}" concatenated as the
  last hyphen-delimited segment, e.g. "chr1.tRNA1-AlaAGC" -> isotype = "Ala"

Any locus whose isotype can't be recovered by either branch is treated as
unmappable (logged, demoted to the default bucket via the existing
UNMAPPABLE_ISOTYPES-style handling) rather than crashing the whole build.

OUTPUT
------
decoding_whitelist.tsv, one row per (isodecoder_id, codon) pair, columns:
  isodecoder_id, isotype, anticodon, position34_base, bucket,
  codon, term_type, gamma_expr, notes

FIXED -- isodecoder_id must be Stage-1's FINAL (post-clustering) ID,
not the raw per-locus anticodon_map ID
-----------------------------------------------------------------------
anticodon_map.tsv's `locus` values are RAW GtRNAdb loci (one row per
gene copy, e.g. "Homo_sapiens_tRNA-Lys-TTT-3-1", "...-3-2", ... "...-3-5").
mim-tRNAseq (Stage 1) merges some of these raw loci into a single
"isodecoder" per cell line, based on that cell line's own observed
coverage/mismatch signal (see e.g. A549_tRNAseq_unsplitClusterInfo.txt:
"insufficient coverage at mismatch X" / "potential mod at mismatch Y") --
producing IDs like "Homo_sapiens_tRNA-Lys-TTT-3/5" (family 3, with the
family-5 raw locus merged in). Every other Stage 2 file that carries an
"isodecoder_id" column (pos34_coverage_matrix.tsv, pos34_mismatch_matrix.tsv,
isodecoder_DESeq2_results.tsv, isodecoder_highconf_intersect.tsv, ...)
already uses these FINAL collapsed IDs, so a whitelist built directly
from anticodon_map's raw locus strings has ZERO overlap with them --
confirmed on a real run: wobble_glm.R found 0/94 I34 isodecoders present
in pos34_coverage_matrix.tsv.

Critically, this clustering is DATA-DRIVEN PER CELL LINE, not a fixed
mapping -- confirmed A549 and THP1 disagree on 41+ loci (e.g. A549 keeps
"Arg-TCT-2", "Arg-TCT-3", "Arg-TCT-5" as three separate isodecoders;
THP1 merges all three into one, "Arg-TCT-3/2/5"). So this whitelist can
no longer be built once, globally -- it is now build per cell_line (see
rule build_decoding_whitelist's {cell_line} wildcard), joining through
that cell line's own mim-tRNAseq annotation output:

  {cell_line}_tRNAseqclusterInfo.txt  -- one row per RAW locus, `parent`
                                          column = raw-locus representative
                                          of whatever cluster it merged into
  Isodecoder_counts.txt               -- one row per FINAL isodecoder,
                                          `parent` column = same raw-locus
                                          representative value, `isodecoder`
                                          column = the final collapsed ID

raw locus --[clusterInfo parent]--> parent raw locus
          --[Isodecoder_counts parent->isodecoder]--> FINAL isodecoder_id

This join (_load_locus_to_isodecoder_map) only relabels the OUTPUT
isodecoder_id -- isotype/anticodon parsing and all bucket/codon/term_type
logic above still operate on the RAW locus string exactly as before,
since that logic is genuinely a structural fact about the individual
gene's anticodon, not something clustering changes. Raw loci with no
entry in the cell line's clustering output (shouldn't normally happen for
in-scope nuclear/cytoplasmic tRNAs) are skipped with a logged warning
rather than crashing. Multiple raw loci that collapse to the same final
isodecoder necessarily produce identical whitelist rows (same anticodon,
same bucket, same codons) -- these are de-duplicated before writing out.
"""

import sys
import re
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

# Sentinel isotype label used when _isotype_from_locus() can't recover an
# amino acid from the locus string at all (neither parsing branch fires).
# Routed through the existing UNMAPPABLE_ISOTYPES handling downstream so
# these loci are excluded from the whitelist with a logged note, rather
# than crashing the build or silently mis-bucketing them.
ISOTYPE_UNPARSEABLE = "UNPARSEABLE"

# Candidate column-name variants to try when reading the Stage-1
# anticodon_map.tsv. CONFIRMED (Stage 1 rule 00d, build_anticodon_map.py)
# to be just `locus` and `anticodon` in the real file -- no isotype column
# exists there at all; isotype is derived separately, see
# _isotype_from_locus() and its call site in build_whitelist().
ISODECODER_ID_CANDIDATES = ["isodecoder_id", "locus", "locus_id", "tRNA_id", "gene_id", "name"]
ANTICODON_CANDIDATES     = ["anticodon", "anticodon_seq", "AC"]


def _detect_columns(df):
    """Auto-detect isodecoder_id / anticodon columns by name. Isotype is
    not a column in the real anticodon_map.tsv -- see _isotype_from_locus()."""
    def pick(candidates, label):
        for c in candidates:
            if c in df.columns:
                return c
        raise ValueError(
            f"Could not find a '{label}' column in anticodon_map among "
            f"candidates {candidates}. Actual columns present: {list(df.columns)}. "
            f"Update ISODECODER_ID_CANDIDATES/ANTICODON_CANDIDATES "
            f"in build_decoding_whitelist.py to match the real file."
        )
    id_col   = pick(ISODECODER_ID_CANDIDATES, "isodecoder_id")
    ac_col   = pick(ANTICODON_CANDIDATES, "anticodon")
    return id_col, ac_col


def _isotype_from_locus(locus, anticodon):
    """
    Recover the isotype (amino acid) from the locus string, since
    anticodon_map.tsv itself has no isotype column. Mirrors the exact two
    parsing branches Stage 1's build_anticodon_map.py used to build this
    same locus/anticodon pairing in the first place (rule 00d):

    New format (GtRNAdb 2.0): "tRNA-{AminoAcid}-{Anticodon}-{Family}-{Copy}"
      e.g. "tRNA-Tyr-GTA-1-1" -> isotype = parts[1] = "Tyr"

    Old format (fallback): last hyphen-delimited segment is the amino acid
      and anticodon concatenated, e.g. "chr1.tRNA1-AlaAGC" -> "Ala"

    Returns ISOTYPE_UNPARSEABLE if neither branch matches.
    """
    parts = locus.split("-")
    if (len(parts) >= 3 and parts[0].endswith("tRNA")
            and len(parts[2]) == 3 and parts[2].upper() == anticodon.upper()):
        return parts[1]

    tail = parts[-1]
    if len(tail) > 3 and tail.upper().endswith(anticodon.upper()):
        return tail[:-3]

    return ISOTYPE_UNPARSEABLE


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
        # REVISED (per supervisor correction, 2026-07-16 meeting + follow-up):
        # native (unmodified) G34 wobble to U-ending codons is real but occurs
        # at reduced efficiency in vivo (Fields et al. 2022, Front. Mol.
        # Biosci. -- "Unmodified G34 tRNAs efficiently decode C-ending codons
        # through traditional Watson-Crick interactions and wobble to
        # U-ending codons with reduced efficiency. Q34 tRNAs appear to read
        # U- and C-ending codons equally well" [Meier et al. 1985]). This
        # whitelist tracks FUNCTIONALLY MEANINGFUL decoding capacity, not
        # every structurally-possible-but-negligible pairing (Crick's 1966
        # rules permit G34->{C,T} unmodified; that structural permission is
        # not in dispute, but is not the standard used here) -- so native
        # G34 is treated as C-ending only, applied uniformly to EVERY G34
        # isodecoder (Q34-eligible or not; see build_whitelist()'s Q34
        # branch and the "default" bucket branch, both of which call this
        # function). U-ending only becomes reachable at meaningful
        # efficiency once Q34-edited -- see TERM TYPES table below.
        pos3_unmod = {"C"}
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


def _strip_copy_suffix(raw_locus):
    """
    Strip the trailing '-<copy_number>' segment from a raw locus string
    (Family-Copy format, e.g. "Homo_sapiens_tRNA-Lys-TTT-3-1" -> family
    "Homo_sapiens_tRNA-Lys-TTT-3").
    """
    return re.sub(r"-\d+$", "", str(raw_locus))


def _load_locus_to_isodecoder_map(unsplit_cluster_info_path, isodecoder_counts_path):
    """
    Build the raw-locus -> final Stage-1 isodecoder-ID map for one cell
    line. See module docstring "FIXED" note for why this join exists and
    why it must be cell-line-specific.

    CORRECTED (second attempt) -- {cell_line}_tRNAseqclusterInfo.txt
    (structural/alignment clustering) is NOT the right source for this
    join: confirmed on real data that clusterInfo.txt groups loci (e.g.
    Lys-TTT-1/2/4 into one "cluster") that Isodecoder_counts.txt does NOT
    actually merge (they remain three separate singleton isodecoders).
    clusterInfo.txt's `parent` column reflects a different, earlier
    clustering pass (likely structural/covariance-model alignment),
    unrelated to isodecoder identity.

    The REAL source of isodecoder merging is
    unsplit_cluster_info_path ({cell_line}_tRNAseq_unsplitClusterInfo.txt):
    one row per DATA-DRIVEN merge decision (columns: Parent, Size, "Unsplit
    transcripts", Reason -- reasons like "insufficient coverage at mismatch
    X" / "potential mod at mismatch Y", confirming this is the per-cell-line
    coverage/mismatch-driven collapsing described in the module docstring).
    `Parent` is a raw locus (Family-Copy); "Unsplit transcripts" is a
    comma-separated list of FAMILY-ONLY names (no copy suffix) that were
    merged into Parent's family.

    Verified against two real, independently-checked cases:
      Parent=Ile-AAT-5-4, Unsplit=Ile-AAT-7,Ile-AAT-3,Ile-AAT-8,Ile-AAT-1
        -> Isodecoder_counts.txt: isodecoder="Ile-AAT-5/1/3/7/8",
           Single_isodecoder=False, parent="Ile-AAT-5"
      Parent=Lys-TTT-3-1, Unsplit=Lys-TTT-5
        -> Isodecoder_counts.txt: isodecoder="Lys-TTT-3/5",
           Single_isodecoder=False, parent="Lys-TTT-3"
    In both cases Isodecoder_counts.txt's `parent` column is reliable ONLY
    for these Single_isodecoder=False (merged) rows -- for ordinary
    Single_isodecoder=True (singleton) rows, `parent` was observed to
    carry unrelated leftover values (e.g. many un-merged Ile-AAT-* rows
    all showing parent="Ile-AAT-5" despite never being merged with it),
    so it is NOT used for singleton lookups here.

    Algorithm:
      1. For each unsplitClusterInfo.txt row, build the family-group =
         {family(Parent)} union {each Unsplit-transcript family}.
      2. Find that group's real final ID: the Isodecoder_counts.txt row
         with Single_isodecoder==False and parent==family(Parent) --
         use its `isodecoder` value as ground truth (not a reconstructed
         string) to avoid depending on exact "/"-join ordering/format.
      3. Any family untouched by any unsplitClusterInfo.txt row defaults
         to its own family-level name (matches ordinary singleton rows).
      4. Defensive cross-check: every computed final ID must actually
         exist as a real `isodecoder` value in Isodecoder_counts.txt --
         anything that doesn't is treated as unresolved (skipped, logged)
         rather than silently written into the whitelist.

    Returns (family_to_final: dict mapping FAMILY-level locus string (not
    raw Family-Copy locus) -> final isodecoder ID, for families touched by
    a merge; valid_isodecoder_ids: the full set of real `isodecoder`
    values in Isodecoder_counts.txt, used by the caller to validate any
    family NOT in family_to_final -- such families default to their own
    family-level name, but that default must still be checked against
    valid_isodecoder_ids to catch anything genuinely unresolvable).
    """
    unsplit = pd.read_csv(unsplit_cluster_info_path, sep="\t")
    required_unsplit_cols = {"Parent", "Unsplit transcripts"}
    if not required_unsplit_cols.issubset(unsplit.columns):
        raise ValueError(
            f"{unsplit_cluster_info_path}: expected columns {required_unsplit_cols}, "
            f"found {list(unsplit.columns)}."
        )

    counts = pd.read_csv(isodecoder_counts_path, sep="\t")
    required_counts_cols = {"isodecoder", "Single_isodecoder", "parent"}
    if not required_counts_cols.issubset(counts.columns):
        raise ValueError(
            f"{isodecoder_counts_path}: expected columns {required_counts_cols}, "
            f"found {list(counts.columns)}."
        )
    merged_rows = counts[counts["Single_isodecoder"] == False]  # noqa: E712 (real bool dtype from True/False strings)
    valid_isodecoder_ids = set(counts["isodecoder"])

    family_to_final = {}
    unresolved_merge_groups = []
    for _, row in unsplit.iterrows():
        parent_family = _strip_copy_suffix(row["Parent"])
        unsplit_families = [
            t.strip() for t in str(row["Unsplit transcripts"]).split(",") if t.strip()
        ]
        group_families = [parent_family] + unsplit_families

        match = merged_rows[merged_rows["parent"] == parent_family]
        if len(match) != 1:
            unresolved_merge_groups.append(
                (parent_family, group_families, len(match))
            )
            continue
        final_id = match.iloc[0]["isodecoder"]
        for fam in group_families:
            family_to_final[fam] = final_id

    if unresolved_merge_groups:
        log.warning(
            f"{len(unresolved_merge_groups)} unsplitClusterInfo.txt merge group(s) "
            f"had != 1 matching Single_isodecoder=False row in Isodecoder_counts.txt "
            f"(expected exactly 1) -- these families will fall through to the default "
            f"(unmerged) lookup, which will likely fail validation below: "
            f"{unresolved_merge_groups[:5]}"
        )

    return family_to_final, valid_isodecoder_ids


def build_whitelist(anticodon_map_path, unsplit_cluster_info_path, isodecoder_counts_path,
                     i34_isotypes, q34_isotypes, out_path, log_path=None):
    df = pd.read_csv(anticodon_map_path, sep="\t")
    id_col, ac_col = _detect_columns(df)
    log.info(f"Using columns: isodecoder_id='{id_col}', anticodon='{ac_col}'")
    log.info("isotype is derived from the locus string, not a column -- see _isotype_from_locus()")

    family_to_final, valid_isodecoder_ids = _load_locus_to_isodecoder_map(
        unsplit_cluster_info_path, isodecoder_counts_path
    )
    log.info(
        f"Loaded {len(family_to_final)} merged-family entries and "
        f"{len(valid_isodecoder_ids)} total valid isodecoder IDs "
        f"({unsplit_cluster_info_path} + {isodecoder_counts_path})"
    )

    i34_set = set(i34_isotypes)
    q34_set = set(q34_isotypes)

    rows = []
    dropped_stop = []
    demoted_warnings = []
    genetic_code_rejections = []
    genetic_code_skipped = []
    unparseable_isotype = []
    raw_locus_not_in_isodecoder_map = []

    for _, r in df.iterrows():
        isodecoder_id = r[id_col]
        anticodon     = str(r[ac_col]).strip().upper()

        if len(anticodon) != 3 or any(b not in "ACGT" for b in anticodon):
            log.warning(f"Skipping {isodecoder_id}: malformed anticodon '{anticodon}'")
            continue

        # FIX: isodecoder_id at this point is still the RAW anticodon_map
        # locus (e.g. "Homo_sapiens_tRNA-Lys-TTT-3-1") -- relabel it to
        # Stage-1's FINAL, cell-line-specific collapsed isodecoder ID
        # (e.g. "Homo_sapiens_tRNA-Lys-TTT-3/5") before it's ever written
        # to an output row. Isotype/anticodon parsing below still uses the
        # RAW locus string (isodecoder_id) unchanged -- only the ID
        # actually stored in `rows` is swapped, via final_isodecoder_id.
        #
        # Resolution: strip the raw locus down to its FAMILY (drop the
        # copy suffix), look it up in family_to_final (families touched
        # by an unsplitClusterInfo.txt merge); if not merged, default to
        # the family-level name itself (ordinary singleton case). Either
        # way, the result is cross-checked against valid_isodecoder_ids --
        # anything that isn't a real Isodecoder_counts.txt value is
        # skipped rather than silently written into the whitelist (this
        # is what caught the previous, wrong join -- see module docstring).
        family = _strip_copy_suffix(isodecoder_id)
        final_isodecoder_id = family_to_final.get(family, family)
        if final_isodecoder_id not in valid_isodecoder_ids:
            raw_locus_not_in_isodecoder_map.append(isodecoder_id)
            log.warning(
                f"Skipping {isodecoder_id}: resolved family-level ID "
                f"'{final_isodecoder_id}' is not a real isodecoder in "
                f"{isodecoder_counts_path} -- likely an organelle/non-nuclear "
                f"locus out of this whitelist's scope, or a genuine data gap."
            )
            continue

        isotype = _isotype_from_locus(str(isodecoder_id), anticodon)
        if isotype == ISOTYPE_UNPARSEABLE:
            unparseable_isotype.append(isodecoder_id)
            log.warning(
                f"Skipping {isodecoder_id}: could not recover isotype from locus "
                f"string (anticodon='{anticodon}') via either GtRNAdb naming "
                f"convention -- neither parsing branch matched."
            )
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
            # I34's wobble range is hardcoded below per the inosine wobble
            # rule rather than re-derived, since inosine is not one of the
            # four standard bases _codon_from_anticodon() switches on).
            #
            # REVISED (per supervisor correction, 2026-07-16 meeting +
            # follow-up): the A-ending (I34:A) candidate is REMOVED. I34:A
            # pairing is structurally real (Crick 1966: I34 -> U,C,A) but
            # occurs at low efficiency in vivo -- Curran, J.F. (1995)
            # "Decoding with the A:I wobble pair is inefficient" (PMC306738)
            # shows, via three independent in vivo assays in E. coli, slow/
            # unstable ribosomal decoding at CGA (an I34:A-read codon), even
            # though the pairing occurs. This whitelist tracks functionally
            # meaningful decoding capacity, not every structurally-possible
            # pairing, so I34:A is dropped rather than included as a third
            # mod_only_I candidate.
            prefix, _pos3_unmod_unused = _codon_from_anticodon(anticodon)
            # T-ending: decoded by whole pool regardless of editing status -> both_I
            candidates = [("T", "both_I", "0",
                           "U-ending codon: decoded by entire pool (edited + unedited); editing is strictly additive for I34."),
                          ("C", "mod_only_I", "1",
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
                    isodecoder_id=final_isodecoder_id, isotype=isotype, anticodon=anticodon,
                    position34_base=n34, bucket=bucket, codon=codon,
                    term_type=term_type, gamma_expr=gamma_expr, notes=note,
                ))

        elif bucket == "Q34":
            # REVISED (per supervisor correction, 2026-07-16 meeting +
            # follow-up), mirroring the I34 term-type structure exactly, as
            # instructed ("revise Q34 to something similar to I34, keep the
            # kappa term for Q34"):
            #   - C-ending: native Watson-Crick pairing for G34, decoded by
            #     the WHOLE pool regardless of Q-editing status (mirrors
            #     both_I) -> term_type "both_Q_C", gamma=0 (was "kappa").
            #   - U-ending: native G34:U wobble is real but low-efficiency
            #     (Fields et al. 2022, Front. Mol. Biosci. -- unmodified G34
            #     reads U-ending "with reduced efficiency"; Q34-edited reads
            #     U- and C-ending "equally well" [Meier et al. 1985]) -- so
            #     U-ending is treated as reachable at meaningful efficiency
            #     ONLY via Q-editing (mirrors mod_only_I), but weighted by
            #     the confidence dial kappa rather than full weight 1, since
            #     Q34 detection is exploratory/low-sensitivity (see rule 11
            #     wobble_glm.R) -> term_type "mod_only_Q", gamma="kappa"
            #     (was term_type "both_Q_U", gamma=0 / no Q-dependence).
            prefix, pos3_unmod = _codon_from_anticodon(anticodon)          # {'C'} (native), U-ending is edit-gated
            candidates = [("C", "both_Q_C", "0",
                           "C-ending codon: native Watson-Crick pairing for G34, decoded by entire pool (edited + unedited) -- mirrors both_I."),
                          ("T", "mod_only_Q", "kappa",
                           "U-ending codon: reachable at meaningful efficiency only via Q34 editing (native G34:U wobble is low-efficiency, Fields et al. 2022) -- mirrors mod_only_I, weighted by kappa (Q34 confidence dial) rather than full weight.")]
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
                    isodecoder_id=final_isodecoder_id, isotype=isotype, anticodon=anticodon,
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
                    isodecoder_id=final_isodecoder_id, isotype=isotype, anticodon=anticodon,
                    position34_base=n34, bucket=bucket, codon=codon,
                    term_type="canonical", gamma_expr="0",
                    notes="No modification pathway modeled for this isodecoder (out of I34/Q34 scope, or non-eligible paralog of a modified isotype).",
                ))

    wl = pd.DataFrame(rows)

    # FIX: multiple RAW loci can now collapse to the same final
    # isodecoder_id (see module docstring "FIXED" note) -- this produces
    # exact duplicate rows (same isodecoder_id/anticodon/bucket/codon/
    # term_type) whenever the merged raw loci share an anticodon, which
    # is the expected case (mim-tRNAseq only merges loci with near-
    # identical mature sequence, so identical anticodon). Drop those
    # exact duplicates first. If merged raw loci ever disagree on
    # term_type/anticodon for the same (isodecoder_id, codon) pair --
    # which would mean mim-tRNAseq clustered together two loci with
    # different anticodons, a genuine data inconsistency rather than
    # redundant rows -- log it loudly and keep only the first (sorted)
    # occurrence, rather than silently emitting duplicate/conflicting
    # keys that would break the one-row-per-(isodecoder_id, codon)
    # assumption compute_delta_c.py depends on.
    if not wl.empty:
        n_before_dedup = len(wl)
        wl = wl.drop_duplicates().reset_index(drop=True)

        key_cols = ["isodecoder_id", "codon"]
        conflict_mask = wl.duplicated(subset=key_cols, keep=False)
        if conflict_mask.any():
            conflicts = wl.loc[conflict_mask].sort_values(key_cols)
            log.warning(
                f"{conflict_mask.sum()} whitelist rows share an (isodecoder_id, codon) "
                f"key with DIFFERING content after merging raw loci into final "
                f"isodecoders -- raw loci mim-tRNAseq clustered together do not agree "
                f"on anticodon/term_type, which should not happen. Keeping only the "
                f"first (sorted) row per key; inspect these manually:\n"
                + conflicts.to_string(index=False)
            )
            wl = wl.drop_duplicates(subset=key_cols, keep="first").reset_index(drop=True)

        n_after_dedup = len(wl)
        if n_after_dedup != n_before_dedup:
            log.info(
                f"Collapsed {n_before_dedup - n_after_dedup} duplicate rows produced by "
                f"raw-locus -> final-isodecoder merging ({n_before_dedup} -> {n_after_dedup})."
            )

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
    if raw_locus_not_in_isodecoder_map:
        summary_lines.append(
            f"Raw loci with no entry in this cell line's isodecoder map "
            f"(excluded from whitelist entirely): {len(raw_locus_not_in_isodecoder_map)} "
            f"(first 10): {raw_locus_not_in_isodecoder_map[:10]}"
            + (" ... (truncated)" if len(raw_locus_not_in_isodecoder_map) > 10 else "")
        )
    if unparseable_isotype:
        summary_lines.append(
            f"Isotype-from-locus parsing FAILED for {len(unparseable_isotype)} loci "
            f"(excluded from whitelist entirely): {unparseable_isotype[:10]}"
            + (" ... (truncated)" if len(unparseable_isotype) > 10 else "")
        )
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
        unsplit_cluster_info_path=snakemake.input.unsplit_cluster_info,
        isodecoder_counts_path=snakemake.input.isodecoder_counts,
        i34_isotypes=snakemake.params.i34_isotypes,
        q34_isotypes=snakemake.params.q34_isotypes,
        out_path=snakemake.output.whitelist,
        log_path=snakemake.log[0] if len(snakemake.log) else None,
    )
