# =============================================================================
# tRNA-seq Analysis Pipeline — STAGE 2
# Project: Characterising tRNA library dynamics in the antiviral response
# Author:  Darren Lim
# =============================================================================
#
# Stage 2 is CHAINED to Stage 1 — it reads Stage 1's already-built outputs
# directly from {scratch}/... and never re-runs alignment, mim-tRNAseq, or
# TRAX. All Stage-2 outputs are written under {scratch}/stage2/ so nothing
# here can collide with or overwrite a Stage-1 output directory.
#
# Pipeline stages
# ---------------
#   09  Reference prep (Stage 2)   (decoding whitelist, codon usage table,
#                                    Watson et al. polysome data, ISG/HK lists)
#   10  Isodecoder/isoacceptor DE  (DESeq2 + edgeR sensitivity check)
#   11  Wobble modification GLM   (I34 confirmatory, Q34 exploratory)
#   12  Pre-tRNA:mature ratio      (per-locus featureCounts + linear model)
#   13  tRF differential abundance (DESeq2 on TRAX per-class counts)
#   14  Per-codon score Delta(c)   (whitelist + FC(i) + f(i) + kappa)
#   15  Gene translational prediction (codon usage . Delta(c))
#   16  Validation                 (Fisher's exact + Spearman vs Watson et al.)
#   17  Cross-cell-line concordance (final gate: A549 vs THP1)
#
# Conda environments
# ------------------
#   envs/r_stats.yaml     — R stack (DESeq2, edgeR, GLM/emmeans, Fisher/Spearman).
#   envs/stage2_python.yaml — Stage-2-owned Python stack (pandas/numpy/scipy/
#   requests, plus gffread and subread/featureCounts) for all Python-only
#   rules (whitelist building, Delta(c), parsing, featureCounts). Kept
#   entirely separate from Stage 1's environment.yaml so that Stage 2's
#   dependencies can never trigger a rebuild/rerun of Stage 1's already-
#   completed pipeline via Snakemake's software-env rerun trigger.
#
# Usage
# -----
#   snakemake -s Snakefile_stage2 -n --use-conda --cores 32
#   snakemake -s Snakefile_stage2 --use-conda --cores 32 --rerun-incomplete \
#             --keep-going --latency-wait 60
#
# =============================================================================

import pandas as pd
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
configfile: "config/config_stage2.yaml"

SCRATCH      = config["scratch"]
STAGE2_ROOT  = config["stage2_root"].format(scratch=SCRATCH)


# ---------------------------------------------------------------------------
# SGE resource helper — identical pattern to Stage 1's Snakefile so that
# fixing an OOM kill or under-provisioning only ever requires editing
# config["resources"][rule_name] in config_stage2.yaml, never a .smk file.
# ---------------------------------------------------------------------------
def sge_extra(rule_name):
    r     = config["resources"][rule_name]
    slots = r.get("slots", 1)
    vmem  = r["vmem_mb"]
    pe    = f" -pe sharedmem {slots}" if slots > 1 else ""
    return f"-V{pe} -l h_vmem={vmem}M"


# ---------------------------------------------------------------------------
# Load sample manifest — Stage 2's OWN copy (sample_manifest.tsv at this
# repo's root), not read from Stage 1's directory. This is a static,
# hand-authored file (not a Stage 1 pipeline output), so Stage 2 keeps its
# own copy rather than reaching across into Stage 1's repo for it. If the
# manifest is ever updated (new samples, corrected metadata), re-copy it
# here manually to keep both stages in sync -- Stage 2 does not read
# Stage 1's copy at runtime.
#
# Resolved relative to this Snakefile's own directory (workflow.basedir)
# rather than the current working directory, so this works regardless of
# where `snakemake` is invoked from.
# ---------------------------------------------------------------------------
_manifest_path = config["manifest"]
if not os.path.isabs(_manifest_path):
    _manifest_path = os.path.join(workflow.basedir, _manifest_path)

manifest = (
    pd.read_csv(_manifest_path, sep="\t", index_col="sample_id")
)

SAMPLES    = manifest.index.tolist()
CELL_LINES = sorted(manifest["cell_line"].unique().tolist())
TIMEPOINTS = sorted(manifest["timepoint"].unique().tolist()) if "timepoint" in manifest.columns else []


def samples_for(cell_line):
    """Return list of sample_ids belonging to a given cell line."""
    return manifest[manifest["cell_line"] == cell_line].index.tolist()


def samples_for_timepoint(cell_line, timepoint):
    """Return list of sample_ids for a given cell line at a given timepoint."""
    sub = manifest[(manifest["cell_line"] == cell_line) & (manifest["timepoint"] == timepoint)]
    return sub.index.tolist()


# ---------------------------------------------------------------------------
# Include rule modules
# ---------------------------------------------------------------------------
include: "workflow/rules/09_reference_prep_stage2.smk"
include: "workflow/rules/10_isodecoder_isoacceptor_de.smk"
include: "workflow/rules/11_wobble_modification_glm.smk"
include: "workflow/rules/12_pretrna_maturation.smk"
include: "workflow/rules/13_trf_differential_abundance.smk"
include: "workflow/rules/14_percodon_score.smk"
include: "workflow/rules/15_gene_translational_prediction.smk"
include: "workflow/rules/16_validation.smk"
include: "workflow/rules/17_cross_cell_line_concordance.smk"

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        # ── Rule 09: reference prep ───────────────────────────────────────
        expand(
            f"{STAGE2_ROOT}/references/{{cell_line}}/decoding_whitelist.tsv",
            cell_line=CELL_LINES,
        ),
        f"{STAGE2_ROOT}/references/gene_codon_usage.tsv",
        f"{STAGE2_ROOT}/references/isg_housekeeping_gene_sets.tsv",
        f"{STAGE2_ROOT}/references/watson_polysome_foldchange.tsv",

        # ── Rule 10: isodecoder/isoacceptor DE ────────────────────────────
        expand(
            f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_DESeq2_results.tsv",
            cell_line=CELL_LINES,
        ),
        expand(
            f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isoacceptor_DESeq2_results.tsv",
            cell_line=CELL_LINES,
        ),
        expand(
            f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_highconf_intersect.tsv",
            cell_line=CELL_LINES,
        ),

        # ── Rule 10 extension: count-based modification score branch ──────
        expand(
            f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/isodecoder_mismatch_table.tsv",
            cell_line=CELL_LINES,
        ),
        expand(
            f"{STAGE2_ROOT}/count_mod_score/{{cell_line}}/anticodon_count_table_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),
        expand(
            f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/anticodon_highconf_intersect_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),
        expand(
            f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_highconf_intersect_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),

        # ── Rule 14: per-codon score, count-based companion (delta_c_v2) ──
        expand(
            f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_v2_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),

        # ── Rule 11: wobble modification GLM ──────────────────────────────
        expand(
            f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/I34_glm_results.tsv",
            cell_line=CELL_LINES,
        ),
        expand(
            f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/Q34_glm_results_exploratory.tsv",
            cell_line=CELL_LINES,
        ),

        # ── Rule 12: pre-tRNA:mature ratio ────────────────────────────────
        expand(
            f"{STAGE2_ROOT}/pretrna_ratio/{{cell_line}}/pretrna_mature_ratio_lm.tsv",
            cell_line=CELL_LINES,
        ),

        # ── Rule 13: tRF differential abundance ───────────────────────────
        expand(
            f"{STAGE2_ROOT}/trf_diff_abundance/{{cell_line}}/trf_DESeq2_results.tsv",
            cell_line=CELL_LINES,
        ),

        # ── Rule 14: per-codon score ───────────────────────────────────────
        expand(
            f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),

        # ── Rule 14: codon-ending stratification (G/C vs A/U hypothesis
        #    test, internal to Delta(c) -- see 14_percodon_score.smk) ──────
        expand(
            f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/codon_ending_wilcoxon_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),

        # ── Rule 15: gene translational prediction ────────────────────────
        expand(
            f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),
        expand(
            f"{STAGE2_ROOT}/gene_prediction/{{cell_line}}/gene_translation_scores_v2_kappa{{kappa}}.tsv",
            cell_line=CELL_LINES,
            kappa=config["wobble_glm"]["kappa_sweep"],
        ),

        # ── Rule 16: validation ────────────────────────────────────────────
        expand(
            f"{STAGE2_ROOT}/validation/{{cell_line}}/validation_summary.tsv",
            cell_line=CELL_LINES,
        ),
        expand(
            f"{STAGE2_ROOT}/validation/{{cell_line}}/validation_summary_v2.tsv",
            cell_line=CELL_LINES,
        ),
        f"{STAGE2_ROOT}/validation/kappa_sweep_summary.tsv",
        f"{STAGE2_ROOT}/validation/kappa_sweep_summary_v2.tsv",

        # ── Rule 17: cross-cell-line concordance (final gate) ─────────────
        f"{STAGE2_ROOT}/concordance/cross_cell_line_concordance_summary.tsv",
        f"{STAGE2_ROOT}/concordance/cross_cell_line_concordance_summary_v2.tsv",

