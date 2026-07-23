# =============================================================================
# workflow/rules/14_percodon_score.smk
# Per-codon score Delta(c) -- the core integration step. Combines:
#   - decoding_whitelist.tsv       (rule 09)
#   - isodecoder DESeq2 FC(i)      (rule 10, high-confidence intersect set)
#   - I34 GLM f_stim/f_ctrl        (rule 11, confirmatory)
#   - Q34 GLM f_stim/f_ctrl        (rule 11, exploratory)
#   - kappa (Q34 confidence dial)
#
# Computed across the full kappa_sweep from config so rule 16 can pick the
# empirically-justified value rather than assuming one. See
# compute_delta_c.py docstring for the full formula and missing-data
# handling -- this script has been unit-tested against every hand-traced
# case from the formula derivation (Ala, Ile 3-way split, Asp, Cys/Ile-TAT
# traps) before being wired in here.
#
# EXTENDED: also computes delta_c_v2, the count/FC-based companion score
# (rule 10's codon_highconf_intersect_kappa{kappa}.tsv, built from real
# anticodon/codon-level counts rather than whitelist-weighted term sums).
# Both versions are computed unconditionally and both feed rule 16/17 --
# v2 does not replace v1. See compute_delta_c_v2.py docstring.
# =============================================================================

rule compute_delta_c:
    input:
        whitelist    = f"{STAGE2_ROOT}/references/{{cell_line}}/decoding_whitelist.tsv",
        isodecoder_fc = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/isodecoder_highconf_intersect.tsv",
        i34_glm      = f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/I34_glm_results.tsv",
        q34_glm      = f"{STAGE2_ROOT}/wobble_glm/{{cell_line}}/Q34_glm_results_exploratory.tsv",
    output:
        delta_c = f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_kappa{{kappa}}.tsv",
    log:
        f"{STAGE2_ROOT}/logs/14_percodon_score/{{cell_line}}_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("percodon_score"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/compute_delta_c.py"


rule compute_delta_c_v2:
    """
    Count/FC-based companion to compute_delta_c above -- reshapes rule
    10's codon-level DESeq2+edgeR high-confidence intersect into the same
    (timepoint, codon, delta_c, ..., kappa) schema so rule 16/17 can
    consume both versions uniformly. No whitelist term-type weighting
    here (already baked into the counts by rule 10's
    build_anticodon_count_matrix / build_codon_count_matrix). See
    compute_delta_c_v2.py docstring for the missing-data handling
    difference from v1 (no FC=1 imputation here).
    """
    input:
        codon_highconf = f"{STAGE2_ROOT}/diff_abundance/{{cell_line}}/codon_highconf_intersect_kappa{{kappa}}.tsv",
    output:
        delta_c_v2 = f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_v2_kappa{{kappa}}.tsv",
    log:
        f"{STAGE2_ROOT}/logs/14_percodon_score/{{cell_line}}_v2_kappa{{kappa}}.log",
    resources:
        sge_extra = sge_extra("percodon_score"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/compute_delta_c_v2.py"


# -----------------------------------------------------------------------
# Codon-ending stratification -- direct internal test of the central
# hypothesis (dissertation Section 1.6: G/C-ending codons favoured,
# A/U-ending disfavoured), using only Delta(c) itself. Deliberately kept
# in rule 14, not rule 16: rule 16 validates Delta(c) against the
# EXTERNAL Watson et al. polysome data (still pending); this rule tests
# whether Delta(c)'s own sign/magnitude already points the predicted
# direction, which requires no external dataset and can be reported today.
# Mann-Whitney U (unpaired rank-sum), NOT a paired Wilcoxon test -- G/C-
# and A/U-ending codons are two independent groups of codons, not paired
# observations. See codon_ending_stratification.py docstring for the full
# rationale, including why this is done per-codon rather than per-
# isodecoder.
#
# NOTE: runs on delta_c (v1) only, not delta_c_v2 -- not extended here
# since it wasn't asked for; delta_c_v2's own directionality can be
# spot-checked the same way later if useful.
# -----------------------------------------------------------------------
rule codon_ending_stratification:
    input:
        delta_c = f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/delta_c_kappa{{kappa}}.tsv",
    output:
        summary = f"{STAGE2_ROOT}/percodon_score/{{cell_line}}/codon_ending_wilcoxon_kappa{{kappa}}.tsv",
    log:
        f"{STAGE2_ROOT}/logs/14_percodon_score/{{cell_line}}_codon_ending_kappa{{kappa}}.log",
    params:
        alpha = config["codon_ending_stratification"]["alpha"],
    resources:
        sge_extra = sge_extra("codon_ending_stratification"),
    conda:
        "../../envs/stage2_python.yaml"
    script:
        "../scripts/codon_ending_stratification.py"
