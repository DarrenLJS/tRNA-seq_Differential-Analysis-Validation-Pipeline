#!/usr/bin/env Rscript
# workflow/scripts/kappa_sweep_summary.R
#
# Sweeps kappa across config wobble_glm.kappa_sweep, recomputing the
# Spearman correlation (predicted translation score vs Watson et al.
# observed log2FC) at each value, across all cell lines and this
# pipeline's own tRNA-seq timepoints. Reports the kappa that maximizes
# mean |rho| (and separately, the kappa that maximizes mean rho with
# correct sign, since a strong NEGATIVE correlation would indicate a
# modelling error, not evidence to prefer that kappa) -- if kappa=0
# comes out best, that is reported as a legitimate finding, not
# suppressed.
#
# SINGLE EXTERNAL TIMEPOINT, NOT A TIMECOURSE -- see
# validate_fisher_spearman.R header and fetch_watson_polysome_data.py
# docstring for the full derivation. watson_fc has no timepoint column
# (Watson et al.'s poly(I:C) stimulation was sequenced at a single 4h
# timepoint only), so it is joined on gene_id alone; each of this
# pipeline's own timepoints is swept against that same single benchmark.

suppressMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

kappa_values <- snakemake@params[["kappa_values"]]
cell_lines   <- snakemake@params[["cell_lines"]]
stage2_root  <- snakemake@params[["stage2_root"]]
# NEW: score_version selects which gene_translation_scores file pattern to
# read -- "v1" (rate/GLM-based Delta(c), rule 14 compute_delta_c) or "v2"
# (count/FC-based delta_c_v2, rule 14 compute_delta_c_v2). Required, no
# default -- every call site (rule 16's validate_kappa_sweep AND its new
# validate_kappa_sweep_v2 sibling) must say explicitly which version it's
# sweeping, rather than one of them silently relying on an implicit default
# that happens to match.
score_version  <- snakemake@params[["score_version"]]
if (!score_version %in% c("v1", "v2")) {
  stop(sprintf("score_version must be 'v1' or 'v2', got: '%s'", score_version))
}
watson_fc_path <- snakemake@input[["watson_fc"]]
out_path     <- snakemake@output[["sweep_summary"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

watson_fc <- read_tsv(watson_fc_path, show_col_types = FALSE)

if (!all(c("gene_id", "log2FC") %in% colnames(watson_fc))) {
  cat("WARNING: watson_fc missing expected columns (gene_id, log2FC) -- see ",
      "validate_fisher_spearman.R / fetch_watson_polysome_data.py notes. ",
      "Writing empty kappa sweep summary.\n")
  write_tsv(data.frame(), out_path)
  quit(save = "no", status = 0)
}

if ("timepoint" %in% colnames(watson_fc)) {
  cat("WARNING: watson_fc unexpectedly has a 'timepoint' column -- this script joins ",
      "on gene_id only, per Watson et al.'s single-4h-timepoint design. If watson_fc ",
      "now legitimately carries multiple external timepoints, update this script to ",
      "join on both keys.\n")
}

results <- list()

for (kappa in kappa_values) {
  for (cl in cell_lines) {
    # FIX (2026-07-19): sprintf("%s", kappa) on a numeric silently drops the
    # trailing ".0" for whole-number kappas (R's as.character(0.0) == "0",
    # as.character(1.0) == "1"), producing paths like ".../kappa0.tsv" and
    # ".../kappa1.tsv" that never exist -- the real files are
    # ".../kappa0.0.tsv" and ".../kappa1.0.tsv" (Python's str(0.0) == "0.0",
    # which is what actually generated them, via Snakemake's own wildcard
    # expansion in rule validate_kappa_sweep / expand()). file.exists() just
    # returned FALSE for those two paths and the loop silently `next`-ed
    # past them with only a `cat()` note, not an error -- so kappa=0.0 and
    # kappa=1.0 were completely missing from every past
    # kappa_sweep_summary.tsv, while 0.1/0.2/0.3/0.5/0.7 (never whole
    # numbers) were unaffected. %.1f formatting matches the actual
    # filenames for every value in config wobble_glm.kappa_sweep.
    scores_path <- file.path(
      stage2_root, "gene_prediction", cl,
      if (score_version == "v2") {
        sprintf("gene_translation_scores_v2_kappa%.1f.tsv", kappa)
      } else {
        sprintf("gene_translation_scores_kappa%.1f.tsv", kappa)
      }
    )
    if (!file.exists(scores_path)) {
      cat(sprintf("Missing gene scores file for kappa=%s, cell_line=%s: %s\n", kappa, cl, scores_path))
      next
    }
    scores <- read_tsv(scores_path, show_col_types = FALSE)

    for (tp in unique(scores$timepoint)) {
      sub_scores <- scores %>% filter(timepoint == tp)
      # Join on gene_id only -- watson_fc has no timepoint dimension.
      sub <- sub_scores %>% inner_join(watson_fc, by = "gene_id")
      if (nrow(sub) < 10) next
      sp <- suppressWarnings(cor.test(sub$predicted_translation_score, sub$log2FC, method = "spearman"))
      results[[length(results) + 1]] <- data.frame(
        kappa = kappa, cell_line = cl, timepoint = tp,
        rho = unname(sp$estimate), pvalue = sp$p.value, n_genes = nrow(sub)
      )
    }
  }
}

final <- bind_rows(results)

if (nrow(final) > 0) {
  summary_by_kappa <- final %>%
    group_by(kappa) %>%
    summarise(
      mean_rho = mean(rho, na.rm = TRUE),
      mean_abs_rho = mean(abs(rho), na.rm = TRUE),
      n_comparisons = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_rho))

  best_kappa_abs <- summary_by_kappa$kappa[1]
  cat(sprintf("Kappa with highest mean |rho|: %s (mean_abs_rho=%.4f)\n",
              best_kappa_abs, summary_by_kappa$mean_abs_rho[1]))

  best_kappa_signed <- summary_by_kappa %>% filter(mean_rho == max(mean_rho, na.rm = TRUE)) %>% pull(kappa)
  cat(sprintf("Kappa with highest mean rho (correct-sign preference): %s\n", paste(best_kappa_signed, collapse=",")))

  if (best_kappa_abs == 0 || (0 %in% best_kappa_signed)) {
    cat("NOTE: kappa=0 (i.e. ignoring Q34 signal entirely) is at or near the top of ",
        "the sweep -- this is a legitimate finding consistent with the assay-",
        "sensitivity discussion (standard-RT misincorporation for queuosine is close ",
        "to background noise), not a failure of the sweep.\n")
  }

  final <- final %>% left_join(summary_by_kappa, by = "kappa")
}

final$watson_benchmark_timepoint <- "4h_single_timepoint"
final$score_version <- score_version
write_tsv(final, out_path)
cat(sprintf("Wrote %d kappa sweep rows -> %s\n", nrow(final), out_path))
cat("NOTE: 'timepoint' is this pipeline's own tRNA-seq timepoint; all rows are swept ",
    "against the SAME single 4h Watson et al. external benchmark.\n")

sink(type = "message")
sink(type = "output")
