#!/usr/bin/env Rscript
# workflow/scripts/kappa_sweep_summary.R
#
# Sweeps kappa across config wobble_glm.kappa_sweep, recomputing the
# Spearman correlation (predicted translation score vs Watson et al.
# observed log2FC) at each value, across all cell lines and timepoints.
# Reports the kappa that maximizes mean |rho| (and separately, the kappa
# that maximizes mean rho with correct sign, since a strong NEGATIVE
# correlation would indicate a modelling error, not evidence to prefer
# that kappa) -- if kappa=0 comes out best, that is reported as a
# legitimate finding, not suppressed.

suppressMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

kappa_values <- snakemake@params[["kappa_values"]]
cell_lines   <- snakemake@params[["cell_lines"]]
stage2_root  <- snakemake@params[["stage2_root"]]
watson_fc_path <- snakemake@input[["watson_fc"]]
out_path     <- snakemake@output[["sweep_summary"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

watson_fc <- read_tsv(watson_fc_path, show_col_types = FALSE)

if (!all(c("gene_id", "timepoint", "log2FC") %in% colnames(watson_fc))) {
  cat("WARNING: watson_fc missing expected columns -- see validate_fisher_spearman.R FIX note. ",
      "Writing empty kappa sweep summary.\n")
  write_tsv(data.frame(), out_path)
  quit(save = "no", status = 0)
}

results <- list()

for (kappa in kappa_values) {
  for (cl in cell_lines) {
    scores_path <- file.path(stage2_root, "gene_prediction", cl, sprintf("gene_translation_scores_kappa%s.tsv", kappa))
    if (!file.exists(scores_path)) {
      cat(sprintf("Missing gene scores file for kappa=%s, cell_line=%s: %s\n", kappa, cl, scores_path))
      next
    }
    scores <- read_tsv(scores_path, show_col_types = FALSE)
    merged <- scores %>% inner_join(watson_fc, by = c("gene_id", "timepoint"))

    for (tp in unique(merged$timepoint)) {
      sub <- merged %>% filter(timepoint == tp)
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

write_tsv(final, out_path)
cat(sprintf("Wrote %d kappa sweep rows -> %s\n", nrow(final), out_path))

sink(type = "message")
sink(type = "output")
