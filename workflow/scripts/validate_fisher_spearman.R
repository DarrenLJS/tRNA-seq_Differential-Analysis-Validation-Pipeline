#!/usr/bin/env Rscript
# workflow/scripts/validate_fisher_spearman.R
#
# Proposal Section 4 validation criteria:
#   - Fisher's exact test: ISG/housekeeping representation among Watson
#     et al. UP/DOWN gene sets, odds ratio + 95% CI.
#   - Spearman correlation: predicted translation score vs observed
#     Watson et al. polysome log2FC.
# Both computed per tRNA-seq timepoint (this pipeline's own 2h/4h/8h
# design, from gene_scores$timepoint).
#
# SINGLE EXTERNAL TIMEPOINT, NOT A TIMECOURSE -- SCOPED 2026-07 (see
# fetch_watson_polysome_data.py docstring for full derivation). Watson et
# al.'s poly(I:C) stimulation was sequenced at a single 4h timepoint only
# -- watson_fc has NO timepoint column (columns: gene_id, log2FC,
# [padj], source). It is therefore NOT joined on timepoint. Instead, each
# of this pipeline's own timepoints (2h/4h/8h) is validated independently
# against that one 4h external benchmark: for every tRNA-seq timepoint tp,
# gene_scores is filtered to tp, then inner-joined to the (timepoint-less)
# watson_fc on gene_id alone. This means the SAME Watson benchmark values
# are compared against three different predicted-score snapshots -- an
# intentional design choice (Delta(c) predictions at 2h/4h/8h are each
# separately tested for agreement with the single 4h ground truth,
# rather than requiring a timepoint match that doesn't exist in the
# source data), not an error. Interpret the 4h comparison as the most
# directly comparable one; 2h/8h comparisons show whether the predicted
# signal is already present earlier or still present later than the
# externally measured snapshot.

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

gene_scores_path <- snakemake@input[["gene_scores"]]
watson_fc_path   <- snakemake@input[["watson_fc"]]
gene_sets_path   <- snakemake@input[["gene_sets"]]
out_path         <- snakemake@output[["summary"]]
fisher_alpha     <- snakemake@params[["fisher_alpha"]]
spearman_alpha   <- snakemake@params[["spearman_alpha"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

gene_scores <- read_tsv(gene_scores_path, show_col_types = FALSE)
watson_fc   <- read_tsv(watson_fc_path, show_col_types = FALSE)
gene_sets   <- read_tsv(gene_sets_path, show_col_types = FALSE)

if (!all(c("gene_id", "log2FC") %in% colnames(watson_fc))) {
  cat("WARNING: watson_fc does not have the expected columns (gene_id, log2FC). ",
      "Found columns: ", paste(colnames(watson_fc), collapse = ", "), "\n",
      "This is expected until fetch_watson_polysome_data.py's COLUMN_MAP is filled ",
      "in against the real NAR supplementary sheet -- validation cannot proceed ",
      "meaningfully until that mapping exists. Writing an empty summary rather than ",
      "crashing.\n")
  write_tsv(data.frame(), out_path)
  quit(save = "no", status = 0)
}

if ("timepoint" %in% colnames(watson_fc)) {
  cat("WARNING: watson_fc unexpectedly has a 'timepoint' column -- Watson et al.'s ",
      "design is a single 4h poly(I:C) timepoint (see fetch_watson_polysome_data.py ",
      "docstring), so this script joins on gene_id only and ignores any timepoint ",
      "column present in watson_fc. If this file now legitimately contains multiple ",
      "external timepoints, this script needs updating to join on both keys again.\n")
}

results <- list()

for (tp in unique(gene_scores$timepoint)) {
  sub_scores <- gene_scores %>% filter(timepoint == tp)

  # Join on gene_id only -- watson_fc has no timepoint dimension (single
  # 4h external benchmark applied to every tRNA-seq timepoint).
  sub <- sub_scores %>%
    inner_join(watson_fc, by = "gene_id")

  if (nrow(sub) < 10) {
    cat(sprintf("Skipping timepoint '%s': too few overlapping genes (n=%d) with Watson et al. data\n", tp, nrow(sub)))
    next
  }

  # ---- Spearman correlation ----
  sp <- suppressWarnings(cor.test(sub$predicted_translation_score, sub$log2FC, method = "spearman"))

  # ---- Fisher's exact: ISG/housekeeping representation among UP/DOWN ----
  sub_labeled <- sub %>%
    mutate(direction = ifelse(log2FC > 0, "UP", "DOWN")) %>%
    left_join(gene_sets, by = "gene_id")

  fisher_rows <- list()
  for (gs in c("ISG", "housekeeping")) {
    tab <- table(
      in_set = sub_labeled$gene_set == gs & !is.na(sub_labeled$gene_set),
      direction = sub_labeled$direction
    )
    if (all(dim(tab) == c(2, 2))) {
      ft <- fisher.test(tab)
      fisher_rows[[gs]] <- data.frame(
        gene_set = gs, timepoint = tp,
        odds_ratio = unname(ft$estimate),
        ci_lower = ft$conf.int[1], ci_upper = ft$conf.int[2],
        pvalue = ft$p.value
      )
    } else {
      cat(sprintf("Skipping Fisher's exact for gene_set '%s', timepoint '%s': degenerate 2x2 table\n", gs, tp))
    }
  }

  fisher_df <- bind_rows(fisher_rows)
  spearman_df <- data.frame(
    timepoint = tp, test = "spearman",
    rho = unname(sp$estimate), pvalue = sp$p.value, n_genes = nrow(sub)
  )

  results[[tp]] <- list(fisher = fisher_df, spearman = spearman_df)
}

fisher_all <- bind_rows(lapply(results, function(x) x$fisher))
spearman_all <- bind_rows(lapply(results, function(x) x$spearman))

fisher_all$significant <- !is.na(fisher_all$pvalue) & fisher_all$pvalue < fisher_alpha
spearman_all$significant <- !is.na(spearman_all$pvalue) & spearman_all$pvalue < spearman_alpha

combined <- bind_rows(
  fisher_all %>% mutate(test = "fisher_exact"),
  spearman_all
)
combined$watson_benchmark_timepoint <- "4h_single_timepoint"
write_tsv(combined, out_path)
cat(sprintf("Wrote %d validation rows -> %s\n", nrow(combined), out_path))
cat("NOTE: every row's 'timepoint' is this pipeline's own tRNA-seq timepoint; all ",
    "rows are validated against the SAME single 4h Watson et al. external benchmark ",
    "(see watson_benchmark_timepoint column) -- not a per-external-timepoint match.\n")

sink(type = "message")
sink(type = "output")
