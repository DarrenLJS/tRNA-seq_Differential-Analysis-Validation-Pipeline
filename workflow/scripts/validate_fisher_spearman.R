#!/usr/bin/env Rscript
# workflow/scripts/validate_fisher_spearman.R
#
# Proposal Section 4 validation criteria:
#   - Fisher's exact test: ISG/housekeeping representation among Watson
#     et al. UP/DOWN gene sets, odds ratio + 95% CI, per timepoint.
#   - Spearman correlation: predicted translation score vs observed
#     Watson et al. polysome log2FC, per timepoint.
#
# FIX-check on first real run: watson_fc's column names/timepoint labels
# need to line up with this pipeline's own timepoint labels (from
# sample_manifest.tsv) -- Watson et al.'s own timepoint scheme almost
# certainly does not match this project's poly(I:C) stimulation timepoints
# 1:1 (different study, different design). A timepoint alignment/mapping
# step will likely be needed here once the real Watson data is parsed
# (see fetch_watson_polysome_data.py FIX notes) -- not implemented yet
# since the source data's actual timepoint structure isn't confirmed.

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

if (!all(c("gene_id", "timepoint", "log2FC") %in% colnames(watson_fc))) {
  cat("WARNING: watson_fc does not have the expected columns (gene_id, timepoint, log2FC). ",
      "Found columns: ", paste(colnames(watson_fc), collapse = ", "), "\n",
      "This is expected until fetch_watson_polysome_data.py's column-mapping FIX is ",
      "completed against the real source file -- validation cannot proceed meaningfully ",
      "until that mapping exists. Writing an empty summary rather than crashing.\n")
  write_tsv(data.frame(), out_path)
  quit(save = "no", status = 0)
}

merged <- gene_scores %>%
  inner_join(watson_fc, by = c("gene_id", "timepoint"))

results <- list()

for (tp in unique(merged$timepoint)) {
  sub <- merged %>% filter(timepoint == tp)
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
write_tsv(combined, out_path)
cat(sprintf("Wrote %d validation rows -> %s\n", nrow(combined), out_path))

sink(type = "message")
sink(type = "output")
