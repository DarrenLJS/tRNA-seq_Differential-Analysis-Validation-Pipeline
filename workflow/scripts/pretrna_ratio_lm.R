#!/usr/bin/env Rscript
# workflow/scripts/pretrna_ratio_lm.R
#
# Linear model of pre-tRNA:mature ratio ~ condition, per locus, per
# timepoint (proposal 3.7). Reported as supplementary/exploratory --
# not gated into rule 14's high-confidence isodecoder set.
#
# Ratio is log-transformed before the LM (log(ratio + pseudocount)) since
# ratios are strictly positive and typically right-skewed -- standard
# practice, flagged explicitly here since it changes the interpretation
# of the model's coefficient (log-fold-change in ratio, not raw
# difference in ratio).

suppressMessages({
  library(dplyr)
  library(broom)
  library(readr)
})

locus_counts_path <- snakemake@input[["locus_counts"]]
coldata_path       <- snakemake@input[["coldata"]]
out_path           <- snakemake@output[["results"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

locus_counts <- read_tsv(locus_counts_path, show_col_types = FALSE)
coldata <- read.table(coldata_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE) %>%
  tibble::rownames_to_column("sample_id")

df <- locus_counts %>%
  filter(!is.na(pretrna_mature_ratio), pretrna_mature_ratio > 0) %>%
  left_join(coldata, by = "sample_id") %>%
  mutate(log_ratio = log(pretrna_mature_ratio + 1e-6))

results <- list()
for (loc in unique(df$locus_id)) {
  for (tp in unique(df$timepoint)) {
    sub <- df %>% filter(locus_id == loc, timepoint == tp)
    if (nrow(sub) < 4 || length(unique(sub$condition)) < 2) next  # need >=2 per group minimum, roughly

    fit <- tryCatch(lm(log_ratio ~ condition, data = sub), error = function(e) NULL)
    if (is.null(fit)) next

    tidy_fit <- tidy(fit) %>% filter(grepl("^condition", term))
    if (nrow(tidy_fit) == 0) next

    results[[length(results) + 1]] <- data.frame(
      locus_id = loc, timepoint = tp,
      log_ratio_effect = tidy_fit$estimate[1],
      std_error = tidy_fit$std.error[1],
      pvalue = tidy_fit$p.value[1],
      n_samples = nrow(sub)
    )
  }
}

final <- bind_rows(results)
if (nrow(final) > 0) {
  final$padj <- p.adjust(final$pvalue, method = "BH")
}
write_tsv(final, out_path)
cat(sprintf("Wrote %d locus x timepoint LM rows -> %s\n", nrow(final), out_path))
cat("NOTE: pre-tRNA:mature ratio results are supplementary/exploratory per proposal 3.7.\n")

sink(type = "message")
sink(type = "output")
