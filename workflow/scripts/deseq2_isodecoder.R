#!/usr/bin/env Rscript
# workflow/scripts/deseq2_isodecoder.R
#
# DESeq2 differential abundance on either the isodecoder or isoacceptor
# count matrix (level switches via snakemake@params$level -- same script
# handles both rules 10's deseq2_isodecoder and deseq2_isoacceptor targets).
#
# Model: pairwise stimulated-vs-control contrast PER TIMEPOINT (proposal
# 3.5), plus a likelihood-ratio-test (LRT) across timepoints using
# diff_abundance.lrt_model from config (default "~ condition * timepoint"
# full model vs "~ condition + timepoint" reduced model, testing for a
# condition:timepoint interaction).
#
# coldata.tsv is expected (from Stage 1's collate_counts.py) to have at
# least: sample_id (rownames), condition, timepoint, replicate. FIX-check
# these exact column names against the real coldata.tsv on first run --
# Stage 1's collate_counts.py builds this from the manifest, so the
# columns should match manifest columns 1:1, but confirm before trusting.

suppressMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
})

counts_path  <- snakemake@input[["counts"]]
coldata_path <- snakemake@input[["coldata"]]
out_results  <- snakemake@output[["results"]]
out_rds      <- snakemake@output[["rds"]]
lrt_model    <- snakemake@params[["lrt_model"]]
fdr          <- snakemake@params[["fdr"]]
level        <- snakemake@params[["level"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat(sprintf("Running DESeq2 at level='%s'\n", level))

counts  <- read.table(counts_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata <- read.table(coldata_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

stopifnot(all(colnames(counts) %in% rownames(coldata)))
coldata <- coldata[colnames(counts), , drop = FALSE]

required_cols <- c("condition", "timepoint")
missing_cols <- setdiff(required_cols, colnames(coldata))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "coldata.tsv is missing required column(s): %s. Found columns: %s. ",
    paste(missing_cols, collapse = ", "), paste(colnames(coldata), collapse = ", ")
  ))
}

coldata$condition <- factor(coldata$condition)
coldata$timepoint <- factor(coldata$timepoint)

# ---- Per-timepoint pairwise contrasts (stimulated vs control) ----
timepoints <- levels(coldata$timepoint)
pairwise_results <- list()

for (tp in timepoints) {
  sub_coldata <- coldata[coldata$timepoint == tp, , drop = FALSE]
  sub_counts  <- counts[, rownames(sub_coldata), drop = FALSE]

  if (length(unique(sub_coldata$condition)) < 2) {
    cat(sprintf("Skipping timepoint '%s': fewer than 2 condition levels present\n", tp))
    next
  }

  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_coldata, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = fdr)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("isodecoder_id") %>%
    mutate(timepoint = tp, level = level, test = "pairwise_stim_vs_ctrl")

  pairwise_results[[tp]] <- res_df
}

all_pairwise <- bind_rows(pairwise_results)

# ---- LRT across timepoints (condition:timepoint interaction) ----
# FIX: lrt_model string from config is currently NOT parsed into an actual
# R formula distinction between full/reduced model -- this needs the
# reduced model spelled out explicitly. Implemented here with the
# straightforward convention (full = condition*timepoint, reduced =
# condition+timepoint) rather than trying to parse an arbitrary formula
# string, since silently eval-ing a config string into a model formula is
# a bigger footgun than hardcoding the one comparison the proposal asks for.
dds_full <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition * timepoint)
dds_lrt  <- DESeq(dds_full, test = "LRT", reduced = ~condition + timepoint)
res_lrt  <- results(dds_lrt, alpha = fdr)
res_lrt_df <- as.data.frame(res_lrt) %>%
  rownames_to_column("isodecoder_id") %>%
  mutate(timepoint = "all_LRT", level = level, test = "LRT_condition_x_timepoint")

final <- bind_rows(all_pairwise, res_lrt_df)
write_tsv(final, out_results)
saveRDS(dds_full, out_rds)

cat(sprintf("Wrote %d rows -> %s\n", nrow(final), out_results))

sink(type = "message")
sink(type = "output")
