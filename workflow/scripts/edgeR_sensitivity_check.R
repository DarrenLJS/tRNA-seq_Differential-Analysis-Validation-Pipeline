#!/usr/bin/env Rscript
# workflow/scripts/edgeR_sensitivity_check.R
#
# Independent edgeR differential abundance run on the same isodecoder
# count matrix DESeq2 used, as a sensitivity check (proposal 3.5). Uses
# edgeR's standard QLF (quasi-likelihood F-test) pipeline, per-timepoint
# pairwise stimulated-vs-control contrasts only (the LRT-across-timepoints
# comparison is DESeq2-only per rule 10's deseq2_isodecoder.R -- edgeR's
# role here is strictly a cross-check on the pairwise calls, not a second
# independent model of the interaction term).

suppressMessages({
  library(edgeR)
  library(dplyr)
  library(tibble)
  library(readr)
})

counts_path  <- snakemake@input[["counts"]]
coldata_path <- snakemake@input[["coldata"]]
out_results  <- snakemake@output[["results"]]
fdr          <- snakemake@params[["fdr"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

counts  <- read.table(counts_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata <- read.table(coldata_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata <- coldata[colnames(counts), , drop = FALSE]
coldata$condition <- factor(coldata$condition)
coldata$timepoint <- factor(coldata$timepoint)

timepoints <- levels(coldata$timepoint)
all_results <- list()

for (tp in timepoints) {
  sub_coldata <- coldata[coldata$timepoint == tp, , drop = FALSE]
  sub_counts  <- counts[, rownames(sub_coldata), drop = FALSE]

  if (length(unique(sub_coldata$condition)) < 2) {
    cat(sprintf("Skipping timepoint '%s': fewer than 2 condition levels present\n", tp))
    next
  }

  y <- DGEList(counts = sub_counts, group = sub_coldata$condition)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~condition, data = sub_coldata)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))  # last coef = condition effect

  res_df <- topTags(qlf, n = Inf)$table %>%
    rownames_to_column("isodecoder_id") %>%
    mutate(timepoint = tp, test = "edgeR_QLF_pairwise")

  all_results[[tp]] <- res_df
}

final <- bind_rows(all_results)
write_tsv(final, out_results)
cat(sprintf("Wrote %d rows -> %s\n", nrow(final), out_results))

sink(type = "message")
sink(type = "output")
