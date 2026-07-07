#!/usr/bin/env Rscript
# workflow/scripts/deseq2_trf.R
#
# DESeq2 differential abundance per tRF class (5'-tRF/3'-tRF/i-tRF/tiRNA),
# stimulated vs control, per timepoint. Reads the per-class matrices
# written by parse_trax_tRF_classes.py. Classes with too few features or
# too few samples are skipped with a logged warning rather than crashing
# the whole rule (tiRNA in particular may have very few features).

suppressMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
})

class_dir   <- snakemake@input[["class_matrices_dir"]]
coldata_path <- snakemake@input[["coldata"]]
out_path    <- snakemake@output[["results"]]
fdr         <- snakemake@params[["fdr"]]
trf_classes <- snakemake@params[["trf_classes"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

coldata <- read.table(coldata_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata$condition <- factor(coldata$condition)
coldata$timepoint <- factor(coldata$timepoint)

all_results <- list()

for (cls in trf_classes) {
  mat_path <- file.path(class_dir, paste0(cls, "_counts_matrix.tsv"))
  if (!file.exists(mat_path)) {
    cat(sprintf("Skipping class '%s': matrix file not found (%s)\n", cls, mat_path))
    next
  }
  counts <- tryCatch(
    read.table(mat_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(counts) || nrow(counts) < 3) {
    cat(sprintf("Skipping class '%s': fewer than 3 features after filtering (n=%s)\n",
                cls, ifelse(is.null(counts), "NA", nrow(counts))))
    next
  }

  sub_coldata <- coldata[colnames(counts), , drop = FALSE]

  for (tp in levels(sub_coldata$timepoint)) {
    tp_coldata <- sub_coldata[sub_coldata$timepoint == tp, , drop = FALSE]
    tp_counts  <- counts[, rownames(tp_coldata), drop = FALSE]
    if (length(unique(tp_coldata$condition)) < 2 || ncol(tp_counts) < 4) {
      cat(sprintf("Skipping class '%s', timepoint '%s': insufficient design (n=%d, conditions=%d)\n",
                  cls, tp, ncol(tp_counts), length(unique(tp_coldata$condition))))
      next
    }

    dds <- tryCatch({
      d <- DESeqDataSetFromMatrix(countData = tp_counts, colData = tp_coldata, design = ~condition)
      DESeq(d)
    }, error = function(e) {
      cat(sprintf("DESeq2 failed for class '%s', timepoint '%s': %s\n", cls, tp, conditionMessage(e)))
      NULL
    })
    if (is.null(dds)) next

    res <- results(dds, alpha = fdr)
    res_df <- as.data.frame(res) %>%
      rownames_to_column("trf_id") %>%
      mutate(trf_class = cls, timepoint = tp)
    all_results[[paste(cls, tp)]] <- res_df
  }
}

final <- bind_rows(all_results)
write_tsv(final, out_path)
cat(sprintf("Wrote %d tRF DESeq2 result rows across %d classes -> %s\n",
            nrow(final), length(unique(final$trf_class)), out_path))
if ("tiRNA" %in% final$trf_class) {
  cat("NOTE: tiRNA class results present -- flagged of particular interest ",
      "(stress-induced translation-initiation inhibitor), see proposal 3.8.\n")
}

sink(type = "message")
sink(type = "output")
