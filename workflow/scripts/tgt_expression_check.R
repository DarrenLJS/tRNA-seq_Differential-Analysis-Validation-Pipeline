#!/usr/bin/env Rscript
# workflow/scripts/tgt_expression_check.R
#
# Simple stimulated-vs-control differential expression check on QTRT1/
# QTRT2 (the TGT enzyme complex), as corroborating (or non-corroborating)
# evidence for Q34 regulation, given that direct Q34 sequencing signal is
# unreliable (see wobble_glm.R). This is NOT a full DESeq2 pipeline --
# just two genes, reported with a simple t-test on normalized counts,
# since building a whole-matrix DESeq2 run for two genes is unnecessary
# overhead.
#
# DOES NOT RUN AS-IS: input$mrna_expression is a placeholder path (see
# 11_wobble_modification_glm.smk rule docstring) -- neither Stage 1 nor
# Stage 2 currently produces a genome-wide mRNA expression matrix. This
# script is written defensively (checks for the genes' presence, fails
# with a clear message rather than a cryptic error) so that once a real
# expression source is identified, only the Snakefile's rule `input:`
# needs updating, not this script.

suppressMessages({
  library(dplyr)
  library(readr)
})

expr_path    <- snakemake@input[["mrna_expression"]]
coldata_path <- snakemake@input[["coldata"]]
out_path     <- snakemake@output[["results"]]

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

if (!file.exists(expr_path)) {
  stop(paste0(
    "mRNA expression matrix not found at: ", expr_path, ". This rule requires a ",
    "genome-wide mRNA expression source that neither Stage 1 nor Stage 2 currently ",
    "produces -- see the rule docstring in 11_wobble_modification_glm.smk. Point ",
    "config/Snakefile at a real total RNA-seq quantification before running this rule."
  ))
}

expr <- read.table(expr_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata <- read.table(coldata_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata <- coldata[colnames(expr), , drop = FALSE]

target_genes <- c("QTRT1", "QTRT2")
present <- intersect(target_genes, rownames(expr))
if (length(present) == 0) {
  stop(paste0(
    "Neither QTRT1 nor QTRT2 found in the expression matrix row names. ",
    "Check gene ID convention (symbol vs Ensembl ID) matches this matrix."
  ))
}

results <- list()
for (gene in present) {
  for (tp in unique(coldata$timepoint)) {
    samples_tp <- rownames(coldata)[coldata$timepoint == tp]
    vals <- as.numeric(expr[gene, samples_tp])
    cond <- coldata[samples_tp, "condition"]
    if (length(unique(cond)) < 2) next
    tt <- tryCatch(t.test(vals ~ cond), error = function(e) NULL)
    if (is.null(tt)) next
    results[[length(results) + 1]] <- data.frame(
      gene = gene, timepoint = tp,
      mean_control = mean(vals[cond == "control"]),
      mean_stimulated = mean(vals[cond == "stimulated"]),
      pvalue = tt$p.value
    )
  }
}

final <- bind_rows(results)
write_tsv(final, out_path)
cat(sprintf("Wrote %d TGT expression check rows -> %s\n", nrow(final), out_path))

sink(type = "message")
sink(type = "output")
