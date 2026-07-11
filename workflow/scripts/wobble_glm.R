#!/usr/bin/env Rscript
# workflow/scripts/wobble_glm.R
#
# Binomial GLM (successes = mismatch count at position 34, trials =
# coverage at position 34) per isodecoder, stimulated vs control, per
# timepoint. Shared between rules wobble_glm_i34 and wobble_glm_q34 --
# params$bucket selects which isodecoders from the whitelist are tested
# (I34 or Q34), params$confidence_tier only affects output labeling/
# framing, not the statistical procedure itself (same GLM either way --
# the DIFFERENCE in trustworthiness is an assay-sensitivity fact discussed
# in the project methods, not something this script should try to encode
# as a different model).
#
# coverage/mismatch matrices are Stage 1 outputs: rows = isodecoder_id,
# cols = sample_id (see 08_data_prep.smk). FIX-check this orientation
# against the real files on first run.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(emmeans)
})

coverage_path  <- snakemake@input[["coverage"]]
mismatch_path  <- snakemake@input[["mismatch"]]
coldata_path   <- snakemake@input[["coldata"]]
whitelist_path <- snakemake@input[["whitelist"]]
out_path       <- snakemake@output[["results"]]
fdr            <- snakemake@params[["fdr"]]
bucket         <- snakemake@params[["bucket"]]           # "I34" or "Q34"
confidence_tier <- snakemake@params[["confidence_tier"]]  # "confirmatory" or "exploratory"

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat(sprintf("Running wobble GLM: bucket=%s, confidence_tier=%s\n", bucket, confidence_tier))

coverage <- read.table(coverage_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
mismatch <- read.table(mismatch_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
coldata  <- read.table(coldata_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
whitelist <- read_tsv(whitelist_path, show_col_types = FALSE)

target_isodecoders <- whitelist %>% filter(bucket == !!bucket) %>% pull(isodecoder_id) %>% unique()
cat(sprintf("Isodecoders in bucket '%s': %d\n", bucket, length(target_isodecoders)))

target_isodecoders <- intersect(target_isodecoders, rownames(coverage))
cat(sprintf("...of which present in coverage matrix: %d\n", length(target_isodecoders)))

if (length(target_isodecoders) == 0) {
  stop(sprintf(
    paste0(
      "No isodecoders from bucket '%s' found in coverage matrix. Check that ",
      "decoding_whitelist.tsv isodecoder_id values match pos34_coverage_matrix.tsv ",
      "row names exactly (naming convention mismatch is the most likely cause)."
    ),
    bucket
  ))
}

coldata <- coldata[colnames(coverage), , drop = FALSE]
coldata$condition <- factor(coldata$condition)
coldata$timepoint <- factor(coldata$timepoint)

results <- list()

for (iso in target_isodecoders) {
  if (!(iso %in% rownames(mismatch))) next
  cov_row <- as.numeric(coverage[iso, ])
  mis_row <- as.numeric(mismatch[iso, ])
  names(cov_row) <- colnames(coverage)
  names(mis_row) <- colnames(mismatch)

  for (tp in levels(coldata$timepoint)) {
    samples_tp <- rownames(coldata)[coldata$timepoint == tp]
    samples_tp <- intersect(samples_tp, names(cov_row))
    if (length(samples_tp) < 2) next

    df <- data.frame(
      sample    = samples_tp,
      cov       = cov_row[samples_tp],
      mis       = mis_row[samples_tp],
      condition = coldata[samples_tp, "condition"]
    )
    df <- df[df$cov > 0, ]
    if (nrow(df) < 2 || length(unique(df$condition)) < 2) next

    fit <- tryCatch(
      glm(cbind(mis, cov - mis) ~ condition, data = df, family = binomial()),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    em <- tryCatch(
      emmeans(fit, ~condition, type = "response"),
      error = function(e) NULL
    )
    if (is.null(em)) next
    em_df <- as.data.frame(em)

    # em_df$prob is the per-condition estimated modification level (see
    # project discussion: this is a point-estimate stoichiometry, not a
    # probability distribution -- treated as such downstream in
    # compute_delta_c.py, which consumes it as f_stim/f_ctrl directly).
    # NOTE: condition labels come from sample_manifest.tsv, which uses
    # "control"/"polyIC" (confirmed against the real Stage 1 manifest), not
    # "stimulated" -- the previous "stimulated" label matched zero rows for
    # EVERY isodecoder/timepoint, so every iteration silently hit the `next`
    # below and the whole rule produced an effectively empty results file
    # rather than erroring loudly.
    stim_row <- em_df[em_df$condition == "polyIC", ]
    ctrl_row <- em_df[em_df$condition == "control", ]
    if (nrow(stim_row) == 0 || nrow(ctrl_row) == 0) next

    # No relevel() is applied anywhere in this pipeline, so glm()/factor()
    # use R's default alphabetical reference level: "control" < "polyIC",
    # so "control" is the reference and the non-reference dummy coefficient
    # is named "conditionpolyIC", not "conditionstimulated".
    pval <- tryCatch({
      s <- summary(fit)
      s$coefficients["conditionpolyIC", "Pr(>|z|)"]
    }, error = function(e) NA_real_)

    results[[length(results) + 1]] <- data.frame(
      isodecoder_id   = iso,
      timepoint       = tp,
      bucket          = bucket,
      confidence_tier = confidence_tier,
      f_ctrl          = ctrl_row$prob,
      f_ctrl_lower    = ctrl_row$asymp.LCL,
      f_ctrl_upper    = ctrl_row$asymp.UCL,
      f_stim          = stim_row$prob,
      f_stim_lower    = stim_row$asymp.LCL,
      f_stim_upper    = stim_row$asymp.UCL,
      pvalue          = pval,
      n_samples       = nrow(df)
    )
  }
}

final <- bind_rows(results)
if (nrow(final) > 0) {
  final$padj <- p.adjust(final$pvalue, method = "BH")
  final$significant <- !is.na(final$padj) & final$padj < fdr
}

write_tsv(final, out_path)
cat(sprintf("Wrote %d GLM result rows -> %s\n", nrow(final), out_path))

if (bucket == "Q34") {
  cat(
    "NOTE: Q34 results are EXPLORATORY. Standard-RT misincorporation for ",
    "queuosine is expected to sit close to background noise (see project ",
    "methods discussion) -- treat f_stim/f_ctrl here as descriptive, not ",
    "confirmatory, regardless of the nominal p-value/CI computed above. ",
    "Interpret only in combination with the TGT expression companion check.\n"
  )
}

sink(type = "message")
sink(type = "output")
