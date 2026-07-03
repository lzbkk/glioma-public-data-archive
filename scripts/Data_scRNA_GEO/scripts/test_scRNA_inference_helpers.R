#!/usr/bin/env Rscript

setwd("/home/lzb/glioma")
source("Data_scRNA_GEO/scripts/helpers/scRNA_inference_helpers.R")

set.seed(20260629)
test_data <- do.call(rbind, lapply(seq_len(12), function(i) {
  data.frame(
    patient = paste0("P", i),
    state = c("AC", "MES", "OPC"),
    x = rnorm(3, mean = i / 5),
    y = rnorm(3, mean = i / 6)
  )
}))

statistic <- function(data) {
  residual_spearman(data, "x", "y", "patient", "state")
}

bootstrap_a <- cluster_bootstrap(
  test_data,
  cluster = "patient",
  statistic = statistic,
  replicates = 200L,
  seed = 42L
)
bootstrap_b <- cluster_bootstrap(
  test_data,
  cluster = "patient",
  statistic = statistic,
  replicates = 200L,
  seed = 42L
)

stopifnot(
  identical(bootstrap_a$bootstrap_values, bootstrap_b$bootstrap_values),
  bootstrap_a$valid_replicates >= 100L,
  is.finite(bootstrap_a$estimate),
  bootstrap_a$ci_low <= bootstrap_a$ci_high
)

lopo <- leave_one_cluster_out(test_data, "patient", statistic)
stopifnot(
  nrow(lopo) == 12L,
  setequal(lopo$omitted_cluster, unique(test_data$patient))
)

fdr_input <- data.frame(
  test = letters[1:4],
  p_value = c(0.01, 0.04, 0.01, 0.04),
  fdr_family = c("primary", "primary", "secondary", "secondary")
)
fdr_output <- adjust_fdr_by_family(fdr_input)
stopifnot(
  identical(round(fdr_output$p_adj_BH, 3), c(0.02, 0.04, 0.02, 0.04))
)

cat("All scRNA inference helper tests passed.\n")
