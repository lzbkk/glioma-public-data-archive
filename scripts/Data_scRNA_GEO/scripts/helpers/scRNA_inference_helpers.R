# Shared patient-level inference helpers for scRNA analyses.

cluster_bootstrap <- function(
    data,
    cluster,
    statistic,
    replicates = 2000L,
    conf = 0.95,
    seed = 20260629L) {
  stopifnot(
    is.data.frame(data),
    cluster %in% names(data),
    is.function(statistic),
    replicates > 0L,
    conf > 0,
    conf < 1
  )

  cluster_values <- unique(as.character(data[[cluster]]))
  cluster_values <- cluster_values[!is.na(cluster_values)]
  if (length(cluster_values) < 2L) {
    stop("cluster_bootstrap requires at least two non-missing clusters")
  }

  observed <- statistic(data)
  if (length(observed) != 1L || !is.numeric(observed)) {
    stop("statistic must return one numeric value")
  }

  set.seed(seed)
  bootstrap_values <- replicate(replicates, {
    sampled <- sample(cluster_values, length(cluster_values), replace = TRUE)
    bootstrap_data <- do.call(rbind, lapply(seq_along(sampled), function(i) {
      block <- data[as.character(data[[cluster]]) == sampled[[i]], , drop = FALSE]
      block$.source_cluster <- sampled[[i]]
      block[[cluster]] <- paste0(sampled[[i]], "__boot_", i)
      block
    }))
    suppressWarnings(statistic(bootstrap_data))
  })

  alpha <- (1 - conf) / 2
  finite_values <- bootstrap_values[is.finite(bootstrap_values)]
  ci <- if (length(finite_values) >= max(20L, ceiling(0.5 * replicates))) {
    stats::quantile(
      finite_values,
      probs = c(alpha, 1 - alpha),
      names = FALSE,
      na.rm = TRUE
    )
  } else {
    c(NA_real_, NA_real_)
  }

  list(
    estimate = observed,
    ci_low = ci[[1]],
    ci_high = ci[[2]],
    replicates = replicates,
    valid_replicates = length(finite_values),
    bootstrap_values = bootstrap_values
  )
}

leave_one_cluster_out <- function(data, cluster, statistic) {
  stopifnot(
    is.data.frame(data),
    cluster %in% names(data),
    is.function(statistic)
  )

  cluster_values <- unique(as.character(data[[cluster]]))
  cluster_values <- cluster_values[!is.na(cluster_values)]
  estimates <- vapply(cluster_values, function(cluster_value) {
    subset_data <- data[
      as.character(data[[cluster]]) != cluster_value,
      ,
      drop = FALSE
    ]
    suppressWarnings(statistic(subset_data))
  }, numeric(1))

  data.frame(
    omitted_cluster = cluster_values,
    estimate = estimates,
    stringsAsFactors = FALSE
  )
}

residual_spearman <- function(data, x, y, patient, state = NULL) {
  required <- c(x, y, patient, state)
  required <- required[!is.na(required) & nzchar(required)]
  stopifnot(all(required %in% names(data)))

  keep <- stats::complete.cases(data[, required, drop = FALSE])
  analysis_data <- data[keep, , drop = FALSE]
  if (
    nrow(analysis_data) < 6L ||
    length(unique(analysis_data[[patient]])) < 3L ||
    length(unique(analysis_data[[x]])) < 3L ||
    length(unique(analysis_data[[y]])) < 3L
  ) {
    return(NA_real_)
  }

  adjustment <- patient
  if (!is.null(state)) {
    adjustment <- c(adjustment, state)
  }
  fit_x <- stats::lm(stats::reformulate(adjustment, response = x), data = analysis_data)
  fit_y <- stats::lm(stats::reformulate(adjustment, response = y), data = analysis_data)
  suppressWarnings(stats::cor(
    stats::residuals(fit_x),
    stats::residuals(fit_y),
    method = "spearman"
  ))
}

adjust_fdr_by_family <- function(data, p_column = "p_value", family_column = "fdr_family") {
  stopifnot(
    is.data.frame(data),
    p_column %in% names(data),
    family_column %in% names(data)
  )

  output <- data
  output$p_adj_BH <- NA_real_
  family_values <- unique(as.character(output[[family_column]]))
  for (family_value in family_values) {
    index <- which(as.character(output[[family_column]]) == family_value)
    output$p_adj_BH[index] <- stats::p.adjust(output[[p_column]][index], method = "BH")
  }
  output
}
