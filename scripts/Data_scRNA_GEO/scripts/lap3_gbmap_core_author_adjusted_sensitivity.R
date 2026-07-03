#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

setwd("/home/lzb/glioma")
data.table::setDTthreads(8)

input_file <- "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/tables/gbmap_core_patient_state_pseudobulk.csv"
out_file <- "Data_scRNA_GEO/GBmap_Core/results/LAP3_CellState_Strict/tables/gbmap_core_author_adjusted_within_state_pathway_associations.csv"

spearman_safe <- function(x, y, min_n = 6L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < min_n || length(unique(x)) < 3L || length(unique(y)) < 3L) {
    return(list(n = length(x), rho = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = length(x), rho = unname(test$estimate), p = test$p.value)
}

ps <- fread(input_file)
pathways <- c(
  "HALLMARK_MTORC1_SIGNALING",
  "LEUCINE_BCAA_CORE",
  "MTORC1_READOUT_CORE",
  "REACTOME_TRANSLATION"
)

rows <- list()
k <- 0L
for (variant in unique(ps$entry_variant)) {
  for (th in sort(unique(ps$threshold))) {
    for (state in sort(unique(ps$author_state))) {
      for (pathway in pathways) {
        d <- ps[
          entry_variant == variant &
            threshold == th &
            author_state == state
        ]
        family <- ifelse(
          pathway %in% c("HALLMARK_MTORC1_SIGNALING", "LEUCINE_BCAA_CORE"),
          "primary",
          "secondary"
        )
        raw <- spearman_safe(d$lap3_mean, d[[pathway]])
        keep <- is.finite(d$lap3_mean) & is.finite(d[[pathway]])
        d2 <- d[keep]

        if (nrow(d2) >= 8L && uniqueN(d2$author) >= 3L && uniqueN(d2$author_donor) >= 6L) {
          lap3_residual <- residuals(lm(lap3_mean ~ author, data = d2))
          pathway_residual <- residuals(lm(reformulate("author", response = pathway), data = d2))
          adjusted <- spearman_safe(lap3_residual, pathway_residual)
        } else {
          adjusted <- list(n = nrow(d2), rho = NA_real_, p = NA_real_)
        }

        k <- k + 1L
        rows[[k]] <- data.table(
          entry_variant = variant,
          threshold = th,
          author_state = state,
          pathway = pathway,
          fdr_family = family,
          n = raw$n,
          raw_rho = raw$rho,
          raw_p = raw$p,
          author_adjusted_rho = adjusted$rho,
          author_adjusted_p = adjusted$p,
          n_authors = uniqueN(d2$author),
          n_donors = uniqueN(d2$author_donor)
        )
      }
    }
  }
}

out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
out[, raw_p_adj := p.adjust(raw_p, method = "BH"), by = fdr_family]
out[, author_adjusted_p_adj := p.adjust(author_adjusted_p, method = "BH"), by = fdr_family]
fwrite(out, out_file)

cat("Wrote:", out_file, "\n")
cat("Rows:", nrow(out), "\n")
