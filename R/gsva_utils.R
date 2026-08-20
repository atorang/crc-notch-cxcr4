run_gsva <- function(exprs, gene_sets, method = "gsva") {
  library(GSVA)
  exprs <- as.matrix(exprs)
  params <- switch(method,
    gsva     = gsvaParam(exprData = exprs, geneSets = gene_sets),
    ssgsea   = ssgseaParam(exprData = exprs, geneSets = gene_sets),
    zscore   = zscoreParam(exprData = exprs, geneSets = gene_sets),
    plage    = plageParam(exprData = exprs, geneSets = gene_sets),
    stop("method must be one of: gsva, ssgsea, zscore, plage")
  )
  gsva(params)
}

#' NOTCH-high/low quartile assignment shared by 08_survival_km.R: high =
#' above the (1-cutoff) quantile, low = below the cutoff quantile, everything
#' in between is NA (excluded from binary survival comparisons).
quartile_label <- function(score, cutoff, high_label = "NOTCH-high", low_label = "NOTCH-low") {
  score <- as.numeric(score)
  hi <- quantile(score, 1 - cutoff, na.rm = TRUE)
  lo <- quantile(score, cutoff, na.rm = TRUE)
  factor(
    ifelse(score > hi, high_label, ifelse(score < lo, low_label, NA)),
    levels = c(low_label, high_label)
  )
}
