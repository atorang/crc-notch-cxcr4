# Mutation-matrix helpers shared by the survival (08) and FBXW7/WNT boxplot
# (09) scripts. 

assert_nonsynonymous_only <- function(consequence, allowed_classes) {
  observed <- unique(consequence)
  unexpected <- setdiff(observed, allowed_classes)
  if (length(unexpected) > 0) {
    stop(sprintf(
      "Unexpected variant consequence classes: %s",
      paste(unexpected, collapse = ", ")
    ))
  }
  invisible(TRUE)
}


collapse_position_to_gene <- function(events, gene_labels) {
  gene_levels <- sort(unique(gene_labels))
  values <- vapply(seq_len(ncol(events)), function(j) {
    tapply(events[, j], gene_labels, function(x) as.integer(any(x == 1, na.rm = TRUE)))[gene_levels]
  }, FUN.VALUE = integer(length(gene_levels)))
  matrix(values, nrow = length(gene_levels), ncol = ncol(events),
         dimnames = list(gene_levels, colnames(events)))
}

#' Classify samples into the 4-level FBXW7 x WNT mutation groupt: 
#' WNT-wt requires no non-synonymous mutation in
#' EITHER APC or CTNNB1 (config: mutation_calling.wnt_wt_definition).
classify_fbxw7_wnt_groups <- function(muts_bin) {
  wnt <- ifelse(muts_bin["APC", ] != 1 & muts_bin["CTNNB1", ] != 1, "WNT-wt", "WNT-mut")
  fbxw7 <- ifelse(muts_bin["FBXW7", ] == 1, "mut", "wt")
  result <- factor(
    paste(paste0("FBXW7-", fbxw7), wnt, sep = "_"),
    levels = c("FBXW7-wt_WNT-mut", "FBXW7-mut_WNT-mut", "FBXW7-wt_WNT-wt", "FBXW7-mut_WNT-wt")
  )
  names(result) <- colnames(muts_bin)
  result
}
