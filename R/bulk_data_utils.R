# Cohort loaders for the survival and mutation-boxplot scripts.

#' Load the TCGA-COAD/READ cohort:
#' stage I-IV primary tumors, sample IDs aligned between expression and the
#' *strict* consensus-damaging mutation matrix, and duplicate gene
#' symbols resolved by keeping the FIRST row for each symbol
load_tcga_survival_cohort <- function(rnaseq_path, alldata_path) {
  e <- new.env()
  load(alldata_path, envir = e)
  muts_bin_full <- e$data.list[["Simple Nucleotide Variation"]]$binary_mutation_matrix_filtered
  colnames(muts_bin_full) <- gsub("-", ".", substr(colnames(muts_bin_full), 1, 16))

  e2 <- new.env()
  load(rnaseq_path, envir = e2)
  d <- e2$data

  d$col.annot$stage_update <- gsub("A|B|C", "", d$col.annot$ajcc_pathologic_stage)
  stage_keep <- d$col.annot$stage_update %in% paste0("Stage ", c("I", "II", "III", "IV"))
  d$col.annot <- d$col.annot[stage_keep, , drop = FALSE]
  for (nm in names(d$data)) d$data[[nm]] <- d$data[[nm]][, stage_keep, drop = FALSE]

  exprs <- d$data$quantile
  colnames(exprs) <- substr(colnames(exprs), 1, 16)
  common <- intersect(colnames(exprs), colnames(muts_bin_full))

  sample_keep <- d$col.annot$sample %in% gsub("\\.", "-", common)
  d$col.annot <- d$col.annot[sample_keep, , drop = FALSE]
  exprs <- d$data$quantile[, sample_keep, drop = FALSE]

  gene <- d$row.annot$gene
  keep_rows <- !is.na(gene) & !duplicated(gene)
  exprs <- exprs[keep_rows, , drop = FALSE]
  rownames(exprs) <- gene[keep_rows]
  colnames(exprs) <- substr(colnames(exprs), 1, 16)

  muts_bin <- muts_bin_full[, colnames(exprs), drop = FALSE]

  list(exprs = exprs, col_annot = d$col.annot, muts_bin = muts_bin)
}

#' Load the Nunes cohort: 
load_nunes_survival_cohort <- function(rnaseq_path, muts_path) {
  e <- new.env()
  load(rnaseq_path, envir = e)
  rna <- e$data

  keep <- rna$col.annot$Tissue.Type == "neoplasm"
  rna$col.annot <- rna$col.annot[keep, , drop = FALSE]
  for (nm in names(rna$data)) rna$data[[nm]] <- rna$data[[nm]][, keep, drop = FALSE]

  rna$col.annot$Recurrence[rna$col.annot$Recurrence == "Not_Applicable"] <- NA
  rna$col.annot$Recurrence[rna$col.annot$Recurrence == "No"] <- 0
  rna$col.annot$Recurrence[rna$col.annot$Recurrence == "Yes"] <- 1

  exprs <- rna$data$quantile
  rownames(exprs) <- rna$row.annot$gene
  stopifnot(!anyDuplicated(rownames(exprs)))

  e2 <- new.env()
  load(muts_path, envir = e2)
  muts <- e2$data
  sample_ids <- gsub("-", ".", colnames(muts$data$events))
  colnames(muts$data$events) <- sample_ids
  events <- muts$data$events[, colnames(exprs), drop = FALSE]

  muts_bin <- collapse_position_to_gene(events, muts$row.annot$gene)

  list(exprs = exprs, col_annot = rna$col.annot, muts_bin = muts_bin)
}
