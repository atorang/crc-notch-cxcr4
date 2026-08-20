# 09_mutation_expression_boxplots.R
#
# Input: data/processed/nunes/{Nunes_Swedish.RNAseq.RData,data.nunes.muts.RData}
# Output: results/figures/HES1_HES5_HES6_by_FBXW7_WNT_Nunes.pdf

library(here)
library(yaml)

cfg <- read_yaml(here("config", "config.yaml"))
source(here("R", "mutation_utils.R"))
source(here("R", "plotting.R"))
source(here("R", "logging.R"))

rnaseq_path <- here(cfg$paths$data_processed, "nunes", "Nunes_Swedish.RNAseq.RData")
muts_path   <- here(cfg$paths$data_processed, "nunes", "data.nunes.muts.RData")
out_dir     <- here(cfg$paths$results_figures)
log_dir     <- here(cfg$paths$results_logs)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load bulk expression, subset to MSS primary tumors --------------------
load(rnaseq_path)         
rna <- data; rm(data)

keep_cols <- rna$col.annot$Tissue.Type == "neoplasm" & rna$col.annot$MSI.Status == "MSS"
exprs <- rna$data$quantile[, keep_cols, drop = FALSE]
rownames(exprs) <- rna$row.annot$gene
stopifnot("gene symbols must be unique in this matrix (see docs/PARAMETERS.md)" = !anyDuplicated(rownames(exprs)))

# ---- load mutation calls, align sample IDs, collapse to gene level --------
load(muts_path)            
muts <- data; rm(data)

assert_nonsynonymous_only(muts$row.annot$consequence, cfg$mutation_calling$non_synonymous_classes_vep)

sample_ids <- gsub("-", ".", colnames(muts$data$events))
colnames(muts$data$events) <- sample_ids
events <- muts$data$events[, colnames(exprs), drop = FALSE]

muts_bin <- collapse_position_to_gene(events, muts$row.annot$gene)
group <- classify_fbxw7_wnt_groups(muts_bin)

# ---- boxplots --------------------------------------------------------------
group_colors <- c("darkblue", "darkred", "#8080C5", "#C58080")
comparisons <- list(
  c("FBXW7-wt_WNT-mut", "FBXW7-mut_WNT-mut"),
  c("FBXW7-wt_WNT-wt", "FBXW7-mut_WNT-wt")
)

out_pdf <- file.path(out_dir, "HES1_HES5_by_FBXW7_WNT_Nunes.pdf")
pdf(out_pdf, width = 12, height = 7)
library(patchwork)
for (add_points in c(TRUE, FALSE)) {
  panels <- lapply(c("HES1", "HES5"), function(gene) {
    boxplot_with_pairwise_pvalues(
      value = exprs[gene, ], group = group, comparisons = comparisons,
      colors = group_colors, test = cfg$boxplot_stats$fbxw7_hes_test,
      title = "Mutation in FBXW7", ylab = paste(gene, "expression"),
      add_points = add_points
    )
  })
  print(panels[[1]] | panels[[2]])
}
dev.off()

cat(sprintf("n samples per group:\n")); print(table(group))
cat(sprintf("\nwrote %s\n", out_pdf))

log_run(
  script_name = "09_mutation_expression_boxplots",
  input_files = c(rnaseq_path, muts_path),
  params = list(
    cohort_filter = "Tissue.Type == neoplasm & MSI.Status == MSS",
    wnt_wt_definition = cfg$mutation_calling$wnt_wt_definition,
    test = cfg$boxplot_stats$fbxw7_hes_test
  ),
  output_files = out_pdf,
  log_dir = log_dir
)
