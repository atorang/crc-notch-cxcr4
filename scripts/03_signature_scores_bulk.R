# 03_signature_scores_bulk.R

library(here)
library(yaml)

cfg <- read_yaml(here("config", "config.yaml"))
source(here("R", "bulk_data_utils.R"))
source(here("R", "mutation_utils.R"))
source(here("R", "gsva_utils.R"))

table_dir <- here(cfg$paths$results_tables)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

pw_sander <- list(Notch_Sander = cfg$signature_scoring$gene_sets$Notch_Sander$genes)

tcga <- load_tcga_survival_cohort(
  here(cfg$paths$data_processed, "tcga", "TCGA.RNAseq.v3.RData"),
  here(cfg$paths$data_processed, "tcga", "TCGA_COAD_READ_AllData_complete.RData")
)
notch_tcga <- run_gsva(tcga$exprs, pw_sander, method = "gsva")[1, ]
write.table(data.frame(sample = names(notch_tcga), Notch_Sander_GSVA = notch_tcga),
            file.path(table_dir, "notch_gsva_tcga.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

nunes <- load_nunes_survival_cohort(
  here(cfg$paths$data_processed, "nunes", "Nunes_Swedish.RNAseq.RData"),
  here(cfg$paths$data_processed, "nunes", "data.nunes.muts.RData")
)
notch_nunes <- run_gsva(nunes$exprs, pw_sander, method = "gsva")[1, ]
write.table(data.frame(sample = names(notch_nunes), Notch_Sander_GSVA = notch_nunes),
            file.path(table_dir, "notch_gsva_nunes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("wrote", file.path(table_dir, "notch_gsva_tcga.tsv"), "\n")
cat("wrote", file.path(table_dir, "notch_gsva_nunes.tsv"), "\n")
