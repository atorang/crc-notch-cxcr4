# 01_preprocess_bulk.R


library(here)
library(yaml)
library(digest)

cfg <- read_yaml(here("config", "config.yaml"))

inputs <- list(
  list(path = here(cfg$paths$data_processed, "tcga", "TCGA_COAD_READ_AllData_complete.RData"),
       note = "Built by TCGA.PanCancer/download_tcga_crc_alldata.R + build_mutation_matrices.R (TCGAbiolinks)."),
  list(path = here(cfg$paths$data_processed, "tcga", "TCGA.RNAseq.v3.RData"),
       note = "Built by TCGA.PanCancer/Data.Handeling.TCGA.RNAseq.v1.v2.v3.R (TCGAbiolinks)."),
  list(path = here(cfg$paths$data_processed, "nunes", "Nunes_Swedish.RNAseq.RData"),
       note = "Built by Nunes_Tobias.Sjoblom/Data.Handelling.Nunes_Swedish.RNAseq.R from ArrayExpress E-MTAB-12862."),
  list(path = here(cfg$paths$data_processed, "nunes", "data.nunes.muts.RData"),
       note = "PINNED, build script not recovered — see docs/PARAMETERS.md and OPEN_QUESTIONS.md #5.")
)

checksum_dir <- here(cfg$paths$data_metadata)
dir.create(checksum_dir, recursive = TRUE, showWarnings = FALSE)
checksum_path <- file.path(checksum_dir, "checksums.md5")

lines <- character(0)
for (item in inputs) {
  if (!file.exists(item$path)) {
    stop(sprintf("Missing required processed input: %s\n  %s", item$path, item$note))
  }
  md5 <- digest(file = item$path, algo = "md5")
  cat(sprintf("%s  %s\n  %s\n", md5, item$path, item$note))
  lines <- c(lines, sprintf("%s  %s", md5, item$path))
}
writeLines(lines, checksum_path)
cat("\nwrote", checksum_path, "\n")
