# 00_download_data.R
#
# Downloads and verifies the bulk data this pipeline needs. Fails loudly
# (stop()).


library(here)
library(yaml)

cfg <- read_yaml(here("config", "config.yaml"))
raw_dir <- here(cfg$paths$data_raw)
dir.create(file.path(raw_dir, "tcga"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(raw_dir, "nunes"), recursive = TRUE, showWarnings = FALSE)

download_tcga <- function() {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks is required (config.yaml: software_versions.tcgabiolinks = ",
         cfg$software_versions$tcgabiolinks, "). install.packages/BiocManager first.")
  }
  library(TCGAbiolinks)
  for (project in c(cfg$data_sources$tcga$project_coad, cfg$data_sources$tcga$project_read)) {
    cat("Querying", project, "Transcriptome Profiling...\n")
    query_expr <- GDCquery(project = project, data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification",
                            workflow.type = cfg$data_sources$tcga$workflow_type_expression)
    GDCdownload(query_expr, directory = file.path(raw_dir, "tcga"))

    cat("Querying", project, "Simple Nucleotide Variation...\n")
    query_maf <- GDCquery(project = project, data.category = "Simple Nucleotide Variation",
                           data.type = "Masked Somatic Mutation",
                           workflow.type = cfg$data_sources$tcga$maf_pipeline)
    GDCdownload(query_maf, directory = file.path(raw_dir, "tcga"))
  }
  cat("TCGA download complete. Run 01_preprocess_bulk.R next ")
}

download_nunes <- function() {
  accession <- cfg$data_sources$nunes$arrayexpress_accession
  base_url <- sprintf("https://www.ebi.ac.uk/biostudies/files/%s", accession)
  files <- c(
    "CRC.SW.mRNA.symbol.count.txt.gz",
    "CRC.SW.mRNA.symbol.TPM.txt.gz",
    "vcfe_nunes1063_vep.txt.gz"
  )
  for (f in files) {
    dest <- file.path(raw_dir, "nunes", f)
    if (file.exists(dest)) { cat("already have", f, "\n"); next }
    url <- paste(base_url, f, sep = "/")
    cat("Downloading", url, "...\n")
    ok <- tryCatch({ download.file(url, dest, mode = "wb"); TRUE },
                    error = function(e) { message(e); FALSE })
    if (!ok) stop(sprintf("Failed to download %s from ArrayExpress. Check the accession/file layout at https://www.ebi.ac.uk/biostudies/arrayexpress/studies/%s", f, accession))
  }
  cat("Nunes download complete.\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || "tcga" %in% args) download_tcga()
if (length(args) == 0 || "nunes" %in% args) download_nunes()
