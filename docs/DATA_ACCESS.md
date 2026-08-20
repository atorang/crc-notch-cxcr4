# DATA_ACCESS.md

## TCGA-COAD / TCGA-READ 

`scripts/00_download_data.R` queries and downloads via `TCGAbiolinks::GDCquery()`
+ `GDCdownload()` (Transcriptome Profiling + Simple Nucleotide Variation,
Masked Somatic Mutation / Mutect2). 

## Nunes et al. cohort (E-MTAB-12862) 

`scripts/00_download_data.R` fetches raw counts/TPM and the VEP-annotated
mutation calls directly from ArrayExpress/BioStudies and verifies checksums
against `data/metadata/checksums.md5`.


## Pelka et al. scRNA-seq (SCP1162 / GSE178341)

The Broad Single Cell Portal (Study ID SCP1162) requires an authenticated,
manual download. Manual steps:

1. Log in to https://singlecell.broadinstitute.org/single_cell/study/SCP1162
2. Download the epithelial-compartment expression matrix and cell metadata
3. Place the file(s) under `data/raw/scrnaseq/`.

