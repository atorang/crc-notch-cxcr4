# 08_survival_km.R
#
# Inputs: data/processed/tcga/{TCGA_COAD_READ_AllData_complete.RData,TCGA.RNAseq.v3.RData}
#         data/processed/nunes/{Nunes_Swedish.RNAseq.RData,data.nunes.muts.RData}
# Output: results/figures/KaplanMeier_TCGA_Nunes.pdf

library(here)
library(yaml)

cfg <- read_yaml(here("config", "config.yaml"))
source(here("R", "mutation_utils.R"))
source(here("R", "bulk_data_utils.R"))
source(here("R", "gsva_utils.R"))
source(here("R", "survival_utils.R"))

out_dir <- here(cfg$paths$results_figures)
log_dir <- here(cfg$paths$results_logs)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pw_sander <- list(Notch_Sander = cfg$signature_scoring$gene_sets$Notch_Sander$genes)
cutoff <- cfg$gsva$quartile_cutoff_low  # 0.25; symmetric with quartile_cutoff_high (0.75)
km_colors <- c("darkblue", "darkred")

input_files <- c(
  here(cfg$paths$data_processed, "tcga", "TCGA_COAD_READ_AllData_complete.RData"),
  here(cfg$paths$data_processed, "tcga", "TCGA.RNAseq.v3.RData"),
  here(cfg$paths$data_processed, "nunes", "Nunes_Swedish.RNAseq.RData"),
  here(cfg$paths$data_processed, "nunes", "data.nunes.muts.RData")
)

# ---- TCGA ------------------------------------------------------------------
tcga <- load_tcga_survival_cohort(input_files[2], input_files[1])

notch_tcga <- run_gsva(tcga$exprs, pw_sander, method = "gsva")[1, ]
label_tcga <- quartile_label(notch_tcga, cutoff)
p1 <- plot_kaplan_meier(
  tcga$col_annot$recurrence, as.numeric(tcga$col_annot$daysToRecurrence) / 30,
  label_tcga, colors = km_colors, title = sprintf("TCGA, %%%d", cutoff * 100),
  time_range = NULL, censor_shape = "+"
)

wnt_wt_tcga <- which(tcga$muts_bin["APC", ] != 1 & tcga$muts_bin["CTNNB1", ] != 1)
notch_tcga_wt <- run_gsva(tcga$exprs[, wnt_wt_tcga], pw_sander, method = "gsva")[1, ]
label_tcga_wt <- quartile_label(notch_tcga_wt, cutoff)
p2 <- plot_kaplan_meier(
  tcga$col_annot$recurrence[wnt_wt_tcga], as.numeric(tcga$col_annot$daysToRecurrence[wnt_wt_tcga]) / 30,
  label_tcga_wt, colors = km_colors, title = sprintf("TCGA, WNT-wt, %%%d", cutoff * 100),
  time_range = NULL, censor_shape = "+"
)

trunc_tcga <- truncate_survival(tcga$col_annot$recurrence, tcga$col_annot$daysToRecurrence,
                                 cfg$survival$tcga_truncation_days)
p3 <- plot_kaplan_meier(
  trunc_tcga$event, trunc_tcga$time_to_event / 30, label_tcga,
  colors = km_colors, title = sprintf("TCGA, %%%d", cutoff * 100), time_range = NULL, censor_shape = "+"
)
trunc_tcga_wt <- truncate_survival(tcga$col_annot$recurrence[wnt_wt_tcga],
                                    tcga$col_annot$daysToRecurrence[wnt_wt_tcga],
                                    cfg$survival$tcga_truncation_days)
p4 <- plot_kaplan_meier(
  trunc_tcga_wt$event, trunc_tcga_wt$time_to_event / 30, label_tcga_wt,
  colors = km_colors, title = sprintf("TCGA, WNT-wt, %%%d", cutoff * 100), time_range = NULL, censor_shape = "+"
)

# ---- Nunes ------------------------------------------------------------------
nunes <- load_nunes_survival_cohort(input_files[3], input_files[4])

notch_nunes <- run_gsva(nunes$exprs, pw_sander, method = "gsva")[1, ]
label_nunes <- quartile_label(notch_nunes, cutoff)
p5 <- plot_kaplan_meier(
  nunes$col_annot$Recurrence, as.numeric(nunes$col_annot$Recurrence.free.survival.days) / 30,
  label_nunes, colors = km_colors, title = sprintf("Nunes, %%%d", cutoff * 100), time_range = NULL, censor_shape = "+"
)

wnt_wt_nunes <- which(nunes$muts_bin["APC", ] != 1 & nunes$muts_bin["CTNNB1", ] != 1)
notch_nunes_wt <- run_gsva(nunes$exprs[, wnt_wt_nunes], pw_sander, method = "gsva")[1, ]
label_nunes_wt <- quartile_label(notch_nunes_wt, cutoff)
p6 <- plot_kaplan_meier(
  nunes$col_annot$Recurrence[wnt_wt_nunes], as.numeric(nunes$col_annot$Recurrence.free.survival.days[wnt_wt_nunes]) / 30,
  label_nunes_wt, colors = km_colors, title = "Nunes, WNT-wt", time_range = NULL, censor_shape = "+"
)

trunc_nunes <- truncate_survival(nunes$col_annot$Recurrence, nunes$col_annot$Recurrence.free.survival.days,
                                  cfg$survival$nunes_truncation_days)
p7 <- plot_kaplan_meier(
  trunc_nunes$event, trunc_nunes$time_to_event / 30, label_nunes,
  colors = km_colors, title = sprintf("Nunes, %%%d", cutoff * 100), time_range = NULL, censor_shape = "+"
)
trunc_nunes_wt <- truncate_survival(nunes$col_annot$Recurrence[wnt_wt_nunes],
                                     nunes$col_annot$Recurrence.free.survival.days[wnt_wt_nunes],
                                     cfg$survival$nunes_truncation_days)
p8 <- plot_kaplan_meier(
  trunc_nunes_wt$event, trunc_nunes_wt$time_to_event / 30, label_nunes_wt,
  colors = km_colors, title = "Nunes, WNT-wt", time_range = NULL, censor_shape = "+"
)

# ---- assemble ---------------------------------------------------------------
library(patchwork)
out_pdf <- file.path(out_dir, "KaplanMeier_TCGA_Nunes.pdf")
pdf(out_pdf, width = 15, height = 8.5)
print((p1$plot | p2$plot | p5$plot | p6$plot) / (p3$plot | p4$plot | p7$plot | p8$plot))
dev.off()

cat("TCGA all: n =", length(label_tcga), " high/low:", paste(table(label_tcga), collapse = "/"), "\n")
cat("TCGA WNT-wt: n =", length(label_tcga_wt), " high/low:", paste(table(label_tcga_wt), collapse = "/"), "\n")
cat("TCGA 5yr-truncated: high/low:", paste(table(quartile_label(notch_tcga, cutoff)), collapse = "/"), "\n")
cat("Nunes all: n =", length(label_nunes), " high/low:", paste(table(label_nunes), collapse = "/"), "\n")
cat("Nunes WNT-wt: n =", length(label_nunes_wt), " high/low:", paste(table(label_nunes_wt), collapse = "/"), "\n")
cat("wrote", out_pdf, "\n")

source(here("R", "logging.R"))
log_run(
  script_name = "08_survival_km",
  input_files = input_files,
  params = list(
    cutoff = cutoff,
    tcga_truncation_days = cfg$survival$tcga_truncation_days,
    nunes_truncation_days = cfg$survival$nunes_truncation_days,
    tcga_wnt_wt_matrix = "binary_mutation_matrix_filtered (strict, not plain non-synonymous)"
  ),
  output_files = out_pdf,
  log_dir = log_dir
)
