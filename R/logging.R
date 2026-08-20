# Shared run-logging helper: every scripts/*.R call log_run() once at the end
# to write a timestamped provenance record (inputs + checksums, parameters,
# package versions, outputs) to results/logs/, per the release requirement
# that every script's provenance is traceable without re-running it.

log_run <- function(script_name, input_files, params, output_files, log_dir) {
  library(digest)
  ts <- format(Sys.time(), "%Y%m%dT%H%M%S")
  log_path <- file.path(log_dir, paste0(script_name, "_", ts, ".log"))

  checksum_block <- function(paths) {
    if (length(paths) == 0) return("  (none)")
    vapply(paths, function(p) {
      if (file.exists(p)) {
        sprintf("  %s  md5=%s  size=%dB", p, digest(file = p, algo = "md5"), file.size(p))
      } else {
        sprintf("  %s  [MISSING]", p)
      }
    }, character(1)) |> paste(collapse = "\n")
  }

  lines <- c(
    sprintf("script: %s", script_name),
    sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "input files:",
    checksum_block(input_files),
    "",
    "parameters:",
    paste0("  ", capture.output(str(params, no.list = TRUE))),
    "",
    "output files:",
    checksum_block(output_files),
    "",
    "session info:",
    paste0("  ", capture.output(sessionInfo()))
  )
  writeLines(lines, log_path)
  invisible(log_path)
}
