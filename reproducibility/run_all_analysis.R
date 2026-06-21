#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) {
    return(default)
  }
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

has_flag <- function(flag) flag %in% args

mode <- get_arg("--mode=", "quick")
if (!mode %in% c("quick", "pilot", "full")) {
  stop("--mode must be one of: quick, pilot, full", call. = FALSE)
}
skip_package_check <- has_flag("--skip-package-check")
python_bin <- get_arg("--python=", Sys.getenv("PYTHON", "python3"))

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "analysis")) &&
        dir.exists(file.path(current, "package", "combodesign")) &&
        file.exists(file.path(current, "package", "combodesign", "DESCRIPTION"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find project root from ", start, call. = FALSE)
    }
    current <- parent
  }
}

root <- find_project_root()
setwd(root)

time_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
repro_dir <- "reproducibility"
log_dir <- file.path(repro_dir, "logs", paste0("run_all_", time_tag))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

records <- data.frame(
  step = character(),
  action = character(),
  mode = character(),
  required = logical(),
  status = character(),
  exit_code = integer(),
  seconds = numeric(),
  log_path = character(),
  note = character(),
  stringsAsFactors = FALSE
)

add_record <- function(step, action, required, status, exit_code = NA_integer_,
                       seconds = NA_real_, log_path = NA_character_,
                       note = NA_character_) {
  records <<- rbind(
    records,
    data.frame(
      step = step,
      action = action,
      mode = mode,
      required = required,
      status = status,
      exit_code = exit_code,
      seconds = seconds,
      log_path = log_path,
      note = note,
      stringsAsFactors = FALSE
    )
  )
}

sanitize_step <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

run_cmd <- function(step, command, command_args = character(),
                    required = TRUE, note = "") {
  log_path <- file.path(log_dir, paste0(sanitize_step(step), ".log"))
  header <- c(
    paste0("Step: ", step),
    paste0("Command: ", paste(c(command, command_args), collapse = " ")),
    paste0("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ""
  )
  writeLines(header, log_path)
  start <- proc.time()[["elapsed"]]
  status <- system2(
    command,
    args = command_args,
    stdout = log_path,
    stderr = log_path,
    wait = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - start
  exit_code <- if (identical(status, 0L)) 0L else as.integer(status)
  ok <- exit_code == 0L
  add_record(
    step = step,
    action = paste(c(command, command_args), collapse = " "),
    required = required,
    status = if (ok) "pass" else "fail",
    exit_code = exit_code,
    seconds = elapsed,
    log_path = log_path,
    note = note
  )
  if (!ok && required) {
    stop("Required step failed: ", step, ". See ", log_path, call. = FALSE)
  }
  invisible(ok)
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "pass")
}

check_file_exists <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing expected files: ", paste(missing, collapse = "; "), call. = FALSE)
  }
  TRUE
}

check_file_nonempty <- function(paths) {
  check_file_exists(paths)
  info <- file.info(paths)
  empty <- paths[is.na(info$size) | info$size <= 0]
  if (length(empty)) {
    stop("Empty expected files: ", paste(empty, collapse = "; "), call. = FALSE)
  }
  TRUE
}

validate_pass_table <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("pass" %in% names(tab)) {
    pass <- truthy(tab$pass)
  } else if ("passed" %in% names(tab)) {
    pass <- truthy(tab$passed)
  } else if ("status" %in% names(tab)) {
    pass <- tolower(trimws(tab$status)) == "pass"
  } else {
    stop(label, " has no pass/passed/status column: ", path, call. = FALSE)
  }
  if (!length(pass)) {
    stop(label, " has no rows to validate: ", path, call. = FALSE)
  }
  if (!all(pass)) {
    failed <- which(!pass)
    stop(label, " has failed rows: ", paste(head(failed, 10), collapse = ", "),
         call. = FALSE)
  }
  TRUE
}

validate_package_check_status <- function() {
  log_path <- file.path("combodesign.Rcheck", "00check.log")
  if (!file.exists(log_path)) {
    stop("Package check log is missing: ", log_path, call. = FALSE)
  }
  lines <- readLines(log_path, warn = FALSE)
  if (!any(grepl("^Status: OK$", lines))) {
    stop("Package check did not end with Status: OK", call. = FALSE)
  }
  TRUE
}

validate_outputs <- function() {
  validate_pass_table(
    "analysis/results/core_redesign_validation_checks_2026-06-19.csv",
    "Core method validation"
  )
  validate_pass_table(
    "simulation/results/participant_failfast_checks_2026-06-19.csv",
    "Fail-fast simulation validation"
  )
  validate_pass_table(
    "simulation/results/participant_ademp_main_2026-06-19_validation_checks.csv",
    "Main simulation validation"
  )
  validate_pass_table(
    "tables/simulation/simulation_display_qa_2026-06-19.csv",
    "Simulation display QA"
  )
  validate_pass_table(
    "tables/simulation/figure_3_combined_qa_2026-06-20.csv",
    "Combined simulation figure QA"
  )
  validate_pass_table(
    "case_study/pdx/pdx_case_study_validation_checks_2026-06-20.csv",
    "PDX case-study validation"
  )
  validate_pass_table(
    "tables/pdx/pdx_display_qa_2026-06-20.csv",
    "PDX display QA"
  )
  check_file_nonempty(c(
    "figures/simulation/figure_3_simulation_error_and_power_2026-06-20.png",
    "figures/simulation/figure_3_simulation_error_and_power_2026-06-20.pdf",
    "figures/pdx/figure_4_pdx_uncertainty_design_2026-06-20.png",
    "figures/pdx/figure_4_pdx_uncertainty_design_2026-06-20.pdf",
    "tables/simulation/table_3_simulation_main_2026-06-20.csv",
    "tables/pdx/table_4_pdx_design_planning_main_2026-06-20.csv"
  ))
  TRUE
}

run_cmd(
  "core_method_validation",
  file.path(R.home("bin"), "Rscript"),
  "analysis/validate_core_method.R",
  required = TRUE,
  note = "Reruns core method checks."
)

if (mode %in% c("pilot", "full")) {
  run_cmd(
    "failfast_simulation",
    file.path(R.home("bin"), "Rscript"),
    "simulation/run_failfast_participant_validation.R",
    required = TRUE,
    note = "Reruns fail-fast participant-level checks."
  )
  sim_mode <- if (mode == "pilot") "--pilot" else "--main"
  run_cmd(
    paste0("participant_simulation_", sub("^--", "", sim_mode)),
    file.path(R.home("bin"), "Rscript"),
    c("simulation/run_main_participant_ademp.R", sim_mode),
    required = TRUE,
    note = "Reruns the participant-level simulation."
  )
  run_cmd(
    "simulation_display_artifacts",
    file.path(R.home("bin"), "Rscript"),
    "simulation/create_simulation_displays.R",
    required = TRUE,
    note = "Rebuilds simulation tables and panel-level figure sources."
  )
  run_cmd(
    "combined_simulation_figure",
    file.path(R.home("bin"), "Rscript"),
    "simulation/create_combined_simulation_figure3.R",
    required = TRUE,
    note = "Rebuilds the combined simulation figure."
  )
}

if (mode == "full") {
  run_cmd(
    "pdx_case_study",
    file.path(R.home("bin"), "Rscript"),
    "case_study/pdx/run_pdx_case_study.R",
    required = TRUE,
    note = "Reruns the PDX case-study analysis."
  )
  run_cmd(
    "pdx_display_artifacts",
    python_bin,
    "case_study/pdx/create_pdx_displays.py",
    required = TRUE,
    note = "Rebuilds PDX tables and figures."
  )
}

if (!skip_package_check) {
  if (file.exists("combodesign_0.1.0.tar.gz")) {
    unlink("combodesign_0.1.0.tar.gz")
  }
  if (dir.exists("combodesign.Rcheck")) {
    unlink("combodesign.Rcheck", recursive = TRUE)
  }
  run_cmd(
    "package_build",
    file.path(R.home("bin"), "R"),
    c("CMD", "build", "package/combodesign"),
    required = TRUE,
    note = "Builds the source package."
  )
  run_cmd(
    "package_check",
    file.path(R.home("bin"), "R"),
    c("CMD", "check", "--no-manual", "combodesign_0.1.0.tar.gz"),
    required = TRUE,
    note = "Checks the built source package."
  )
  validate_package_check_status()
}

validate_outputs()
add_record(
  step = "output_integrity_validation",
  action = "validate_outputs()",
  required = TRUE,
  status = "pass",
  exit_code = 0L,
  seconds = 0,
  log_path = NA_character_,
  note = "Validated core, simulation, PDX, and display output files."
)

step_log_path <- file.path(repro_dir, paste0("run_all_step_log_", time_tag, ".csv"))
write.csv(records, step_log_path, row.names = FALSE)

report_path <- file.path(repro_dir, paste0("RUN_ALL_REPORT_", time_tag, ".md"))
report_lines <- c(
  "# Run-all analysis report",
  "",
  paste0("Mode: `", mode, "`"),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Step log CSV: `", step_log_path, "`"),
  paste0("Command log directory: `", log_dir, "`"),
  "",
  "## Step Summary",
  "",
  paste0("- ", records$step, ": ", records$status, " (", round(records$seconds, 2), " seconds)"),
  "",
  "## Overall Status",
  "",
  if (all(records$status == "pass")) "PASS" else "FAIL"
)
writeLines(report_lines, report_path)

message("Run-all analysis complete: ", report_path)
