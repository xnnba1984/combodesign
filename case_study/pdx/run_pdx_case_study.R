#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

cmd_args <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)][1] %||% NA_character_)
root_dir <- if (is.na(file_arg)) {
  normalizePath(".", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(normalizePath(file_arg, mustWork = TRUE)), "..", ".."), mustWork = TRUE)
}
setwd(root_dir)

source("analysis/contribution_sample_size_allocation.R")

set.seed(20260620)

out_dir <- file.path("case_study", "pdx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_data <- file.path("data", "combination.csv")
stopifnot(file.exists(source_data))

alpha <- 0.025
target_power <- 0.80
min_allocation <- 0.10
n_boot <- 2000L
arms <- c("A", "B", "AB")
covariance_model <- "independent_arms"

candidate_plan <- tibble::tribble(
  ~case_role, ~component_a, ~component_b, ~combination,
  "primary_continuity_case", "BYL719", "binimetinib", "BYL719 + binimetinib",
  "strong_comparison_case", "LEE011", "encorafenib", "LEE011 + encorafenib",
  "strong_comparison_case", "BYL719", "LEE011", "BYL719 + LEE011",
  "cautionary_sensitivity_case", "BKM120", "binimetinib", "BKM120 + binimetinib"
)

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, prob, names = FALSE, type = 8))
}

pooled_sd_common <- function(...) {
  vectors <- list(...)
  n <- vapply(vectors, function(x) sum(is.finite(x)), integer(1))
  vars <- vapply(vectors, function(x) stats::var(x[is.finite(x)]), numeric(1))
  if (any(n < 2) || any(!is.finite(vars))) {
    return(NA_real_)
  }
  denom <- sum(n - 1)
  if (denom <= 0) {
    return(NA_real_)
  }
  out <- sqrt(sum((n - 1) * vars) / denom)
  ifelse(is.finite(out) && out > 0, out, NA_real_)
}

pooled_sd_pair <- function(x, y) {
  pooled_sd_common(x, y)
}

extract_triplet <- function(model_means, component_a, component_b, combination) {
  needed <- c(component_a, component_b, combination, "untreated")
  wide <- model_means %>%
    filter(Treatment %in% needed) %>%
    select(Model, Treatment, favorable_response) %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = favorable_response)

  required <- c(component_a, component_b, combination)
  missing_required <- setdiff(required, names(wide))
  if (length(missing_required)) {
    stop("Missing treatment columns for ", combination, ": ",
         paste(missing_required, collapse = ", "), call. = FALSE)
  }

  wide %>%
    filter(
      is.finite(.data[[component_a]]),
      is.finite(.data[[component_b]]),
      is.finite(.data[[combination]])
    )
}

summarise_triplet <- function(triplet, component_a, component_b, combination) {
  a <- triplet[[component_a]]
  b <- triplet[[component_b]]
  ab <- triplet[[combination]]
  sd_common <- pooled_sd_common(a, b, ab)
  sd_ab_a_pair <- pooled_sd_pair(ab, a)
  sd_ab_b_pair <- pooled_sd_pair(ab, b)

  additive_excess <- if ("untreated" %in% names(triplet) &&
                         sum(is.finite(triplet[["untreated"]])) >= 3) {
    mean(ab - a - b + triplet[["untreated"]], na.rm = TRUE)
  } else {
    NA_real_
  }

  tibble(
    n_triplet = nrow(triplet),
    mean_favorable_a = mean(a, na.rm = TRUE),
    mean_favorable_b = mean(b, na.rm = TRUE),
    mean_favorable_ab = mean(ab, na.rm = TRUE),
    common_pooled_sd = sd_common,
    delta_ab_minus_a = (mean(ab, na.rm = TRUE) - mean(a, na.rm = TRUE)) / sd_common,
    delta_ab_minus_b = (mean(ab, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sd_common,
    pairwise_delta_ab_minus_a = (mean(ab, na.rm = TRUE) - mean(a, na.rm = TRUE)) / sd_ab_a_pair,
    pairwise_delta_ab_minus_b = (mean(ab, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sd_ab_b_pair,
    additive_excess = additive_excess
  )
}

bootstrap_triplet <- function(triplet, component_a, component_b, combination,
                              n_boot = 2000L) {
  replicate_mat <- replicate(n_boot, {
    idx <- sample.int(nrow(triplet), size = nrow(triplet), replace = TRUE)
    samp <- triplet[idx, , drop = FALSE]
    a <- samp[[component_a]]
    b <- samp[[component_b]]
    ab <- samp[[combination]]
    sd_common <- pooled_sd_common(a, b, ab)
    delta_a <- (mean(ab, na.rm = TRUE) - mean(a, na.rm = TRUE)) / sd_common
    delta_b <- (mean(ab, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sd_common
    c(
      delta_ab_minus_a = delta_a,
      delta_ab_minus_b = delta_b,
      min_contribution = min(delta_a, delta_b)
    )
  })
  as_tibble(t(replicate_mat)) %>%
    mutate(
      valid_effect_input = is.finite(delta_ab_minus_a) &
        is.finite(delta_ab_minus_b) &
        is.finite(min_contribution)
    )
}

make_input_scenarios <- function(point_row, boot_df) {
  valid_boot <- boot_df %>% filter(valid_effect_input)

  point <- tibble(
    input_scenario = "plug_in",
    delta_ab_minus_a = point_row$delta_ab_minus_a,
    delta_ab_minus_b = point_row$delta_ab_minus_b
  )

  lower_bound <- tibble(
    input_scenario = "bootstrap_lower_95_effects",
    delta_ab_minus_a = safe_quantile(valid_boot$delta_ab_minus_a, 0.025),
    delta_ab_minus_b = safe_quantile(valid_boot$delta_ab_minus_b, 0.025)
  )

  bind_rows(point, lower_bound) %>%
    mutate(
      alpha = alpha,
      target_joint_power = target_power,
      clinical_covariance_model = covariance_model,
      valid_design_input = is.finite(delta_ab_minus_a) &
        is.finite(delta_ab_minus_b) &
        delta_ab_minus_a > 0 &
        delta_ab_minus_b > 0
    )
}

fit_design_methods <- function(delta, candidate_label, input_scenario) {
  equal <- sample_size_joint_success(
    delta = delta,
    allocation = equal_allocation(arms),
    rho = NULL,
    sigma2 = 1,
    target_power = target_power,
    alpha = alpha,
    arms = arms,
    covariance_model = covariance_model,
    N_upper = 1e6
  )
  joint <- optimize_allocation_for_sample_size(
    delta = delta,
    rho = NULL,
    sigma2 = 1,
    target_power = target_power,
    alpha = alpha,
    arms = arms,
    covariance_model = covariance_model,
    min_allocation = min_allocation,
    N_upper = 1e6,
    control = list(maxit = 140, reltol = 1e-7)
  )

  bind_rows(
    allocation_summary_row("equal allocation", equal),
    allocation_summary_row("joint-decision optimized allocation", joint)
  ) %>%
    mutate(
      candidate = candidate_label,
      input_scenario = input_scenario,
      alpha = alpha,
      clinical_covariance_model = covariance_model,
      .before = method
    )
}

design_from_scenario <- function(input_row, candidate_label) {
  if (!isTRUE(input_row$valid_design_input)) {
    return(tibble(
      candidate = candidate_label,
      input_scenario = input_row$input_scenario,
      method = "no finite design under selected lower-bound effects",
      alpha = alpha,
      clinical_covariance_model = covariance_model,
      target_joint_success_probability = target_power,
      N_continuous = NA_real_,
      N = NA_integer_,
      achieved_joint_success_probability = NA_real_,
      allocation_A = NA_real_,
      allocation_B = NA_real_,
      allocation_AB = NA_real_,
      marginal_power_AB_minus_A = NA_real_,
      marginal_power_AB_minus_B = NA_real_,
      reason = "At least one selected contribution-effect input is not positive."
    ))
  }

  delta <- c(
    AB_minus_A = input_row$delta_ab_minus_a,
    AB_minus_B = input_row$delta_ab_minus_b
  )

  out <- fit_design_methods(delta, candidate_label, input_row$input_scenario)
  out$reason <- NA_character_
  out
}

evaluate_design_under_input <- function(design_row, input_row) {
  if (!is.finite(design_row$N) || !isTRUE(input_row$valid_design_input)) {
    return(NA_real_)
  }
  allocation <- c(
    A = design_row$allocation_A,
    B = design_row$allocation_B,
    AB = design_row$allocation_AB
  )
  delta <- c(
    AB_minus_A = input_row$delta_ab_minus_a,
    AB_minus_B = input_row$delta_ab_minus_b
  )
  joint_success_probability(
    N = design_row$N,
    allocation = allocation,
    delta = delta,
    rho = NULL,
    sigma2 = 1,
    alpha = alpha,
    arms = arms,
    covariance_model = covariance_model
  )
}

variance_sensitivity_design <- function(input_row, clinical_sd_multiplier) {
  if (!isTRUE(input_row$valid_design_input)) {
    return(tibble(
      candidate = input_row$candidate,
      input_scenario = input_row$input_scenario,
      clinical_sd_multiplier = clinical_sd_multiplier,
      optimized_total_N = NA_integer_,
      allocation_A = NA_real_,
      allocation_B = NA_real_,
      allocation_AB = NA_real_,
      achieved_joint_power = NA_real_,
      reason = "At least one selected contribution-effect input is not positive."
    ))
  }
  delta <- c(
    AB_minus_A = input_row$delta_ab_minus_a,
    AB_minus_B = input_row$delta_ab_minus_b
  )
  sigma2 <- rep(clinical_sd_multiplier^2, 3)
  names(sigma2) <- arms
  design <- optimize_allocation_for_sample_size(
    delta = delta,
    rho = NULL,
    sigma2 = sigma2,
    target_power = target_power,
    alpha = alpha,
    arms = arms,
    covariance_model = covariance_model,
    min_allocation = min_allocation,
    N_upper = 1e6,
    control = list(maxit = 140, reltol = 1e-7)
  )
  tibble(
    candidate = input_row$candidate,
    input_scenario = input_row$input_scenario,
    clinical_sd_multiplier = clinical_sd_multiplier,
    optimized_total_N = design$N,
    allocation_A = design$allocation[["A"]],
    allocation_B = design$allocation[["B"]],
    allocation_AB = design$allocation[["AB"]],
    achieved_joint_power = design$achieved_power,
    reason = NA_character_
  )
}

round_numeric_df <- function(x, digits = 3) {
  is_num <- vapply(x, is.numeric, logical(1))
  x[is_num] <- lapply(x[is_num], function(v) round(v, digits))
  x
}

input_df <- readr::read_csv(source_data, show_col_types = FALSE) %>%
  mutate(
    Model = as.character(Model),
    Treatment = stringr::str_squish(gsub('"', "", Treatment)),
    favorable_response = -BestAvgResponse
  )

required_cols <- c("Model", "Treatment", "Treatment type", "BestAvgResponse", "favorable_response")
missing_cols <- setdiff(required_cols, names(input_df))
if (length(missing_cols)) {
  stop("Missing required columns in combination.csv: ", paste(missing_cols, collapse = ", "))
}

model_means <- input_df %>%
  group_by(Model, Treatment) %>%
  summarise(
    favorable_response = mean(favorable_response, na.rm = TRUE),
    .groups = "drop"
  )

candidate_results <- pmap(candidate_plan, function(case_role, component_a, component_b, combination) {
  triplet <- extract_triplet(model_means, component_a, component_b, combination)
  if (nrow(triplet) < 20) {
    stop("Candidate has fewer than 20 complete triplets: ", combination, call. = FALSE)
  }
  point <- summarise_triplet(triplet, component_a, component_b, combination) %>%
    mutate(
      case_role = case_role,
      component_a = component_a,
      component_b = component_b,
      combination = combination,
      candidate = combination,
      .before = n_triplet
    )
  boot <- bootstrap_triplet(triplet, component_a, component_b, combination, n_boot = n_boot)
  scenarios <- make_input_scenarios(point, boot) %>%
    mutate(
      case_role = case_role,
      component_a = component_a,
      component_b = component_b,
      combination = combination,
      candidate = combination,
      .before = input_scenario
    )
  designs <- pmap_dfr(
    list(split(scenarios, seq_len(nrow(scenarios))), scenarios$candidate),
    function(input_row, candidate) design_from_scenario(input_row, candidate)
  ) %>%
    left_join(
      scenarios %>%
        select(candidate, input_scenario, case_role, component_a, component_b,
               combination, delta_ab_minus_a, delta_ab_minus_b,
               valid_design_input),
      by = c("candidate", "input_scenario")
    )

  list(point = point, boot = boot, scenarios = scenarios, designs = designs)
})

input_table <- bind_rows(lapply(candidate_results, `[[`, "point")) %>%
  mutate(
    bootstrap_reps = n_boot,
    endpoint_scale = "favorable_response = -BestAvgResponse",
    standardization = "common pooled SD across A, B, and AB within complete PDX triplets",
    clinical_covariance_model = "independent randomized clinical-arm means"
  ) %>%
  select(
    case_role, component_a, component_b, combination, n_triplet, endpoint_scale,
    standardization, mean_favorable_a, mean_favorable_b, mean_favorable_ab,
    common_pooled_sd, delta_ab_minus_a, delta_ab_minus_b,
    pairwise_delta_ab_minus_a, pairwise_delta_ab_minus_b, additive_excess,
    clinical_covariance_model, bootstrap_reps
  )

bootstrap_table <- bind_rows(lapply(seq_along(candidate_results), function(i) {
  result <- candidate_results[[i]]
  meta <- candidate_plan[i, ]
  valid <- result$boot %>% filter(valid_effect_input)
  tibble(
    case_role = meta$case_role,
    component_a = meta$component_a,
    component_b = meta$component_b,
    combination = meta$combination,
    n_boot_requested = n_boot,
    n_boot_valid = nrow(valid),
    valid_boot_fraction = nrow(valid) / n_boot,
    delta_ab_minus_a_median = safe_quantile(valid$delta_ab_minus_a, 0.50),
    delta_ab_minus_a_lcl = safe_quantile(valid$delta_ab_minus_a, 0.025),
    delta_ab_minus_a_ucl = safe_quantile(valid$delta_ab_minus_a, 0.975),
    delta_ab_minus_b_median = safe_quantile(valid$delta_ab_minus_b, 0.50),
    delta_ab_minus_b_lcl = safe_quantile(valid$delta_ab_minus_b, 0.025),
    delta_ab_minus_b_ucl = safe_quantile(valid$delta_ab_minus_b, 0.975),
    min_contribution_median = safe_quantile(valid$min_contribution, 0.50),
    min_contribution_lcl = safe_quantile(valid$min_contribution, 0.025),
    min_contribution_ucl = safe_quantile(valid$min_contribution, 0.975),
    prob_both_contributions_positive = mean(
      valid$delta_ab_minus_a > 0 & valid$delta_ab_minus_b > 0,
      na.rm = TRUE
    ),
    prob_min_contribution_gt_0_10 = mean(valid$min_contribution > 0.10, na.rm = TRUE)
  )
}))

scenario_input_table <- bind_rows(lapply(candidate_results, `[[`, "scenarios"))

design_table <- bind_rows(lapply(candidate_results, `[[`, "designs")) %>%
  arrange(candidate, factor(input_scenario, levels = c("plug_in", "bootstrap_lower_95_effects")),
          factor(method, levels = c("equal allocation", "joint-decision optimized allocation",
                                    "no finite design under selected lower-bound effects")))

plugin_joint <- design_table %>%
  filter(input_scenario == "plug_in", method == "joint-decision optimized allocation") %>%
  select(candidate, plugin_joint_N = N, plugin_allocation_A = allocation_A,
         plugin_allocation_B = allocation_B, plugin_allocation_AB = allocation_AB,
         plugin_joint_power = achieved_joint_success_probability)

lower_bound_inputs <- scenario_input_table %>%
  filter(input_scenario == "bootstrap_lower_95_effects")

uncertainty_table <- plugin_joint %>%
  left_join(
    design_table %>%
      filter(input_scenario == "bootstrap_lower_95_effects",
             method == "joint-decision optimized allocation") %>%
      select(candidate, lower_bound_joint_N = N,
             lower_bound_allocation_A = allocation_A,
             lower_bound_allocation_B = allocation_B,
             lower_bound_allocation_AB = allocation_AB,
             lower_bound_joint_power = achieved_joint_success_probability,
             lower_bound_reason = reason),
    by = "candidate"
  ) %>%
  left_join(
    lower_bound_inputs %>%
      select(candidate, case_role, component_a, component_b, combination,
             lower_bound_delta_ab_minus_a = delta_ab_minus_a,
             lower_bound_delta_ab_minus_b = delta_ab_minus_b,
             valid_lower_bound_input = valid_design_input),
    by = "candidate"
  ) %>%
  mutate(
    plugin_design_joint_power_under_lower_bound_effects =
      pmap_dbl(
        list(
          plugin_joint_N,
          plugin_allocation_A,
          plugin_allocation_B,
          plugin_allocation_AB,
          lower_bound_delta_ab_minus_a,
          lower_bound_delta_ab_minus_b,
          valid_lower_bound_input
        ),
        function(N, allocation_A, allocation_B, allocation_AB,
                 delta_a, delta_b, valid_input) {
          evaluate_design_under_input(
            tibble(
              N = N,
              allocation_A = allocation_A,
              allocation_B = allocation_B,
              allocation_AB = allocation_AB
            ),
            tibble(
              input_scenario = "bootstrap_lower_95_effects",
              delta_ab_minus_a = delta_a,
              delta_ab_minus_b = delta_b,
              valid_design_input = valid_input
            )
          )
        }
      ),
    uncertainty_sample_size_inflation =
      ifelse(is.finite(lower_bound_joint_N) & is.finite(plugin_joint_N),
             lower_bound_joint_N / plugin_joint_N,
             NA_real_),
    planning_interpretation = dplyr::case_when(
      !valid_lower_bound_input ~
        "lower-bound effects do not support both contribution inputs",
      is.finite(plugin_design_joint_power_under_lower_bound_effects) &
        plugin_design_joint_power_under_lower_bound_effects < 0.50 ~
        "plug-in design has low joint power under lower-bound effects",
      is.finite(uncertainty_sample_size_inflation) &
        uncertainty_sample_size_inflation > 1.50 ~
        "lower-bound effects materially increase the sample size",
      TRUE ~
        "lower-bound effects preserve the planning target"
    ),
    alpha = alpha,
    target_joint_power = target_power,
    clinical_covariance_model = "independent randomized clinical-arm means"
  )

clinical_variance_sensitivity <- scenario_input_table %>%
  filter(candidate == "BYL719 + binimetinib") %>%
  tidyr::crossing(clinical_sd_multiplier = c(0.75, 1.00, 1.25)) %>%
  split(seq_len(nrow(.))) %>%
  purrr::map_dfr(function(row) {
    variance_sensitivity_design(row, row$clinical_sd_multiplier)
  }) %>%
  mutate(
    alpha = alpha,
    target_joint_power = target_power,
    clinical_covariance_model = "independent randomized clinical-arm means"
  )

validation_checks <- list()
record_check <- function(check, pass, observed = NA_character_,
                         expected = NA_character_, detail = NA_character_) {
  tibble(
    check = check,
    pass = isTRUE(pass),
    observed = as.character(observed),
    expected = as.character(expected),
    detail = as.character(detail)
  )
}

validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "source data has expected dimensions from PDX contract",
  nrow(input_df) == 4758 && ncol(input_df) == 12,
  paste(nrow(input_df), ncol(input_df), sep = " x "),
  "4758 x 12 after adding favorable_response",
  "Raw file has 11 columns; script adds one derived endpoint column."
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "all planned PDX candidates have at least 20 complete triplets",
  all(input_table$n_triplet >= 20),
  paste(input_table$combination, input_table$n_triplet, sep = "=", collapse = "; "),
  "all n_triplet >= 20"
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "bootstrap effect validity fraction is high",
  all(bootstrap_table$valid_boot_fraction >= 0.95),
  paste(bootstrap_table$combination, round(bootstrap_table$valid_boot_fraction, 3), sep = "=", collapse = "; "),
  "all valid fractions >= 0.95"
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "clinical covariance model is independent randomized arms",
  all(design_table$clinical_covariance_model == covariance_model),
  paste(unique(design_table$clinical_covariance_model), collapse = "; "),
  covariance_model,
  "PDX-derived arm correlations are not used in clinical design calculations."
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "primary one-sided alpha is 0.025",
  all(design_table$alpha == alpha) && identical(alpha, 0.025),
  alpha,
  "0.025"
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "plug-in joint optimized designs reach target joint power",
  all(plugin_joint$plugin_joint_power >= target_power),
  paste(plugin_joint$candidate, round(plugin_joint$plugin_joint_power, 3), sep = "=", collapse = "; "),
  paste0("all >= ", target_power)
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "plug-in optimized allocations respect lower bound",
  all(plugin_joint$plugin_allocation_A >= min_allocation - 1e-8 &
        plugin_joint$plugin_allocation_B >= min_allocation - 1e-8 &
        plugin_joint$plugin_allocation_AB >= min_allocation - 1e-8),
  paste(plugin_joint$candidate,
        paste(round(plugin_joint$plugin_allocation_A, 3),
              round(plugin_joint$plugin_allocation_B, 3),
              round(plugin_joint$plugin_allocation_AB, 3), sep = "/"),
        sep = "=", collapse = "; "),
  paste0("all arms >= ", min_allocation)
)
validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "invalid lower-bound effect case is explicitly flagged rather than forced",
  any(!uncertainty_table$valid_lower_bound_input) ||
    all(is.finite(uncertainty_table$lower_bound_joint_N)),
  paste(uncertainty_table$combination,
        uncertainty_table$valid_lower_bound_input,
        sep = "=", collapse = "; "),
  "invalid lower-bound effects flagged, valid lower-bound effects produce designs"
)

design_comparison_check <- design_table %>%
  filter(method %in% c("equal allocation", "joint-decision optimized allocation")) %>%
  select(candidate, input_scenario, method, N) %>%
  tidyr::pivot_wider(names_from = method, values_from = N)

validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "optimized allocation does not require more total N than equal allocation",
  all(design_comparison_check[["joint-decision optimized allocation"]] <=
        design_comparison_check[["equal allocation"]], na.rm = TRUE),
  paste(
    design_comparison_check$candidate,
    design_comparison_check$input_scenario,
    paste0(
      design_comparison_check[["equal allocation"]],
      "->",
      design_comparison_check[["joint-decision optimized allocation"]]
    ),
    sep = "=",
    collapse = "; "
  ),
  "joint optimized N <= equal allocation N for all valid scenarios"
)

primary_stress <- uncertainty_table %>%
  filter(combination == "BYL719 + binimetinib")

validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "primary plug-in design is evaluated under lower-bound contribution effects",
  nrow(primary_stress) == 1 &&
    is.finite(primary_stress$plugin_design_joint_power_under_lower_bound_effects),
  round(primary_stress$plugin_design_joint_power_under_lower_bound_effects, 3),
  "finite lower-bound joint power value"
)

validation_checks[[length(validation_checks) + 1L]] <- record_check(
  "clinical variance sensitivity generated for primary PDX case",
  nrow(clinical_variance_sensitivity) == 6 &&
    all(clinical_variance_sensitivity$clinical_sd_multiplier %in% c(0.75, 1, 1.25)),
  paste(nrow(clinical_variance_sensitivity), "rows"),
  "plug-in and lower-bound inputs at three clinical SD multipliers"
)

validation_checks <- bind_rows(validation_checks)

input_table_out <- round_numeric_df(input_table)
bootstrap_table_out <- round_numeric_df(bootstrap_table)
scenario_input_table_out <- round_numeric_df(scenario_input_table)
design_table_out <- round_numeric_df(design_table)
uncertainty_table_out <- round_numeric_df(uncertainty_table)
clinical_variance_sensitivity_out <- round_numeric_df(clinical_variance_sensitivity)

write_csv(input_table_out, file.path(out_dir, "pdx_case_study_inputs_2026-06-20.csv"))
write_csv(bootstrap_table_out, file.path(out_dir, "pdx_case_study_bootstrap_input_summary_2026-06-20.csv"))
write_csv(scenario_input_table_out, file.path(out_dir, "pdx_case_study_input_scenarios_2026-06-20.csv"))
write_csv(design_table_out, file.path(out_dir, "pdx_case_study_design_outputs_2026-06-20.csv"))
write_csv(uncertainty_table_out, file.path(out_dir, "pdx_case_study_uncertainty_planning_2026-06-20.csv"))
write_csv(clinical_variance_sensitivity_out, file.path(out_dir, "pdx_case_study_clinical_variance_sensitivity_2026-06-20.csv"))
write_csv(validation_checks, file.path(out_dir, "pdx_case_study_validation_checks_2026-06-20.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "pdx_case_study_session_info_2026-06-20.txt"))

checks_passed <- sum(validation_checks$pass)
checks_total <- nrow(validation_checks)
cat("PDX case-study analysis complete. Validation checks passed: ",
    checks_passed, "/", checks_total, "
", sep = "")
