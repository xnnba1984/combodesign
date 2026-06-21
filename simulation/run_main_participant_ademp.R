#!/usr/bin/env Rscript

# Main participant-level ADEMP simulation runner for the SBR
# contribution-of-components method.

source("analysis/contribution_operating_characteristics.R")
source("analysis/contribution_sample_size_allocation.R")

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--main" %in% args) {
  "main"
} else if ("--pilot" %in% args || length(args) == 0L) {
  "pilot"
} else {
  stop("Use --pilot or --main", call. = FALSE)
}

stamp <- "2026-06-19"
tag <- paste0("participant_ademp_", mode, "_", stamp)
out_dir <- file.path("simulation", "results")
log_dir <- file.path("simulation", "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

arms <- c("A", "B", "AB")
alpha_levels <- c(0.025, 0.050)
nrep <- if (mode == "main") 50000L else 3000L
seed_base <- if (mode == "main") 202606190L else 202606191L
min_allocation <- 0.08
analytic_match_tolerance <- if (mode == "main") 0.020 else 0.060
error_tolerance <- if (mode == "main") 0.008 else 0.020

round_numeric <- function(x, digits = 6) {
  if (is.numeric(x)) round(x, digits) else x
}

mc_se <- function(p, nrep) {
  if (!is.finite(p)) return(NA_real_)
  sqrt(p * (1 - p) / nrep)
}

mc_ci <- function(p, nrep, z = 1.96) {
  se <- mc_se(p, nrep)
  if (!is.finite(se)) return(c(lower = NA_real_, upper = NA_real_))
  c(lower = max(0, p - z * se), upper = min(1, p + z * se))
}

format_alpha <- function(alpha) {
  gsub("\\.", "p", sprintf("%.3f", alpha))
}

contribution_delta_from_additive_excess <- function(theta_A, theta_B, gamma) {
  additive_excess_to_contribution(theta_A, theta_B, gamma)
}

scenario_row <- function(
    scenario_label,
    scenario_family,
    purpose,
    effect_source,
    truth_delta,
    plan_delta,
    margin = c(AB_minus_A = 0, AB_minus_B = 0),
    sigma2 = c(A = 1, B = 1, AB = 1),
    N = 240,
    variance_method = "common",
    distribution = "normal",
    theta_A = NA_real_,
    theta_B = NA_real_,
    gamma = NA_real_,
    uncertainty_width = 0.10,
    expected_analytic_match = TRUE) {
  data.frame(
    scenario_label = scenario_label,
    scenario_family = scenario_family,
    purpose = purpose,
    effect_source = effect_source,
    theta_A = theta_A,
    theta_B = theta_B,
    gamma = gamma,
    truth_delta_AB_minus_A = unname(truth_delta["AB_minus_A"]),
    truth_delta_AB_minus_B = unname(truth_delta["AB_minus_B"]),
    plan_delta_AB_minus_A = unname(plan_delta["AB_minus_A"]),
    plan_delta_AB_minus_B = unname(plan_delta["AB_minus_B"]),
    margin_AB_minus_A = unname(margin["AB_minus_A"]),
    margin_AB_minus_B = unname(margin["AB_minus_B"]),
    sigma2_A = unname(sigma2["A"]),
    sigma2_B = unname(sigma2["B"]),
    sigma2_AB = unname(sigma2["AB"]),
    N = N,
    variance_method = variance_method,
    distribution = distribution,
    uncertainty_width = uncertainty_width,
    expected_analytic_match = expected_analytic_match,
    stringsAsFactors = FALSE
  )
}

scenario_rows <- list(
  scenario_row(
    scenario_label = "global_null",
    scenario_family = "global_null",
    purpose = "False full-contribution conclusion when neither component contribution is present",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.00, AB_minus_B = 0.00),
    plan_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    N = 240
  ),
  scenario_row(
    scenario_label = "partial_null_AB_minus_A",
    scenario_family = "partial_null",
    purpose = "False full-contribution conclusion when AB does not beat A",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.00, AB_minus_B = 0.45),
    plan_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    N = 240
  ),
  scenario_row(
    scenario_label = "partial_null_AB_minus_B",
    scenario_family = "partial_null",
    purpose = "False full-contribution conclusion when AB does not beat B",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.00),
    plan_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    N = 240
  ),
  scenario_row(
    scenario_label = "balanced_positive",
    scenario_family = "alternative",
    purpose = "Balanced contribution effects under correctly planned effects",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    plan_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    N = 240
  ),
  scenario_row(
    scenario_label = "limiting_AB_minus_A",
    scenario_family = "alternative",
    purpose = "The first contribution contrast is harder and the planned vector is correct",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.28, AB_minus_B = 0.55),
    plan_delta = c(AB_minus_A = 0.28, AB_minus_B = 0.55),
    N = 260
  ),
  scenario_row(
    scenario_label = "limiting_AB_minus_B",
    scenario_family = "alternative",
    purpose = "The second contribution contrast is harder and the planned vector is correct",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.55, AB_minus_B = 0.28),
    plan_delta = c(AB_minus_A = 0.55, AB_minus_B = 0.28),
    N = 260
  ),
  scenario_row(
    scenario_label = "positive_margin",
    scenario_family = "margin_sensitivity",
    purpose = "Both contributions must exceed a positive clinical margin",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.55, AB_minus_B = 0.55),
    plan_delta = c(AB_minus_A = 0.55, AB_minus_B = 0.55),
    margin = c(AB_minus_A = 0.10, AB_minus_B = 0.10),
    N = 280
  ),
  scenario_row(
    scenario_label = "unequal_variance_positive",
    scenario_family = "variance_sensitivity",
    purpose = "Correct effects with unequal arm variances and Welch standard errors",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    plan_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    sigma2 = c(A = 1.40, B = 1.20, AB = 0.90),
    N = 270,
    variance_method = "welch"
  ),
  scenario_row(
    scenario_label = "heavy_tail_global_null",
    scenario_family = "distribution_sensitivity",
    purpose = "False full-contribution conclusion under modest nonnormality",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.00, AB_minus_B = 0.00),
    plan_delta = c(AB_minus_A = 0.45, AB_minus_B = 0.45),
    N = 300,
    distribution = "t5",
    expected_analytic_match = FALSE
  ),
  scenario_row(
    scenario_label = "overoptimistic_planning",
    scenario_family = "effect_misspecification",
    purpose = "The planned effects are stronger than the true effects",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.32, AB_minus_B = 0.32),
    plan_delta = c(AB_minus_A = 0.55, AB_minus_B = 0.55),
    N = 240
  ),
  scenario_row(
    scenario_label = "reversed_limiting_component",
    scenario_family = "allocation_fragility",
    purpose = "The planned harder contribution is the opposite of the true harder contribution",
    effect_source = "direct_delta",
    truth_delta = c(AB_minus_A = 0.25, AB_minus_B = 0.55),
    plan_delta = c(AB_minus_A = 0.55, AB_minus_B = 0.25),
    N = 260
  ),
  {
    delta <- contribution_delta_from_additive_excess(
      theta_A = 0.15,
      theta_B = 0.20,
      gamma = 0.15
    )
    scenario_row(
      scenario_label = "additive_positive_excess",
      scenario_family = "synergy_motivation",
      purpose = "Additive excess parameter translates into both required contribution contrasts",
      effect_source = "additive_excess",
      truth_delta = delta,
      plan_delta = delta,
      N = 320,
      theta_A = 0.15,
      theta_B = 0.20,
      gamma = 0.15
    )
  }
)

base_scenarios <- do.call(rbind, scenario_rows)
scenario_grid <- do.call(rbind, lapply(alpha_levels, function(alpha) {
  out <- base_scenarios
  out$alpha <- alpha
  out$scenario_id <- paste(out$scenario_label, paste0("alpha_", format_alpha(alpha)), sep = "__")
  out
}))
scenario_grid$nrep <- nrep
scenario_grid$seed <- seed_base + seq_len(nrow(scenario_grid)) * 1000L
scenario_grid <- scenario_grid[
  c(
    "scenario_id", "scenario_label", "scenario_family", "purpose",
    "effect_source", "theta_A", "theta_B", "gamma",
    "truth_delta_AB_minus_A", "truth_delta_AB_minus_B",
    "plan_delta_AB_minus_A", "plan_delta_AB_minus_B",
    "margin_AB_minus_A", "margin_AB_minus_B",
    "sigma2_A", "sigma2_B", "sigma2_AB", "N", "alpha",
    "variance_method", "distribution", "uncertainty_width",
    "expected_analytic_match", "nrep", "seed"
  )
]

method_specs <- data.frame(
  method_id = c(
    "equal_iut_full",
    "shared_arm_iut_full",
    "plugin_iut_full",
    "conservative_iut_full",
    "marginal_iut_full",
    "equal_maxT_separate"
  ),
  claim_regime = c(
    "full_contribution",
    "full_contribution",
    "full_contribution",
    "full_contribution",
    "full_contribution",
    "separate_component_claims"
  ),
  adjustment_method = c("iut", "iut", "iut", "iut", "iut", "maxT"),
  allocation_strategy = c(
    "equal",
    "shared_arm_benchmark",
    "plugin_joint_success",
    "conservative_joint_success",
    "marginal_power_comparator",
    "equal"
  ),
  manuscript_role = c(
    "transparent equal-allocation full-contribution baseline",
    "shared-arm benchmark for two contribution comparisons",
    "planned-effect allocation for the full-contribution IUT decision",
    "conservative planned-effect allocation for uncertain contribution inputs",
    "diagnostic comparator that balances the weaker marginal contribution",
    "optional separate component-claim family-wise analysis"
  ),
  stringsAsFactors = FALSE
)

scenario_path <- file.path(out_dir, paste0(tag, "_scenario_grid.csv"))
method_path <- file.path(out_dir, paste0(tag, "_method_specs.csv"))
write.csv(scenario_grid, scenario_path, row.names = FALSE)
write.csv(method_specs, method_path, row.names = FALSE)

arm_means_from_delta <- function(delta) {
  c(
    A = 0,
    B = unname(delta["AB_minus_A"] - delta["AB_minus_B"]),
    AB = unname(delta["AB_minus_A"])
  )
}

draw_arm_matrix <- function(nrep, n, mean, sd, distribution) {
  if (distribution == "normal") {
    return(matrix(stats::rnorm(nrep * n, mean = mean, sd = sd), nrow = nrep))
  }
  if (distribution == "t5") {
    z <- stats::rt(nrep * n, df = 5) / sqrt(5 / 3)
    return(matrix(mean + sd * z, nrow = nrep))
  }
  stop("Unknown distribution: ", distribution, call. = FALSE)
}

simulate_trial_summaries <- function(nrep, n, mu, sigma, distribution, seed) {
  set.seed(seed)
  means <- matrix(NA_real_, nrow = nrep, ncol = length(n))
  vars <- matrix(NA_real_, nrow = nrep, ncol = length(n))
  colnames(means) <- names(n)
  colnames(vars) <- names(n)

  for (arm in names(n)) {
    y <- draw_arm_matrix(
      nrep = nrep,
      n = n[[arm]],
      mean = mu[[arm]],
      sd = sigma[[arm]],
      distribution = distribution
    )
    means[, arm] <- rowMeans(y)
    vars[, arm] <- apply(y, 1L, stats::var)
  }

  list(means = means, vars = vars)
}

estimated_statistics <- function(means, vars, n, margin, variance_method) {
  d1 <- means[, "AB"] - means[, "A"]
  d2 <- means[, "AB"] - means[, "B"]
  u1 <- d1 - unname(margin["AB_minus_A"])
  u2 <- d2 - unname(margin["AB_minus_B"])

  if (variance_method == "common") {
    pooled <- ((n["A"] - 1) * vars[, "A"] +
                 (n["B"] - 1) * vars[, "B"] +
                 (n["AB"] - 1) * vars[, "AB"]) /
      (sum(n) - length(n))
    v_ab <- pooled / n["AB"]
    se1 <- sqrt(pooled / n["AB"] + pooled / n["A"])
    se2 <- sqrt(pooled / n["AB"] + pooled / n["B"])
  } else if (variance_method == "welch") {
    v_ab <- vars[, "AB"] / n["AB"]
    se1 <- sqrt(vars[, "AB"] / n["AB"] + vars[, "A"] / n["A"])
    se2 <- sqrt(vars[, "AB"] / n["AB"] + vars[, "B"] / n["B"])
  } else {
    stop("Unknown variance_method: ", variance_method, call. = FALSE)
  }

  statistic <- cbind(
    AB_minus_A = u1 / se1,
    AB_minus_B = u2 / se2
  )
  list(
    statistic = statistic,
    rho = pmin(pmax(v_ab / (se1 * se2), -0.999), 0.999)
  )
}

maxT_cutoff_cache <- new.env(parent = emptyenv())

maxT_cutoffs_from_rho <- function(rho, alpha) {
  rounded <- round(pmin(pmax(rho, -0.999), 0.999), 5)
  unique_rho <- sort(unique(rounded))
  cutoffs <- vapply(unique_rho, function(r) {
    key <- paste(alpha, r, sep = "|")
    if (!exists(key, envir = maxT_cutoff_cache, inherits = FALSE)) {
      R <- matrix(c(1, r, r, 1), nrow = 2)
      dimnames(R) <- list(
        c("AB_minus_A", "AB_minus_B"),
        c("AB_minus_A", "AB_minus_B")
      )
      assign(
        key,
        critical_value_maxT_one_sided(R, alpha = alpha),
        envir = maxT_cutoff_cache
      )
    }
    get(key, envir = maxT_cutoff_cache, inherits = FALSE)
  }, numeric(1))
  cutoffs[match(rounded, unique_rho)]
}

analyze_rejections <- function(statistic, rho, alpha, claim_regime) {
  if (claim_regime == "full_contribution") {
    cutoff <- critical_value_iut(alpha)
    reject <- statistic > cutoff
    return(list(
      cutoff = rep(cutoff, nrow(statistic)),
      reject = reject,
      joint_success = rowSums(reject) == ncol(reject)
    ))
  }

  if (claim_regime == "separate_component_claims") {
    cutoff <- maxT_cutoffs_from_rho(rho, alpha = alpha)
    reject <- cbind(
      AB_minus_A = statistic[, "AB_minus_A"] > cutoff,
      AB_minus_B = statistic[, "AB_minus_B"] > cutoff
    )
    return(list(
      cutoff = cutoff,
      reject = reject,
      joint_success = rowSums(reject) == ncol(reject)
    ))
  }

  stop("Unknown claim_regime: ", claim_regime, call. = FALSE)
}

conservative_delta <- function(plan_delta, margin, width) {
  out <- plan_delta - width
  lower <- margin + 0.02
  pmax(out, lower)
}

shared_arm_allocation <- function() {
  shared_weight <- sqrt(2)
  out <- c(A = 1, B = 1, AB = shared_weight)
  out / sum(out)
}

allocation_cache <- new.env(parent = emptyenv())

allocation_cache_key <- function(row, spec, planning_delta) {
  paste(
    spec$allocation_strategy,
    spec$claim_regime,
    spec$adjustment_method,
    row$N,
    row$alpha,
    row$margin_AB_minus_A,
    row$margin_AB_minus_B,
    row$sigma2_A,
    row$sigma2_B,
    row$sigma2_AB,
    sprintf("%.6f", unname(planning_delta["AB_minus_A"])),
    sprintf("%.6f", unname(planning_delta["AB_minus_B"])),
    sep = "|"
  )
}

allocation_for_method <- function(row, spec) {
  allocation_strategy <- spec$allocation_strategy
  margin <- c(
    AB_minus_A = row$margin_AB_minus_A,
    AB_minus_B = row$margin_AB_minus_B
  )
  plan_delta <- c(
    AB_minus_A = row$plan_delta_AB_minus_A,
    AB_minus_B = row$plan_delta_AB_minus_B
  )
  planning_delta <- if (allocation_strategy == "conservative_joint_success") {
    conservative_delta(plan_delta, margin, row$uncertainty_width)
  } else {
    plan_delta
  }
  sigma2 <- c(A = row$sigma2_A, B = row$sigma2_B, AB = row$sigma2_AB)
  key <- allocation_cache_key(row, spec, planning_delta)
  if (exists(key, envir = allocation_cache, inherits = FALSE)) {
    return(get(key, envir = allocation_cache, inherits = FALSE))
  }

  allocation <- if (allocation_strategy == "equal") {
    equal_allocation(arms)
  } else if (allocation_strategy == "shared_arm_benchmark") {
    shared_arm_allocation()
  } else if (allocation_strategy %in% c(
    "plugin_joint_success",
    "conservative_joint_success"
  )) {
    optimize_allocation_for_joint_success(
      N = row$N,
      delta = planning_delta,
      sigma2 = sigma2,
      alpha = row$alpha,
      margin = margin,
      claim_regime = spec$claim_regime,
      adjustment_method = spec$adjustment_method,
      covariance_model = "independent_arms",
      arms = arms,
      min_allocation = min_allocation,
      control = list(maxit = 120)
    )$allocation
  } else if (allocation_strategy == "marginal_power_comparator") {
    optimize_allocation_for_min_marginal_power(
      N = row$N,
      delta = planning_delta,
      sigma2 = sigma2,
      alpha = row$alpha,
      margin = margin,
      claim_regime = spec$claim_regime,
      adjustment_method = spec$adjustment_method,
      covariance_model = "independent_arms",
      arms = arms,
      min_allocation = min_allocation,
      control = list(maxit = 120)
    )$allocation
  } else {
    stop("Unknown allocation_strategy: ", allocation_strategy, call. = FALSE)
  }
  assign(key, allocation, envir = allocation_cache)
  allocation
}

analytic_design <- function(row, spec, allocation, delta) {
  design_from_allocation(
    N = row$N,
    allocation = allocation,
    delta = delta,
    sigma2 = c(A = row$sigma2_A, B = row$sigma2_B, AB = row$sigma2_AB),
    alpha = row$alpha,
    margin = c(
      AB_minus_A = row$margin_AB_minus_A,
      AB_minus_B = row$margin_AB_minus_B
    ),
    claim_regime = spec$claim_regime,
    adjustment_method = spec$adjustment_method,
    covariance_model = "independent_arms",
    integer_allocation = TRUE,
    arms = arms
  )
}

run_one <- function(row, spec, row_index, method_index) {
  truth_delta <- c(
    AB_minus_A = row$truth_delta_AB_minus_A,
    AB_minus_B = row$truth_delta_AB_minus_B
  )
  plan_delta <- c(
    AB_minus_A = row$plan_delta_AB_minus_A,
    AB_minus_B = row$plan_delta_AB_minus_B
  )
  margin <- c(
    AB_minus_A = row$margin_AB_minus_A,
    AB_minus_B = row$margin_AB_minus_B
  )
  sigma2 <- c(A = row$sigma2_A, B = row$sigma2_B, AB = row$sigma2_AB)
  allocation <- allocation_for_method(row, spec)
  n <- allocation_to_integer_sample_sizes(row$N, allocation, arms = arms)
  mu <- arm_means_from_delta(truth_delta)
  seed <- row$seed + method_index

  sim <- simulate_trial_summaries(
    nrep = row$nrep,
    n = n,
    mu = mu,
    sigma = sqrt(sigma2),
    distribution = row$distribution,
    seed = seed
  )
  est <- estimated_statistics(
    means = sim$means,
    vars = sim$vars,
    n = n,
    margin = margin,
    variance_method = row$variance_method
  )
  rej <- analyze_rejections(
    statistic = est$statistic,
    rho = est$rho,
    alpha = row$alpha,
    claim_regime = spec$claim_regime
  )

  true_null <- truth_delta <= margin + 1e-12
  empirical_joint_success <- mean(rej$joint_success)
  empirical_marginal <- colMeans(rej$reject)
  empirical_any_true_reject <- if (any(true_null)) {
    mean(rowSums(rej$reject[, true_null, drop = FALSE]) > 0)
  } else {
    NA_real_
  }
  empirical_false_full_contribution <- if (any(true_null)) {
    empirical_joint_success
  } else {
    NA_real_
  }

  analytic_truth <- analytic_design(row, spec, allocation, truth_delta)
  analytic_plan <- analytic_design(row, spec, allocation, plan_delta)
  analytic_any_true_reject <- if (any(true_null)) {
    prob_any_true_null_rejected(
      lambda = analytic_truth$lambda,
      R = analytic_truth$R_Delta,
      c_alpha = analytic_truth$c_alpha,
      true_null = true_null
    )
  } else {
    NA_real_
  }

  joint_ci <- mc_ci(empirical_joint_success, row$nrep)
  any_true_ci <- mc_ci(empirical_any_true_reject, row$nrep)
  false_full_ci <- mc_ci(empirical_false_full_contribution, row$nrep)

  data.frame(
    scenario_id = row$scenario_id,
    scenario_label = row$scenario_label,
    scenario_family = row$scenario_family,
    purpose = row$purpose,
    effect_source = row$effect_source,
    theta_A = row$theta_A,
    theta_B = row$theta_B,
    gamma = row$gamma,
    method_id = spec$method_id,
    claim_regime = spec$claim_regime,
    adjustment_method = spec$adjustment_method,
    allocation_strategy = spec$allocation_strategy,
    manuscript_role = spec$manuscript_role,
    N = row$N,
    alpha = row$alpha,
    nrep = row$nrep,
    variance_method = row$variance_method,
    distribution = row$distribution,
    expected_analytic_match = row$expected_analytic_match,
    truth_delta_AB_minus_A = unname(truth_delta["AB_minus_A"]),
    truth_delta_AB_minus_B = unname(truth_delta["AB_minus_B"]),
    plan_delta_AB_minus_A = unname(plan_delta["AB_minus_A"]),
    plan_delta_AB_minus_B = unname(plan_delta["AB_minus_B"]),
    margin_AB_minus_A = unname(margin["AB_minus_A"]),
    margin_AB_minus_B = unname(margin["AB_minus_B"]),
    sigma2_A = unname(sigma2["A"]),
    sigma2_B = unname(sigma2["B"]),
    sigma2_AB = unname(sigma2["AB"]),
    allocation_A = unname(allocation["A"]),
    allocation_B = unname(allocation["B"]),
    allocation_AB = unname(allocation["AB"]),
    n_A = unname(n["A"]),
    n_B = unname(n["B"]),
    n_AB = unname(n["AB"]),
    true_null_AB_minus_A = unname(true_null["AB_minus_A"]),
    true_null_AB_minus_B = unname(true_null["AB_minus_B"]),
    cutoff_mean = mean(rej$cutoff),
    cutoff_min = min(rej$cutoff),
    cutoff_max = max(rej$cutoff),
    empirical_joint_success = empirical_joint_success,
    empirical_joint_success_lcl = joint_ci["lower"],
    empirical_joint_success_ucl = joint_ci["upper"],
    analytic_joint_success = analytic_truth$joint_success_probability,
    planned_analytic_joint_success = analytic_plan$joint_success_probability,
    joint_success_absolute_difference =
      abs(empirical_joint_success - analytic_truth$joint_success_probability),
    empirical_any_true_reject = empirical_any_true_reject,
    empirical_any_true_reject_lcl = any_true_ci["lower"],
    empirical_any_true_reject_ucl = any_true_ci["upper"],
    analytic_any_true_reject = analytic_any_true_reject,
    empirical_false_full_contribution = empirical_false_full_contribution,
    empirical_false_full_contribution_lcl = false_full_ci["lower"],
    empirical_false_full_contribution_ucl = false_full_ci["upper"],
    analytic_false_full_contribution =
      if (any(true_null)) analytic_truth$joint_success_probability else NA_real_,
    empirical_marginal_power_AB_minus_A =
      unname(empirical_marginal["AB_minus_A"]),
    empirical_marginal_power_AB_minus_B =
      unname(empirical_marginal["AB_minus_B"]),
    analytic_marginal_power_AB_minus_A =
      unname(analytic_truth$marginal_power["AB_minus_A"]),
    analytic_marginal_power_AB_minus_B =
      unname(analytic_truth$marginal_power["AB_minus_B"]),
    seed = seed,
    stringsAsFactors = FALSE
  )
}

result_rows <- list()
allocation_rows <- list()
row_counter <- 0L
total_rows <- nrow(scenario_grid) * nrow(method_specs)

for (i in seq_len(nrow(scenario_grid))) {
  row <- scenario_grid[i, ]
  for (j in seq_len(nrow(method_specs))) {
    spec <- method_specs[j, ]
    row_counter <- row_counter + 1L
    result <- run_one(row, spec, i, j)
    result_rows[[row_counter]] <- result
    allocation_rows[[row_counter]] <- result[
      c(
        "scenario_id", "scenario_label", "scenario_family", "method_id",
        "claim_regime", "adjustment_method", "allocation_strategy",
        "alpha", "N", "truth_delta_AB_minus_A", "truth_delta_AB_minus_B",
        "plan_delta_AB_minus_A", "plan_delta_AB_minus_B",
        "margin_AB_minus_A", "margin_AB_minus_B",
        "allocation_A", "allocation_B", "allocation_AB",
        "n_A", "n_B", "n_AB",
        "planned_analytic_joint_success"
      )
    ]
    message(sprintf(
      "[%03d/%03d] %s | %s | empirical %.4f | analytic %.4f",
      row_counter,
      total_rows,
      row$scenario_id,
      spec$method_id,
      result$empirical_joint_success,
      result$analytic_joint_success
    ))
  }
}

results <- do.call(rbind, result_rows)
allocations <- do.call(rbind, allocation_rows)
results <- as.data.frame(lapply(results, round_numeric), stringsAsFactors = FALSE)
allocations <- as.data.frame(lapply(allocations, round_numeric), stringsAsFactors = FALSE)

result_path <- file.path(out_dir, paste0(tag, "_operating_characteristics.csv"))
allocation_path <- file.path(out_dir, paste0(tag, "_method_allocations.csv"))
write.csv(results, result_path, row.names = FALSE)
write.csv(allocations, allocation_path, row.names = FALSE)

is_iut <- results$claim_regime == "full_contribution"
is_maxT <- results$claim_regime == "separate_component_claims"
is_primary_alpha <- abs(results$alpha - 0.025) < 1e-12
has_true_null <- results$true_null_AB_minus_A | results$true_null_AB_minus_B
matched_normal <- results$expected_analytic_match &
  results$distribution == "normal"

check_rows <- list()
error_bound_label <- if (mode == "main") "main bound" else "pilot bound"
add_check <- function(check, passed, observed, expected) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

add_check(
  "complete method by scenario grid",
  nrow(results) == nrow(scenario_grid) * nrow(method_specs),
  nrow(results),
  nrow(scenario_grid) * nrow(method_specs)
)
add_check(
  "all integer arm sizes positive and sum to total N",
  all(results$n_A > 0 & results$n_B > 0 & results$n_AB > 0 &
        results$n_A + results$n_B + results$n_AB == results$N),
  "checked all rows",
  "positive arm sizes and n_A + n_B + n_AB = N"
)
add_check(
  "reported allocations sum to one within CSV rounding tolerance",
  max(abs(results$allocation_A + results$allocation_B +
            results$allocation_AB - 1)) < 1e-5,
  max(abs(results$allocation_A + results$allocation_B +
            results$allocation_AB - 1)),
  "< 1e-5"
)
primary_iut_null <- is_iut & is_primary_alpha & has_true_null &
  results$scenario_family %in% c("global_null", "partial_null")
add_check(
  paste("primary IUT false full-contribution error not above", error_bound_label),
  max(results$empirical_false_full_contribution[primary_iut_null], na.rm = TRUE) <=
    0.025 + error_tolerance,
  max(results$empirical_false_full_contribution[primary_iut_null], na.rm = TRUE),
  paste0("<= ", 0.025 + error_tolerance)
)
primary_maxT_global <- is_maxT & is_primary_alpha &
  results$scenario_family == "global_null"
add_check(
  paste("primary maxT separate-claim family-wise error not above", error_bound_label),
  max(results$empirical_any_true_reject[primary_maxT_global], na.rm = TRUE) <=
    0.025 + error_tolerance,
  max(results$empirical_any_true_reject[primary_maxT_global], na.rm = TRUE),
  paste0("<= ", 0.025 + error_tolerance)
)
matched_rows_for_diff <- matched_normal & is_iut
add_check(
  "matched normal IUT empirical joint power near analytic reference",
  max(results$joint_success_absolute_difference[matched_rows_for_diff], na.rm = TRUE) <=
    analytic_match_tolerance,
  max(results$joint_success_absolute_difference[matched_rows_for_diff], na.rm = TRUE),
  paste0("<= ", analytic_match_tolerance)
)
positive_iut <- is_iut & is_primary_alpha &
  results$scenario_label == "balanced_positive"
global_iut <- is_iut & is_primary_alpha &
  results$scenario_label == "global_null"
add_check(
  "balanced positive scenario has higher joint power than global null",
  min(results$empirical_joint_success[positive_iut], na.rm = TRUE) >
    max(results$empirical_joint_success[global_iut], na.rm = TRUE),
  paste(
    min(results$empirical_joint_success[positive_iut], na.rm = TRUE),
    max(results$empirical_joint_success[global_iut], na.rm = TRUE),
    sep = " > "
  ),
  "minimum positive-scenario value exceeds maximum global-null value"
)
add_check(
  "alpha 0.05 analytic joint power at least alpha 0.025 for same design rows",
  {
    keys <- paste(
      results$scenario_label,
      results$method_id,
      results$claim_regime,
      results$allocation_strategy,
      sep = "|"
    )
    pass <- TRUE
    for (key in unique(keys)) {
      idx025 <- which(keys == key & abs(results$alpha - 0.025) < 1e-12)
      idx050 <- which(keys == key & abs(results$alpha - 0.050) < 1e-12)
      if (length(idx025) == 1L && length(idx050) == 1L) {
        pass <- pass && results$analytic_joint_success[idx050] + 1e-12 >=
          results$analytic_joint_success[idx025]
      }
    }
    pass
  },
  "checked paired alpha rows",
  "alpha 0.05 >= alpha 0.025 analytically"
)
nonmisspec_alt <- is_iut & is_primary_alpha &
  results$scenario_family %in% c("alternative", "margin_sensitivity",
                                 "variance_sensitivity", "synergy_motivation")
plugin_rows <- nonmisspec_alt & results$method_id == "plugin_iut_full"
equal_rows <- nonmisspec_alt & results$method_id == "equal_iut_full"
plugin_compare <- merge(
  results[plugin_rows, c("scenario_label", "analytic_joint_success")],
  results[equal_rows, c("scenario_label", "analytic_joint_success")],
  by = "scenario_label",
  suffixes = c("_plugin", "_equal")
)
add_check(
  "plug-in IUT allocation is not analytically worse than equal allocation in correctly planned alternatives",
  all(plugin_compare$analytic_joint_success_plugin + 0.005 >=
        plugin_compare$analytic_joint_success_equal),
  paste(round(plugin_compare$analytic_joint_success_plugin -
                plugin_compare$analytic_joint_success_equal, 6),
        collapse = "; "),
  "plugin minus equal >= -0.005"
)

checks <- do.call(rbind, check_rows)
check_path <- file.path(out_dir, paste0(tag, "_validation_checks.csv"))
write.csv(checks, check_path, row.names = FALSE)

error_control_table <- subset(
  results,
  has_true_null & alpha == 0.025,
  select = c(
    "scenario_label", "scenario_family", "method_id", "claim_regime",
    "allocation_strategy", "alpha", "N", "nrep",
    "truth_delta_AB_minus_A", "truth_delta_AB_minus_B",
    "margin_AB_minus_A", "margin_AB_minus_B",
    "empirical_false_full_contribution",
    "empirical_false_full_contribution_lcl",
    "empirical_false_full_contribution_ucl",
    "analytic_false_full_contribution",
    "empirical_any_true_reject",
    "empirical_any_true_reject_lcl",
    "empirical_any_true_reject_ucl",
    "analytic_any_true_reject"
  )
)
joint_success_table <- subset(
  results,
  !has_true_null & alpha == 0.025,
  select = c(
    "scenario_label", "scenario_family", "method_id", "claim_regime",
    "allocation_strategy", "alpha", "N", "nrep",
    "truth_delta_AB_minus_A", "truth_delta_AB_minus_B",
    "plan_delta_AB_minus_A", "plan_delta_AB_minus_B",
    "margin_AB_minus_A", "margin_AB_minus_B",
    "allocation_A", "allocation_B", "allocation_AB",
    "empirical_joint_success", "empirical_joint_success_lcl",
    "empirical_joint_success_ucl", "analytic_joint_success",
    "planned_analytic_joint_success",
    "empirical_marginal_power_AB_minus_A",
    "empirical_marginal_power_AB_minus_B"
  )
)
fragility_table <- subset(
  results,
  scenario_family %in% c("effect_misspecification", "allocation_fragility") &
    claim_regime == "full_contribution" & alpha == 0.025,
  select = c(
    "scenario_label", "method_id", "allocation_strategy", "N", "nrep",
    "truth_delta_AB_minus_A", "truth_delta_AB_minus_B",
    "plan_delta_AB_minus_A", "plan_delta_AB_minus_B",
    "allocation_A", "allocation_B", "allocation_AB",
    "empirical_joint_success", "analytic_joint_success",
    "planned_analytic_joint_success"
  )
)
robustness_table <- subset(
  results,
  scenario_family %in% c(
    "variance_sensitivity", "distribution_sensitivity", "margin_sensitivity"
  ) & alpha == 0.025,
  select = c(
    "scenario_label", "scenario_family", "method_id", "claim_regime",
    "variance_method", "distribution", "alpha", "N", "nrep",
    "empirical_joint_success", "empirical_joint_success_lcl",
    "empirical_joint_success_ucl", "analytic_joint_success",
    "empirical_false_full_contribution",
    "empirical_any_true_reject"
  )
)

table_error_path <- file.path(out_dir, paste0(tag, "_table_error_control.csv"))
table_joint_path <- file.path(out_dir, paste0(tag, "_table_joint_success.csv"))
table_fragility_path <- file.path(out_dir, paste0(tag, "_table_allocation_fragility.csv"))
table_robustness_path <- file.path(out_dir, paste0(tag, "_table_robustness.csv"))
write.csv(error_control_table, table_error_path, row.names = FALSE)
write.csv(joint_success_table, table_joint_path, row.names = FALSE)
write.csv(fragility_table, table_fragility_path, row.names = FALSE)
write.csv(robustness_table, table_robustness_path, row.names = FALSE)

get_value <- function(filter, column, fun = identity) {
  x <- results[filter, column]
  if (length(x) == 0L) return(NA_real_)
  fun(x)
}

max_iut_false_full <- get_value(
  is_iut & is_primary_alpha & has_true_null &
    results$scenario_family %in% c("global_null", "partial_null"),
  "empirical_false_full_contribution",
  max
)
max_maxT_fwer <- get_value(
  is_maxT & is_primary_alpha & results$scenario_family == "global_null",
  "empirical_any_true_reject",
  max
)
balanced_equal <- get_value(
  is_iut & is_primary_alpha & results$scenario_label == "balanced_positive" &
    results$method_id == "equal_iut_full",
  "empirical_joint_success"
)
balanced_plugin <- get_value(
  is_iut & is_primary_alpha & results$scenario_label == "balanced_positive" &
    results$method_id == "plugin_iut_full",
  "empirical_joint_success"
)
positive_margin_plugin <- get_value(
  is_iut & is_primary_alpha & results$scenario_label == "positive_margin" &
    results$method_id == "plugin_iut_full",
  "empirical_joint_success"
)
misspec_plugin <- get_value(
  is_iut & is_primary_alpha &
    results$scenario_label == "reversed_limiting_component" &
    results$method_id == "plugin_iut_full",
  "empirical_joint_success"
)
misspec_conservative <- get_value(
  is_iut & is_primary_alpha &
    results$scenario_label == "reversed_limiting_component" &
    results$method_id == "conservative_iut_full",
  "empirical_joint_success"
)

finding_heading <- if (mode == "main") "## Main Findings" else "## Pilot Findings"
interpretation_lines <- if (mode == "main") {
  c(
    "- This main run is the stable simulation evidence source for selecting manuscript tables and figures.",
    "- The primary checks focus on false full-contribution conclusions for IUT and family-wise error for optional separate claims under maxT.",
    "- Positive-margin, variance, distribution, and planning-misspecification scenarios support both strengths and limitations.",
    "- The next step is to regenerate simulation displays from the main-mode results."
  )
} else {
  c(
    "- This pilot is a computation and design gate, not final manuscript evidence.",
    "- The primary checks focus on false full-contribution conclusions for IUT and family-wise error for optional separate claims under maxT.",
    "- Positive-margin, variance, distribution, and planning-misspecification scenarios are included so the full simulation can support both strengths and limitations.",
    "- If the pilot checks pass, run the main mode with more replicates."
  )
}

if (!all(checks$passed)) {
  failed <- checks$check[!checks$passed]
  stop(
    "Participant-level ADEMP ", mode, " run failed checks: ",
    paste(failed, collapse = "; "),
    call. = FALSE
  )
}

message(sprintf(
  "Participant-level ADEMP %s run passed: %d / %d checks.",
  mode,
  sum(checks$passed),
  nrow(checks)
))
