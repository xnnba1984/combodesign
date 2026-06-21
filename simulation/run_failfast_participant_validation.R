#!/usr/bin/env Rscript

# Participant-level fail-fast runner for the SBR simulation.

source("analysis/contribution_operating_characteristics.R")
source("analysis/contribution_sample_size_allocation.R")

stamp <- "2026-06-19"
out_dir <- file.path("simulation", "results")
log_dir <- file.path("simulation", "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

round_numeric <- function(x, digits = 6) {
  if (is.numeric(x)) round(x, digits) else x
}

mc_se <- function(p, nrep) {
  if (!is.finite(p)) return(NA_real_)
  sqrt(p * (1 - p) / nrep)
}

mc_ci <- function(p, nrep, z = 1.96) {
  se <- mc_se(p, nrep)
  c(lower = max(0, p - z * se), upper = min(1, p + z * se))
}

arm_means_from_delta <- function(delta) {
  required <- c("AB_minus_A", "AB_minus_B")
  if (!all(required %in% names(delta))) {
    stop("delta must contain AB_minus_A and AB_minus_B", call. = FALSE)
  }
  c(
    A = 0,
    B = unname(delta["AB_minus_A"] - delta["AB_minus_B"]),
    AB = unname(delta["AB_minus_A"])
  )
}

draw_arm_matrix <- function(nrep, n, mean, sd, distribution) {
  if (distribution == "normal") {
    matrix(stats::rnorm(nrep * n, mean = mean, sd = sd), nrow = nrep)
  } else if (distribution == "t5") {
    z <- stats::rt(nrep * n, df = 5) / sqrt(5 / 3)
    matrix(mean + sd * z, nrow = nrep)
  } else {
    stop("Unknown distribution: ", distribution, call. = FALSE)
  }
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

  rho <- v_ab / (se1 * se2)
  z1 <- u1 / se1
  z2 <- u2 / se2
  statistic <- cbind(AB_minus_A = z1, AB_minus_B = z2)

  list(statistic = statistic, rho = rho)
}

maxT_cutoffs_from_rho <- function(rho, alpha) {
  rounded <- round(pmin(pmax(rho, -0.999), 0.999), 5)
  unique_rho <- sort(unique(rounded))
  cutoffs <- vapply(unique_rho, function(r) {
    R <- matrix(c(1, r, r, 1), nrow = 2)
    dimnames(R) <- list(c("AB_minus_A", "AB_minus_B"),
                        c("AB_minus_A", "AB_minus_B"))
    critical_value_maxT_one_sided(R, alpha = alpha)
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

analytic_for_scenario <- function(scenario, n) {
  delta <- c(
    AB_minus_A = scenario$delta_AB_A,
    AB_minus_B = scenario$delta_AB_B
  )
  margin <- c(
    AB_minus_A = scenario$margin_AB_A,
    AB_minus_B = scenario$margin_AB_B
  )
  sigma2 <- c(
    A = scenario$sigma2_A,
    B = scenario$sigma2_B,
    AB = scenario$sigma2_AB
  )
  design_from_allocation(
    N = sum(n),
    allocation = n / sum(n),
    delta = delta,
    sigma2 = sigma2,
    alpha = scenario$alpha,
    margin = margin,
    claim_regime = scenario$claim_regime,
    adjustment_method = if (scenario$claim_regime == "full_contribution") {
      NULL
    } else {
      "maxT"
    },
    covariance_model = "independent_arms",
    integer_allocation = TRUE
  )
}

scenario_rows <- list(
  data.frame(
    scenario_id = "global_null_iut",
    purpose = "IUT full contribution error under global null",
    N = 240,
    allocation_A = 1 / 3,
    allocation_B = 1 / 3,
    allocation_AB = 1 / 3,
    delta_AB_A = 0,
    delta_AB_B = 0,
    margin_AB_A = 0,
    margin_AB_B = 0,
    sigma2_A = 1,
    sigma2_B = 1,
    sigma2_AB = 1,
    alpha = 0.025,
    claim_regime = "full_contribution",
    variance_method = "common",
    distribution = "normal",
    nrep = 8000,
    seed = 19101,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario_id = "partial_null_iut",
    purpose = "IUT full contribution error under partial null",
    N = 240,
    allocation_A = 1 / 3,
    allocation_B = 1 / 3,
    allocation_AB = 1 / 3,
    delta_AB_A = 0,
    delta_AB_B = 0.50,
    margin_AB_A = 0,
    margin_AB_B = 0,
    sigma2_A = 1,
    sigma2_B = 1,
    sigma2_AB = 1,
    alpha = 0.025,
    claim_regime = "full_contribution",
    variance_method = "common",
    distribution = "normal",
    nrep = 8000,
    seed = 19102,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario_id = "positive_alt_iut",
    purpose = "IUT participant-level power versus analytic approximation",
    N = 240,
    allocation_A = 1 / 3,
    allocation_B = 1 / 3,
    allocation_AB = 1 / 3,
    delta_AB_A = 0.45,
    delta_AB_B = 0.55,
    margin_AB_A = 0,
    margin_AB_B = 0,
    sigma2_A = 1,
    sigma2_B = 1,
    sigma2_AB = 1,
    alpha = 0.025,
    claim_regime = "full_contribution",
    variance_method = "common",
    distribution = "normal",
    nrep = 8000,
    seed = 19103,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario_id = "separate_claim_global_maxT",
    purpose = "maxT family-wise error for separate claims under global null",
    N = 240,
    allocation_A = 1 / 3,
    allocation_B = 1 / 3,
    allocation_AB = 1 / 3,
    delta_AB_A = 0,
    delta_AB_B = 0,
    margin_AB_A = 0,
    margin_AB_B = 0,
    sigma2_A = 1,
    sigma2_B = 1,
    sigma2_AB = 1,
    alpha = 0.025,
    claim_regime = "separate_component_claims",
    variance_method = "common",
    distribution = "normal",
    nrep = 8000,
    seed = 19104,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario_id = "positive_margin_iut",
    purpose = "positive margins reduce full contribution power",
    N = 240,
    allocation_A = 1 / 3,
    allocation_B = 1 / 3,
    allocation_AB = 1 / 3,
    delta_AB_A = 0.45,
    delta_AB_B = 0.55,
    margin_AB_A = 0.10,
    margin_AB_B = 0.10,
    sigma2_A = 1,
    sigma2_B = 1,
    sigma2_AB = 1,
    alpha = 0.025,
    claim_regime = "full_contribution",
    variance_method = "common",
    distribution = "normal",
    nrep = 8000,
    seed = 19105,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario_id = "unequal_variance_alt_iut",
    purpose = "IUT under unequal variances with Welch standard errors",
    N = 260,
    allocation_A = 0.30,
    allocation_B = 0.30,
    allocation_AB = 0.40,
    delta_AB_A = 0.45,
    delta_AB_B = 0.55,
    margin_AB_A = 0,
    margin_AB_B = 0,
    sigma2_A = 1.20,
    sigma2_B = 1.60,
    sigma2_AB = 1.00,
    alpha = 0.025,
    claim_regime = "full_contribution",
    variance_method = "welch",
    distribution = "normal",
    nrep = 8000,
    seed = 19106,
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario_id = "heavy_tail_null_iut",
    purpose = "IUT robustness warning under heavy-tailed null outcomes",
    N = 300,
    allocation_A = 1 / 3,
    allocation_B = 1 / 3,
    allocation_AB = 1 / 3,
    delta_AB_A = 0,
    delta_AB_B = 0,
    margin_AB_A = 0,
    margin_AB_B = 0,
    sigma2_A = 1,
    sigma2_B = 1,
    sigma2_AB = 1,
    alpha = 0.025,
    claim_regime = "full_contribution",
    variance_method = "common",
    distribution = "t5",
    nrep = 8000,
    seed = 19107,
    stringsAsFactors = FALSE
  )
)

scenarios <- do.call(rbind, scenario_rows)
scenario_path <- file.path(out_dir, paste0("participant_failfast_scenarios_", stamp, ".csv"))
write.csv(scenarios, scenario_path, row.names = FALSE)

run_one_scenario <- function(row) {
  allocation <- c(
    A = row$allocation_A,
    B = row$allocation_B,
    AB = row$allocation_AB
  )
  n <- allocation_to_integer_sample_sizes(row$N, allocation)
  delta <- c(
    AB_minus_A = row$delta_AB_A,
    AB_minus_B = row$delta_AB_B
  )
  margin <- c(
    AB_minus_A = row$margin_AB_A,
    AB_minus_B = row$margin_AB_B
  )
  sigma2 <- c(
    A = row$sigma2_A,
    B = row$sigma2_B,
    AB = row$sigma2_AB
  )
  sigma <- sqrt(sigma2)
  mu <- arm_means_from_delta(delta)

  sim <- simulate_trial_summaries(
    nrep = row$nrep,
    n = n,
    mu = mu,
    sigma = sigma,
    distribution = row$distribution,
    seed = row$seed
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
    claim_regime = row$claim_regime
  )

  true_null <- delta <= margin + 1e-12
  joint_success <- mean(rej$joint_success)
  marginal_reject <- colMeans(rej$reject)
  any_true_reject <- if (any(true_null)) {
    mean(rowSums(rej$reject[, true_null, drop = FALSE]) > 0)
  } else {
    NA_real_
  }
  full_type_i <- if (any(true_null)) joint_success else NA_real_

  analytic <- analytic_for_scenario(row, n)
  if (row$claim_regime == "separate_component_claims") {
    analytic_error <- evaluate_null_states(
      delta_alt = delta,
      Sigma_D = analytic$Sigma_D,
      alpha = row$alpha,
      margin = margin,
      claim_regime = row$claim_regime,
      adjustment_method = "maxT"
    )
    analytic_target <- analytic_error$familywise_error_any_true_null[
      analytic_error$state == "global_null"
    ][1]
  } else {
    analytic_target <- analytic$joint_success_probability
  }
  ci <- mc_ci(joint_success, row$nrep)
  fwer_ci <- mc_ci(any_true_reject, row$nrep)

  data.frame(
    scenario_id = row$scenario_id,
    purpose = row$purpose,
    claim_regime = row$claim_regime,
    variance_method = row$variance_method,
    distribution = row$distribution,
    alpha = row$alpha,
    nrep = row$nrep,
    n_A = unname(n["A"]),
    n_B = unname(n["B"]),
    n_AB = unname(n["AB"]),
    true_null_AB_minus_A = unname(true_null["AB_minus_A"]),
    true_null_AB_minus_B = unname(true_null["AB_minus_B"]),
    empirical_joint_success = joint_success,
    empirical_joint_success_lcl = ci["lower"],
    empirical_joint_success_ucl = ci["upper"],
    empirical_any_true_reject = any_true_reject,
    empirical_any_true_reject_lcl = fwer_ci["lower"],
    empirical_any_true_reject_ucl = fwer_ci["upper"],
    empirical_full_type_i = full_type_i,
    analytic_reference = analytic_target,
    absolute_difference_from_analytic =
      ifelse(is.finite(analytic_target), abs(joint_success - analytic_target), NA_real_),
    marginal_reject_AB_minus_A = unname(marginal_reject["AB_minus_A"]),
    marginal_reject_AB_minus_B = unname(marginal_reject["AB_minus_B"]),
    mean_cutoff = mean(rej$cutoff),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  run_one_scenario(scenarios[i, ])
}))
results <- as.data.frame(lapply(results, round_numeric), stringsAsFactors = FALSE)
result_path <- file.path(out_dir, paste0("participant_failfast_results_", stamp, ".csv"))
write.csv(results, result_path, row.names = FALSE)

get_result <- function(id, col) {
  results[results$scenario_id == id, col][[1]]
}

check_rows <- list()
add_check <- function(check, passed, observed, expected) {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

alpha_main <- 0.025
tolerance <- 0.012
add_check(
  "IUT global-null full contribution error controlled",
  get_result("global_null_iut", "empirical_full_type_i") <= alpha_main + tolerance,
  get_result("global_null_iut", "empirical_full_type_i"),
  paste0("<= ", alpha_main + tolerance)
)
add_check(
  "IUT partial-null full contribution error controlled",
  get_result("partial_null_iut", "empirical_full_type_i") <= alpha_main + tolerance,
  get_result("partial_null_iut", "empirical_full_type_i"),
  paste0("<= ", alpha_main + tolerance)
)
add_check(
  "maxT separate-claim family-wise error controlled under global null",
  get_result("separate_claim_global_maxT", "empirical_any_true_reject") <= alpha_main + tolerance,
  get_result("separate_claim_global_maxT", "empirical_any_true_reject"),
  paste0("<= ", alpha_main + tolerance)
)
add_check(
  "positive alternative agrees with analytic approximation",
  get_result("positive_alt_iut", "absolute_difference_from_analytic") <= 0.035,
  get_result("positive_alt_iut", "absolute_difference_from_analytic"),
  "<= 0.035"
)
add_check(
  "positive margins reduce joint power",
  get_result("positive_margin_iut", "empirical_joint_success") <
    get_result("positive_alt_iut", "empirical_joint_success"),
  paste(
    get_result("positive_margin_iut", "empirical_joint_success"),
    get_result("positive_alt_iut", "empirical_joint_success"),
    sep = " < "
  ),
  "margin scenario lower than zero-margin scenario"
)
add_check(
  "unequal variance alternative produces finite power",
  is.finite(get_result("unequal_variance_alt_iut", "empirical_joint_success")) &&
    get_result("unequal_variance_alt_iut", "empirical_joint_success") > 0,
  get_result("unequal_variance_alt_iut", "empirical_joint_success"),
  "finite and positive"
)
add_check(
  "heavy-tail null full contribution error not above fail-fast bound",
  get_result("heavy_tail_null_iut", "empirical_full_type_i") <= alpha_main + 0.020,
  get_result("heavy_tail_null_iut", "empirical_full_type_i"),
  paste0("<= ", alpha_main + 0.020)
)

checks <- do.call(rbind, check_rows)
check_path <- file.path(out_dir, paste0("participant_failfast_checks_", stamp, ".csv"))
write.csv(checks, check_path, row.names = FALSE)

print(results)
print(checks)

if (!all(checks$passed)) {
  stop("Participant-level fail-fast simulation checks failed", call. = FALSE)
}

cat("Participant-level fail-fast simulation passed: ",
    sum(checks$passed), " / ", nrow(checks), "\n", sep = "")
