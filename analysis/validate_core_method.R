#!/usr/bin/env Rscript

# Fail-fast validation for the Section 3 core code.

source("analysis/contribution_operating_characteristics.R")
source("analysis/contribution_sample_size_allocation.R")

stamp <- "2026-06-19"
result_dir <- file.path("analysis", "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, digits = 10) {
  paste(signif(as.numeric(x), digits), collapse = "; ")
}

check_rows <- list()
add_check <- function(check, passed, observed, expected = "") {
  check_rows[[length(check_rows) + 1L]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

expect_close <- function(x, y, tol = 1e-8) {
  isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
}

alpha <- 0.025
arms <- c("A", "B", "AB")
allocation <- c(A = 0.30, B = 0.30, AB = 0.40)
n <- c(A = 40, B = 50, AB = 60)
sigma2_vec <- c(A = 1.20, B = 1.50, AB = 2.00)
C0 <- make_contribution_contrast_matrix(arms)

Sigma_y <- make_arm_mean_covariance(
  n = n,
  sigma2 = sigma2_vec,
  covariance_model = "independent_arms"
)
covariance <- contrast_covariance(C0, Sigma_y)
manual_Sigma_D <- matrix(
  c(
    sigma2_vec["AB"] / n["AB"] + sigma2_vec["A"] / n["A"],
    sigma2_vec["AB"] / n["AB"],
    sigma2_vec["AB"] / n["AB"],
    sigma2_vec["AB"] / n["AB"] + sigma2_vec["B"] / n["B"]
  ),
  nrow = 2,
  dimnames = list(c("AB_minus_A", "AB_minus_B"),
                  c("AB_minus_A", "AB_minus_B"))
)
add_check(
  "independent-arm contrast covariance matches manual formula",
  expect_close(covariance$Sigma_D, manual_Sigma_D, tol = 1e-12),
  fmt(covariance$Sigma_D),
  fmt(manual_Sigma_D)
)

rho_test <- diag(3)
dimnames(rho_test) <- list(arms, arms)
rho_block <- try(
  make_arm_mean_covariance(n = n, sigma2 = 1, rho = rho_test),
  silent = TRUE
)
add_check(
  "rho is blocked under independent-arm default",
  inherits(rho_block, "try-error"),
  class(rho_block)[1],
  "try-error"
)

Sigma_y_general <- make_arm_mean_covariance(
  n = n,
  sigma2 = 1,
  rho = rho_test,
  covariance_model = "general_correlated"
)
add_check(
  "rho is allowed only in explicit general-correlated extension",
  expect_close(Sigma_y_general, diag(1 / n, 3), tol = 1e-12),
  fmt(diag(Sigma_y_general)),
  fmt(1 / n)
)

iut_cutoff <- critical_value_for_regime(
  covariance$R_Delta,
  alpha = alpha,
  claim_regime = "full_contribution"
)
maxT_cutoff <- critical_value_for_regime(
  covariance$R_Delta,
  alpha = alpha,
  claim_regime = "separate_component_claims",
  adjustment_method = "maxT"
)
add_check(
  "IUT cutoff equals marginal one-sided cutoff",
  expect_close(iut_cutoff, stats::qnorm(1 - alpha), tol = 1e-12),
  fmt(iut_cutoff),
  fmt(stats::qnorm(1 - alpha))
)
add_check(
  "maxT cutoff is not smaller than IUT cutoff",
  maxT_cutoff >= iut_cutoff - 1e-10,
  fmt(maxT_cutoff),
  paste0(">= ", fmt(iut_cutoff))
)

delta_alt <- c(AB_minus_A = 0.45, AB_minus_B = 0.55)
margin <- c(AB_minus_A = 0.05, AB_minus_B = 0.05)
full_oc <- operating_characteristics(
  delta = delta_alt,
  Sigma_D = covariance$Sigma_D,
  alpha = alpha,
  margin = margin,
  claim_regime = "full_contribution"
)
separate_oc <- operating_characteristics(
  delta = delta_alt,
  Sigma_D = covariance$Sigma_D,
  alpha = alpha,
  margin = margin,
  claim_regime = "separate_component_claims",
  adjustment_method = "maxT"
)
add_check(
  "full contribution default uses IUT",
  identical(full_oc$adjustment_method, "iut"),
  full_oc$adjustment_method,
  "iut"
)
add_check(
  "separate component claim mode uses maxT when requested",
  identical(separate_oc$adjustment_method, "maxT"),
  separate_oc$adjustment_method,
  "maxT"
)
add_check(
  "full contribution joint power is at least maxT all-claims power for same alpha",
  full_oc$joint_success >= separate_oc$joint_success - 1e-12,
  fmt(c(full_oc$joint_success, separate_oc$joint_success)),
  "IUT joint power >= maxT all-claims power"
)

full_null_states <- evaluate_null_states(
  delta_alt = delta_alt,
  Sigma_D = covariance$Sigma_D,
  alpha = alpha,
  margin = margin,
  claim_regime = "full_contribution"
)
type_i_rows <- full_null_states[full_null_states$true_null_count > 0, ]
add_check(
  "IUT controls full contribution type I error in checked null states",
  all(type_i_rows$full_contribution_type_i_error <= alpha + 1e-10),
  fmt(type_i_rows$full_contribution_type_i_error),
  paste0("<= ", alpha)
)

separate_null_states <- evaluate_null_states(
  delta_alt = delta_alt,
  Sigma_D = covariance$Sigma_D,
  alpha = alpha,
  margin = margin,
  claim_regime = "separate_component_claims",
  adjustment_method = "maxT"
)
sep_rows <- separate_null_states[separate_null_states$true_null_count > 0, ]
add_check(
  "maxT controls separate-claim family-wise error in checked null states",
  all(sep_rows$familywise_error_any_true_null <= alpha + 1e-10),
  fmt(sep_rows$familywise_error_any_true_null),
  paste0("<= ", alpha)
)

balanced_allocation <- c(A = 1 / 3, B = 1 / 3, AB = 1 / 3)
size_zero_margin <- sample_size_joint_success(
  delta = c(AB_minus_A = 0.50, AB_minus_B = 0.55),
  allocation = balanced_allocation,
  target_power = 0.60,
  alpha = alpha,
  margin = c(AB_minus_A = 0, AB_minus_B = 0),
  N_upper = 5000
)
size_positive_margin <- sample_size_joint_success(
  delta = c(AB_minus_A = 0.50, AB_minus_B = 0.55),
  allocation = balanced_allocation,
  target_power = 0.60,
  alpha = alpha,
  margin = c(AB_minus_A = 0.08, AB_minus_B = 0.08),
  N_upper = 5000
)
add_check(
  "positive contribution margins increase required sample size",
  size_positive_margin$N > size_zero_margin$N,
  paste(size_zero_margin$N, size_positive_margin$N, sep = " -> "),
  "positive-margin N > zero-margin N"
)

bad_margin_input <- try(
  sample_size_joint_success(
    delta = c(AB_minus_A = 0.04, AB_minus_B = 0.55),
    allocation = balanced_allocation,
    target_power = 0.60,
    alpha = alpha,
    margin = c(AB_minus_A = 0.05, AB_minus_B = 0.05),
    N_upper = 5000
  ),
  silent = TRUE
)
add_check(
  "nonpositive margin-adjusted contribution input fails",
  inherits(bad_margin_input, "try-error"),
  class(bad_margin_input)[1],
  "try-error"
)

integer_n <- allocation_to_integer_sample_sizes(
  N = 101,
  allocation = c(A = 0.20, B = 0.30, AB = 0.50)
)
add_check(
  "integer allocation sums to total N",
  sum(integer_n) == 101,
  sum(integer_n),
  101
)

draws <- rbind(
  c(AB_minus_A = 0.55, AB_minus_B = 0.60),
  c(AB_minus_A = 0.45, AB_minus_B = 0.50),
  c(AB_minus_A = 0.25, AB_minus_B = 0.30)
)
assurance <- assurance_probability(
  N = 180,
  allocation = balanced_allocation,
  delta_draws = draws,
  target_power = 0.30,
  alpha = alpha
)
manual_assurance <- mean(assurance$power_by_draw >= 0.30)
add_check(
  "assurance calculation matches draw-level powers",
  expect_close(assurance$assurance, manual_assurance, tol = 1e-12),
  fmt(assurance$assurance),
  fmt(manual_assurance)
)

checks <- do.call(rbind, check_rows)
check_path <- file.path(result_dir, paste0("core_redesign_validation_checks_", stamp, ".csv"))
write.csv(checks, check_path, row.names = FALSE)

if (!all(checks$passed)) {
  print(checks)
  stop("Core redesign validation failed", call. = FALSE)
}

print(checks)
cat("Core redesign validation passed: ", sum(checks$passed), " / ",
    nrow(checks), "\n", sep = "")
