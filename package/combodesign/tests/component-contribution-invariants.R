library(combodesign)

expect_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

expect_close <- function(observed, expected, tolerance, message) {
  ok <- isTRUE(all.equal(
    observed,
    expected,
    tolerance = tolerance,
    check.attributes = FALSE
  ))
  if (!ok) {
    stop(
      message,
      "\nObserved: ", paste(signif(as.vector(observed), 10), collapse = ", "),
      "\nExpected: ", paste(signif(as.vector(expected), 10), collapse = ", "),
      call. = FALSE
    )
  }
}

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("mvtnorm is required for invariant tests", call. = FALSE)
}

allocation <- c(A = 1 / 3, B = 1 / 3, AB = 1 / 3)
delta <- c(AB_minus_A = 0.55, AB_minus_B = 0.65)

fixed <- contribution_operating_characteristics(
  N = 120,
  allocation = allocation,
  delta = delta,
  sigma2 = 1,
  alpha = 0.025
)

expected_C <- rbind(
  AB_minus_A = c(A = -1, B = 0, AB = 1),
  AB_minus_B = c(A = 0, B = -1, AB = 1)
)
expect_close(
  fixed$C,
  expected_C,
  tolerance = 1e-12,
  message = "Contribution contrast matrix must encode AB - A and AB - B."
)

expected_Sigma_D <- matrix(
  c(0.05, 0.025, 0.025, 0.05),
  nrow = 2,
  dimnames = list(c("AB_minus_A", "AB_minus_B"),
                  c("AB_minus_A", "AB_minus_B"))
)
expect_close(
  fixed$Sigma_D,
  expected_Sigma_D,
  tolerance = 1e-12,
  message = "Contrast covariance is not the analytic AB - A / AB - B covariance."
)
expect_close(
  fixed$R_Delta,
  matrix(c(1, 0.5, 0.5, 1), nrow = 2),
  tolerance = 1e-12,
  message = "Induced contribution correlation should be 0.5 in the balanced independent case."
)

expect_true(
  identical(fixed$claim_regime, "full_contribution"),
  "Default claim regime should be full contribution."
)
expect_true(
  identical(fixed$adjustment_method, "iut"),
  "Default full contribution method should be IUT."
)
expect_close(
  fixed$c_alpha,
  stats::qnorm(0.975),
  tolerance = 1e-12,
  message = "IUT critical value should equal the marginal one-sided cutoff."
)

global_full_error <- as.numeric(mvtnorm::pmvnorm(
  lower = rep(fixed$c_alpha, 2),
  upper = c(Inf, Inf),
  mean = c(0, 0),
  sigma = fixed$R_Delta
))
expect_true(
  global_full_error <= 0.025 + 1e-10,
  "IUT must control the global-null full contribution error."
)

separate <- contribution_operating_characteristics(
  N = 120,
  allocation = allocation,
  delta = delta,
  sigma2 = 1,
  alpha = 0.025,
  claim_regime = "separate_component_claims",
  adjustment_method = "maxT"
)
expect_true(
  separate$c_alpha >= fixed$c_alpha,
  "maxT cutoff should not be smaller than the IUT cutoff for the same alpha."
)
separate_global_fwer <- 1 - as.numeric(mvtnorm::pmvnorm(
  lower = c(-Inf, -Inf),
  upper = rep(separate$c_alpha, 2),
  mean = c(0, 0),
  sigma = separate$R_Delta
))
expect_close(
  separate_global_fwer,
  0.025,
  tolerance = 5e-6,
  message = "maxT must control the two-contrast global-null family-wise error."
)

rho_bad <- diag(3)
dimnames(rho_bad) <- list(c("A", "B", "AB"), c("A", "B", "AB"))
rho_block <- try(
  contribution_operating_characteristics(
    N = 120,
    allocation = allocation,
    delta = delta,
    rho = rho_bad
  ),
  silent = TRUE
)
expect_true(
  inherits(rho_block, "try-error"),
  "rho should be blocked unless the general-correlated covariance model is explicit."
)

low_N <- contribution_operating_characteristics(
  N = 90,
  allocation = allocation,
  delta = delta
)
high_N <- contribution_operating_characteristics(
  N = 180,
  allocation = allocation,
  delta = delta
)
stronger_delta <- contribution_operating_characteristics(
  N = 90,
  allocation = allocation,
  delta = delta + c(AB_minus_A = 0.10, AB_minus_B = 0.10)
)
expect_true(
  high_N$joint_success_probability > low_N$joint_success_probability,
  "Joint power should increase when total sample size increases."
)
expect_true(
  stronger_delta$joint_success_probability > low_N$joint_success_probability,
  "Joint power should increase when both contribution effects increase."
)

partial_null <- contribution_operating_characteristics(
  N = 120,
  allocation = allocation,
  delta = c(AB_minus_A = 0, AB_minus_B = 0.65),
  alpha = 0.025
)
expect_true(
  unname(partial_null$marginal_power["AB_minus_A"]) <= 0.025 + 1e-10,
  "A true contribution null should not have marginal rejection probability above alpha under IUT."
)
expect_true(
  partial_null$joint_success_probability <=
    unname(partial_null$marginal_power["AB_minus_A"]) + 1e-12,
  "False contribution selection probability should not exceed the rejection probability of the true null."
)

fit <- component_contribution_design(
  delta = c(AB_minus_A = 0.65, AB_minus_B = 0.45),
  target_power = 0.65,
  min_allocation = 0.10,
  alpha = 0.025,
  N_upper = 5000
)
summary <- fit$summary
expect_true(nrow(summary) == 3, "Design comparison should return three method rows.")
for (i in seq_len(nrow(summary))) {
  allocation_row <- as.numeric(summary[i, c("allocation_A", "allocation_B", "allocation_AB")])
  expect_close(
    sum(allocation_row),
    1,
    tolerance = 1e-8,
    message = "Allocation proportions must sum to one."
  )
  expect_true(
    all(allocation_row >= 0.10 - 1e-8),
    "Allocation proportions must respect the minimum allocation constraint."
  )
}
expect_true(
  all(summary$achieved_joint_success_probability >= 0.65),
  "All reported sample-size designs must reach the target joint power."
)
expect_true(
  fit$joint_optimized$N <= fit$equal$N,
  "Joint-decision optimized allocation should not require more total N than equal allocation in this fitted comparison."
)

invalid_lower_input <- try(
  component_contribution_design(
    delta = c(AB_minus_A = 0.40, AB_minus_B = -0.05),
    target_power = 0.65,
    min_allocation = 0.10,
    alpha = 0.025,
    N_upper = 5000
  ),
  silent = TRUE
)
expect_true(
  inherits(invalid_lower_input, "try-error"),
  "A nonpositive required contribution input should fail rather than produce a misleading design."
)
