#!/usr/bin/env Rscript

# Sample-size and allocation functions for the SBR
# component-contribution design. Full conjunctive contribution with IUT is the
# default claim regime. Separate component claims must be requested explicitly.

if (!exists("operating_characteristics", mode = "function")) {
  source("analysis/contribution_operating_characteristics.R")
}

assert_positive_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(arg_name, " must be a positive finite scalar", call. = FALSE)
  }
  invisible(TRUE)
}

assert_nonnegative_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0) {
    stop(arg_name, " must be a nonnegative finite scalar", call. = FALSE)
  }
  invisible(TRUE)
}

validate_min_allocation <- function(min_allocation, arms) {
  assert_nonnegative_scalar(min_allocation, "min_allocation")
  if (min_allocation >= 1 / length(arms)) {
    stop("min_allocation must be less than 1 / number of arms", call. = FALSE)
  }
  invisible(TRUE)
}

validate_allocation <- function(allocation, arms = c("A", "B", "AB"),
                                min_allocation = 0, tol = 1e-8) {
  assert_named_numeric(allocation, "allocation")
  if (anyDuplicated(names(allocation))) {
    stop("allocation names must be unique", call. = FALSE)
  }
  if (length(arms) < 2 || anyDuplicated(arms)) {
    stop("arms must contain at least two unique arm labels", call. = FALSE)
  }
  validate_min_allocation(min_allocation, arms)

  missing_arms <- setdiff(arms, names(allocation))
  if (length(missing_arms) > 0) {
    stop("allocation is missing arms: ", paste(missing_arms, collapse = ", "),
         call. = FALSE)
  }

  allocation <- allocation[arms]
  if (any(!is.finite(allocation))) {
    stop("allocation must contain finite values", call. = FALSE)
  }
  if (any(allocation <= 0)) {
    stop("allocation must assign positive mass to every modeled arm",
         call. = FALSE)
  }
  if (any(allocation < min_allocation - tol)) {
    stop("allocation violates min_allocation", call. = FALSE)
  }
  if (abs(sum(allocation) - 1) > tol) {
    stop("allocation must sum to 1", call. = FALSE)
  }
  allocation / sum(allocation)
}

equal_allocation <- function(arms = c("A", "B", "AB")) {
  allocation <- rep(1 / length(arms), length(arms))
  names(allocation) <- arms
  allocation
}

allocation_to_integer_sample_sizes <- function(N, allocation,
                                               arms = c("A", "B", "AB")) {
  assert_positive_scalar(N, "N")
  if (abs(N - round(N)) > 1e-8) {
    stop("N must be an integer when integer_allocation = TRUE", call. = FALSE)
  }
  allocation <- validate_allocation(allocation, arms = arms)
  raw <- as.numeric(N) * allocation
  n <- floor(raw)
  remainder <- as.integer(round(N - sum(n)))
  if (remainder > 0L) {
    order_index <- order(raw - n, decreasing = TRUE)
    n[order_index[seq_len(remainder)]] <- n[order_index[seq_len(remainder)]] + 1L
  }
  names(n) <- arms
  if (any(n <= 0)) {
    stop("integer allocation produced an empty arm; increase N or change allocation",
         call. = FALSE)
  }
  n
}

allocation_to_sample_sizes <- function(N, allocation,
                                       arms = c("A", "B", "AB"),
                                       integer_allocation = FALSE) {
  assert_positive_scalar(N, "N")
  allocation <- validate_allocation(allocation, arms = arms)
  if (isTRUE(integer_allocation)) {
    return(allocation_to_integer_sample_sizes(N, allocation, arms = arms))
  }
  n <- N * allocation
  names(n) <- arms
  n
}

align_delta_to_contrast <- function(delta, C) {
  assert_named_numeric(delta, "delta")
  missing_delta <- setdiff(rownames(C), names(delta))
  if (length(missing_delta) > 0) {
    stop("delta is missing contribution contrasts: ",
         paste(missing_delta, collapse = ", "), call. = FALSE)
  }
  delta[rownames(C)]
}

design_from_allocation <- function(
    N,
    allocation,
    delta,
    rho = NULL,
    sigma2 = 1,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    integer_allocation = TRUE) {
  assert_positive_scalar(N, "N")
  assert_probability(alpha, "alpha")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  covariance_model <- match.arg(covariance_model)
  allocation <- validate_allocation(allocation, arms = arms)

  C <- make_contribution_contrast_matrix(
    arms = arms,
    contribution_map = contribution_map
  )
  delta <- align_delta_to_contrast(delta, C)
  margin <- align_margin(margin, names(delta))
  n <- allocation_to_sample_sizes(
    N,
    allocation,
    arms = arms,
    integer_allocation = integer_allocation
  )
  Sigma_y <- make_arm_mean_covariance(
    n = n,
    sigma2 = sigma2,
    rho = rho,
    covariance_model = covariance_model
  )
  covariance <- contrast_covariance(C, Sigma_y)
  operating <- operating_characteristics(
    delta = delta,
    Sigma_D = covariance$Sigma_D,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method
  )

  list(
    N = N,
    allocation = allocation,
    n = n,
    arms = arms,
    C = C,
    delta = delta,
    margin = margin,
    adjusted_delta = operating$adjusted_delta,
    Sigma_y = Sigma_y,
    Sigma_D = covariance$Sigma_D,
    R_Delta = covariance$R_Delta,
    alpha = alpha,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    c_alpha = operating$c_alpha,
    lambda = operating$lambda,
    joint_success_probability = operating$joint_success,
    marginal_power = operating$marginal_power
  )
}

joint_success_probability <- function(
    N,
    allocation,
    delta,
    rho = NULL,
    sigma2 = 1,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    integer_allocation = TRUE) {
  design_from_allocation(
    N = N,
    allocation = allocation,
    delta = delta,
    rho = rho,
    sigma2 = sigma2,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    integer_allocation = integer_allocation
  )$joint_success_probability
}

sample_size_joint_success <- function(
    delta,
    allocation,
    rho = NULL,
    sigma2 = 1,
    target_power = 0.80,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    N_lower = NULL,
    N_upper = 1e6,
    integer_total = TRUE,
    tol = 1e-7,
    max_expansions = 40L) {
  assert_probability(target_power, "target_power")
  assert_positive_scalar(N_upper, "N_upper")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  covariance_model <- match.arg(covariance_model)
  allocation <- validate_allocation(allocation, arms = arms)
  C <- make_contribution_contrast_matrix(
    arms = arms,
    contribution_map = contribution_map
  )
  delta <- align_delta_to_contrast(delta, C)
  margin <- align_margin(margin, names(delta))
  if (any(delta - margin <= 0)) {
    stop("All planned contribution effects must exceed their margins",
         call. = FALSE)
  }

  if (is.null(N_lower)) {
    N_lower <- ceiling(1 / min(allocation))
  }
  assert_positive_scalar(N_lower, "N_lower")
  if (N_lower >= N_upper) {
    stop("N_lower must be smaller than N_upper", call. = FALSE)
  }

  power_at <- function(N) {
    joint_success_probability(
      N = N,
      allocation = allocation,
      delta = delta,
      rho = rho,
      sigma2 = sigma2,
      alpha = alpha,
      margin = margin,
      claim_regime = claim_regime,
      adjustment_method = adjustment_method,
      covariance_model = covariance_model,
      arms = arms,
      contribution_map = contribution_map,
      integer_allocation = FALSE
    )
  }

  lower_power <- power_at(N_lower)
  if (lower_power >= target_power) {
    N_continuous <- N_lower
  } else {
    upper <- min(max(N_lower * 2, N_lower + 1), N_upper)
    upper_power <- power_at(upper)
    expansions <- 0L
    while (upper_power < target_power && upper < N_upper &&
           expansions < max_expansions) {
      upper <- min(upper * 2, N_upper)
      upper_power <- power_at(upper)
      expansions <- expansions + 1L
    }
    if (upper_power < target_power) {
      stop(
        "target_power is not reached before N_upper; check effect sizes, ",
        "allocation, or target_power",
        call. = FALSE
      )
    }
    root <- stats::uniroot(
      function(N) power_at(N) - target_power,
      interval = c(N_lower, upper),
      tol = tol
    )
    N_continuous <- root$root
  }

  N_report <- if (isTRUE(integer_total)) ceiling(N_continuous) else N_continuous
  achieved <- design_from_allocation(
    N = N_report,
    allocation = allocation,
    delta = delta,
    rho = rho,
    sigma2 = sigma2,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    integer_allocation = isTRUE(integer_total)
  )
  if (isTRUE(integer_total)) {
    while (achieved$joint_success_probability < target_power) {
      N_report <- N_report + 1L
      achieved <- design_from_allocation(
        N = N_report,
        allocation = allocation,
        delta = delta,
        rho = rho,
        sigma2 = sigma2,
        alpha = alpha,
        margin = margin,
        claim_regime = claim_regime,
        adjustment_method = adjustment_method,
        covariance_model = covariance_model,
        arms = arms,
        contribution_map = contribution_map,
        integer_allocation = TRUE
      )
    }
  }

  list(
    target_power = target_power,
    N_continuous = N_continuous,
    N = N_report,
    achieved_joint_success_probability =
      achieved$joint_success_probability,
    allocation = achieved$allocation,
    n = achieved$n,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    marginal_power = achieved$marginal_power,
    c_alpha = achieved$c_alpha,
    design = achieved
  )
}

allocation_to_theta <- function(allocation, arms = c("A", "B", "AB"),
                                min_allocation = 0) {
  allocation <- validate_allocation(
    allocation,
    arms = arms,
    min_allocation = min_allocation
  )
  validate_min_allocation(min_allocation, arms)

  residual <- 1 - length(arms) * min_allocation
  raw <- (allocation - min_allocation) / residual
  if (any(raw <= 0)) {
    stop("allocation must be strictly above min_allocation for optimization",
         call. = FALSE)
  }
  log(raw[-length(raw)] / raw[length(raw)])
}

theta_to_allocation <- function(theta, arms = c("A", "B", "AB"),
                                min_allocation = 0) {
  if (!is.numeric(theta) || length(theta) != length(arms) - 1 ||
      any(!is.finite(theta))) {
    stop("theta must be a finite numeric vector with length number of arms - 1",
         call. = FALSE)
  }
  validate_min_allocation(min_allocation, arms)

  eta <- c(theta, 0)
  exp_eta <- exp(eta - max(eta))
  raw <- exp_eta / sum(exp_eta)
  allocation <- min_allocation +
    (1 - length(arms) * min_allocation) * raw
  names(allocation) <- arms
  allocation
}

default_start_allocations <- function(arms = c("A", "B", "AB"),
                                      min_allocation = 0.02) {
  validate_min_allocation(min_allocation, arms)
  m <- length(arms)
  starts <- list(equal_allocation(arms))
  high <- min(0.70, 1 - (m - 1) * min_allocation - 1e-4)
  if (high > 1 / m) {
    for (i in seq_len(m)) {
      p <- rep((1 - high) / (m - 1), m)
      p[i] <- high
      names(p) <- arms
      starts[[length(starts) + 1L]] <- p
    }
  }
  starts
}

optimize_from_starts <- function(objective, starts, arms, min_allocation,
                                 method = "BFGS", control = list(maxit = 300)) {
  results <- lapply(starts, function(start) {
    theta_start <- allocation_to_theta(
      start,
      arms = arms,
      min_allocation = min_allocation
    )
    fit <- try(
      stats::optim(
        par = theta_start,
        fn = objective,
        method = method,
        control = control
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      return(NULL)
    }
    fit
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    stop("allocation optimization failed from all starts", call. = FALSE)
  }
  values <- vapply(results, function(x) x$value, numeric(1))
  finite_values <- is.finite(values)
  if (!any(finite_values)) {
    stop("allocation optimization returned no finite objective values",
         call. = FALSE)
  }
  results <- results[finite_values]
  values <- values[finite_values]
  results[[which.min(values)]]
}

optimize_allocation_for_joint_success <- function(
    N,
    delta,
    rho = NULL,
    sigma2 = 1,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    min_allocation = 0.02,
    start_allocations = NULL,
    control = list(maxit = 300)) {
  assert_positive_scalar(N, "N")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  covariance_model <- match.arg(covariance_model)
  validate_min_allocation(min_allocation, arms)
  if (is.null(start_allocations)) {
    start_allocations <- default_start_allocations(arms, min_allocation)
  }

  objective <- function(theta) {
    allocation <- theta_to_allocation(
      theta,
      arms = arms,
      min_allocation = min_allocation
    )
    value <- try(
      joint_success_probability(
        N = N,
        allocation = allocation,
        delta = delta,
        rho = rho,
        sigma2 = sigma2,
        alpha = alpha,
        margin = margin,
        claim_regime = claim_regime,
        adjustment_method = adjustment_method,
        covariance_model = covariance_model,
        arms = arms,
        contribution_map = contribution_map
      ),
      silent = TRUE
    )
    if (inherits(value, "try-error") || !is.finite(value)) {
      return(Inf)
    }
    -value
  }

  fit <- optimize_from_starts(
    objective = objective,
    starts = start_allocations,
    arms = arms,
    min_allocation = min_allocation,
    control = control
  )
  allocation <- theta_to_allocation(
    fit$par,
    arms = arms,
    min_allocation = min_allocation
  )
  design <- design_from_allocation(
    N = N,
    allocation = allocation,
    delta = delta,
    rho = rho,
    sigma2 = sigma2,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map
  )

  list(
    objective = "maximize joint power at fixed N",
    convergence = fit$convergence,
    value = -fit$value,
    allocation = allocation,
    N = N,
    joint_success_probability = design$joint_success_probability,
    marginal_power = design$marginal_power,
    design = design,
    optim = fit
  )
}

optimize_allocation_for_min_marginal_power <- function(
    N,
    delta,
    rho = NULL,
    sigma2 = 1,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    min_allocation = 0.02,
    start_allocations = NULL,
    control = list(maxit = 300)) {
  assert_positive_scalar(N, "N")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  covariance_model <- match.arg(covariance_model)
  validate_min_allocation(min_allocation, arms)
  if (is.null(start_allocations)) {
    start_allocations <- default_start_allocations(arms, min_allocation)
  }

  objective <- function(theta) {
    allocation <- theta_to_allocation(
      theta,
      arms = arms,
      min_allocation = min_allocation
    )
    design <- try(
      design_from_allocation(
        N = N,
        allocation = allocation,
        delta = delta,
        rho = rho,
        sigma2 = sigma2,
        alpha = alpha,
        margin = margin,
        claim_regime = claim_regime,
        adjustment_method = adjustment_method,
        covariance_model = covariance_model,
        arms = arms,
        contribution_map = contribution_map
      ),
      silent = TRUE
    )
    if (inherits(design, "try-error")) {
      return(Inf)
    }
    -min(design$marginal_power)
  }

  fit <- optimize_from_starts(
    objective = objective,
    starts = start_allocations,
    arms = arms,
    min_allocation = min_allocation,
    control = control
  )
  allocation <- theta_to_allocation(
    fit$par,
    arms = arms,
    min_allocation = min_allocation
  )
  design <- design_from_allocation(
    N = N,
    allocation = allocation,
    delta = delta,
    rho = rho,
    sigma2 = sigma2,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map
  )

  list(
    objective = paste(
      "maximize weaker marginal contribution power at fixed N;",
      "diagnostic comparator, not proposed design target"
    ),
    convergence = fit$convergence,
    value = -fit$value,
    allocation = allocation,
    N = N,
    joint_success_probability = design$joint_success_probability,
    min_marginal_power = min(design$marginal_power),
    marginal_power = design$marginal_power,
    design = design,
    optim = fit
  )
}

optimize_allocation_for_sample_size <- function(
    delta,
    rho = NULL,
    sigma2 = 1,
    target_power = 0.80,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    min_allocation = 0.02,
    start_allocations = NULL,
    N_upper = 1e6,
    control = list(maxit = 160)) {
  assert_probability(target_power, "target_power")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  covariance_model <- match.arg(covariance_model)
  validate_min_allocation(min_allocation, arms)
  if (is.null(start_allocations)) {
    start_allocations <- default_start_allocations(arms, min_allocation)
  }

  objective <- function(theta) {
    allocation <- theta_to_allocation(
      theta,
      arms = arms,
      min_allocation = min_allocation
    )
    size <- try(
      sample_size_joint_success(
        delta = delta,
        allocation = allocation,
        rho = rho,
        sigma2 = sigma2,
        target_power = target_power,
        alpha = alpha,
        margin = margin,
        claim_regime = claim_regime,
        adjustment_method = adjustment_method,
        covariance_model = covariance_model,
        arms = arms,
        contribution_map = contribution_map,
        N_upper = N_upper,
        integer_total = FALSE
      ),
      silent = TRUE
    )
    if (inherits(size, "try-error") || !is.finite(size$N_continuous)) {
      return(Inf)
    }
    size$N_continuous
  }

  fit <- optimize_from_starts(
    objective = objective,
    starts = start_allocations,
    arms = arms,
    min_allocation = min_allocation,
    method = "Nelder-Mead",
    control = control
  )
  allocation <- theta_to_allocation(
    fit$par,
    arms = arms,
    min_allocation = min_allocation
  )
  size <- sample_size_joint_success(
    delta = delta,
    allocation = allocation,
    rho = rho,
    sigma2 = sigma2,
    target_power = target_power,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    N_upper = N_upper,
    integer_total = TRUE
  )

  list(
    objective = "minimize total N for target joint power",
    convergence = fit$convergence,
    allocation = allocation,
    target_power = target_power,
    N_continuous = fit$value,
    N = size$N,
    achieved_joint_success_probability =
      size$achieved_joint_success_probability,
    marginal_power = size$marginal_power,
    sample_size = size,
    optim = fit
  )
}

allocation_summary_row <- function(label, size_fit) {
  allocation <- size_fit$allocation
  marginal <- size_fit$marginal_power
  data.frame(
    method = label,
    target_joint_success_probability = size_fit$target_power,
    N_continuous = size_fit$N_continuous,
    N = size_fit$N,
    achieved_joint_success_probability =
      size_fit$achieved_joint_success_probability,
    allocation_A = unname(allocation["A"]),
    allocation_B = unname(allocation["B"]),
    allocation_AB = unname(allocation["AB"]),
    marginal_power_AB_minus_A = unname(marginal["AB_minus_A"]),
    marginal_power_AB_minus_B = unname(marginal["AB_minus_B"]),
    stringsAsFactors = FALSE
  )
}

compare_allocation_strategies <- function(
    delta,
    rho = NULL,
    sigma2 = 1,
    target_power = 0.80,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL,
    min_allocation = 0.02,
    N_upper = 1e6) {
  if (!identical(arms, c("A", "B", "AB"))) {
    stop("compare_allocation_strategies currently reports the A, B, AB case",
         call. = FALSE)
  }

  equal <- sample_size_joint_success(
    delta = delta,
    allocation = equal_allocation(arms),
    rho = rho,
    sigma2 = sigma2,
    target_power = target_power,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    N_upper = N_upper
  )
  joint <- optimize_allocation_for_sample_size(
    delta = delta,
    rho = rho,
    sigma2 = sigma2,
    target_power = target_power,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    min_allocation = min_allocation,
    N_upper = N_upper
  )
  marginal_comparator <- optimize_allocation_for_min_marginal_power(
    N = equal$N,
    delta = delta,
    rho = rho,
    sigma2 = sigma2,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    min_allocation = min_allocation
  )
  marginal_size <- sample_size_joint_success(
    delta = delta,
    allocation = marginal_comparator$allocation,
    rho = rho,
    sigma2 = sigma2,
    target_power = target_power,
    alpha = alpha,
    margin = margin,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    covariance_model = covariance_model,
    arms = arms,
    contribution_map = contribution_map,
    N_upper = N_upper
  )

  summary <- rbind(
    allocation_summary_row("equal allocation", equal),
    allocation_summary_row("joint-decision optimized allocation", joint),
    allocation_summary_row(
      "minimum marginal-power comparator",
      marginal_size
    )
  )
  rownames(summary) <- NULL

  list(
    summary = summary,
    equal = equal,
    joint_optimized = joint,
    marginal_power_comparator = marginal_comparator,
    marginal_power_comparator_sample_size = marginal_size
  )
}

assurance_probability <- function(
    N,
    allocation,
    delta_draws,
    target_power = 0.80,
    rho = NULL,
    sigma2 = 1,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    covariance_model = c("independent_arms", "general_correlated"),
    arms = c("A", "B", "AB"),
    contribution_map = NULL) {
  assert_positive_scalar(N, "N")
  assert_probability(target_power, "target_power")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  covariance_model <- match.arg(covariance_model)
  if (is.data.frame(delta_draws)) {
    delta_draws <- as.matrix(delta_draws)
  }
  if (!is.matrix(delta_draws) || !is.numeric(delta_draws) ||
      is.null(colnames(delta_draws))) {
    stop("delta_draws must be a numeric matrix or data frame with column names",
         call. = FALSE)
  }

  powers <- apply(delta_draws, 1L, function(delta_row) {
    delta <- as.numeric(delta_row)
    names(delta) <- colnames(delta_draws)
    design_from_allocation(
      N = N,
      allocation = allocation,
      delta = delta,
      rho = rho,
      sigma2 = sigma2,
      alpha = alpha,
      margin = margin,
      claim_regime = claim_regime,
      adjustment_method = adjustment_method,
      covariance_model = covariance_model,
      arms = arms,
      contribution_map = contribution_map,
      integer_allocation = TRUE
    )$joint_success_probability
  })

  list(
    N = N,
    allocation = validate_allocation(allocation, arms = arms),
    target_power = target_power,
    assurance = mean(powers >= target_power),
    power_by_draw = powers,
    min_power = min(powers),
    lower_quartile_power = unname(stats::quantile(powers, 0.25)),
    median_power = unname(stats::quantile(powers, 0.50)),
    lower_05_power = unname(stats::quantile(powers, 0.05))
  )
}
