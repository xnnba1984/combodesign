#!/usr/bin/env Rscript

# Operating-characteristic functions for the SBR component-contribution
# framework. The default claim regime is full
# conjunctive contribution, analyzed by intersection union testing. maxT and
# other family-wise procedures are retained for optional separate component
# claims only.

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is not installed: ", package, call. = FALSE)
  }
}

require_namespace("mvtnorm")

assert_named_numeric <- function(x, arg_name) {
  if (!is.numeric(x) || is.null(names(x)) || any(names(x) == "")) {
    stop(arg_name, " must be a named numeric vector", call. = FALSE)
  }
  invisible(TRUE)
}

assert_probability <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0 || x >= 1) {
    stop(arg_name, " must be a single number in (0, 1)", call. = FALSE)
  }
  invisible(TRUE)
}

validate_claim_regime <- function(claim_regime) {
  match.arg(claim_regime, c("full_contribution", "separate_component_claims"))
}

validate_adjustment_method <- function(adjustment_method) {
  match.arg(adjustment_method, c("iut", "maxT", "bonferroni", "holm"))
}

resolve_adjustment_method <- function(claim_regime, adjustment_method = NULL) {
  claim_regime <- validate_claim_regime(claim_regime)
  if (is.null(adjustment_method)) {
    return(if (claim_regime == "full_contribution") "iut" else "maxT")
  }
  adjustment_method <- validate_adjustment_method(adjustment_method)
  if (claim_regime == "full_contribution" && adjustment_method != "iut") {
    stop(
      "Full contribution uses IUT. Use claim_regime = ",
      "'separate_component_claims' for maxT, Holm, or Bonferroni.",
      call. = FALSE
    )
  }
  adjustment_method
}

align_margin <- function(margin, contrast_names) {
  if (is.null(margin)) {
    out <- rep(0, length(contrast_names))
    names(out) <- contrast_names
    return(out)
  }
  if (!is.numeric(margin) || any(!is.finite(margin))) {
    stop("margin must be a finite numeric vector", call. = FALSE)
  }
  if (length(margin) == 1L) {
    out <- rep(margin, length(contrast_names))
    names(out) <- contrast_names
  } else {
    if (is.null(names(margin)) || any(names(margin) == "")) {
      stop("margin must be named unless it has length one", call. = FALSE)
    }
    missing_margin <- setdiff(contrast_names, names(margin))
    if (length(missing_margin) > 0) {
      stop("margin is missing contribution contrasts: ",
           paste(missing_margin, collapse = ", "), call. = FALSE)
    }
    out <- margin[contrast_names]
  }
  if (any(out < 0)) {
    stop("contribution margins must be nonnegative", call. = FALSE)
  }
  out
}

assert_square_named_matrix <- function(x, arg_name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
    stop(arg_name, " must be a numeric square matrix", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(arg_name, " must have row and column names", call. = FALSE)
  }
  if (!identical(rownames(x), colnames(x))) {
    stop(arg_name, " row names and column names must be identical", call. = FALSE)
  }
  invisible(TRUE)
}

assert_psd <- function(x, arg_name, tol = 1e-8) {
  eig <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < -tol) {
    stop(arg_name, " must be positive semi-definite; minimum eigenvalue is ",
         signif(min(eig), 4), call. = FALSE)
  }
  invisible(TRUE)
}

as_correlation_matrix <- function(R, arg_name = "R") {
  assert_square_named_matrix(R, arg_name)
  if (max(abs(R - t(R))) > 1e-10) {
    stop(arg_name, " must be symmetric", call. = FALSE)
  }
  if (max(abs(diag(R) - 1)) > 1e-10) {
    stop(arg_name, " must have unit diagonal", call. = FALSE)
  }
  if (any(abs(R) > 1 + 1e-10)) {
    stop(arg_name, " entries must be between -1 and 1", call. = FALSE)
  }
  assert_psd(R, arg_name)
  (R + t(R)) / 2
}

make_default_contribution_map <- function(arms) {
  if (!all(c("A", "B", "AB") %in% arms)) {
    stop("Default contribution map requires arms A, B, and AB", call. = FALSE)
  }
  data.frame(
    combination = "AB",
    base = "A",
    component = "B",
    stringsAsFactors = FALSE
  )
}

make_contribution_contrast_matrix <- function(
    arms = c("A", "B", "AB"),
    contribution_map = NULL) {
  if (is.null(contribution_map)) {
    contribution_map <- make_default_contribution_map(arms)
  }

  required_columns <- c("combination", "base", "component")
  if (!all(required_columns %in% names(contribution_map))) {
    stop("contribution_map must contain columns: ",
         paste(required_columns, collapse = ", "), call. = FALSE)
  }

  all_map_arms <- unique(unlist(contribution_map[required_columns], use.names = FALSE))
  missing_arms <- setdiff(all_map_arms, arms)
  if (length(missing_arms) > 0) {
    stop("contribution_map contains arms not present in arms: ",
         paste(missing_arms, collapse = ", "), call. = FALSE)
  }

  row_labels <- as.vector(rbind(
    paste0(contribution_map$combination, "_minus_", contribution_map$base),
    paste0(contribution_map$combination, "_minus_", contribution_map$component)
  ))
  if (anyDuplicated(row_labels)) {
    stop("Contribution contrast row labels are duplicated", call. = FALSE)
  }

  C <- matrix(
    0,
    nrow = length(row_labels),
    ncol = length(arms),
    dimnames = list(row_labels, arms)
  )

  row_index <- 1L
  for (i in seq_len(nrow(contribution_map))) {
    combination <- contribution_map$combination[i]
    base <- contribution_map$base[i]
    component <- contribution_map$component[i]

    C[row_index, combination] <- 1
    C[row_index, base] <- -1
    row_index <- row_index + 1L

    C[row_index, combination] <- 1
    C[row_index, component] <- -1
    row_index <- row_index + 1L
  }

  C
}

make_arm_mean_covariance <- function(
    n,
    sigma2 = 1,
    rho = NULL,
    covariance_model = c("independent_arms", "general_correlated")) {
  assert_named_numeric(n, "n")
  if (any(!is.finite(n)) || any(n <= 0)) {
    stop("n must contain positive finite sample sizes", call. = FALSE)
  }

  arms <- names(n)
  covariance_model <- match.arg(covariance_model)
  if (length(sigma2) == 1L) {
    if (!is.finite(sigma2) || sigma2 <= 0) {
      stop("sigma2 must be positive and finite", call. = FALSE)
    }
    sigma2_vec <- rep(sigma2, length(arms))
    names(sigma2_vec) <- arms
  } else {
    assert_named_numeric(sigma2, "sigma2")
    missing_sigma2 <- setdiff(arms, names(sigma2))
    if (length(missing_sigma2) > 0) {
      stop("sigma2 is missing arms: ", paste(missing_sigma2, collapse = ", "),
           call. = FALSE)
    }
    sigma2_vec <- sigma2[arms]
    if (any(!is.finite(sigma2_vec)) || any(sigma2_vec <= 0)) {
      stop("sigma2 must contain positive finite variances", call. = FALSE)
    }
  }

  if (covariance_model == "independent_arms") {
    if (!is.null(rho)) {
      stop(
        "rho is not allowed when covariance_model = 'independent_arms'. ",
        "Use covariance_model = 'general_correlated' only with a justified ",
        "correlated-estimator design.",
        call. = FALSE
      )
    }
    Sigma_y <- diag(sigma2_vec / n, length(arms))
    dimnames(Sigma_y) <- list(arms, arms)
    assert_psd(Sigma_y, "Sigma_y")
    return(Sigma_y)
  }

  if (is.null(rho)) {
    rho <- diag(1, length(arms))
    dimnames(rho) <- list(arms, arms)
  }
  rho <- as_correlation_matrix(rho, "rho")
  if (!all(arms %in% rownames(rho))) {
    stop("rho must contain all arms named in n", call. = FALSE)
  }
  rho <- rho[arms, arms, drop = FALSE]

  scale <- sqrt(sigma2_vec / n)
  Sigma_y <- diag(scale, length(scale)) %*% rho %*% diag(scale, length(scale))
  dimnames(Sigma_y) <- list(arms, arms)
  assert_psd(Sigma_y, "Sigma_y")
  Sigma_y
}

contrast_covariance <- function(C, Sigma_y) {
  if (!is.matrix(C) || !is.numeric(C) || is.null(rownames(C)) || is.null(colnames(C))) {
    stop("C must be a numeric matrix with row and column names", call. = FALSE)
  }
  assert_square_named_matrix(Sigma_y, "Sigma_y")
  missing_arms <- setdiff(colnames(C), rownames(Sigma_y))
  if (length(missing_arms) > 0) {
    stop("Sigma_y is missing contrast columns: ",
         paste(missing_arms, collapse = ", "), call. = FALSE)
  }
  Sigma_y <- Sigma_y[colnames(C), colnames(C), drop = FALSE]
  Sigma_D <- C %*% Sigma_y %*% t(C)
  Sigma_D <- (Sigma_D + t(Sigma_D)) / 2
  if (any(diag(Sigma_D) <= 0)) {
    stop("All contribution contrast variances must be positive", call. = FALSE)
  }
  assert_psd(Sigma_D, "Sigma_D")
  R_Delta <- stats::cov2cor(Sigma_D)
  list(
    C = C,
    Sigma_D = Sigma_D,
    R_Delta = R_Delta,
    variance = diag(Sigma_D),
    sd = sqrt(diag(Sigma_D))
  )
}

pmvnorm_prob <- function(lower, upper, mean, sigma) {
  as.numeric(mvtnorm::pmvnorm(
    lower = lower,
    upper = upper,
    mean = mean,
    sigma = sigma
  ))
}

critical_value_iut <- function(alpha = 0.025) {
  assert_probability(alpha, "alpha")
  stats::qnorm(1 - alpha)
}

critical_value_bonferroni <- function(m, alpha = 0.025) {
  assert_probability(alpha, "alpha")
  if (!is.numeric(m) || length(m) != 1 || !is.finite(m) || m < 1) {
    stop("m must be a positive finite scalar", call. = FALSE)
  }
  stats::qnorm(1 - alpha / m)
}

critical_value_maxT_one_sided <- function(R, alpha = 0.025, tol = 1e-10) {
  assert_probability(alpha, "alpha")
  R <- as_correlation_matrix(R, "R")
  m <- nrow(R)

  if (m == 1) {
    return(stats::qnorm(1 - alpha))
  }
  if (max(abs(R - matrix(1, m, m))) < 1e-10) {
    return(stats::qnorm(1 - alpha))
  }

  target <- 1 - alpha
  f <- function(c_value) {
    pmvnorm_prob(
      lower = rep(-Inf, m),
      upper = rep(c_value, m),
      mean = rep(0, m),
      sigma = R
    ) - target
  }

  lower <- -10
  upper <- 10
  if (f(lower) > 0 || f(upper) < 0) {
    stop("Could not bracket one-sided critical value", call. = FALSE)
  }
  stats::uniroot(f, interval = c(lower, upper), tol = tol)$root
}

critical_value_one_sided <- function(...) {
  stop(
    "critical_value_one_sided is retired because it is ambiguous. ",
    "Use critical_value_iut() for the primary full-contribution IUT rule ",
    "or critical_value_maxT_one_sided() for optional separate component claims.",
    call. = FALSE
  )
}

critical_value_for_regime <- function(
    R,
    alpha = 0.025,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL) {
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  R <- as_correlation_matrix(R, "R")
  if (adjustment_method == "iut") {
    return(critical_value_iut(alpha))
  }
  if (adjustment_method == "maxT") {
    return(critical_value_maxT_one_sided(R, alpha = alpha))
  }
  if (adjustment_method == "bonferroni") {
    return(critical_value_bonferroni(nrow(R), alpha = alpha))
  }
  stop("Holm does not have a single common analytic cutoff; use p-value ",
       "adjustment in the separate-claim simulation layer.", call. = FALSE)
}

operating_characteristics <- function(
    delta,
    Sigma_D,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    c_alpha = NULL) {
  assert_named_numeric(delta, "delta")
  assert_probability(alpha, "alpha")
  assert_square_named_matrix(Sigma_D, "Sigma_D")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  if (!all(names(delta) %in% rownames(Sigma_D))) {
    stop("delta names must be present in Sigma_D", call. = FALSE)
  }
  Sigma_D <- Sigma_D[names(delta), names(delta), drop = FALSE]
  R_Delta <- stats::cov2cor(Sigma_D)
  R_Delta <- as_correlation_matrix(R_Delta, "R_Delta")
  margin <- align_margin(margin, names(delta))
  adjusted_delta <- delta - margin

  if (is.null(c_alpha)) {
    c_alpha <- critical_value_for_regime(
      R_Delta,
      alpha = alpha,
      claim_regime = claim_regime,
      adjustment_method = adjustment_method
    )
  }
  lambda <- adjusted_delta / sqrt(diag(Sigma_D))
  names(lambda) <- names(delta)
  m <- length(lambda)

  joint_success <- pmvnorm_prob(
    lower = rep(c_alpha, m),
    upper = rep(Inf, m),
    mean = lambda,
    sigma = R_Delta
  )
  marginal_power <- 1 - stats::pnorm(c_alpha - lambda)
  names(marginal_power) <- names(delta)

  list(
    alpha = alpha,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method,
    c_alpha = c_alpha,
    delta = delta,
    margin = margin,
    adjusted_delta = adjusted_delta,
    lambda = lambda,
    R_Delta = R_Delta,
    joint_success = joint_success,
    marginal_power = marginal_power
  )
}

prob_any_true_null_rejected <- function(lambda, R, c_alpha, true_null) {
  if (!is.logical(true_null) || length(true_null) != length(lambda)) {
    stop("true_null must be a logical vector with the same length as lambda", call. = FALSE)
  }
  true_index <- which(true_null)
  if (length(true_index) == 0) {
    return(NA_real_)
  }
  if (length(true_index) == 1) {
    return(1 - stats::pnorm(c_alpha - lambda[true_index]))
  }
  R_true <- R[true_index, true_index, drop = FALSE]
  lambda_true <- lambda[true_index]
  1 - pmvnorm_prob(
    lower = rep(-Inf, length(true_index)),
    upper = rep(c_alpha, length(true_index)),
    mean = lambda_true,
    sigma = R_true
  )
}

default_null_state_grid <- function(delta_alt) {
  m <- length(delta_alt)
  states <- list(global_null = rep(TRUE, m))
  if (m == 2) {
    states$partial_AB_A_null <- c(TRUE, FALSE)
    states$partial_AB_B_null <- c(FALSE, TRUE)
  } else {
    for (i in seq_len(m)) {
      state <- rep(FALSE, m)
      state[i] <- TRUE
      states[[paste0("partial_", names(delta_alt)[i], "_null")]] <- state
    }
  }
  states$all_alternatives <- rep(FALSE, m)
  states
}

evaluate_null_states <- function(
    delta_alt,
    Sigma_D,
    alpha = 0.025,
    margin = NULL,
    claim_regime = c("full_contribution", "separate_component_claims"),
    adjustment_method = NULL,
    null_state_grid = NULL) {
  assert_named_numeric(delta_alt, "delta_alt")
  assert_square_named_matrix(Sigma_D, "Sigma_D")
  claim_regime <- validate_claim_regime(claim_regime)
  adjustment_method <- resolve_adjustment_method(claim_regime, adjustment_method)
  if (!all(names(delta_alt) %in% rownames(Sigma_D))) {
    stop("delta_alt names must be present in Sigma_D", call. = FALSE)
  }

  Sigma_D <- Sigma_D[names(delta_alt), names(delta_alt), drop = FALSE]
  R_Delta <- stats::cov2cor(Sigma_D)
  c_alpha <- critical_value_for_regime(
    R_Delta,
    alpha = alpha,
    claim_regime = claim_regime,
    adjustment_method = adjustment_method
  )
  margin <- align_margin(margin, names(delta_alt))
  if (is.null(null_state_grid)) {
    null_state_grid <- default_null_state_grid(delta_alt)
  }

  rows <- lapply(names(null_state_grid), function(state_name) {
    true_null <- null_state_grid[[state_name]]
    delta <- delta_alt
    delta[true_null] <- margin[true_null]
    oc <- operating_characteristics(
      delta,
      Sigma_D,
      alpha = alpha,
      margin = margin,
      claim_regime = claim_regime,
      adjustment_method = adjustment_method,
      c_alpha = c_alpha
    )
    fwer <- prob_any_true_null_rejected(oc$lambda, oc$R_Delta, c_alpha, true_null)
    data.frame(
      state = state_name,
      t(as.data.frame(delta, optional = TRUE)),
      true_null_count = sum(true_null),
      true_null_labels = paste(names(delta_alt)[true_null], collapse = ";"),
      familywise_error_any_true_null = fwer,
      full_contribution_type_i_error =
        if (any(true_null)) oc$joint_success else NA_real_,
      false_contribution_selection_probability =
        if (any(true_null)) oc$joint_success else NA_real_,
      joint_success_probability = oc$joint_success,
      min_marginal_power = min(oc$marginal_power),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

additive_excess_to_contribution <- function(theta_A, theta_B, gamma) {
  if (!all(vapply(list(theta_A, theta_B, gamma), function(x) {
    is.numeric(x) && length(x) == 1 && is.finite(x)
  }, logical(1)))) {
    stop("theta_A, theta_B, and gamma must be finite numeric scalars", call. = FALSE)
  }
  c(
    AB_minus_A = theta_B + gamma,
    AB_minus_B = theta_A + gamma
  )
}

ratio_sensitivity_summary <- function(delta, theta_A, theta_B) {
  assert_named_numeric(delta, "delta")
  required <- c("AB_minus_A", "AB_minus_B")
  if (!all(required %in% names(delta))) {
    stop("delta must contain AB_minus_A and AB_minus_B", call. = FALSE)
  }
  c(
    r_A = unname(delta["AB_minus_B"]) / theta_A,
    r_B = unname(delta["AB_minus_A"]) / theta_B
  )
}
