#' Fixed-design operating characteristics for contribution contrasts
#'
#' Computes operating characteristics for a fixed total sample size and
#' allocation. The default claim regime is full conjunctive contribution, which
#' uses intersection union testing. The optional separate component-claim
#' regime uses maxT by default. The default arm labels are `A`, `B`, and `AB`,
#' with contribution contrasts `AB_minus_A` and `AB_minus_B`.
#'
#' @param N Total sample size.
#' @param allocation Named allocation proportions for the trial arms.
#' @param delta Named standardized contribution effects. For the default
#'   three-arm setting, use names `AB_minus_A` and `AB_minus_B`.
#' @param rho Optional arm-level correlation matrix with row and column names
#'   matching `arms`. This is allowed only when `covariance_model` is
#'   `"general_correlated"`.
#' @param sigma2 Common endpoint variance.
#' @param alpha One-sided type I error rate. The default is `0.025`.
#' @param margin Optional contribution margins. Use names `AB_minus_A` and
#'   `AB_minus_B` in the default setting. The default is zero margins.
#' @param claim_regime Either `"full_contribution"` for the primary
#'   intersection union rule or `"separate_component_claims"` for optional
#'   separately reportable component claims.
#' @param adjustment_method Optional adjustment method. Defaults to `"iut"` for
#'   full contribution and `"maxT"` for separate component claims.
#' @param covariance_model `"independent_arms"` for the primary randomized
#'   clinical trial model or `"general_correlated"` for a justified
#'   correlated-estimator extension.
#' @param arms Character vector of arm names.
#' @param contribution_map Optional data frame with columns `combination`,
#'   `base`, and `component`.
#' @param integer_allocation Whether deterministic integer arm sizes should be
#'   used for the fixed total sample size.
#'
#' @return A list containing the contribution contrast matrix, contrast
#'   covariance, one-sided critical value, joint power, marginal
#'   powers, allocation, and sample sizes.
#' @examples
#' delta <- c(AB_minus_A = 0.45, AB_minus_B = 0.55)
#' allocation <- c(A = 1/3, B = 1/3, AB = 1/3)
#' contribution_operating_characteristics(
#'   N = 120,
#'   allocation = allocation,
#'   delta = delta
#' )
#' @export
contribution_operating_characteristics <- function(
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
  )
}

#' Joint-power design for a component-contribution trial
#'
#' Compares equal allocation, joint-decision optimized allocation, and a
#' diagnostic marginal-power allocation under the component-contribution
#' framework. The default primary design target is full contribution joint
#' power under IUT.
#'
#' @param delta Named standardized contribution effects. For the default
#'   three-arm setting, use names `AB_minus_A` and `AB_minus_B`.
#' @param rho Optional arm-level correlation matrix with row and column names
#'   matching `arms`. This is allowed only when `covariance_model` is
#'   `"general_correlated"`.
#' @param sigma2 Common endpoint variance.
#' @param target_power Target joint power.
#' @param alpha One-sided type I error rate. The default is `0.025`.
#' @param margin Optional contribution margins. The default is zero margins.
#' @param claim_regime Either `"full_contribution"` or
#'   `"separate_component_claims"`.
#' @param adjustment_method Optional adjustment method. Defaults to `"iut"` for
#'   full contribution and `"maxT"` for separate component claims.
#' @param covariance_model Primary `"independent_arms"` model or
#'   `"general_correlated"` extension.
#' @param arms Character vector of arm names.
#' @param contribution_map Optional data frame with columns `combination`,
#'   `base`, and `component`.
#' @param min_allocation Minimum allowed arm allocation proportion.
#' @param N_upper Upper bound used by the sample-size search.
#'
#' @return A list with a comparison summary and fitted design objects.
#' @examples
#' delta <- c(AB_minus_A = 0.55, AB_minus_B = 0.65)
#' fit <- component_contribution_design(
#'   delta = delta,
#'   target_power = 0.70,
#'   min_allocation = 0.10
#' )
#' fit$summary
#' @export
component_contribution_design <- function(
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
  compare_allocation_strategies(
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
}

#' Optimal design wrapper for the contribution framework
#'
#' This function is retained as a compatibility name for earlier package users.
#' It calls [component_contribution_design()] and uses one-sided `AB_minus_A`
#' and `AB_minus_B` contribution effects. Older argument sets are not supported
#' by this wrapper.
#'
#' @inheritParams component_contribution_design
#'
#' @return A list with a comparison summary and fitted design objects.
#' @examples
#' delta <- c(AB_minus_A = 0.55, AB_minus_B = 0.65)
#' optimal_design(delta = delta, target_power = 0.70, min_allocation = 0.10)
#' @export
optimal_design <- function(
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
  component_contribution_design(
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
}
