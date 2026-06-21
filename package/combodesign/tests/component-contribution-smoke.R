library(combodesign)

delta <- c(AB_minus_A = 0.55, AB_minus_B = 0.65)
allocation <- c(A = 1 / 3, B = 1 / 3, AB = 1 / 3)

fixed <- contribution_operating_characteristics(
  N = 120,
  allocation = allocation,
  delta = delta
)

stopifnot(isTRUE(fixed$joint_success_probability > 0))
stopifnot(isTRUE(fixed$joint_success_probability < 1))
stopifnot(identical(names(fixed$delta), c("AB_minus_A", "AB_minus_B")))

fit <- component_contribution_design(
  delta = delta,
  target_power = 0.70,
  min_allocation = 0.10
)

stopifnot(is.data.frame(fit$summary))
stopifnot(nrow(fit$summary) == 3)
stopifnot(all(fit$summary$achieved_joint_success_probability >= 0.70))

alias_fit <- optimal_design(
  delta = delta,
  target_power = 0.70,
  min_allocation = 0.10
)

stopifnot(identical(names(alias_fit), names(fit)))
