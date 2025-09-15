# Smooth approximation of the indicator function I[x > t] with smoothness gamma
# Note that, by construction, s(x) <= I[x > t].
s_gt <- function(x, t, gamma) ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 1, 1 - exp(1) * exp(1 / (((x - t) / gamma)^2 - 1))))

# Derivative of smooth approximation of indicator function I[X > t] with smoothness gamma
s_gt_dot <- function(x, t, gamma) {
  ifelse(x - t <= 0, 0, ifelse(x - t >= 0 + gamma, 0, 2 * gamma^2 * exp(1 / ((t - x)^2 / gamma^2 - 1) + 1) * (x - t) / (gamma^2 - (t - x)^2)^2))
}

# Smooth approximation of the indicator function I[x < t] with smoothness gamma
# Note that, by construction, s(x) <= I[x < t]
s_lt <- function(x, t, gamma) ifelse(x <= t - gamma, 1, ifelse(x >= t, 0,  1 - exp(1 / (((x - t) / gamma)^2 - 1))/ exp(-1)))

# Derivative of smooth approximation of indicator function I[X > t] with smoothness gamma
s_lt_dot <- function(x, t, gamma) {
  ifelse(x <= t - gamma, 0, ifelse(x >= t, 0, 2 * gamma^2 * exp(1 / ((t - x)^2 / gamma^2 - 1) + 1) * (x - t) / (gamma^2 - (t - x)^2)^2))
}

#' Compute one-step estimator for one-dimensional target parameter
#'
#' @param plugin plugin point estimate
#' @param eif EIF vector
#' @noRd
onestep <- function(plugin, eif) plugin + mean(eif)

#' Compute CI from influence function estimate
#'
#' @param estimate point estimtae
#' @param eif estimate of EIF (vector)
#' @param N sample size
#' @param alpha significance level
#'
#' @noRd
ci <- function(estimate, eif, N, alpha = 0.05) {
  estimate + stats::qnorm(c(alpha / 2, 1 - alpha / 2)) * stats::sd(eif) / sqrt(N)
}

#' Apply multiplier bootstrap for lower and upper effect bounds
#'
#' @param lower point estimates of lower bound for set of propensity score thresholds (vector)
#' @param upper point estimates of upper bound for set of propensity score thresholds (vector)
#' @param eif_lower efficient influence function estimates for lower bound (matrix)
#' @param eif_upper efficient influence function estimates for upper bound (matrix)
#' @param draws number of bootstrap values
#' @param alpha significance level, leading to (1 - alpha)% uniform confidence set
#'
#' @noRd
multiplier_bootstrap <- function(lower, upper, eif_lower, eif_upper, draws = 1e3, alpha = 0.05) {
  N <- nrow(eif_lower)
  K <- ncol(eif_lower)

  se_lower <- apply(eif_lower, 2, stats::sd) |> matrix(ncol = K, nrow = N, byrow = TRUE)
  se_upper <- apply(eif_upper, 2, stats::sd) |> matrix(ncol = K, nrow = N, byrow = TRUE)
  mean_lower <- apply(eif_lower, 2, mean) |> matrix(ncol = K, nrow = N, byrow = TRUE)
  mean_upper <- apply(eif_upper, 2, mean) |> matrix(ncol = K, nrow = N, byrow = TRUE)
  eif_lower_scaled <- (eif_lower - mean_lower) / se_lower
  eif_upper_scaled <- (eif_upper - mean_upper) / se_upper

  zs <- matrix(2 * stats::rbinom(draws * N, 1, 0.5) - 1, nrow = N, ncol = draws)

  T_lower <- lapply(1:draws, \(draw) colSums(zs[, draw] * eif_lower_scaled) / sqrt(N))
  T_upper <- lapply(1:draws, \(draw) colSums(zs[, draw] * eif_upper_scaled) / sqrt(N))

  T_max <- lapply(1:draws, \(draw) pmax(T_lower[[draw]], -T_upper[[draw]]))

  max_max <- unlist(lapply(T_max, max, na.rm = TRUE))

  # Critical value
  calpha <- stats::quantile(max_max, 1 - alpha)

  list(
    ci = matrix(c(
      lower - calpha * apply(eif_lower, 2, stats::sd) / sqrt(N),
      upper + calpha * apply(eif_upper, 2, stats::sd) / sqrt(N)
    ), ncol = 2, nrow = K, byrow = FALSE),
    critical_value = calpha
  )
}
