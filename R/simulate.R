#' Simulate from a simple illustrative data-generating process designed
#' to illustrate ATE non-overlap bounds
#'
#' @param seed random number seed
#' @param N number of observations
#' @param alpha strength of non-overlap violations (higher alpha implies more extreme non-overlap violation)
#' @param beta proportion of population subject to non-overlap violation
#' @param gamma controls magnitude of the effect
#' @param binary_outcome specifies whether to generate a binary (default) or continuous outcome
#'
#' @return data frame with columns X1, X2, A, and Y
#'
#' @export
simulate_ate_example <- function(seed = 1, N = 1e3, alpha = 1, beta = 0.05, gamma = 1, binary_outcome = TRUE) {
  set.seed(seed)
  X1 <- stats::runif(N, -1, 1)
  x  <- stats::runif(N)
  X2 <- ifelse(x < beta / 2, -1, ifelse(x > 1 - beta / 2, 1, 0))
  A  <- stats::rbinom(N, 1, stats::plogis(X1 + alpha * X2))

  if(binary_outcome == TRUE) {
    Y  <- stats::rbinom(N, 1, stats::plogis(X1 + (A - 0.5) * gamma))
  }
  else {
    Y  <- stats::rnorm(N, X1 + (A - 0.5) * gamma, 0.1)
  }

  data.frame(X1 = X1, X2 = X2, A = A, Y = Y)
}

#' Simulate from a simple illustrative data-generating process designed
#' to illustrate survival non-overlap bounds
#'
#' @param seed random number seed
#' @param N number of observations
#' @param alpha strength of non-overlap violations (higher alpha implies more extreme non-overlap violation)
#' @param beta proportion of population subject to non-overlap violation
#' @param gamma controls magnitude of the effect
#' @param max_time maximum observed event time
#'
#' @return data frame with columns X1, X2, A, C, and Y
#'
#' @export
simulate_survival_example <- function(seed = 1, N = 1e3, alpha = 1, beta = 0.05, gamma = 1, max_time = 5) {
  set.seed(seed)
  X1 <- stats::runif(N, -1, 1)
  x  <- stats::runif(N)
  X2 <- ifelse(x < beta / 2, -1, ifelse(x > 1 - beta / 2, 1, 0))
  pi <- stats::plogis(X1 + alpha * X2)
  A  <- stats::rbinom(N, 1, stats::plogis(X1 + alpha * X2))

  event_time <- stats::rexp(N, 1 + A * gamma)
  censoring_time <- pmin(stats::rexp(N, 0.01), max_time)

  Y <- pmin(event_time, censoring_time)
  C <- event_time < censoring_time

  data.frame(X1 = X1, X2 = X2, A = A, C = C, Y = Y, pi = pi)
}
