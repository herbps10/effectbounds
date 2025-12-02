bound <- function(x, lower = 0, upper = 1) pmin(upper, pmax(lower, x))

cdrf_kernel <- function(a0, bw) {
  k <- \(a, a0, bw) 1 / (sqrt(2 * pi) * bw) * exp(-(a - a0)^2 / (2 * bw^2))
  n <- integrate(\(a) k(a, a0, bw), -Inf, Inf)$value
  \(a) k(a, a0, bw) / n
}

cdrf_onestep <- function(A, Y, trt_grid, a_grid, nuisance, bw) {
  N <- length(Y)
  eif <- matrix(ncol = length(trt_grid), nrow = length(Y))

  for(index in seq_along(trt_grid)) {
    k <- cdrf_kernel(trt_grid[index], bw)
    eif[, index] <- eif_cdrf(A, Y, nuisance$mu_hat, nuisance$mu_a_hat, nuisance$pi_hat, nuisance$pi_a_hat, a_grid, k)
  }

  cdrf <- apply(eif, 2, \(x) mean(x[!is.infinite(x)]))
  se <- apply(eif, 2, \(x) sd(x[!is.infinite(x)]))
  lower <- cdrf + qnorm(0.025) * se / sqrt(N)
  upper <- cdrf + qnorm(0.975) * se / sqrt(N)

  list(
    cdrf = bound(cdrf, 0, 1),
    lower = bound(lower, 0, 1),
    upper = bound(upper, 0, 1)
  )
}

estimate_cdrf_nuisance <- function(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, a_grid, outcome_type) {
  N <- nrow(data)
  mu_a_hat <- pi_a_hat <-numeric(N)
  pi_hat <- mu_hat <- matrix(nrow = N, ncol = length(a_grid))

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)

  outcome_family <- stats::gaussian()
  if(all(data[[Y]] %in% c(0, 1))) outcome_family <- stats::binomial()

  if(outer_folds > 1) {
    for(fold in seq_along(cv)) {
      training   <- cv[[fold]]$training_set
      validation <- cv[[fold]]$validation_set

      a_model <- SuperLearner::SuperLearner(
        Y = data[[A]][training],
        X = data[training, X, drop = FALSE],
        SL.library = learners_trt,
        family = "gaussian",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      squared_residuals <- (data[[A]][training] - a_model$SL.predict)^2

      a2_model <- SuperLearner::SuperLearner(
        Y = squared_residuals,
        X = data[training, X, drop = FALSE],
        SL.library = learners_trt,
        family = "gaussian",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      a_std <- (data[[A]][training] - a_model$SL.predict) / sqrt(a2_model$SL.predict)

      pi_dens <- density(a_std[!is.na(a_std)])
      pi_fun  <- approxfun(pi_dens$x, pi_dens$y, yleft = 0, yright = 0)
      pi_hat_mat <- (matrix(pi_fun(((rep(a_grid, length(training)) - rep(a_model$SL.predict, each = length(a_grid)))) / rep(sqrt(a2_model$SL.predict), each = length(a_grid))), ncol = length(a_grid), nrow = length(training), byrow = TRUE))
      var_pi_fun <- approxfun(a_grid, colMeans(pi_hat_mat), rule = 2)

      a_pred  <- SuperLearner::predict.SuperLearner(a_model, newdata = data[validation, X, drop = FALSE])$pred[, 1]
      a2_pred <- SuperLearner::predict.SuperLearner(a2_model, newdata = data[validation, X, drop = FALSE])$pred[, 1]

      pi_a_hat[validation] <- pi_fun((data[[A]][validation] - a_pred) / sqrt(a2_pred)) #/ var_pi_fun(data[[A]][validation])
      a_std_validation <- ((rep(a_grid, length(validation)) - rep(a_pred, each = length(a_grid)))) / rep(sqrt(a2_pred), each = length(a_grid))
      pi_hat[validation, ] <- (matrix(pi_fun(a_std_validation) / var_pi_fun(data[[A]][validation]), ncol = length(a_grid), nrow = length(validation), byrow = TRUE))

      mu_model <- SuperLearner::SuperLearner(
        Y = data[[Y]][training],
        X = data[training, c(X, A), drop = FALSE],
        SL.library = learners_outcome,
        family = outcome_family,
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      mu_a_hat[validation] <- SuperLearner::predict.SuperLearner(mu_model, newdata = data[validation, c(X, A)], onlySL = TRUE)$pred

      dataA <- data[rep(validation, times = length(a_grid)), X, drop = FALSE]
      dataA[[A]] <- rep(a_grid, each = length(validation))
      mu_hat[validation, ] <- matrix(SuperLearner::predict.SuperLearner(mu_model, newdata = dataA, onlySL = TRUE)$pred, ncol = length(a_grid), nrow = length(validation))
    }
  }
  else {
    pi_model <- SuperLearner::SuperLearner(
      Y = data[[A]],
      X = data[, X, drop = FALSE],
      SL.library = learners_trt,
      cvControl = cv_control,
      family = "binomial",
      env = environment(SuperLearner::SuperLearner)
    )

    mu_model <- SuperLearner::SuperLearner(
      Y = data[[Y]],
      X = data[, c(X, A), drop = FALSE],
      SL.library = learners_outcome,
      family = outcome_family,
      cvControl = cv_control,
      env = environment(SuperLearner::SuperLearner)
    )

    pi_hat  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data, onlySL = TRUE)$pred
    mu0_hat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data0, onlySL = TRUE)$pred
    mu1_hat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data1, onlySL = TRUE)$pred
  }

  #eps <- 1e-8
  #pi_hat[pi_hat == 0] <- eps
  #pi_hat[pi_hat == 1] <- 1 - eps
  #pi_a_hat[pi_a_hat == 0] <- eps
  #pi_a_hat[pi_a_hat == 1] <- 1 - eps

  list(
    pi_hat = pi_hat,
    pi_a_hat = pi_a_hat,
    mu_hat = mu_hat,
    mu_a_hat = mu_a_hat
  )
}

# One-step algorithm for non-overlap causal dose-response function
onestep_smooth_cdrf <- function(A, Y, mu, mu_a, pi, pi_a, trt_grid, a_grid, bw, threshold, smoothness, parameter = "trimmed") {
  N <- length(Y)

  eif <- matrix(nrow = N, ncol = length(trt_grid))
  for(index in seq_along(trt_grid)) {
    k <- cdrf_kernel(trt_grid[index], bw)
    if(parameter == "upper") {
      eif[, index] <- eif_cdrf_upper(A, Y, mu, mu_a, pi, pi_a, a_grid, k, threshold, smoothness)
    }
    else {
      eif[, index] <- eif_cdrf_lower(A, Y, mu, mu_a, pi, pi_a, a_grid, k, threshold, smoothness)
    }
  }
  psi <- colMeans(eif)

  ci <- matrix(NA, ncol = length(trt_grid), nrow = 2)
  ci[1, ] <- psi + qnorm(0.025) * apply(eif, 2, sd) / sqrt(N)
  ci[2, ] <- psi + qnorm(0.975) * apply(eif, 2, sd) / sqrt(N)

  list(
    psi = psi,
    eif = eif,
    ci = ci
  )
}

#' Estimate non-overlap bounds for the Causal Dose-Response Function (CDRF)
#'
#' @param data data frame containing data estimating CDRF bounds
#' @param X vector of covariate column names
#' @param A name of column containing binary treatment indicator
#' @param Y name of column containing outcome variable (bounded between zero and one)
#' @param learners_trt SuperLearner learners for estimating propensity
#' score model.
#' @param learners_outcome SuperLearner learners for estimtaing outcome
#' model.
#' @param trt_grid grid of treatment values
#' @param thresholds vector of propensity score thresholds
#' @param bandwidth kernel bandwidth
#' @param smoothness tuning parameter controlling smoothness of the indicator
#' unction approximations (smaller values imply a less smooth approximation)
#' @param alpha significance level of pointwise and
#' uniform confidence intervals (default 5%)
#' @param outer_folds Number of folds in outer cross-fitting loop
#' @param inner_folds Number of folds used by SuperLearner within
#' each outer cross-fitting loop
#' @param bootstrap whether to use multiplier bootstrap to compute
#'  uniform confidence sets
#' @param bootstrap_draws number of multiplier bootstrap draws
#' @param nuisance list containing estimated nuisance parameters (optional)
#'
#' @details
#' This function estimates non-overlap bounds for the Causal Dose-Response
#' Function (CDRF) on a user-specified grid of exposure values, propensity score
#' thresholds and smoothness, and tuning parameter values.
#'
#' By default, 95% uniform confidence sets are formed. You can configure the
#' significance level with the \code{alpha} argument; the uniform confidence
#' sets are designed to be valid with probability \eqn{(1 - \alpha) \times 100\%}.
#'
#' An important tuning parameter is the \code{smoothness} argument, which
#' controls the smoothness of an inner approximation of certain indicator
#' functions that arise in the definition of the non-overlap bounds. In practice,
#' a small value (like \code{10e-2}) can typically be used. We
#' recommend trying several small values, like \code{10e-3}, \code{10e-2},
#' and \code{10e-1} in a sensitivity analysis.
#'
#' @return A list of class \code{cdrfbounds} containing the following elements:
#' \describe{
#'  \item{bounds}{List containing estimated bounds for each smoothness parameter.}
#'  \item{smoothness}{Vector of smoothness parameters.}
#'  \item{thresholds}{Vector of propensity score thresholds.}
#'  \item{onestep}{One-step point estimate and confidence interval for CDRF.}
#'  \item{alpha}{Significance level.}
#'  \item{N}{Number of observations.}
#'  \item{K}{Number of propensity score thresholds.}
#'  \item{nuisance}{propensity score and conditional mean outcome predictions.}
#' }
#'
#' @seealso [summary.cdrfbounds]
#' @seealso [plot.cdrfbounds]
#'
#' @name cdrf_bounds
#'
#' @examples
#' dat <- simulate_cdrf_example(
#'   seed = 1,
#'   N = 5e2,
#'   alpha = 3,
#'   beta = 0.1,
#'   gamma = 1
#' )
#'
#' bounds <- cdrf_bounds(
#'   dat,
#'   X = c("X1", "X2"), A = "A", Y = "Y",
#'   thresholds = c(10^seq(-3, -0.5, 0.1)),
#'   smoothness = c(0.005)
#' )
#'
#' @export
cdrf_bounds <- function(data, X, A, Y, learners_trt = c("SL.glm"), learners_outcome = c("SL.glm"), trt_grid = seq(0, 1, 0.1), thresholds = c(10^seq(-4, -1, 0.05)), smoothness = 1e-2, bw = 0.1, alpha = 0.05, outer_folds = 5, inner_folds = 5, bootstrap = TRUE, bootstrap_draws = 1e3, nuisance = NULL) {
  #assert_ate_data(data, X, A, Y)
  #assert_folds(outer_folds)
  #assert_folds(inner_folds)
  #assert_thresholds(thresholds)
  #assert_smoothness(smoothness)
  #assert_bootstrap(bootstrap, bootstrap_draws)
  #assert_outcome_bounds(data, Y)

  K <- length(thresholds)
  N <- nrow(data)

  a_grid <- seq(min(data[[A]]) - bw * 10, max(data[[A]]) + bw * 10, length.out = 150)

  # Cross-fitted nuisance models
  if(!is.null(nuisance)) {
    assert_ate_nuisance(nuisance, N)
    nuisance$mu_hat <- ifelse(data[[A]] == 1, nuisance$mu1_hat, nuisance$mu0_hat)
  }
  else {
    nuisance <- estimate_cdrf_nuisance(data, X, A, Y, learners_trt, learners_outcome, outer_folds, inner_folds, a_grid)
  }

  onestep <- cdrf_onestep(data[[A]], data[[Y]], trt_grid, a_grid, nuisance, bw)

  results <- lapply(smoothness, \(smoothness) {
    # Set up output
    lower     <- upper     <- matrix(nrow = K, ncol = length(trt_grid))
    lower_ci  <- upper_ci  <- array(dim = c(2, K, length(trt_grid)))
    lower_eif <- upper_eif <- array(dim = c(N, K, length(trt_grid)))

    # One-step
    for(index in seq_along(thresholds)) {
      threshold <- thresholds[index]

      onestep_lower   <- onestep_smooth_cdrf(data[[A]], data[[Y]], nuisance$mu_hat, nuisance$mu_a_hat, nuisance$pi_hat, nuisance$pi_a_hat, trt_grid, a_grid, bw, threshold, smoothness, parameter = "lower")
      onestep_upper   <- onestep_smooth_cdrf(data[[A]], data[[Y]], nuisance$mu_hat, nuisance$mu_a_hat, nuisance$pi_hat, nuisance$pi_a_hat, trt_grid, a_grid, bw, threshold, smoothness, parameter = "upper")

      lower[index, ]   <- onestep_lower$psi
      upper[index, ]   <- onestep_upper$psi

      lower_eif[, index, ]   <- onestep_lower$eif
      upper_eif[, index, ]   <- onestep_upper$eif

      lower_ci[1, index, ] <- bound(onestep_lower$ci[1, ], 0, 1)
      lower_ci[2, index, ] <- bound(onestep_lower$ci[2, ], 0, 1)
      upper_ci[1, index, ] <- bound(onestep_upper$ci[1, ], 0, 1)
      upper_ci[2, index, ] <- bound(onestep_upper$ci[2, ], 0, 1)
    }

    list(
      thresholds = thresholds,
      lower = lower,
      upper = upper,
      lower_eif = lower_eif,
      upper_eif = upper_eif,
      lower_pointwise = lower_ci[1, , ],
      upper_pointwise = upper_ci[2, , ],
      lower_uniform = matrix(ncol = length(trt_grid), nrow = K),
      upper_uniform = matrix(ncol = length(trt_grid), nrow = K)
    )
  })

  uniform_critical_value <- NA
  if(bootstrap == TRUE) {
    # Multiplier bootstrap
    uniform_ci <- matrix(NA, K * length(smoothness) * length(trt_grid), 2)

    # Combine point estimates and EIFs from all smoothness options into combined vectors/matrices
    lower <- numeric(K * length(trt_grid) * length(smoothness))
    upper <- numeric(K * length(trt_grid) * length(smoothness))

    lower_eif <- matrix(nrow = N, ncol = K * length(trt_grid) * length(smoothness))
    upper_eif <- matrix(nrow = N, ncol = K * length(trt_grid) * length(smoothness))

    lower <- unlist(lapply(results, \(x) t(x$lower)))
    upper <- unlist(lapply(results, \(x) t(x$upper)))
    lower_eif <- Reduce(cbind, lapply(1:length(smoothness), \(smoothness_index) Reduce(cbind, lapply(1:K, \(threshold_index) results[[smoothness_index]]$lower_eif[,threshold_index,]))))
    upper_eif <- Reduce(cbind, lapply(1:length(smoothness), \(smoothness_index) Reduce(cbind, lapply(1:K, \(threshold_index) results[[smoothness_index]]$upper_eif[,threshold_index,]))))

    uniform_ci <- multiplier_bootstrap(lower, upper, lower_eif, upper_eif, draws = bootstrap_draws, alpha = alpha)

    uniform_critical_value <- uniform_ci$critical_value

    for(smoothness_index in seq_along(smoothness)) {
      for(threshold_index in seq_along(thresholds)) {
        ri <- ((smoothness_index - 1) * K * length(trt_grid) + (threshold_index - 1) * length(trt_grid) + 1):((smoothness_index - 1) * K * length(trt_grid) + threshold_index * length(trt_grid))
        results[[smoothness_index]]$lower_uniform[threshold_index, ] <- bound(uniform_ci$ci[ri, 1], 0, 1)
        results[[smoothness_index]]$upper_uniform[threshold_index, ] <- bound(uniform_ci$ci[ri, 2], 0, 1)
      }
    }

    tightest_bounds <- matrix(nrow = 2, ncol = length(trt_grid))
    for(tindex in seq_along(trt_grid)) {
      tightest_bounds[1, tindex] <- max(unlist(lapply(results, \(x) max(x$lower_uniform[, tindex]))))
      tightest_bounds[2, tindex] <- min(unlist(lapply(results, \(x) min(x$upper_uniform[, tindex]))))
    }
    tightest_bounds[1, ] <- tightest_bounds[1, ]
    tightest_bounds[2, ] <- tightest_bounds[2, ]
  }

  out <- list(
    trt = trt_grid,
    bounds = results,
    smoothness = smoothness,
    thresholds = thresholds,
    onestep = onestep,
    alpha = alpha,
    tightest_bounds = tightest_bounds,
    uniform_critical_value = uniform_critical_value,
    N = N,
    K = K,
    nuisance = nuisance
  )


  class(out) <- "cdrfbounds"
  out
}
