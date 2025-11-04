isotonize <- function(s) 1 - isoreg(1:length(s), 1 - s)$yf

bound <- function(x, lower = 0, upper = 1) pmin(upper, pmax(lower, x))

survival_onestep <- function(A, C, Y, times, nuisance, bootstrap, draws, alpha) {
  N <- length(Y)

  h0 <- matrix((A == 0) / (1 - nuisance$pi_hat), nrow = N, ncol = length(times))
  h1 <- matrix((A == 1) / (nuisance$pi_hat), nrow = N, ncol = length(times))

  int0 <- nuisance$hazard0_hat / (nuisance$surv0_hat * nuisance$cens_hat)
  int1 <- nuisance$hazard1_hat / (nuisance$surv1_hat * nuisance$cens_hat)

  # Indicator for discrete integral
  tind <- matrix(Y, nrow = N, ncol = length(times), byrow = FALSE) >= matrix(times, nrow = N, ncol = length(times), byrow = TRUE)

  # 1 / (S(y | a, w) * G(y | a, w)
  yind <- matrix(times, nrow = N, ncol = length(times), byrow = TRUE) == matrix(Y, nrow = N, ncol = length(times), byrow = FALSE)
  yinv0 <- matrix(1 / rowSums(nuisance$surv0_hat * nuisance$cens_hat * yind), nrow = N, ncol = length(times), byrow = FALSE)
  yinv1 <- matrix(1 / rowSums(nuisance$surv1_hat * nuisance$cens_hat * yind), nrow = N, ncol = length(times), byrow = FALSE)

  cens_lt_mat <- (matrix(C, nrow = length(Y), ncol = length(times), byrow = FALSE) == 1) & (matrix(Y, nrow = length(Y), ncol = length(times), byrow = FALSE) <= matrix(times, nrow = length(Y), ncol = length(times), byrow = TRUE))

  cumsums0 <- t(apply(tind * int0, 1, cumsum))
  cumsums1 <- t(apply(tind * int1, 1, cumsum))

  eif0 <- nuisance$surv0_hat * (1 - h0 * (cens_lt_mat * yinv0 - cumsums0))
  eif1 <- nuisance$surv1_hat * (1 - h1 * (cens_lt_mat * yinv1 - cumsums1))

  curve0 <- isotonize(bound(colMeans(eif0)))
  curve0_lower <- isotonize(bound(curve0 + stats::qnorm(0.025) * apply(eif0, 2, stats::sd) / sqrt(N)))
  curve0_upper <- isotonize(bound(curve0 + stats::qnorm(0.975) * apply(eif0, 2, stats::sd) / sqrt(N)))

  curve1 <- isotonize(bound(colMeans(eif1)))
  curve1_lower <- isotonize(bound(curve1 + stats::qnorm(0.025) * apply(eif1, 2, stats::sd) / sqrt(N)))
  curve1_upper <- isotonize(bound(curve1 + stats::qnorm(0.975) * apply(eif1, 2, stats::sd) / sqrt(N)))

  effect <- curve1 - curve0
  effect_lower <- effect + stats::qnorm(0.025) * apply(eif1 - eif0, 2, stats::sd) / sqrt(N)
  effect_upper <- effect + stats::qnorm(0.975) * apply(eif1 - eif0, 2, stats::sd) / sqrt(N)

  effect_lower_uniform <- effect_upper_uniform <- rep(NA, length(times))
  calpha <- NA
  if(bootstrap == TRUE) {
    zs <- matrix(2 * stats::rbinom(draws * N, 1, 0.5) - 1, nrow = N, ncol = draws)
    eif_scaled <- (eif1 - eif0 - (curve1 - curve0))
    eif_scaled <- eif_scaled / matrix(apply(eif_scaled, 2, stats::sd), ncol = length(times), nrow = N, byrow = TRUE)
    T_effect <- unlist(lapply(1:draws, \(draw) max(colSums(zs[, draw] * eif_scaled) / sqrt(N))))
    calpha <- stats::quantile(T_effect, 1 - alpha)
    effect_lower_uniform <- effect - calpha * apply(eif1 - eif0, 2, stats::sd) / sqrt(N)
    effect_upper_uniform <- effect + calpha * apply(eif1 - eif0, 2, stats::sd) / sqrt(N)
  }

  list(
    effect = effect,
    effect_lower = effect_lower,
    effect_upper = effect_upper,
    effect_lower_uniform = effect_lower_uniform,
    effect_upper_uniform = effect_upper_uniform,
    calpha = calpha,

    curve1 = curve1,
    curve1_lower = curve1_lower,
    curve1_upper = curve1_upper,

    curve0 = curve0,
    curve0_lower = curve0_lower,
    curve0_upper = curve0_upper
  )
}

estimate_survival_nuisance <- function(data, X, A, C, Y, learners_trt, learners_event, learners_cens, times, outer_folds, inner_folds, outcome_type) {
  N <- nrow(data)
  data0 <- data1 <- data
  data0[[A]] <- 0
  data1[[A]] <- 1
  pi_hat <- numeric(N)
  surv0_hat <- surv1_hat <- cens0_hat <- cens1_hat <- hazard0_hat <- hazard1_hat <- matrix(nrow = N, ncol = length(times))

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)

  outcome_family <- stats::gaussian()
  if(all(data[[Y]] %in% c(0, 1))) outcome_family <- stats::binomial()

  if(outer_folds > 1) {
    for(fold in seq_along(cv)) {
      training   <- cv[[fold]]$training_set
      validation <- cv[[fold]]$validation_set

      pi_model <- SuperLearner::SuperLearner(
        Y = data[[A]][training],
        X = data[training, X, drop = FALSE],
        SL.library = learners_trt,
        family = "binomial",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      newX <- rbind(data0[validation, ], data1[validation, ])

      survival_model <- survSuperLearner::survSuperLearner(
        time = data[[Y]][training],
        event = data[[C]][training],
        X = data[training, c(X, A), drop = FALSE],
        newX = newX[, c(X, A), drop = FALSE],
        new.times = times,
        event.SL.library = learners_event,
        cens.SL.library = learners_cens,
        cvControl = cv_control
      )

      part0 <- 1:length(validation)
      part1 <- (length(validation) + 1):(2*length(validation))

      pi_hat[validation]  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data[validation, ], onlySL = TRUE)$pred
      surv0_hat[validation, ] <- survival_model$event.SL.predict[part0, ]
      surv1_hat[validation, ] <- survival_model$event.SL.predict[part1, ]
      cens0_hat[validation, ] <- cbind(1, survival_model$cens.SL.predict[part0, -length(times)])
      cens1_hat[validation, ] <- cbind(1, survival_model$cens.SL.predict[part1, -length(times)])
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

    newX <- rbind(data0, data1)

    survival_model <- survSuperLearner::survSuperLearner(
      time = data[[Y]],
      event = data[[C]],
      X = data[, c(X, A), drop = FALSE],
      newX = newX[c(X, A), drop = FALSE],
      new.times = times,
      event.SL.library = learners_event,
      cens.SL.library = learners_cens,
      cvControl = cv_control
    )

    pi_hat  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data, onlySL = TRUE)$pred
    surv0_hat[validation, ] <- survival_model$event.SL.predict[1:nrow(data), ]
    surv1_hat[validation, ] <- survival_model$event.SL.predict[(nrow(data) + 1):(2 * nrow(data)), ]
    cens0_hat[validation, ] <- survival_model$cens.SL.predict[1:nrow(data), ]
    cens1_hat[validation, ] <- survival_model$cens.SL.predict[(nrow(data) + 1):(2 * nrow(data)), ]
  }
  surv_hat <- ifelse(matrix(data[[A]] == 1, nrow = N, ncol = length(times), byrow = FALSE), surv1_hat, surv0_hat)
  cens_hat <- ifelse(matrix(data[[A]] == 1, nrow = N, ncol = length(times), byrow = FALSE), cens1_hat, cens0_hat)

  surv0_hat <- surv0_hat
  surv1_hat <- surv1_hat
  surv_hat <- surv_hat

  hazard0_hat <- cbind(1 - surv0_hat[, 1], 1 - surv0_hat[, 2:ncol(surv0_hat)] / surv0_hat[, 1:(ncol(surv0_hat) - 1)])
  hazard1_hat <- cbind(1 - surv1_hat[, 1], 1 - surv1_hat[, 2:ncol(surv0_hat)] / surv1_hat[, 1:(ncol(surv1_hat) - 1)])
  hazard_hat  <- cbind(1 - surv_hat[, 1], 1 - surv_hat[, 2:ncol(surv0_hat)] / surv_hat[, 1:(ncol(surv_hat) - 1)])

  list(
    pi_hat      = pi_hat,
    surv0_hat   = surv0_hat,
    surv1_hat   = surv1_hat,
    surv_hat    = surv_hat,
    cens0_hat   = cens0_hat,
    cens1_hat   = cens1_hat,
    cens_hat    = cens_hat,
    hazard0_hat = hazard0_hat,
    hazard1_hat = hazard1_hat,
    hazard_hat  = hazard_hat
  )
}

# One-step algorithm for non-overlap survival bound parameters
onestep_smooth_survival <- function(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness, parameter = "trimmed") {
  N <- length(Y)

  psi_trimmed <- colMeans(
    matrix(s_gt(pi, threshold, smoothness), nrow = N, ncol = length(times)) * surv1 -
    matrix(s_lt(pi, 1 - threshold, smoothness), nrow = N, ncol = length(times)) * surv0
  )

  if(parameter == "trimmed") {
    psi <- psi_trimmed
    eif <- eif_survival_trimmed(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness)
  }
  else if(parameter == "upper") {
    omega <- matrix(1 - pi, nrow = N, ncol = length(times), byrow = FALSE) * cens
    psi <- psi_trimmed + (1 - colMeans(s_gt(omega, threshold, smoothness)))
    eif <- eif_survival_upper(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness)
  }
  else if(parameter == "lower") {
    omega <- matrix(pi, nrow = N, ncol = length(times), byrow = FALSE) * cens
    psi <- psi_trimmed + (1 - colMeans(s_gt(omega, threshold, smoothness)))
    eif <- eif_survival_lower(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness)
  }

  psi <- psi + colMeans(eif)

  ci <- matrix(NA, ncol = length(times), nrow = 2)
  ci[1, ] <- psi + qnorm(0.025) * apply(eif, 2, sd) / sqrt(N)
  ci[2, ] <- psi + qnorm(0.975) * apply(eif, 2, sd) / sqrt(N)

  list(
    psi = psi,
    eif = eif,
    ci = ci
  )
}

#' Estimate non-overlap bounds for the survival curves
#'
#' @param data data frame containing data estimating survival curve bounds
#' @param X vector of covariate column names
#' @param A name of column containing binary treatment indicator
#' @param C name of column containing censoring indicator (C = 1 indicates observed, C = 0 indicates censored)
#' @param Y name of column containing time-to-event variable
#' @param times vector of times to estimate survival curve
#' @param learners_trt survSuperLearner learners for estimating survival model
#' @param learners_event survSuperLearner learners for estimating survival model
#' @param learners_cens survSuperLearner learners for estimating censoring model
#' @param thresholds vector of propensity score thresholds
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
#' This function estimates non-overlap bounds for the ATE on a
#' user-specified grid of propensity score thresholds and smoothness
#' tuning parameter values.
#'
#' The basic idea of non-overlap bounds is to divide the population into
#' two parts: a subpopulation in which overlap is satisfied, and a subpopulation
#' in which overlap is violated. The propensity score thresholds defines how the
#' division is done: all units with estimated propensity score below the
#' threshold are considered to be in the non-overlap subpopulation, and the rest
#' of the units are considered to be in the overlap subpopulation.
#'
#' A-priori, we don't typically know what propensity score threshold will lead
#' to the tightest bounds on the ATE. Therefore, we try a range of values
#' (using the \code{threshold} argument), and the function uses a multiplier
#' bootstrap technique to form a uniform confidence set. You can configure the
#' multiplier bootstrap through the \code{bootstrap_draws}
#' argument, which sets how many bootstrap replications to use.
#'
#' By default, 95% uniform confidence sets are formed. You can configure the
#' significance level with the \code{alpha} argument; the uniform confidence
#' sets are designed to be valid with probability \eqn{(1 - \alpha) \times 100\%}.
#'
#' The propensity scores and outcome regression models needed to estimate the
#' bounds can be estimated using any method, including logistic regressions or
#' flexible machine-learning algorithms. To control overfitting, we use
#' sample-splitting methods. There is an outer layer of cross-fitting,
#' in which the sample is split into folds (you can set how many using the
#' \code{outer_folds} argument), models are fit on the training sample for each
#' fold, and propensity scores and conditional mean outcomes are predicted on
#' the validation sample. Within each of these outer folds, \code{SuperLearner}
#' is used to form an ensemble of learners to optimally predict the propensity
#' scores and conditional mean outcomes. \code{SuperLearner} itself uses
#' cross-validation to estimate the performance of each learner, and you can
#' configure the number of cross-validation folds using the \code{inner_folds}
#' argument.
#'
#' You can use the function \code{SuperLearner::listWrappers} for a list of the
#' algorithms available for inclusion in the \code{learners_trt} and
#' \code{learners_outcome} arguments.
#'
#' An important tuning parameter is the \code{smoothness} argument, which
#' controls the smoothness of an inner approximation of certain indicator
#' functions that arise in the definition of the non-overlap bounds. In practice,
#' a small value (like \code{10e-2}) can typically be used. We
#' recommend trying several small values, like \code{10e-3}, \code{10e-2},
#' and \code{10e-1} in a sensitivity analysis.
#'
#' @return A list of class \code{surivalbounds} containing the following elements:
#' \describe{
#'  \item{bounds}{List containing estimated bounds for each smoothness parameter.}
#'  \item{smoothness}{Vector of smoothness parameters.}
#'  \item{thresholds}{Vector of propensity score thresholds.}
#'  \item{onestep}{Doubly-robust one-step point estimate and confidence interval for ATE.}
#'  \item{alpha}{Significance level.}
#'  \item{N}{Number of observations.}
#'  \item{K}{Number of propensity score thresholds.}
#'  \item{nuisance}{propensity score and conditional mean outcome predictions.}
#' }
#'
#' @seealso [summary.survivalbounds]
#' @seealso [plot.survivalbounds]
#'
#' @name survival_bounds
#'
#' @examples
#' dat <- simulate_survival_example(
#'   seed = 1,
#'   N = 5e2,
#'   alpha = 3,
#'   beta = 0.1,
#'   gamma = 1
#' )
#'
#' bounds <- survival_bounds(
#'   dat,
#'   X = c("X1", "X2"), A = "A", Y = "Y",
#'   thresholds = c(10^seq(-3, -0.5, 0.1)),
#'   smoothness = c(0.005)
#' )
#'
#' @export
survival_bounds <- function(data, X, A, C, Y, times = NULL, learners_trt = c("SL.glm"), learners_event = c("survSL.coxph"), learners_cens = c("survSL.coxph"), thresholds = c(10^seq(-4, -1, 0.05)), smoothness = 1e-2, alpha = 0.05, outer_folds = 5, inner_folds = 5, bootstrap = TRUE, bootstrap_draws = 1e3, nuisance = NULL) {
  assert_survival_data(data, X, A, C, Y)
  assert_folds(outer_folds)
  assert_folds(inner_folds)
  assert_thresholds(thresholds)
  assert_smoothness(smoothness)
  assert_bootstrap(bootstrap, bootstrap_draws)
  assert_outcome_nonnegative(data, Y)

  if(is.null(times)) {
    times <- sort(unique(data[[Y]]))
  }

  K <- length(thresholds)
  N <- nrow(data)

  # Cross-fitted nuisance models
  if(!is.null(nuisance)) {
    assert_survival_nuisance(nuisance, N)
    nuisance$surv_hat <- ifelse(matrix(data[[A]] == 1, nrow = N, ncol = length(times), byrow = FALSE), nuisance$surv1_hat, nuisance$surv0_hat)
    nuisance$cens_hat <- ifelse(matrix(data[[A]] == 1, nrow = N, ncol = length(times), byrow = FALSE), nuisance$cens1_hat, nuisance$cens0_hat)
  }
  else {
    nuisance <- estimate_survival_nuisance(data, X, A, C, Y, learners_trt, learners_event, learners_cens, times, outer_folds, inner_folds)
  }

  results <- lapply(smoothness, \(smoothness) {
    # Set up output
    trimmed     <- lower     <- upper     <- matrix(nrow = K, ncol = length(times))
    trimmed_ci  <- lower_ci  <- upper_ci  <- array(dim = c(2, K, length(times)))
    trimmed_eif <- lower_eif <- upper_eif <- array(dim = c(N, K, length(times)))

    for(index in seq_along(thresholds)) {
      threshold <- thresholds[index]

      onestep_lower   <- onestep_smooth_survival(data[[A]], data[[C]], data[[Y]], nuisance$surv0_hat, nuisance$surv1_hat, nuisance$hazard0_hat, nuisance$hazard1_hat, nuisance$cens_hat, nuisance$pi_hat, times, threshold, smoothness, parameter = "lower")
      onestep_upper   <- onestep_smooth_survival(data[[A]], data[[C]], data[[Y]], nuisance$surv0_hat, nuisance$surv1_hat, nuisance$hazard0_hat, nuisance$hazard1_hat, nuisance$cens_hat, nuisance$pi_hat, times, threshold, smoothness, parameter = "upper")
      onestep_trimmed <- onestep_smooth_survival(data[[A]], data[[C]], data[[Y]], nuisance$surv0_hat, nuisance$surv1_hat, nuisance$hazard0_hat, nuisance$hazard1_hat, nuisance$cens_hat, nuisance$pi_hat, times, threshold, smoothness, parameter = "trimmed")

      trimmed[index, ] <- onestep_trimmed$psi
      lower[index, ]   <- onestep_lower$psi
      upper[index, ]   <- onestep_upper$psi

      trimmed_eif[, index, ] <- onestep_trimmed$eif
      lower_eif[, index, ]   <- onestep_lower$eif
      upper_eif[, index, ]   <- onestep_upper$eif

      trimmed_ci[1, index, ] <- onestep_trimmed$ci[1,]
      trimmed_ci[2, index, ] <- onestep_trimmed$ci[2,]
      lower_ci[1, index, ]   <- bound(onestep_lower$ci[1,], -1, 1)
      lower_ci[2, index, ]   <- bound(onestep_lower$ci[2,], -1, 1)
      upper_ci[1, index, ]   <- bound(onestep_upper$ci[1,], -1, 1)
      upper_ci[2, index, ]   <- bound(onestep_upper$ci[2,], -1, 1)
    }

    list(
      thresholds = thresholds,
      lower = lower,
      upper = upper,
      lower_eif = lower_eif,
      upper_eif = upper_eif,
      lower_pointwise = lower_ci[1, , ],
      upper_pointwise = upper_ci[2, , ],
      lower_uniform = matrix(ncol = length(times), nrow = K),
      upper_uniform = matrix(ncol = length(times), nrow = K)
    )
  })

  uniform_critical_value <- NA
  tightest_bounds <- NA
  if(bootstrap == TRUE) {
    # Multiplier bootstrap

    uniform_ci <- matrix(NA, K * length(smoothness) * length(times), 2)

    # Combine point estimates and EIFs from all smoothness options into combined vectors/matrices
    lower <- numeric(K * length(times) * length(smoothness))
    upper <- numeric(K * length(times) * length(smoothness))
    lower_eif <- matrix(nrow = N, ncol = K * length(times) * length(smoothness))
    upper_eif <- matrix(nrow = N, ncol = K * length(times) * length(smoothness))

    lower <- unlist(lapply(results, \(x) t(x$lower)))
    upper <- unlist(lapply(results, \(x) t(x$upper)))
    lower_eif <- Reduce(cbind, lapply(1:length(smoothness), \(smoothness_index) Reduce(cbind, lapply(1:K, \(threshold_index) results[[smoothness_index]]$lower_eif[,threshold_index,]))))
    upper_eif <- Reduce(cbind, lapply(1:length(smoothness), \(smoothness_index) Reduce(cbind, lapply(1:K, \(threshold_index) results[[smoothness_index]]$upper_eif[,threshold_index,]))))

    uniform_ci <- multiplier_bootstrap(lower, upper, lower_eif, upper_eif, draws = bootstrap_draws, alpha = alpha)

    uniform_critical_value <- uniform_ci$critical_value

    for(smoothness_index in seq_along(smoothness)) {
      for(threshold_index in seq_along(thresholds)) {
        ri <- ((smoothness_index - 1) * K * length(times) + (threshold_index - 1) * length(times) + 1):((smoothness_index - 1) * K * length(times) + threshold_index * length(times))
        results[[smoothness_index]]$lower_uniform[threshold_index, ] <- bound(uniform_ci$ci[ri, 1], -1, 1)
        results[[smoothness_index]]$upper_uniform[threshold_index, ] <- bound(uniform_ci$ci[ri, 2], -1, 1)
      }
    }

    tightest_bounds <- matrix(nrow = 2, ncol = length(times))
    for(tindex in seq_along(times)) {
      tightest_bounds[1, tindex] <- max(unlist(lapply(results, \(x) max(x$lower_uniform[, tindex]))))
      tightest_bounds[2, tindex] <- min(unlist(lapply(results, \(x) min(x$upper_uniform[, tindex]))))
    }
    tightest_bounds[1, ] <- tightest_bounds[1, ]
    tightest_bounds[2, ] <- tightest_bounds[2, ]
  }

  onestep_results <- survival_onestep(data[[A]], data[[C]], data[[Y]], times, nuisance, bootstrap, bootstrap_draws, alpha)

  out <- list(
    bounds = results,
    tightest_bounds = tightest_bounds,
    smoothness = smoothness,
    thresholds = thresholds,
    times = times,
    onestep = onestep_results,
    alpha = alpha,
    uniform_critical_value = uniform_critical_value,
    N = N,
    K = K,
    nuisance = nuisance
  )

  class(out) <- "atebounds"
  out
}
