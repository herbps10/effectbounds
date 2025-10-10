transport_onestep <- function(S, A, Y, nuisance) {
  plugin <- mean((nuisance$mu1_hat - nuisance$mu0_hat)[S == 0])
  eif <- with(nuisance, 1 / mean(S == 0) * (
    (1 - phi_hat) / phi_hat * (S == 1) * (A / pi_hat - (1 - A) / (1 - pi_hat)) * (Y - mu_hat) + (S == 0) * (mu1_hat - mu0_hat)
  )) - plugin

  effect <- plugin + mean(eif)

  ci <- effect + stats::qnorm(c(0.025, 0.975)) * stats::sd(eif) / sqrt(length(Y))

  list(
    effect = effect,
    lower = ci[1],
    upper = ci[2],
    test = ci[2] < 0 || ci[1] > 0
  )
}

estimate_transport_nuisance <- function(data, X, S, A, Y, learners_trt, learners_source, learners_outcome, outer_folds, inner_folds, outcome_type) {
  N <- nrow(data)
  data0 <- data1 <- data
  data0[[A]] <- 0
  data1[[A]] <- 1
  phi_hat <- pi_hat <- mu0_hat <- mu1_hat <- numeric(N)

  cv <- origami::make_folds(nrow(data), origami::folds_vfold, V = outer_folds)
  cv_control <- SuperLearner::SuperLearner.CV.control(V = inner_folds)

  outcome_family <- stats::gaussian()
  if(all(data[[Y]] %in% c(0, 1))) outcome_family <- stats::binomial()

  if(outer_folds > 1) {
    for(fold in seq_along(cv)) {
      training   <- cv[[fold]]$training_set
      validation <- cv[[fold]]$validation_set

      training_source <- training[which(data[[S]][training] == 1)]

      phi_model <- SuperLearner::SuperLearner(
        Y = data[[S]][training],
        X = data[training, X, drop = FALSE],
        SL.library = learners_trt,
        family = "binomial",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      pi_model <- SuperLearner::SuperLearner(
        Y = data[[A]][training_source],
        X = data[training_source, X, drop = FALSE],
        SL.library = learners_trt,
        family = "binomial",
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      mu_model <- SuperLearner::SuperLearner(
        Y = data[[Y]][training_source],
        X = data[training_source, c(X, A), drop = FALSE],
        SL.library = learners_outcome,
        family = outcome_family,
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )

      phi_hat[validation] <- SuperLearner::predict.SuperLearner(phi_model, newdata = data[validation, ], onlySL = TRUE)$pred
      pi_hat[validation]  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data[validation, ], onlySL = TRUE)$pred
      mu0_hat[validation] <- SuperLearner::predict.SuperLearner(mu_model, newdata = data0[validation, ], onlySL = TRUE)$pred
      mu1_hat[validation] <- SuperLearner::predict.SuperLearner(mu_model, newdata = data1[validation, ], onlySL = TRUE)$pred
    }
  }
  else {
    phi_model <- SuperLearner::SuperLearner(
      Y = data[[S]],
      X = data[, X, drop = FALSE],
      SL.library = learners_trt,
      cvControl = cv_control,
      family = "binomial",
      env = environment(SuperLearner::SuperLearner)
    )

    pi_model <- SuperLearner::SuperLearner(
      Y = data[[A]][data[[S]] == 1],
      X = data[data[[S]] == 1, X, drop = FALSE],
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

    phi_hat <- SuperLearner::predict.SuperLearner(phi_model, newdata = data, onlySL = TRUE)$pred
    pi_hat  <- SuperLearner::predict.SuperLearner(pi_model, newdata = data, onlySL = TRUE)$pred
    mu0_hat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data0, onlySL = TRUE)$pred
    mu1_hat <- SuperLearner::predict.SuperLearner(mu_model, newdata = data1, onlySL = TRUE)$pred
  }
  mu_hat <- ifelse(data[[A]] == 1, mu1_hat, mu0_hat)

  list(
    phi_hat = phi_hat,
    pi_hat = pi_hat,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    mu_hat = mu_hat
  )
}

# TMLE algorithm for non-overlap transport bound parameters
tmle_smooth_transport <- function(S, A, Y, mu0, mu1, phi, pi, threshold, smoothness, parameter = "trimmed", maxiter = 25, verbose = FALSE) {
  fluctuation <- \(epsilon, mu0, mu1, phi, pi) {
    cleverA <- rep(0, length(pi))
    cleverS <- rep(0, length(phi))

    w0 <- s_gt((1 - pi) * phi, threshold, smoothness)
    w1 <- s_gt(pi * phi, threshold, smoothness)

    clevermu0 <- -(S == 1) / mean(S == 0) * (1 - phi) / phi * (1 - A) * w0 / (1 - pi)
    clevermu1 <-  (S == 1) / mean(S == 0) * (1 - phi) / phi * A * w1 / pi

    if(smoothness > 0) {
      w0_dot <- s_gt_dot((1 - pi) * phi, threshold, smoothness)
      w1_dot <- s_gt_dot(pi * phi, threshold, smoothness)

      cleverA <- (S == 1) / mean(S == 0) * (1 - phi) * (mu1 * w1_dot - mu0 * w0_dot)
      cleverS <- 1 / mean(S == 0) * (1 - phi) * (mu1 * w1_dot * pi - mu0 * w0_dot * (1 - pi))

      if(parameter == "lower") {
        cleverA <- cleverA - (S == 1) / mean(S == 0) * w0_dot * (1 - phi)
        cleverS <- cleverS + 1 / mean(S == 0) * w0_dot * (1 - pi) * (1 - phi)
      }
      if(parameter == "upper") {
        cleverA <- cleverA - (S == 1) / mean(S == 0) * w1_dot * (1 - phi)
        cleverS <- cleverS - 1 / mean(S == 0) * w1_dot * pi * (1 - phi)
      }
    }

    list(
      mu0 = stats::plogis(stats::qlogis(mu0) + epsilon * clevermu0),
      mu1 = stats::plogis(stats::qlogis(mu1) + epsilon * clevermu1),
      pi  = stats::plogis(stats::qlogis(pi)  + epsilon * cleverA),
      phi = stats::plogis(stats::qlogis(phi) + epsilon * cleverS)
    )
  }

  loss <- \(params, mu0, mu1, phi, pi) {
    f <- fluctuation(params, mu0, mu1, phi, pi)
    x <- mean(
      (S == 1) * (-A * log(f$pi) - (1 - A) * log(1 - f$pi) + ifelse(A == 1, -Y * log(f$mu1) - (1 - Y) * log(1 - f$mu1), -Y * log(f$mu0) - (1 - Y) * log(1 - f$mu0))) #- S * log(f$phi) - (1 - S) * log(1 - f$phi)
    )
    if(is.infinite(x) || is.nan(x)) return(Inf)
    x
  }

  # Start at initial estimators
  mu0_star <- mu0
  mu1_star <- mu1
  pi_star  <- pi
  phi_star <- phi
  converged <- FALSE
  for(iter in 1:maxiter) {
    # evaluate loss at bounds
    left_bound <- c(-1)
    right_bound <- 1
    if(is.infinite(loss(0, mu0_star, mu1_star, phi_star, pi_star))) stop("Infinite loss")
    while(is.infinite(loss(left_bound, mu0_star, mu1_star, phi_star, pi_star))) left_bound <- left_bound / 2
    while(is.infinite(loss(right_bound, mu0_star, mu1_star, phi_star, pi_star))) right_bound <- right_bound / 2

    epsilon_star <- stats::optimize(loss, interval = c(left_bound, right_bound), mu0 = mu0_star, mu1 = mu1_star, phi = phi_star, pi = pi_star)$minimum
    f <- fluctuation(epsilon_star, mu0_star, mu1_star, phi_star, pi_star)

    mu0_star <- f$mu0
    mu1_star <- f$mu1
    pi_star  <- f$pi
    phi_star <- f$phi

    if(verbose) cat(glue::glue("Iter: {iter} epsilon_star: {epsilon_star} smoothness: {smoothness} threshold: {threshold} bounds: ({left_bound}, {right_bound}) \n\n"))

    if(abs(epsilon_star) < 1e-2) {
      converged <- TRUE
      break
    }
  }
  if(converged == FALSE) {
    warning("TMLE Failed to converge")
    return(list(
      psi = NA,
      phi = NA,
      eif = NA,
      mu0 = NA,
      mu1 = NA,
      pi = NA,
      ci = NA
    ))
  }

  w0 <- s_gt((1 - pi_star) * phi_star, threshold, smoothness)
  w1 <- s_gt(pi_star * phi_star, threshold, smoothness)
  psi_trimmed <- mean((mu1_star * w1 - mu0_star * w0) * (S == 0)) / mean(S == 0)

  if(parameter == "trimmed") {
    psi <- psi_trimmed
    eif <- eif_transport_trimmed(S, A, Y, mu0_star, mu1_star, phi_star, pi_star, threshold, smoothness)
  }
  else if(parameter == "upper") {
    psi <- psi_trimmed + 1 / mean(S == 0) * mean((S == 0) * (1 - s_gt(pi_star * phi_star, threshold, smoothness)))
    eif <- eif_transport_upper(S, A, Y, mu0_star, mu1_star, phi_star, pi_star, threshold, smoothness)
  }
  else if(parameter == "lower") {
    psi <- psi_trimmed - 1 / mean(S == 0) * mean((S == 0) * (1 - s_gt((1 - pi_star) * phi_star, threshold, smoothness)))
    eif <- eif_transport_lower(S, A, Y, mu0_star, mu1_star, phi_star, pi_star, threshold, smoothness)
  }

  list(
    psi = psi,
    eif = eif,
    mu0 = mu0_star,
    mu1 = mu1_star,
    pi = pi_star,
    phi = phi_star,
    ci = psi + stats::qnorm(c(0.025, 0.975)) * stats::sd(eif) / sqrt(length(Y))
  )
}

#' Estimate non-overlap bounds for the Transported Average Treatment Effect
#'
#' @param data data frame containing data estimating transport effect bounds
#' @param X vector of covariate column names
#' @param S name of column containing binary source indicator
#' @param A name of column containing binary treatment indicator
#' @param Y name of column containing outcome variable (bounded between zero and one)
#' @param learners_trt SuperLearner learners for estimating propensity
#' score model.
#' @param learners_source SuperLearner learners for estimating source model
#' @param learners_outcome SuperLearner learners for estimtaing outcome
#' model.
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
#' This function estimates non-overlap bounds for the transport effect on a
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
#' @return A list of class \code{transportbounds} containing the following elements:
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
#' @seealso [summary.transportbounds]
#' @seealso [plot.transportbounds]
#'
#' @name transport_bounds
#'
#' @examples
#' dat <- simulate_transport_example(
#'   seed = 1,
#'   N = 5e2,
#'   alpha = 3,
#'   beta = 0.1,
#'   gamma = 1
#' )
#'
#' bounds <- transport_bounds(
#'   dat,
#'   X = c("X1", "X2"), S = "S", A = "A", Y = "Y",
#'   thresholds = c(10^seq(-3, -0.5, 0.1)),
#'   smoothness = c(0.005)
#' )
#'
#' @export
transport_bounds <- function(data, X, S, A, Y, learners_trt = c("SL.glm"), learners_source = c("SL.glm"), learners_outcome = c("SL.glm"), thresholds = c(10^seq(-4, -1, 0.05)), smoothness = 1e-2, alpha = 0.05, outer_folds = 5, inner_folds = 5, bootstrap = TRUE, bootstrap_draws = 1e3, nuisance = NULL) {
  assert_ate_data(data, X, A, Y)
  assert_folds(outer_folds)
  assert_folds(inner_folds)
  assert_thresholds(thresholds)
  assert_smoothness(smoothness)
  assert_bootstrap(bootstrap, bootstrap_draws)
  assert_outcome_bounds(data, Y)

  K <- length(thresholds)
  N <- nrow(data)

  # Cross-fitted nuisance models
  if(!is.null(nuisance)) {
    assert_ate_nuisance(nuisance, N)
    nuisance$mu_hat <- ifelse(data[[A]] == 1, nuisance$mu1_hat, nuisance$mu0_hat)
  }
  else {
    nuisance <- estimate_transport_nuisance(data, X, S, A, Y, learners_trt, learners_source, learners_outcome, outer_folds, inner_folds)
  }

  results <- lapply(smoothness, \(smoothness) {
    # Set up output
    trimmed     <- lower     <- upper     <- numeric(K)
    trimmed_ci  <- lower_ci  <- upper_ci  <- matrix(ncol = 2, nrow = K)
    trimmed_eif <- lower_eif <- upper_eif <- matrix(nrow = N, ncol = K)

    # TMLE
    for(index in seq_along(thresholds)) {
      threshold <- thresholds[index]

      tmle_lower   <- tmle_smooth_transport(data[[S]], data[[A]], data[[Y]], nuisance$mu0_hat, nuisance$mu1_hat, nuisance$phi_hat, nuisance$pi_hat, threshold, smoothness, parameter = "lower")
      tmle_upper   <- tmle_smooth_transport(data[[S]], data[[A]], data[[Y]], nuisance$mu0_hat, nuisance$mu1_hat, nuisance$phi_hat, nuisance$pi_hat, threshold, smoothness, parameter = "upper")
      tmle_trimmed <- tmle_smooth_transport(data[[S]], data[[A]], data[[Y]], nuisance$mu0_hat, nuisance$mu1_hat, nuisance$phi_hat, nuisance$pi_hat, threshold, smoothness, parameter = "trimmed")

      trimmed[index] <- tmle_trimmed$psi
      lower[index]   <- tmle_lower$psi
      upper[index]   <- tmle_upper$psi

      trimmed_eif[, index] <- tmle_trimmed$eif
      lower_eif[, index]   <- tmle_lower$eif
      upper_eif[, index]   <- tmle_upper$eif

      trimmed_ci[index, 1] <- tmle_trimmed$ci[1]
      trimmed_ci[index, 2] <- tmle_trimmed$ci[2]
      lower_ci[index, 1]   <- tmle_lower$ci[1]
      lower_ci[index, 2]   <- tmle_lower$ci[2]
      upper_ci[index, 1]   <- tmle_upper$ci[1]
      upper_ci[index, 2]   <- tmle_upper$ci[2]
    }

    list(
      thresholds = thresholds,
      lower = lower,
      upper = upper,
      lower_eif = lower_eif,
      upper_eif = upper_eif,
      lower_pointwise = lower_ci[, 1],
      upper_pointwise = upper_ci[, 2]
    )
  })

  uniform_critical_value <- NA
  if(bootstrap == TRUE) {
    # Multiplier bootstrap

    uniform_ci <- matrix(NA, K * length(smoothness), 2)

    # Combine point estimates and EIFs from all smoothness options into combined vectors/matrices
    lower <- unlist(lapply(results, `[[`, "lower"))
    upper <- unlist(lapply(results, `[[`, "upper"))
    lower_eif <- do.call(cbind, lapply(results, `[[`, "lower_eif"))
    upper_eif <- do.call(cbind, lapply(results, `[[`, "upper_eif"))

    uniform_ci <- multiplier_bootstrap(lower, upper, lower_eif, upper_eif, draws = bootstrap_draws, alpha = alpha)

    uniform_critical_value <- uniform_ci$critical_value

    for(index in seq_along(smoothness)) {
      ri <- ((index - 1) * K + 1):(index * K)
      results[[index]]$lower_uniform <- pmax(-1, uniform_ci$ci[ri, 1])
      results[[index]]$upper_uniform <- pmin(1, uniform_ci$ci[ri, 2])
    }
  }

  out <- list(
    bounds = results,
    smoothness = smoothness,
    thresholds = thresholds,
    onestep = transport_onestep(data[[S]], data[[A]], data[[Y]], nuisance),
    alpha = alpha,
    uniform_critical_value = uniform_critical_value,
    N = N,
    K = K,
    nuisance = nuisance
  )

  class(out) <- "transportbounds"
  out
}
