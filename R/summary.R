###### ATE #####

#' Summarize ATE non-overlap bounds.
#'
#' @param object list with class \code{atebounds}
#' @param ... additional arguments (not used)
#'
#' @return list with class \code{summary.atebounds}
#' @name summary.atebounds
#'
#' @export
summary.atebounds <- function(object, ...) {
  shortest_intervals <- data.frame(
    smoothness = object$smoothness,
    threshold = rep(NA, length(object$smoothness)),
    lower_uniform = rep(NA, length(object$smoothness)),
    upper_uniform = rep(NA, length(object$smoothness))
  )

  for(index in seq_along(object$smoothness)) {
    shortest_intervals$lower_uniform[index] <- max(object$bounds[[index]]$lower_uniform)
    shortest_intervals$upper_uniform[index] <- min(object$bounds[[index]]$upper_uniform)
  }

  out <- list(
    propensity_score_range = range(object$nuisance$pi_hat),
    alpha = object$alpha,
    onestep = object$onestep,
    shortest_intervals = shortest_intervals
  )
  class(out) <- "summary.atebounds"
  out
}

#' Print summary of non-overlap bounds for ATE.
#' @param x list with class \code{summary.atebounds}
#' @param ... additional arguments (not used)
#'
#' @return printed summary of ATE non-overlap bounds
#' @name summary.atebounds
#'
#' @export
print.summary.atebounds <- function(x, ...) {
  cli::cli_text("{.strong Summary: ATE non-overlap bounds}")

  l <- format(x$propensity_score_range[1], digits = 3)
  u <- format(x$propensity_score_range[2], digits = 3)
  cat(glue::glue("\n\n Range of estimated propensity scores: {l} - {u}"))

  l <- format(1 / x$propensity_score_range[2], digits = 3)
  u <- format(1 / x$propensity_score_range[1], digits = 3)
  cat(glue::glue("\n\n Range of estimated inverse propensity scores: {l} - {u}"))

  point <- format(x$onestep$ate, digits = 3)
  l <- format(x$onestep$lower, digits = 3)
  u <- format(x$onestep$upper, digits = 3)
  cat(glue::glue("\n\nATE point estimate: {point}. 95% CI: ({l}, {u})\n\n"))

  cat(glue::glue("\n\nSmoothness \t Shortest Uniform {(1-x$alpha)*100}% Bounds \n\n"))
  for(index in seq_along(1:nrow(x$shortest_intervals))) {
    s <- x$shortest_intervals$smoothness[index]
    l <- format(x$shortest_intervals$lower_uniform[index], digits = 3)
    u <- format(x$shortest_intervals$upper_uniform[index], digits = 3)
    cat(glue::glue("{s} \t\t ({l}, {u}) \n\n"))
  }
}

##### Transport ATE #####

#' Summarize non-overlap bounds for transport ATE.
#'
#' @param object list with class \code{transportbounds}
#' @param ... additional arguments (not used)
#'
#' @return list with class \code{summary.transportbounds}
#' @name summary.atebounds
#'
#' @export
summary.transportbounds <- function(object, ...) {
  shortest_intervals <- data.frame(
    smoothness = object$smoothness,
    threshold = rep(NA, length(object$smoothness)),
    lower_uniform = rep(NA, length(object$smoothness)),
    upper_uniform = rep(NA, length(object$smoothness))
  )

  for(index in seq_along(object$smoothness)) {
    shortest_intervals$lower_uniform[index] <- max(object$bounds[[index]]$lower_uniform)
    shortest_intervals$upper_uniform[index] <- min(object$bounds[[index]]$upper_uniform)
  }

  out <- list(
    propensity_score_range = range(object$nuisance$pi_hat),
    alpha = object$alpha,
    onestep = object$onestep,
    shortest_intervals = shortest_intervals
  )
  class(out) <- "summary.transportbounds"
  out
}

#' Print summary of non-overlap bounds for ATE.
#' @param x list with class \code{summary.atebounds}
#' @param ... additional arguments (not used)
#'
#' @return printed summary of ATE non-overlap bounds
#' @name summary.atebounds
#'
#' @export
print.summary.transportbounds <- function(x, ...) {
  cli::cli_text("{.strong Summary: Transport ATE non-overlap bounds}")

  l <- format(x$propensity_score_range[1], digits = 3)
  u <- format(x$propensity_score_range[2], digits = 3)
  cat(glue::glue("\n\n Range of estimated propensity scores: {l} - {u}"))

  l <- format(1 / x$propensity_score_range[2], digits = 3)
  u <- format(1 / x$propensity_score_range[1], digits = 3)
  cat(glue::glue("\n\n Range of estimated inverse propensity scores: {l} - {u}"))

  point <- format(x$onestep$ate, digits = 3)
  l <- format(x$onestep$lower, digits = 3)
  u <- format(x$onestep$upper, digits = 3)
  cat(glue::glue("\n\nTransport ATE point estimate: {point}. 95% CI: ({l}, {u})\n\n"))

  cat(glue::glue("\n\nSmoothness \t Shortest Uniform {(1-x$alpha)*100}% Bounds \n\n"))
  for(index in seq_along(1:nrow(x$shortest_intervals))) {
    s <- x$shortest_intervals$smoothness[index]
    l <- format(x$shortest_intervals$lower_uniform[index], digits = 3)
    u <- format(x$shortest_intervals$upper_uniform[index], digits = 3)
    cat(glue::glue("{s} \t\t ({l}, {u}) \n\n"))
  }
}
