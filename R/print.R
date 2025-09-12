#' Print estimated non-overlap bounds for ATE
#'
#' @param x list with class \code{atebounds}
#' @param ... additional arguments for print.atebounds (not used)
#'
#' @name ate_bounds
#'
#' @export
print.atebounds <- function(x, ...) {
  cli::cli_text("{.strong ATE non-overlap bounds}")
  cli::cli_text("{x$N} observations, {x$K} propensity score thresholds")
}
