#' Print estimated non-overlap bounds for ATE
#'
#' @param x list with class \code{atebounds}
#' @param ... additional arguments (not used)
#'
#' @return printed \code{atebounds} object.
#'
#' @name ate_bounds
#'
#' @examples
#' @export
print.atebounds <- function(x, ...) {
  cli::cli_text("{.strong ATE non-overlap bounds}")
  cli::cli_text("{x$N} observations, {x$K} propensity score thresholds")
}
