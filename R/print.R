#' Print estimated ATE non-overlap bounds
#'
#' @param x object of class "atebounds"
#' @param ... further arguments (not used)
#'
#' @export
print.atebounds <- function(x, ...) {
  cli::cli_text("{.strong ATE non-overlap bounds}")
  cli::cli_text("{x$N} observations, {x$K} propensity score thresholds")
}
