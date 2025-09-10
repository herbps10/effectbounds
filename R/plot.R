#' Plot estimated ATE non-overlap bounds
#'
#' @param x object of type "atebounds"
#' @param smoothness which smoothness tuning parameter bounds to plot; set to NA to print all
#' @param point_estimate whether to plot point estimate and 95% confidence interval
#' @param ylim limits of y axis
#' @param xlab x axis label
#' @param ylab y axis label
#' @param legend_position where to show legend (set to "none" to hide legend)
#' @param bounds_color color of point estimate and 95% confidence interval
#' @param point_estimate_color color of point estimate and 95% confidence interval
#' @param ... additional arguments passed to plot.default
#' @export
plot.atebounds <- function(x, smoothness = x$smoothness[1], point_estimate = FALSE, ylim = NA, xlab = "Propensity Score Threshold", ylab = "ATE", legend_position = "bottomright", bounds_color = "black", point_estimate_color = "blue", ...) {
  if(is.na(smoothness)) {
    indexes <- seq_along(x$smoothness)
  }
  else {
    if(!any(x$smoothness == smoothness)) stop(glue::glue("Smoothness {smoothness} not found in ATE bounds object."))
    indexes <- which(x$smoothness == smoothness)
  }

  if(any(is.na(ylim))) {
    ylim <- range(unlist(lapply(x$bounds, \(bounds) range(c(bounds$lower_uniform, bounds$upper_uniform)))))
    if(point_estimate == TRUE) ylim <- range(c(x$onestep$lower, x$onestep$upper, ylim))
  }

  plot(1, type = "n", xlim = range(c(0, x$thresholds)), ylim = ylim, xlab = xlab, ylab = ylab)
  abline(h = 0, col = "gray")
  for(index in indexes) {
    points(x = x$thresholds, y = x$bounds[[index]]$lower_uniform, pch = 20, col = bounds_color)
    points(x = x$thresholds, y = x$bounds[[index]]$upper_uniform, pch = 20, col = bounds_color)

    lines(x = x$thresholds, y = x$bounds[[index]]$lower_uniform, col = bounds_color)
    lines(x = x$thresholds, y = x$bounds[[index]]$upper_uniform, col = bounds_color)
  }

  bound_title <- glue::glue("Non-overlap {(1 - x$alpha) * 100}% bounds")

  if(point_estimate == TRUE) {
    points(x = 0, y = x$onestep$ate, col = point_estimate_color, pch = 20)
    lines(x = c(0, 0), y = c(x$onestep$lower, x$onestep$upper), col = point_estimate_color)

    if(legend_position != "none") legend(legend_position, c(bound_title, "ATE point estimate and 95% CI"), fill = c(bounds_color, point_estimate_color))
  }
  else {
    if(legend_position != "none") legend(legend_position, c(bound_title), fill = c(bounds_color))
  }
}
