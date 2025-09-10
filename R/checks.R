check_ate_data <- function(data, X, A, Y) {
  if(TRUE != checkmate::check_data_frame(data, min.rows = 1, min.cols = 2)) return("data must be a data frame")

  if(TRUE != checkmate::check_character(X)) return("X must be a vector of column names")
  if(TRUE != checkmate::check_string(A)) return("A must be a column name")
  if(TRUE != checkmate::check_string(Y)) return("Y must be a column name")

  if(!all(X %in% names(data))) return("All X columns must be in input data")
  if(!(A %in% names(data))) return(glue::glue("treatment column '{A}' must be in input data"))
  if(!(Y %in% names(data))) return(glue::glue("outcome column '{Y}' must be in input data"))
  if(any(sort(unique(data[[Y]])) != c(0, 1))) return("Outcome must be binary")

  TRUE
}
assert_ate_data <- checkmate::makeAssertionFunction(check_ate_data)

check_folds <- function(folds) {
  if(TRUE != checkmate::check_number(folds, lower = 1)) return("folds must be a number >= 1")
  TRUE
}
assert_folds <- checkmate::makeAssertionFunction(check_folds)

check_thresholds <- function(thresholds) {
  if(TRUE != checkmate::check_numeric(thresholds, lower = 0, upper = 1, min.len = 1)) return("thresholds must be a vector of numbers, all between 0 and 1")
  TRUE
}
assert_thresholds <- checkmate::makeAssertionFunction(check_thresholds)

check_smoothness <- function(smoothness) {
  if(TRUE != checkmate::check_numeric(smoothness, lower = 0)) return("smoothness must be a vector, all >= 0")
  TRUE
}
assert_smoothness <- checkmate::makeAssertionFunction(check_smoothness)

check_bootstrap <- function(bootstrap, bootstrap_draws) {
  if(TRUE != checkmate::check_logical(bootstrap)) return("bootstrap must be TRUE or FALSE")
  if(TRUE != checkmate::check_number(bootstrap_draws, lower = 1)) return("bootstrap_draws must be a number >= 1")
  TRUE
}
assert_bootstrap <- checkmate::makeAssertionFunction(check_bootstrap)

check_ate_nuisance <- function(nuisance, N) {
  if(TRUE != checkmate::assert_list(nuisance)) return("nuisance must be a list")
  if(TRUE != checkmate::assert_numeric(nuisance$pi_hat, lower = 0, upper = 1, len = N)) return("nuisance$pi_hat must be a vector")
  if(TRUE != checkmate::assert_numeric(nuisance$mu0_hat, lower = 0, upper = 1, len = N)) return("nuisance$mu0_hat must be a vector")
  if(TRUE != checkmate::assert_numeric(nuisance$mu1_hat, lower = 0, upper = 1, len = N)) return("nuisance$mu1_hat must be a vector")

  TRUE
}
assert_ate_nuisance <- checkmate::makeAssertionFunction(check_ate_nuisance)
