##### ATE #####

# Efficient Influence Function (EIF) for the Average Treatment Effect (ATE) parameter
#
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
#
# @return vector of efficient influence function values
eif_ate <- function(A, Y, mu0, mu1, pi) {
  mu <- ifelse(A == 1, mu1, mu0)
  mu1 - mu0 - mean(mu1 - mu0) + (A / pi - (1 - A) / (1 - pi)) * (Y - mu)
}

# Efficient Influence Function (EIF) for the smooth trimmed ATE parameter
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_trimmed <- function(A, Y, mu0, mu1, pi, threshold, smoothness) {
  w1 <- s_gt(pi, threshold, smoothness)
  w1_dot <- s_gt_dot(pi, threshold, smoothness)
  w0 <- s_lt(pi, 1 - threshold, smoothness)
  w0_dot <- s_lt_dot(pi, 1 - threshold, smoothness)

  A / pi * w1 * (Y - mu1) - (1 - A) / (1 - pi) * w0 * (Y - mu0) +
    (mu1 * w1_dot - mu0 * w0_dot) * (A - pi) +
    mu1 * w1 - mu0 * w0 - mean(mu1 * w1 - mu0 * w0)
}

# Efficient Influence Function (EIF) for the smooth lower ATE bound
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_lower <- function(A, Y, mu0, mu1, pi, threshold, smoothness) {
  eif_trimmed(A, Y, mu0, mu1, pi, threshold, smoothness) +
    s_lt_dot(pi, 1 - threshold, smoothness) * (A - pi) + s_lt(pi, 1 - threshold, smoothness) - mean(s_lt(pi, 1 - threshold, smoothness))
}

# Efficient Influence Function (EIF) for the smooth lower ATE bound
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_upper <- function(A, Y, mu0, mu1, pi, threshold, smoothness) {
  eif_trimmed(A, Y, mu0, mu1, pi, threshold, smoothness) -
    s_gt_dot(pi, threshold, smoothness) * (A - pi) - s_gt(pi, threshold, smoothness) + mean(s_gt(pi, threshold, smoothness))
}

##### Survival #####

# Efficient Influence Function (EIF) for the smooth survival effect parameter
# @param A binary treatment indicator
# @param C censoring indicator
# @param Y binary outcome
# @param surv0 survival given A = 0
# @param surv1 survival given A = 1
# @param pi propensity score P(A = 1 | X)
# @param times
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_survival_trimmed <- function(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness) {
  N <- length(Y)

  eif0 <- matrix(nrow = N, ncol = length(times))
  eif1 <- matrix(nrow = N, ncol = length(times))

  h0 <- matrix((A == 0) / (1 - pi), nrow = N, ncol = length(times), byrow = FALSE)
  h1 <- matrix((A == 1) / (pi), nrow = N, ncol = length(times), byrow = FALSE)

  int0 <- hazard0 / (surv0 * cens)
  int1 <- hazard1 / (surv1 * cens)

  # Indicator for discrete integral
  # 1 / (S(y | a, w) * G(y | a, w)
  yind <- matrix(times, nrow = N, ncol = length(times), byrow = TRUE) == matrix(Y, nrow = N, ncol = length(times), byrow = FALSE)
  yinv0 <- matrix(1 / rowSums(surv0 * cens * yind), nrow = N, ncol = length(times), byrow = FALSE)
  yinv1 <- matrix(1 / rowSums(surv1 * cens * yind), nrow = N, ncol = length(times), byrow = FALSE)

  #browser()

  tind <- matrix(Y, nrow = N, ncol = length(times), byrow = FALSE) >= matrix(times, nrow = N, ncol = length(times), byrow = TRUE)

  omega0 <- matrix((1 - pi), nrow = N, ncol = length(times), byrow = FALSE) * cens
  omega1 <- matrix(pi, nrow = N, ncol = length(times), byrow = FALSE) * cens

  w0     <- s_gt(omega0, threshold, smoothness)
  w0_dot <- s_gt_dot(omega0, threshold, smoothness)
  w1     <- s_gt(omega1, threshold, smoothness)
  w1_dot <- s_gt_dot(omega1, threshold, smoothness)

  cens_mat <- (matrix(C, nrow = length(Y), ncol = length(times), byrow = FALSE) == 1) & (matrix(Y, nrow = length(Y), ncol = length(times), byrow = FALSE) >= matrix(times, nrow = length(Y), ncol = length(times), byrow = TRUE))
  cens_lt_mat <- (matrix(C, nrow = length(Y), ncol = length(times), byrow = FALSE) == 1) & (matrix(Y, nrow = length(Y), ncol = length(times), byrow = FALSE) <= matrix(times, nrow = length(Y), ncol = length(times), byrow = TRUE))

  cumsums0 <- t(apply(tind * int0, 1, cumsum))
  cumsums1 <- t(apply(tind * int1, 1, cumsum))

  eif0 <- w0 * surv0 * (1 - h0 * (cens_lt_mat * yinv0 - cumsums0)) +
    w0_dot * surv0 * matrix(A == 0, nrow = N, ncol = length(times), byrow = FALSE) * (cens_mat - cens) +
    w0_dot * surv0 * matrix((A == 0) - (1 - pi), nrow = N, ncol = length(times), byrow = FALSE) -
    matrix(colMeans(w0 * surv0), nrow = N, ncol = length(times), byrow = TRUE)

  eif1 <- w1 * surv1 * (1 - h1 * (cens_lt_mat * yinv1 - cumsums1)) +
    w1_dot * surv1 * matrix(A == 1, nrow = N, ncol = length(times), byrow = FALSE) * (cens_mat - cens) +
    w1_dot * surv1 * matrix((A == 1) - (pi), nrow = N, ncol = length(times), byrow = FALSE) -
    matrix(colMeans(w1 * surv1), nrow = N, ncol = length(times), byrow = TRUE)

  eif1 - eif0
}

# Efficient Influence Function (EIF) for the smooth survival effect lower bound
# @param A binary treatment indicator
# @param C binary outcome
# @param Y binary outcome
# @param surv0 survival given A = 0
# @param surv1 survival given A = 1
# @param pi propensity score P(A = 1 | X)
# @param times
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_survival_lower <- function(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness) {
  x <- eif_survival_trimmed(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness)

  omega <- matrix(pi, nrow = length(Y), ncol = length(times), byrow = FALSE) * cens
  w <- s_gt(omega, threshold, smoothness)
  w_dot <- s_gt_dot(omega, threshold, smoothness)

  cens_mat <- (matrix(C, nrow = length(Y), ncol = length(times), byrow = FALSE) == 1) & (matrix(Y, nrow = length(Y), ncol = length(times), byrow = FALSE) >= matrix(times, nrow = length(Y), ncol = length(times), byrow = TRUE))

  x + w_dot * matrix(((A == 1) - pi), nrow = length(Y), ncol = length(times), byrow = FALSE) +
    w_dot * (cens_mat - cens) +
    w - colMeans(w)
}

# Efficient Influence Function (EIF) for the smooth survival effect upper bound
# @param A binary treatment indicator
# @param C censoring indicator
# @param Y binary outcome
# @param surv0 survival given A = 0
# @param surv1 survival given A = 1
# @param pi propensity score P(A = 1 | X)
# @param times
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_survival_upper <- function(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness) {
  x <- eif_survival_trimmed(A, C, Y, surv0, surv1, hazard0, hazard1, cens, pi, times, threshold, smoothness)

  omega <- matrix(1 - pi, nrow = length(Y), ncol = length(times), byrow = FALSE) * cens
  w <- s_gt(omega, threshold, smoothness)
  w_dot <- s_gt_dot(omega, threshold, smoothness)

  cens_mat <- (matrix(C, nrow = length(Y), ncol = length(times), byrow = FALSE) == 1) & (matrix(Y, nrow = length(Y), ncol = length(times), byrow = FALSE) >= matrix(times, nrow = length(Y), ncol = length(times), byrow = TRUE))

  x + w_dot * matrix(((A == 0) - (1 - pi)), nrow = length(Y), ncol = length(times), byrow = FALSE) +
    w_dot * (cens_mat - cens) +
    w - colMeans(w)

  x
}
