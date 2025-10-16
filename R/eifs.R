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

  w1     <- s_gt(pi, threshold, smoothness)
  w1_dot <- s_gt_dot(pi, threshold, smoothness)
  w0     <- s_lt(pi, 1 - threshold, smoothness)
  w0_dot <- s_lt_dot(pi, 1 - threshold, smoothness)

  eif0 <- matrix(nrow = N, ncol = length(times))
  eif1 <- matrix(nrow = N, ncol = length(times))

  h0 <- (A == 0) / (1 - pi)
  h1 <- (A == 1) / (pi)

  int0 <- hazard0 / (surv0 * cens)
  int1 <- hazard1 / (surv1 * cens)

  # Indicator for discrete integral
  tind <- matrix(0, nrow = N, ncol = length(times))

  # 1 / (S(y | a, w) * G(y | a, w)
  yind <- matrix(times, nrow = N, ncol = length(times), byrow = TRUE) == matrix(Y, nrow = N, ncol = length(times), byrow = FALSE)
  yinv0 <- 1 / rowSums(surv0 * cens * yind)
  yinv1 <- 1 / rowSums(surv1 * cens * yind)

  for(tindex in 1:length(times)) {
    t0 <- times[tindex]
    tind[, tindex] <- Y >= t0

    eif0[, tindex] <- w0 * surv0[, tindex] * (1 - h0 * ((Y <= t0 & C == 1) * yinv0 - rowSums(tind * int0))) + (w0_dot * surv0[, tindex]) * (A - pi) - mean(w0 * surv0[, tindex])
    eif1[, tindex] <- w1 * surv1[, tindex] * (1 - h1 * ((Y <= t0 & C == 1) * yinv1 - rowSums(tind * int1))) + (w1_dot * surv1[, tindex]) * (A - pi) - mean(w1 * surv1[, tindex])
  }

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
  x <- x + s_lt_dot(pi, 1 - threshold, smoothness) * (A - pi) +
    s_lt(pi, 1 - threshold, smoothness) - mean(s_lt(pi, 1 - threshold, smoothness))
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
  x <- x - s_gt_dot(pi, threshold, smoothness) * (A - pi) -
    s_gt(pi, threshold, smoothness) + mean(s_gt(pi, threshold, smoothness))
}
