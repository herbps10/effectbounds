###### EIFs for ATE parameter #####

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


###### EIFs for transport parameter #####

# Efficient Influence Function (EIF) for the Transport Average Treatment Effect (TATE) parameter
#
# @param S binary source indicator
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
#
# @return vector of efficient influence function values
eif_transport <- function(S, A, Y, mu0, mu1, phi, pi) {
  mu <- ifelse(A == 1, mu1, mu0)
  1 / mean(S == 0) * (
    (S == 0) * (mu1 - mu0) - mean((mu1 - mu0) * (S == 0)) + (S == 1) / phi *  (A / pi - (1 - A) / (1 - pi)) * (Y - mu)
  )
}

# Efficient Influence Function (EIF) for the smooth trimmed transport parameter
# @param A binary source indicator
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_transport_trimmed <- function(S, A, Y, mu0, mu1, phi, pi, threshold, smoothness) {
  w1 <- s_gt(pi * phi, threshold, smoothness)
  w1_dot <- s_gt_dot(pi * phi, threshold, smoothness)
  w0 <- s_gt((1 - pi) * phi, threshold, smoothness)
  w0_dot <- s_gt_dot((1 - pi) * phi, threshold, smoothness)

  mu <- ifelse(A == 1, mu1, mu0)

  s1_weight <- (S == 1) / mean(S == 0)
  s0_weight <- (S == 0) / mean(S == 0)
  phi_weight <- (1 - phi) / phi

  plugin <- (S == 0) / mean(S == 0) * (mu1 * w1 - mu0 * w0)

  eif <- s1_weight * phi_weight * (A / pi * w1  - (1 - A) / (1 - pi) * w0) * (Y - mu)
  eif <- eif + s1_weight * (1 - phi) * (mu1 * w1_dot + mu0 * w0_dot ) * (A - pi)
  eif <- eif + 1 / mean(S == 0) * (mu1 * w1_dot * pi - mu0 * w0_dot * (1 - pi)) * (S - phi)
  eif <- eif + s0_weight * (mu1 * w1 - mu0 * w0 - plugin)
  eif
}

# Efficient Influence Function (EIF) for the smooth lower transport bound
# @param A binary source indicator
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_transport_lower <- function(S, A, Y, mu0, mu1, phi, pi, threshold, smoothness) {
  w0     <- s_gt((1 - pi) * phi, threshold, smoothness)
  w0_dot <- s_gt_dot((1 - pi) * phi, threshold, smoothness)
  plugin <- mean(w0 * (S == 0)) / mean(S == 0)

  eif <- eif_transport_trimmed(S, A, Y, mu0, mu1, phi, pi, threshold, smoothness)

  eif <- eif - (S == 1) / mean(S == 0) * w0_dot * (1 - phi) * (A - pi)
  eif <- eif + 1 / mean(S == 0) * w0_dot * (1 - pi) * (1 - phi) * (S - phi)
  eif <- eif + (S == 0) / mean(S == 0) * (w0 - plugin)
  eif
}

# Efficient Influence Function (EIF) for the smooth lower transport bound
# @param A binary source indicator
# @param A binary treatment indicator
# @param Y binary outcome
# @param mu0 conditional mean of Y given A = 0
# @param mu1 conditional mean of Y given A = 0
# @param pi propensity score P(A = 1 | X)
# @param threshold propensity score threshold
# @param smoothness smooth approximation tuning parameter
#
# @return vector of efficient influence function values
eif_transport_upper <- function(S, A, Y, mu0, mu1, phi, pi, threshold, smoothness) {
  w1     <- s_gt(pi * phi, threshold, smoothness)
  w1_dot <- s_gt_dot(pi * phi, threshold, smoothness)
  plugin <- mean(w1 * (S == 0)) / mean(S == 0)

  eif <- eif_transport_trimmed(S, A, Y, mu0, mu1, phi, pi, threshold, smoothness)

  eif <- eif - (S == 1) / mean(S == 0) * w1_dot * (1 - phi) * (A - pi)
  eif <- eif - 1 / mean(S == 0) * w1_dot * pi * (1 - phi) * (S - phi)
  eif <- eif - (S == 0) / mean(S == 0) * (w1 - plugin)
  eif
}
