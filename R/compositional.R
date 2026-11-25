# Compositional correction: Bayesian version of maaslin3's median comparison.
#
# With relative abundances (e.g. TSS), when some features truly increase, all
# others appear to decrease. maaslin3 therefore tests each coefficient against
# the median of that coefficient across features instead of against zero.
#
# Here the comparison is done on the posterior. For coefficient k (e.g.
# diagnosisCD) and posterior draw s:
#   m_k^(s)      = median over features j of beta_jk^(s)
#   delta_jk^(s) = beta_jk^(s) - m_k^(s)
# Features are fitted separately, so their posteriors are independent and
# the draws can be combined directly. The posterior of delta_jk accounts for
# the uncertainty of the median, including the fact that feature j itself
# contributes to it (maaslin3 approximates this with normal simulations and
# a covariance adjustment).
#
# The local false sign rate is then computed for delta instead of beta. With
# a threshold t > 0 an effect only counts if it is beyond the median by more
# than t:  lfsr = min(P(delta <= t), P(delta >= -t)).

.bm_median_comparison <- function(res, draws, part, threshold = 0,
                                  subtract_median = FALSE, min_features = 3) {
  if (!"null_hypothesis" %in% names(res)) res$null_hypothesis <- 0
  rows_part <- which(res$model == part)
  for (nm in unique(res$name[rows_part])) {
    idx <- rows_part[res$name[rows_part] == nm]
    has_draws <- !is.na(draws[idx, 1])
    use <- idx[has_draws & is.na(res$error[idx])]
    if (length(use) < min_features) next

    m <- apply(draws[use, , drop = FALSE], 2, stats::median)
    res$null_hypothesis[idx] <- mean(m)

    for (i in use) {
      d <- draws[i, ] - m
      res$prob_positive[i] <- mean(d > 0)
      res$lfsr[i] <- min(mean(d <= threshold), mean(d >= -threshold))
      if (subtract_median) {
        res$coef[i] <- mean(d)
        res$stderr[i] <- stats::sd(d)
        res$lower_95[i] <- stats::quantile(d, 0.025, names = FALSE)
        res$upper_95[i] <- stats::quantile(d, 0.975, names = FALSE)
      }
    }
  }
  res
}
