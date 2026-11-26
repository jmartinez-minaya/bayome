#' Priors for bayome models
#'
#' Collects the priors used by bayome. Precisions of random effects and of
#' `ordered()` random walks get penalised-complexity (PC) priors, set so that
#' `P(sd > u) = alpha`. Fixed effects get Gaussian priors with mean 0.
#'
#' @param random PC prior for random-effect standard deviations, as
#'   `c(u, alpha)`.
#' @param ordered PC prior for the standard deviation of the `rw1` effect
#'   used for `ordered()` variables, as `c(u, alpha)`.
#' @param strata_log_prec Fixed log-precision for `strata()` effects. A very
#'   small precision (large variance) makes each stratum act as its own free
#'   intercept, which mimics conditioning on the strata.
#' @param fixed_sd_prevalence Prior standard deviation of the fixed effects
#'   (log-odds scale) in the prevalence model. The default, 2.5, is weakly
#'   informative for standardised covariates and keeps estimates finite under
#'   complete separation (for example a feature present in every case).
#' @param fixed_sd_abundance Prior standard deviation of the fixed effects
#'   (log2 abundance scale) in the abundance model.
#' @param adjustment_sd Prior standard deviation of the effect of each
#'   adjustment covariate on each feature (the covariates that are not the
#'   `exposure` of [bayome()] with `shrinkage = "joint"`), in both parts:
#'   log-odds for prevalence, log2 for abundance. These effects are not
#'   shared across features, so the prior is fixed; the default, 2.5, is
#'   weakly informative for standardised covariates.
#'
#' @return An object of class `bm_prior`: a list with `hyper` (the `hyper`
#'   lists passed to INLA's `f()`) and `fixed` (prior precisions of the
#'   fixed effects for the `prevalence` and `abundance` parts, and of the
#'   per-feature `adjustment` effects).
#' @examples
#' bm_prior(random = c(0.5, 0.01))
#' @export
bm_prior <- function(random = c(1, 0.01), ordered = c(1, 0.01),
                     strata_log_prec = -10, fixed_sd_prevalence = 2.5,
                     fixed_sd_abundance = 10, adjustment_sd = 2.5) {
  check_pc <- function(p, name) {
    if (!is.numeric(p) || length(p) != 2 || p[1] <= 0 || p[2] <= 0 || p[2] >= 1) {
      stop(sprintf("`%s` must be c(u, alpha) with u > 0 and 0 < alpha < 1.", name),
           call. = FALSE)
    }
  }
  check_pc(random, "random")
  check_pc(ordered, "ordered")
  if (!is.numeric(fixed_sd_prevalence) || length(fixed_sd_prevalence) != 1 ||
      fixed_sd_prevalence <= 0 || !is.numeric(fixed_sd_abundance) ||
      length(fixed_sd_abundance) != 1 || fixed_sd_abundance <= 0 ||
      !is.numeric(adjustment_sd) || length(adjustment_sd) != 1 ||
      adjustment_sd <= 0) {
    stop("Fixed-effect prior standard deviations must be positive numbers.",
         call. = FALSE)
  }

  structure(
    list(
      hyper = list(
        random = list(prec = list(prior = "pc.prec", param = random)),
        ordered = list(prec = list(prior = "pc.prec", param = ordered)),
        strata = list(prec = list(initial = strata_log_prec, fixed = TRUE))
      ),
      fixed = list(prevalence = 1 / fixed_sd_prevalence^2,
                   abundance = 1 / fixed_sd_abundance^2,
                   adjustment = 1 / adjustment_sd^2)
    ),
    class = "bm_prior"
  )
}
