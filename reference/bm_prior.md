# Priors for bayome models

Collects the priors used by bayome. Precisions of random effects and of
[`ordered()`](https://rdrr.io/r/base/factor.html) random walks get
penalised-complexity (PC) priors, set so that `P(sd > u) = alpha`. Fixed
effects get Gaussian priors with mean 0.

## Usage

``` r
bm_prior(
  random = c(1, 0.01),
  ordered = c(1, 0.01),
  strata_log_prec = -10,
  fixed_sd_prevalence = 2.5,
  fixed_sd_abundance = 10,
  adjustment_sd = 2.5
)
```

## Arguments

- random:

  PC prior for random-effect standard deviations, as `c(u, alpha)`.

- ordered:

  PC prior for the standard deviation of the `rw1` effect used for
  [`ordered()`](https://rdrr.io/r/base/factor.html) variables, as
  `c(u, alpha)`.

- strata_log_prec:

  Fixed log-precision for `strata()` effects. A very small precision
  (large variance) makes each stratum act as its own free intercept,
  which mimics conditioning on the strata.

- fixed_sd_prevalence:

  Prior standard deviation of the fixed effects (log-odds scale) in the
  prevalence model. The default, 2.5, is weakly informative for
  standardised covariates and keeps estimates finite under complete
  separation (for example a feature present in every case).

- fixed_sd_abundance:

  Prior standard deviation of the fixed effects (log2 abundance scale)
  in the abundance model.

- adjustment_sd:

  Prior standard deviation of the effect of each adjustment covariate on
  each feature (the covariates that are not the `exposure` of
  [`bayome()`](https://jmartinez-minaya.github.io/bayome/reference/bayome.md)
  with `shrinkage = "joint"`), in both parts: log-odds for prevalence,
  log2 for abundance. These effects are not shared across features, so
  the prior is fixed; the default, 2.5, is weakly informative for
  standardised covariates.

## Value

An object of class `bm_prior`: a list with `hyper` (the `hyper` lists
passed to INLA's `f()`) and `fixed` (prior precisions of the fixed
effects for the `prevalence` and `abundance` parts, and of the
per-feature `adjustment` effects).

## Examples

``` r
bm_prior(random = c(0.5, 0.01))
#> $hyper
#> $hyper$random
#> $hyper$random$prec
#> $hyper$random$prec$prior
#> [1] "pc.prec"
#> 
#> $hyper$random$prec$param
#> [1] 0.50 0.01
#> 
#> 
#> 
#> $hyper$ordered
#> $hyper$ordered$prec
#> $hyper$ordered$prec$prior
#> [1] "pc.prec"
#> 
#> $hyper$ordered$prec$param
#> [1] 1.00 0.01
#> 
#> 
#> 
#> $hyper$strata
#> $hyper$strata$prec
#> $hyper$strata$prec$initial
#> [1] -10
#> 
#> $hyper$strata$prec$fixed
#> [1] TRUE
#> 
#> 
#> 
#> 
#> $fixed
#> $fixed$prevalence
#> [1] 0.16
#> 
#> $fixed$abundance
#> [1] 0.01
#> 
#> $fixed$adjustment
#> [1] 0.16
#> 
#> 
#> attr(,"class")
#> [1] "bm_prior"
```
