# Translate a parsed formula into an INLA model

Creates the INLA formula and the index columns it needs. Random effects,
[`ordered()`](https://rdrr.io/r/base/factor.html) and `strata()` terms
become `f()` terms with integer indices added to `data`.

## Usage

``` r
bm_inla_formula(formula, data, response = "y", prior = bm_prior())
```

## Arguments

- formula:

  A formula, a character string, or a `bm_formula` object from
  [`bm_parse_formula()`](https://jmartinez-minaya.github.io/bayome/reference/bm_parse_formula.md).

- data:

  Metadata data frame (samples in rows).

- response:

  Name of the response variable in the INLA formula.

- prior:

  A list of priors, as returned by
  [`bm_prior()`](https://jmartinez-minaya.github.io/bayome/reference/bm_prior.md).

## Value

An object of class `bm_inla_model`: a list with `formula` (the INLA
formula, whose environment holds the hyperprior lists), `data` (the
metadata with the added index columns), `parsed` (the `bm_formula`) and
`terms` (a data frame describing each term).

## Examples

``` r
meta <- data.frame(diagnosis = factor(c("CD", "UC", "nonIBD", "CD")),
                   age = c(20, 35, 50, 41),
                   subject = c("a", "a", "b", "b"))
m <- bm_inla_formula(~ diagnosis + age + (1 | subject), meta)
m$formula
#> y ~ 1 + diagnosis + age + f(.bm_re_subject, model = "iid", hyper = .bm_hyper$random)
#> <environment: 0x5cd476332448>
```
