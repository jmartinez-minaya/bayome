# Parse a MaAsLin3-style formula

Splits a one-sided formula (as used by `maaslin3()`) into fixed effects,
the MaAsLin3 specials `group()`,
[`ordered()`](https://rdrr.io/r/base/factor.html) and `strata()`, and
lme4-style random-effect terms.

## Usage

``` r
bm_parse_formula(formula)
```

## Arguments

- formula:

  A formula or a character string such as
  `"~ diagnosis + age + (1 | subject)"`. A left-hand side, if present,
  is ignored.

## Value

An object of class `bm_formula`: a list with elements `fixed` (character
vector of fixed-effect term labels), `intercept` (logical), `group`,
`ordered`, `strata` (character vectors of variable names) and `random`
(a data frame with one row per random-effect term, with columns `group`
and `slope`; `slope` is `NA` for a random intercept).

## Examples

``` r
bm_parse_formula(~ diagnosis + age + ordered(stage) + (1 | subject))
#> <bm_formula>
#>   intercept: TRUE 
#>   fixed:    diagnosis, age
#>   group:    -
#>   ordered:  stage
#>   strata:   -
#>   random:   (1 | subject)
```
