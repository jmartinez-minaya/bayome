# Bayesian q-values from local false sign rates

For each result, the q-value is the expected proportion of sign errors
among all results with a local false sign rate (lfsr) at most as large
(the running mean of the sorted lfsr values). Selecting the results with
`qval <= 0.1` controls the Bayesian false sign rate at 10%.

## Usage

``` r
bm_qvalue(lfsr)
```

## Arguments

- lfsr:

  Numeric vector of local false sign rates; `NA` values get `NA`.

## Value

Numeric vector of q-values, same length as `lfsr`.

## Examples

``` r
bm_qvalue(c(0.001, 0.2, 0.01, NA))
#> [1] 0.00100000 0.07033333 0.00550000         NA
```
