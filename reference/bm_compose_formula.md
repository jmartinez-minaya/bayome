# Build a MaAsLin3-style formula from separate effect arguments

Mirrors how `maaslin3()` combines `fixed_effects`, `random_effects`,
`group_effects`, `ordered_effects` and `strata_effects` when no
`formula` is given.

## Usage

``` r
bm_compose_formula(
  fixed_effects = NULL,
  random_effects = NULL,
  group_effects = NULL,
  ordered_effects = NULL,
  strata_effects = NULL
)
```

## Arguments

- fixed_effects, random_effects, group_effects, ordered_effects,
  strata_effects:

  Character vectors of metadata column names, or `NULL`.

## Value

A one-sided formula.

## Examples

``` r
bm_compose_formula(fixed_effects = c("diagnosis", "age"),
                   random_effects = "subject")
#> ~diagnosis + age + (1 | subject)
#> <environment: 0x5cd477cb2c80>
```
