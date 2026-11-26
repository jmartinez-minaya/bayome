# Bayesian occurrence-abundance microbiome estimation (MaAsLin3-compatible)

Main entry point of bayome. The arguments follow `maaslin3::maaslin3()`
so that an existing MaAsLin3 call can be reused by changing only the
function name. Bayesian settings are optional extras.

## Usage

``` r
bayome(
  input_data,
  input_metadata = NULL,
  output = NULL,
  formula = NULL,
  fixed_effects = NULL,
  reference = NULL,
  random_effects = NULL,
  group_effects = NULL,
  ordered_effects = NULL,
  strata_effects = NULL,
  min_abundance = 0,
  min_prevalence = 0,
  max_prevalence = 1.01,
  zero_threshold = 0,
  min_variance = 0,
  max_significance = 0.1,
  normalization = c("TSS", "CLR", "NONE"),
  transform = c("LOG", "PLOG", "NONE"),
  standardize = TRUE,
  evaluate_only = NULL,
  prior = bm_prior(),
  link = c("none", "shared"),
  shrinkage = c("joint", "none", "hierarchical"),
  exposure = NULL,
  hierarchical_reference = c("mean", "median"),
  hierarchical_variance = c("vb", "plugin", "estimated"),
  int_strategy = c("auto", "ccd", "grid", "eb"),
  median_comparison_abundance = TRUE,
  median_comparison_prevalence = FALSE,
  median_comparison_abundance_threshold = 0,
  median_comparison_prevalence_threshold = 0,
  subtract_median = FALSE,
  n_draws = 4000,
  selection = c("bh_lfsr", "rope", "lfsr"),
  rope_delta = 0.5,
  save_models = FALSE,
  cores = 1,
  verbosity = "FINEST",
  ...
)
```

## Arguments

- input_data:

  Feature table (data frame or matrix). Samples may be in rows or
  columns; they are matched to `input_metadata` by name.

- input_metadata:

  Metadata data frame with samples in rows.

- output:

  Output directory for `all_results.tsv` and `significant_results.tsv`,
  or `NULL` to write nothing.

- formula:

  One-sided formula or string in MaAsLin3 notation, e.g.
  `~ diagnosis + age + ordered(stage) + (1 | subject)`.

- fixed_effects, random_effects, group_effects, ordered_effects,
  strata_effects:

  Alternative to `formula`, as in `maaslin3()`. Ignored if `formula` is
  given.

- reference:

  Reference levels, as in `maaslin3()`, e.g. `"diagnosis,nonIBD"` or
  `c("diagnosis,nonIBD", "site,A")`.

- min_abundance, min_prevalence, max_prevalence, zero_threshold,
  min_variance:

  Feature filtering thresholds, as in `maaslin3()`.

- max_significance:

  Threshold on `qval_individual` / `qval_joint` for
  `significant_results.tsv`.

- normalization, transform:

  Pre-processing, as in `maaslin3()`.

- standardize:

  Z-score continuous metadata used as fixed effects.

- evaluate_only:

  `NULL` (both parts), `"abundance"` or `"prevalence"`.

- prior:

  Priors, see
  [`bm_prior()`](https://jmartinez-minaya.github.io/bayome/reference/bm_prior.md).

- link:

  `"none"` (separate random effects per part) or `"shared"` (abundance
  random intercepts are a scaled copy of the prevalence ones).

- shrinkage:

  How information is shared across features and between the two parts:

  - `"joint"` (default): one INLA model with both parts and all
    features. Prevalence: one effect per feature, `u_P ~ N(0, 2.5^2)`
    (fixed, no pooling). Abundance: `u_A = lambda * u_P + e`, where
    `lambda` (one per covariate, `N(0, 1)` prior) links abundance to
    prevalence and `e ~ N(0, 1/tau)` is shared across features (PC prior
    on `tau`). This structure is used for the `exposure` only; the other
    covariates (adjustment covariates) get one effect per feature and
    part with a fixed `N(0, 2.5^2)` prior (see
    `bm_prior(adjustment_sd = )`): no pooling, no link and no extra
    hyperparameters. Formula random effects get one copy per part,
    replicated per feature. Residual variances by `"vb"`; compositional
    reference: the median over features, on joint posterior samples.
    Features without both presences and absences only enter the
    abundance part. With `evaluate_only`, `"hierarchical"` is used
    instead. Does not support
    [`ordered()`](https://rdrr.io/r/base/factor.html) yet.

  - `"none"`: each feature on its own, with the median comparison on
    posterior draws for compositionality.

  - `"hierarchical"`: one INLA model per part with all features. Without
    median comparison, `beta_jk ~ N(0, tau_k^2)`; with it,
    `beta_jk = mu_k + u_jk`, `u_jk ~ N(0, tau_k^2)`, with a PC prior on
    `tau_k`; see `hierarchical_reference`. Random effects are replicated
    per feature with shared variances. The Gaussian prior has no point
    mass at the reference, so `lfsr` controls sign errors, not point
    nulls. Does not support
    [`ordered()`](https://rdrr.io/r/base/factor.html) or
    `link = "shared"` yet; thresholds are not used.

- exposure:

  With `shrinkage = "joint"`, the variable(s) of interest: names of
  fixed-effect variables of the formula (for example `"diagnosis"`).
  Their coefficients get the full joint structure (link `lambda` to
  prevalence and information shared across features, two hyperparameters
  per coefficient); the other covariates are adjustment covariates with
  one effect per feature. Default: the first variable of the formula. To
  give every covariate the full structure, list them all. Not used by
  the other `shrinkage` options.

- hierarchical_reference:

  Compositional reference used by `shrinkage = "hierarchical"` under
  median comparison:

  - `"mean"` (default): the latent common centre `mu_k`. The corrected
    effect is `u_jk = beta_jk - mu_k`, whose posterior marginal INLA
    gives directly.

  - `"median"`: the median over features of `beta_jk`, computed on joint
    posterior samples. More robust when the true effects go mostly in
    one direction, which pulls `mu_k` away from the features with no
    effect.

- hierarchical_variance:

  Residual variances of the abundance part with
  `shrinkage = "hierarchical"`, one per feature:

  - `"vb"` (default): mean-field variational updates alternating with
    INLA. INLA fits the model given the current variances; each variance
    then gets a conjugate update from the INLA residuals (including
    their uncertainty and the random effects), with a common prior on
    the variances estimated from all features as in limma. Repeats until
    the precisions change by less than 1% (at most 10 INLA fits). No
    [`lm()`](https://rdrr.io/r/stats/lm.html).

  - `"plugin"`: estimated with [`lm()`](https://rdrr.io/r/stats/lm.html)
    per feature (without the random effects) and passed to INLA as fixed
    weights; INLA only estimates a common multiplier. Fastest, but with
    random effects the relative variances are less accurate.

  - `"estimated"`: one Gaussian likelihood per feature, so INLA
    estimates every precision with its uncertainty. Adds one
    hyperparameter per feature, which makes the fit much slower as the
    number of features grows; needs `int_strategy = "eb"` with many
    features.

- int_strategy:

  How INLA integrates over the hyperparameters (INLA's `int.strategy`):
  `"auto"`, `"ccd"`, `"grid"` or `"eb"` (empirical Bayes: only the
  posterior mode of the hyperparameters, much faster with many of them).
  With `shrinkage = "joint"`, `"auto"` means `"eb"`: the model has
  several hyperparameters per coefficient, well estimated from all
  features, and integrating over them (`"ccd"`) changes the results
  little but makes the final fit several times slower. Otherwise it
  applies to the hierarchical fits and to
  `hierarchical_variance = "estimated"`.

- median_comparison_abundance, median_comparison_prevalence:

  Compare each coefficient with the median of that coefficient across
  features instead of with zero, as in `maaslin3()`. This corrects for
  compositionality: with relative abundances, when some features
  increase the others appear to decrease. On by default for abundance,
  off for prevalence (as in MaAsLin3). The comparison is done on
  posterior draws, so `lfsr`, `prob_positive` and the q-values refer to
  `coef - median`; see Details.

- median_comparison_abundance_threshold,
  median_comparison_prevalence_threshold:

  Minimum distance from the median for an effect to count, on the
  coefficient scale. With threshold `t`,
  `lfsr = min(P(delta <= t), P(delta >= -t))` where
  `delta = coef - median`.

- subtract_median:

  If `TRUE`, `coef`, `stderr` and the credible interval describe
  `coef - median` instead of the raw coefficient.

- n_draws:

  Number of posterior draws per coefficient used for the median
  comparison (at most 1000 joint samples for `"joint"` and
  `"hierarchical"`).

- selection:

  Rule turning posteriors into `qval_individual` and `qval_joint` (and
  `significant_results.tsv`):

  - `"bh_lfsr"` (default): Benjamini-Hochberg on `2 * lfsr` within each
    model part; joint as MaAsLin3, `pbeta(min(p), 1, 2)` then BH. It
    kept the false discovery rate at the nominal level in simulations.

  - `"rope"`: probability that the corrected effect is practically null,
    `P(|effect| <= rope_delta)`, from a normal approximation of the
    posterior; q-value = running mean of the sorted probabilities
    (Bayesian FDR for relevant effects). Too liberal for features whose
    effects are not shrunk.

  - `"lfsr"`: running mean of the sorted `lfsr` (Bayesian false sign
    rate). Too liberal when many features have no effect.

- rope_delta:

  Half-width of the region of practical equivalence for
  `selection = "rope"`, on the coefficient scale (0.5 on log2 is a
  1.4-fold change).

- save_models:

  Keep the INLA fit of every feature in the result.

- cores:

  Number of threads used by INLA.

- verbosity:

  `"FINEST"`, `"DEBUG"` or `"INFO"` print progress (with
  `shrinkage = "joint"`, also each `"vb"` iteration with its time, the
  final fit and the posterior sampling); any other value, such as
  `"WARN"`, is quiet.

- ...:

  Other `maaslin3()` arguments. They are accepted for compatibility but
  not used yet; a warning lists them.

## Value

An object of class `bayome`: a list with `results` (one row per feature,
coefficient and model part), `fits` (INLA fits if `save_models = TRUE`),
`lambda` (the estimated prevalence-abundance links, one per coefficient,
with `shrinkage = "joint"`; otherwise `NULL`), the processed `data` and
`metadata`, the INLA `model` specification and the `settings`. Results
columns follow MaAsLin3 (`feature`, `metadata`, `value`, `name`, `coef`,
`stderr`, `qval_individual`, `qval_joint`, `error`, `model`, `N`,
`N_not_zero`): `coef` and `stderr` are the posterior mean and standard
deviation, and the q-values follow `selection` (Benjamini-Hochberg on
`2 * lfsr` by default). With `shrinkage = "joint"` and
`save_models = TRUE`, the fit is in `fits$joint`, with the estimated
links in `attr(fits$joint, "bm_lambda")`. Extra columns: `lower_95`,
`upper_95` (95% credible interval), `prob_positive` (posterior
probability that the effect is positive), `lfsr` and `lfsr_joint` (local
false sign rates) and `null_hypothesis` (the posterior mean of the
median used as reference; 0 when the coefficient was compared with
zero).

## Details

bayome fits a two-part (hurdle) model:

- **prevalence**: presence/absence, binomial likelihood with logit link,
  using all samples;

- **abundance**: transformed abundance (log2 of the normalised value by
  default), Gaussian likelihood, using only the samples where the
  feature is present.

By default (`shrinkage = "joint"`) all features and both parts are
fitted in one INLA model. Prevalence effects are estimated feature by
feature; abundance effects are shrunk towards what each feature's
prevalence effect suggests, `u_A = lambda * u_P + e`, with one `lambda`
per covariate estimated from all features and `e` shared across
features.

As in MaAsLin3, the pipeline is normalise, filter, transform, then fit.

**Median comparison.** For each coefficient (for example `diagnosisCD`)
and posterior draw, the median across features is computed, and each
feature's draw is compared with it: `delta = beta - median`. Features
are fitted separately, so their posteriors are independent and the draws
can be combined directly; the posterior of `delta` therefore includes
the uncertainty of the median, and the fact that each feature
contributes to it. Rows with an `error`, and coefficients present in
fewer than three features, do not take part. Draws are random: call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible results.

## References

Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA,
Thompson KN, Nearing JT, Huttenhower C (2026). MaAsLin 3: refining and
extending generalized multivariable linear models for meta-omic
association discovery. Nature Methods, 23, 554-564.
[doi:10.1038/s41592-025-02923-9](https://doi.org/10.1038/s41592-025-02923-9)
