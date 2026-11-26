# bayome

<!-- badges: start -->
[![R-CMD-check](https://github.com/jmartinez-minaya/bayome/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jmartinez-minaya/bayome/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/jmartinez-minaya/bayome/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/jmartinez-minaya/bayome/actions/workflows/test-coverage.yaml)
[![lint](https://github.com/jmartinez-minaya/bayome/actions/workflows/lint.yaml/badge.svg)](https://github.com/jmartinez-minaya/bayome/actions/workflows/lint.yaml)
[![codecov](https://codecov.io/gh/jmartinez-minaya/bayome/graph/badge.svg)](https://app.codecov.io/gh/jmartinez-minaya/bayome)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**BAY**esian **O**ccurrence–abundance **M**icrobiome **E**stimation

bayome finds associations between microbial features (taxa, genes,
pathways) and sample metadata. For each feature it models the
**occurrence** (presence/absence) and the **abundance** together, links the
abundance effects to the occurrence effects, and **shares information
across features**, so that features with few observations still get stable
estimates and honest intervals. Every model is fitted with the integrated
nested Laplace approximation ([R-INLA](https://www.r-inla.org/)).

The arguments, formulas and output columns follow
[MaAsLin3](https://github.com/biobakery/maaslin3), so an existing MaAsLin3
analysis runs by changing only the function name.

Website and worked example: <https://jmartinez-minaya.github.io/bayome/>

> **Status:** early development (version 0.0.0.9000). The default model
> (`shrinkage = "joint"`) and the default selection rule
> (`selection = "bh_lfsr"`) work and are tested on simulated data, not yet
> on real data.

## Installation

INLA is not on CRAN. Install it first:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

Then install bayome from GitHub:

```r
# install.packages("remotes")
remotes::install_github("jmartinez-minaya/bayome")
```

## Quick start

```r
library(bayome)

fit <- bayome(
  input_data     = taxa,      # samples x features (or features x samples)
  input_metadata = meta,      # samples in rows
  output         = "bayome_output",
  formula        = ~ diagnosis + age + (1 | subject),
  exposure       = "diagnosis", # variable of interest; age is an adjustment
  reference      = "diagnosis,nonIBD"
)

fit
head(fit$results)
```

`output` gets `all_results.tsv` and `significant_results.tsv`. Use
`output = NULL` to keep the results in R only.

Plots (ggplot2 objects):

```r
bm_plot_summary(fit)                              # top features, all coefficients
bm_plot_effects(fit, "diagnosisCD")               # effects with credible intervals
bm_plot_link(fit, "diagnosisCD")                  # prevalence vs abundance effects
bm_plot_association(fit, "Faecalibacterium_prausnitzii", "diagnosis")  # the data
```

The vignette ("Getting started with bayome") runs the whole workflow on the
HMP2 inflammatory bowel disease data included in the package.

## The model

For each feature *j* and sample *i*, with *Y*<sub>ij</sub> the normalised
abundance:

- **Prevalence** (all samples):
  logit P(*Y*<sub>ij</sub> > 0) = **x**<sub>i</sub>' **β**<sub>j</sub><sup>prev</sup> + random effects
- **Abundance** (samples where the feature is present):
  log<sub>2</sub> *Y*<sub>ij</sub> = **x**<sub>i</sub>' **β**<sub>j</sub><sup>abun</sup> + random effects + ε<sub>ij</sub>

Both parts are fitted together in one INLA call, with a binomial and a
Gaussian likelihood. Pre-processing: normalise (`TSS`, `CLR` or `NONE`),
filter, then transform (`LOG`, `PLOG` or `NONE`). Continuous covariates are
standardised by default.

### The default model (`shrinkage = "joint"`)

All features and both parts go into **one** INLA model. The model
distinguishes the **exposure** (`exposure = "diagnosis"`: the variable of
interest; by default the first variable of the formula) from the
**adjustment covariates** (all the others). For an exposure coefficient *k*
(for example `diagnosisCD`), the effect on feature *j* is split as:

- **Prevalence:** β<sub>jk</sub><sup>P</sup> = μ<sub>k</sub><sup>P</sup> + *u*<sub>jk</sub><sup>P</sup>,
  with *u*<sub>jk</sub><sup>P</sup> ~ N(0, 2.5²) fixed: each feature is
  estimated from its own presence data, without pooling.
- **Abundance:** β<sub>jk</sub><sup>A</sup> = μ<sub>k</sub><sup>A</sup> + **λ<sub>k</sub> · *u*<sub>jk</sub><sup>P</sup>** + *e*<sub>jk</sub>,
  with *e*<sub>jk</sub> ~ N(0, 1/τ<sub>k</sub>) shared across features.

**λ<sub>k</sub>** is one number per exposure coefficient, estimated from all features:
how much abundance tends to follow prevalence in these data. A feature with
few non-zero samples has its abundance effect shrunk towards what its
prevalence effect suggests, instead of towards zero. When prevalence and
abundance are unrelated, λ<sub>k</sub> ≈ 0 and the model reduces to an
independent hierarchical abundance model.

Each **adjustment covariate** *l* (for example age or antibiotic use) gets
one effect per feature and part, γ<sub>jl</sub> ~ N(0, 2.5²) with a fixed
prior: no pooling, no link and no extra hyperparameters. So the model has two
hyperparameters per exposure coefficient (λ<sub>k</sub> and τ<sub>k</sub>)
plus the variances of the random effects, however many adjustment
covariates there are. To give every covariate the full structure, list them
all in `exposure`.

Random effects of the formula (for example `(1 | subject)`) get one copy
per part, replicated per feature with a shared variance. Residual variances
are one per feature, estimated with `"vb"` (variational updates with a
common prior). The abundance effects are compared with their median over
features (compositional correction). The estimated λ values are printed
during the fit and stored in `fit$lambda`.

The alternatives are `shrinkage = "hierarchical"` (each part separately,
with a common prior for each) and `shrinkage = "none"` (each feature on its
own).

Priors are set with `bm_prior()`:

| Component | Default prior |
|---|---|
| Fixed effects, prevalence | Normal(0, 2.5²) on the log-odds scale |
| Fixed effects, abundance | Normal(0, 10²) on the log<sub>2</sub> scale |
| Random-effect SDs | PC prior, P(σ > 1) = 0.01 |
| `ordered()` random walk | PC prior, P(σ > 1) = 0.01 |
| Between-feature SD (`"joint"`: abundance *e*; `"hierarchical"`: both parts) | PC prior, P(σ > 1) = 0.01 (set with `bm_prior(random = )`) |
| Prevalence effect per feature *u*<sup>P</sup> (`"joint"`) | Normal(0, 2.5²), fixed (no pooling; set with `bm_prior(fixed_sd_prevalence = )`) |
| Prevalence–abundance link λ (`"joint"`) | Normal(0, 1) |
| Adjustment covariate effect per feature γ (`"joint"`) | Normal(0, 2.5²), fixed, both parts (set with `bm_prior(adjustment_sd = )`) |
| `strata()` | iid effect with fixed, very small precision |

### Compositional correction

With relative abundances (TSS), when some features truly increase, all the
others appear to decrease. Each abundance coefficient is therefore compared
with the **median of that coefficient across features** instead of with
zero (`median_comparison_abundance = TRUE`, the default).

bayome does this on the posterior. For each coefficient and posterior
draw *s*:

- *m*<sup>(s)</sup> = median over features *j* of β<sub>j</sub><sup>(s)</sup>
- δ<sub>j</sub><sup>(s)</sup> = β<sub>j</sub><sup>(s)</sup> − *m*<sup>(s)</sup>

`lfsr`, `prob_positive` and the q-values then refer to δ, whose posterior
includes the uncertainty of the median. `coef` keeps the raw coefficient
unless `subtract_median = TRUE`. `null_hypothesis` reports the median used.
Draws are random: call `set.seed()` first for reproducible results.

### Sharing information across features

| `shrinkage` | Prevalence | Abundance | Link between parts |
|---|---|---|---|
| `"joint"` (default) | Each feature on its own (fixed N(0, 2.5²) prior) | Common prior across features, centred on λ·*u*<sup>P</sup> | λ per covariate |
| `"hierarchical"` | Common prior across features, centred at 0 | Common prior across features (`hierarchical_reference`: μ or median) | None |
| `"none"` | Each feature on its own | Each feature on its own | None (or shared random intercepts with `link = "shared"`) |

With `"hierarchical"`, `hierarchical_variance` sets the residual variances
(`"vb"` by default, `"plugin"` from an `lm()` per feature, or `"estimated"`
inside INLA, which is much slower).

## Formula notation

Formulas are written in lme4-style notation and translated into INLA
internally. `bm_inla_formula()` shows the result of the translation.

| You write | Meaning | INLA term |
|---|---|---|
| `diagnosis + age` | Fixed effects | unchanged |
| `group(site)` | Categorical variable tested as a whole | fixed factor (joint inference: planned) |
| `ordered(stage)` | Ordered categorical; each level is compared with the previous one | `f(stage, model = "rw1")` (only with `shrinkage = "none"`) |
| `strata(pair)` | Matched designs (only one allowed) | `f(pair, model = "iid")` with fixed small precision |
| `(1 \| subject)` | Random intercept | `f(subject, model = "iid")` |
| `(age \|\| subject)` | Uncorrelated random intercept and slope | two `iid` terms |
| `(age \| subject)` | Correlated random intercept and slope | not supported yet |

The arguments `fixed_effects`, `random_effects`, `group_effects`,
`ordered_effects` and `strata_effects` also work when `formula` is not given.

```r
meta <- data.frame(diagnosis = factor(c("CD", "UC", "nonIBD", "CD")),
                   age = c(20, 35, 50, 41),
                   subject = c("a", "a", "b", "b"))
bm_inla_formula(~ diagnosis + age + (1 | subject), meta)
#> <bm_inla_model>
#> INLA formula:
#>   y ~ 1 + diagnosis + age + f(.bm_re_subject, model = "iid", hyper = .bm_hyper$random)
```

## Output

One row per feature, covariate level and model part:

| Column | Meaning |
|---|---|
| `feature`, `metadata`, `value`, `name` | Feature and coefficient |
| `coef` | Posterior mean |
| `stderr` | Posterior standard deviation |
| `lower_95`, `upper_95` | 95% credible interval |
| `prob_positive` | Posterior probability that the effect is positive |
| `lfsr` | Local false sign rate: probability that the sign of `coef` is wrong |
| `qval_individual` | Q-value within each coefficient and model part, from the `selection` rule |
| `lfsr_joint`, `qval_joint` | Association in at least one of the two parts |
| `null_hypothesis` | Reference the coefficient is compared with: the median across features under median comparison, otherwise 0 |
| `error` | Why a part was not fitted, or a warning about the estimate |
| `model` | `"prevalence"` or `"abundance"` |
| `N`, `N_not_zero` | Number of samples, and of samples where the feature is present |

There are no p-values.

### Selecting associations

`selection` turns the posteriors into `qval_individual` and `qval_joint`
(and `significant_results.tsv`):

| `selection` | Rule |
|---|---|
| `"bh_lfsr"` (default) | Benjamini–Hochberg on 2·lfsr within each coefficient and part; joint q-value from `pbeta(min(p), 1, 2)` then BH |
| `"rope"` | Probability that the corrected effect is practically null, P(\|effect\| ≤ `rope_delta`); q-value = running mean of the sorted probabilities |
| `"lfsr"` | Running mean of the sorted lfsr (Bayesian false sign rate) |

The lfsr is a sign-error probability. With continuous priors the posterior
puts no mass on "no effect", so for null features it behaves like a
one-sided p-value. That is why the default applies Benjamini–Hochberg to
2·lfsr, which controls the false discovery rate, while the running-mean
rules can let null features through when many features have no effect.

Every rule is applied **separately for each coefficient** (for example
`diagnosisCD`): the multiple-testing family is "which features are
associated with this coefficient". MaAsLin3 corrects all coefficients
together; then a covariate with many associations (such as sequencing depth)
lets null features of the exposure through.

## Main arguments

Data and model: `input_data`, `input_metadata`, `output`, `formula`,
`fixed_effects`, `random_effects`, `group_effects`, `ordered_effects`,
`strata_effects`, `reference`.

Filtering and pre-processing: `min_abundance`, `min_prevalence`,
`max_prevalence`, `zero_threshold`, `min_variance`, `normalization`,
`transform`, `standardize`.

Other: `max_significance`, `evaluate_only`, `save_models`, `cores`,
`verbosity`.

Bayesian settings:

- `shrinkage`: `"joint"` (default), `"hierarchical"` or `"none"`; see
  [Sharing information across features](#sharing-information-across-features).
  With `evaluate_only`, `"hierarchical"` is used instead of `"joint"`.
- `exposure`: the variable(s) of interest for `"joint"` (default: the first
  variable of the formula); the other covariates are adjustment covariates.
- `selection`: `"bh_lfsr"` (default), `"rope"` or `"lfsr"`, and
  `rope_delta` (default 0.5); see [Selecting associations](#selecting-associations).
- `prior`: priors, from `bm_prior()`.
- `hierarchical_reference`, `hierarchical_variance`, `int_strategy`: options
  of `shrinkage = "hierarchical"`.
- `link`: `"none"` (default) or `"shared"` (abundance random intercepts are
  a scaled copy of the prevalence ones; only with `shrinkage = "none"`).
- `n_draws`: posterior draws for the median comparison (default 4000; at
  most 1000 joint samples for `"joint"` and `"hierarchical"`).

Median comparison: `median_comparison_abundance`,
`median_comparison_prevalence`, `median_comparison_abundance_threshold`,
`median_comparison_prevalence_threshold` and `subtract_median`.

Unsupported arguments, such as plotting or `augment`, are accepted but
ignored, with a warning.

## Roadmap

- [x] Familiar interface and formula translation to INLA
- [x] Two-part hurdle model in a single INLA fit
- [x] Compositional correction: median comparison on posterior draws
- [x] Information sharing across features, with per-feature residual
      variances estimated by `"vb"`
- [x] Joint prevalence–abundance model with an estimated link λ
      (`"joint"`, default)
- [x] Selection rule that controls the FDR (BH on 2·lfsr, default)
- [ ] Fully Bayesian selection with a spike-and-slab prior (posterior
      probability of no effect, less bias on large effects)
- [ ] Faster fits with subject effects
- [ ] `ordered()` with `"joint"` and `"hierarchical"`
- [ ] Evaluation on real data
- [ ] Joint inference for `group()` variables
- [ ] Correlated random effects `(x | g)` (`iid2d` / `iidkd`)
- [ ] Count likelihoods with sequencing depth as an offset
- [x] Plots
- [x] Vignette with the HMP2 example data

## Acknowledgements

The interface, the pre-processing steps, the two-part structure, the median
comparison and the output format follow
[MaAsLin3](https://github.com/biobakery/maaslin3) (Nickols, Nearing et al.,
Huttenhower lab), distributed under the MIT license. bayome is an
independent package, not part of MaAsLin3 or bioBakery. If you use bayome,
please also cite MaAsLin3:

> Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA, Thompson
> KN, Nearing JT, Huttenhower C. MaAsLin 3: refining and extending
> generalized multivariable linear models for meta-omic association
> discovery. *Nature Methods* 23, 554–564 (2026). doi:
> [10.1038/s41592-025-02923-9](https://doi.org/10.1038/s41592-025-02923-9)

## Author

Joaquín Martínez-Minaya

## License

MIT. See [LICENSE.md](LICENSE.md).
