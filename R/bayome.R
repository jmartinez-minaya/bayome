#' Bayesian occurrence-abundance microbiome estimation (MaAsLin3-compatible)
#'
#' Main entry point of bayome. The arguments follow `maaslin3::maaslin3()`
#' so that an existing MaAsLin3 call can be reused by changing only the
#' function name. Bayesian settings are optional extras.
#'
#' bayome fits a two-part (hurdle) model:
#' * **prevalence**: presence/absence, binomial likelihood with logit link,
#'   using all samples;
#' * **abundance**: transformed abundance (log2 of the normalised value by
#'   default), Gaussian likelihood, using only the samples where the feature
#'   is present.
#'
#' By default (`shrinkage = "joint"`) all features and both parts are fitted
#' in one INLA model. Prevalence effects are estimated feature by feature;
#' abundance effects are shrunk towards what each feature's prevalence effect
#' suggests, `u_A = lambda * u_P + e`, with one `lambda` per covariate
#' estimated from all features and `e` shared across features.
#'
#' As in MaAsLin3, the pipeline is normalise, filter, transform, then fit.
#'
#' @param input_data Feature table (data frame or matrix). Samples may be in
#'   rows or columns; they are matched to `input_metadata` by name.
#' @param input_metadata Metadata data frame with samples in rows.
#' @param output Output directory for `all_results.tsv` and
#'   `significant_results.tsv`, or `NULL` to write nothing.
#' @param formula One-sided formula or string in MaAsLin3 notation, e.g.
#'   `~ diagnosis + age + ordered(stage) + (1 | subject)`.
#' @param fixed_effects,random_effects,group_effects,ordered_effects,strata_effects
#'   Alternative to `formula`, as in `maaslin3()`. Ignored if `formula` is
#'   given.
#' @param reference Reference levels, as in `maaslin3()`, e.g.
#'   `"diagnosis,nonIBD"` or `c("diagnosis,nonIBD", "site,A")`.
#' @param min_abundance,min_prevalence,max_prevalence,zero_threshold,min_variance
#'   Feature filtering thresholds, as in `maaslin3()`.
#' @param max_significance Threshold on `qval_individual` / `qval_joint` for
#'   `significant_results.tsv`.
#' @param normalization,transform Pre-processing, as in `maaslin3()`.
#' @param standardize Z-score continuous metadata used as fixed effects.
#' @param evaluate_only `NULL` (both parts), `"abundance"` or
#'   `"prevalence"`.
#' @param prior Priors, see [bm_prior()].
#' @param link `"none"` (separate random effects per part) or `"shared"`
#'   (abundance random intercepts are a scaled copy of the prevalence ones).
#' @param shrinkage How information is shared across features and between
#'   the two parts:
#'   * `"joint"` (default): one INLA model with both parts and all features.
#'     Prevalence: one effect per feature, `u_P ~ N(0, 2.5^2)` (fixed, no
#'     pooling). Abundance: `u_A = lambda * u_P + e`, where `lambda` (one per
#'     covariate, `N(0, 1)` prior) links abundance to prevalence and
#'     `e ~ N(0, 1/tau)` is shared across features (PC prior on `tau`).
#'     Formula random effects get one copy per part, replicated per feature.
#'     Residual variances by `"vb"`; compositional reference: the median over
#'     features, on joint posterior samples. Features without both presences
#'     and absences only enter the abundance part. With `evaluate_only`,
#'     `"hierarchical"` is used instead. Does not support `ordered()` yet.
#'   * `"none"`: each feature on its own, with the median comparison on
#'     posterior draws for compositionality.
#'   * `"hierarchical"`: one INLA model per part with all features. Without
#'     median comparison, `beta_jk ~ N(0, tau_k^2)`; with it,
#'     `beta_jk = mu_k + u_jk`, `u_jk ~ N(0, tau_k^2)`, with a PC prior on
#'     `tau_k`; see `hierarchical_reference`. Random effects are replicated
#'     per feature with shared variances. The Gaussian prior has no point
#'     mass at the reference, so `lfsr` controls sign errors, not point
#'     nulls. Does not support `ordered()` or `link = "shared"` yet;
#'     thresholds are not used.
#' @param hierarchical_reference Compositional reference used by
#'   `shrinkage = "hierarchical"` under median comparison:
#'   * `"mean"` (default): the latent common centre `mu_k`. The corrected
#'     effect is `u_jk = beta_jk - mu_k`, whose posterior marginal INLA
#'     gives directly.
#'   * `"median"`: the median over features of `beta_jk`, computed on joint
#'     posterior samples. More robust when the true effects go mostly in one
#'     direction, which pulls `mu_k` away from the features with no effect.
#' @param hierarchical_variance Residual variances of the abundance part with
#'   `shrinkage = "hierarchical"`, one per feature:
#'   * `"vb"` (default): mean-field variational updates alternating with
#'     INLA. INLA fits the model given the current variances; each variance
#'     then gets a conjugate update from the INLA residuals (including their
#'     uncertainty and the random effects), with a common prior on the
#'     variances estimated from all features as in limma. Repeats until the
#'     precisions change by less than 1% (at most 10 INLA fits). No `lm()`.
#'   * `"plugin"`: estimated with `lm()` per feature (without the random
#'     effects) and passed to INLA as fixed weights; INLA only estimates a
#'     common multiplier. Fastest, but with random effects the relative
#'     variances are less accurate.
#'   * `"estimated"`: one Gaussian likelihood per feature, so INLA estimates
#'     every precision with its uncertainty. Adds one hyperparameter per
#'     feature, which makes the fit much slower as the number of features
#'     grows; needs `int_strategy = "eb"` with many features.
#' @param int_strategy How INLA integrates over the hyperparameters when
#'   `hierarchical_variance = "estimated"` (INLA's `int.strategy`):
#'   `"auto"`, `"ccd"`, `"grid"` or `"eb"` (empirical Bayes: only the
#'   posterior mode of the hyperparameters, much faster with many of them).
#' @param median_comparison_abundance,median_comparison_prevalence Compare
#'   each coefficient with the median of that coefficient across features
#'   instead of with zero, as in `maaslin3()`. This corrects for
#'   compositionality: with relative abundances, when some features increase
#'   the others appear to decrease. On by default for abundance, off for
#'   prevalence (as in MaAsLin3). The comparison is done on posterior draws,
#'   so `lfsr`, `prob_positive` and the q-values refer to
#'   `coef - median`; see Details.
#' @param median_comparison_abundance_threshold,median_comparison_prevalence_threshold
#'   Minimum distance from the median for an effect to count, on the
#'   coefficient scale. With threshold `t`, `lfsr = min(P(delta <= t),
#'   P(delta >= -t))` where `delta = coef - median`.
#' @param subtract_median If `TRUE`, `coef`, `stderr` and the credible
#'   interval describe `coef - median` instead of the raw coefficient.
#' @param n_draws Number of posterior draws per coefficient used for the
#'   median comparison (at most 1000 joint samples for `"joint"` and
#'   `"hierarchical"`).
#' @param selection Rule turning posteriors into `qval_individual` and
#'   `qval_joint` (and `significant_results.tsv`):
#'   * `"bh_lfsr"` (default): Benjamini-Hochberg on `2 * lfsr` within each
#'     model part; joint as MaAsLin3, `pbeta(min(p), 1, 2)` then BH. It kept
#'     the false discovery rate at the nominal level in simulations.
#'   * `"rope"`: probability that the corrected effect is practically null,
#'     `P(|effect| <= rope_delta)`, from a normal approximation of the
#'     posterior; q-value = running mean of the sorted probabilities
#'     (Bayesian FDR for relevant effects). Too liberal for features whose
#'     effects are not shrunk.
#'   * `"lfsr"`: running mean of the sorted `lfsr` (Bayesian false sign
#'     rate). Too liberal when many features have no effect.
#' @param rope_delta Half-width of the region of practical equivalence for
#'   `selection = "rope"`, on the coefficient scale (0.5 on log2 is a
#'   1.4-fold change).
#' @param save_models Keep the INLA fit of every feature in the result.
#' @param cores Number of threads used by INLA.
#' @param verbosity `"FINEST"`, `"DEBUG"` or `"INFO"` print progress; any
#'   other value is quiet.
#' @param ... Other `maaslin3()` arguments. They are accepted for
#'   compatibility but not used yet; a warning lists them.
#'
#' @return An object of class `bayome`: a list with `results` (one row per
#'   feature, coefficient and model part), `fits` (INLA fits if
#'   `save_models = TRUE`), the processed `data` and `metadata`, the INLA
#'   `model` specification and the `settings`. Results columns follow
#'   MaAsLin3 (`feature`, `metadata`, `value`, `name`, `coef`, `stderr`,
#'   `qval_individual`, `qval_joint`, `error`, `model`, `N`, `N_not_zero`):
#'   `coef` and `stderr` are the posterior mean and standard deviation, and
#'   the q-values follow `selection` (Benjamini-Hochberg on `2 * lfsr` by
#'   default). With `shrinkage = "joint"` and `save_models = TRUE`, the fit is
#'   in `fits$joint`, with the estimated links in `attr(fits$joint, "bm_lambda")`. Extra
#'   columns: `lower_95`, `upper_95` (95% credible interval),
#'   `prob_positive` (posterior probability that the effect is positive),
#'   `lfsr` and `lfsr_joint` (local false sign rates) and `null_hypothesis`
#'   (the posterior mean of the median used as reference; 0 when the
#'   coefficient was compared with zero).
#'
#' @details
#' **Median comparison.** For each coefficient (for example `diagnosisCD`)
#' and posterior draw, the median across features is computed, and each
#' feature's draw is compared with it: `delta = beta - median`. Features are
#' fitted separately, so their posteriors are independent and the draws can
#' be combined directly; the posterior of `delta` therefore includes the
#' uncertainty of the median, and the fact that each feature contributes to
#' it. Rows with an `error`, and coefficients present in fewer than three
#' features, do not take part. Draws are random: call `set.seed()` first for
#' reproducible results.
#' @references
#' Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA,
#' Thompson KN, Nearing JT, Huttenhower C (2026). MaAsLin 3: refining and
#' extending generalized multivariable linear models for meta-omic
#' association discovery. Nature Methods, 23, 554-564.
#' \doi{10.1038/s41592-025-02923-9}
#' @export
bayome <- function(input_data,
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
                     ...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("bayome needs the INLA package: install.packages(\"INLA\", ",
         "repos = c(getOption(\"repos\"), ",
         "INLA = \"https://inla.r-inla-download.org/R/stable\"))",
         call. = FALSE)
  }
  normalization <- match.arg(normalization)
  transform <- match.arg(transform)
  link <- match.arg(link)
  shrinkage <- match.arg(shrinkage)
  selection <- match.arg(selection)
  hierarchical_reference <- match.arg(hierarchical_reference)
  hierarchical_variance <- match.arg(hierarchical_variance)
  int_strategy <- match.arg(int_strategy)
  parts <- if (is.null(evaluate_only)) .bm_parts else
    match.arg(evaluate_only, .bm_parts)
  if (shrinkage == "joint" && !is.null(evaluate_only)) {
    message("shrinkage = \"joint\" models both parts together; with evaluate_only, ",
            "shrinkage = \"hierarchical\" is used instead.")
    shrinkage <- "hierarchical"
  }
  if (!inherits(prior, "bm_prior")) {
    stop("`prior` must be created with bm_prior().", call. = FALSE)
  }
  extra <- list(...)
  if (length(extra)) {
    warning(sprintf("These maaslin3() arguments are not supported yet and were ignored: %s",
                    paste(names(extra), collapse = ", ")), call. = FALSE)
  }
  if (is.null(input_metadata)) {
    stop("`input_metadata` is required.", call. = FALSE)
  }
  say <- function(...) {
    if (verbosity %in% c("FINEST", "DEBUG", "INFO")) message(...)
  }

  if (is.null(formula)) {
    formula <- bm_compose_formula(fixed_effects, random_effects, group_effects,
                                  ordered_effects, strata_effects)
  }
  parsed <- bm_parse_formula(formula)

  aligned <- .bm_align_samples(input_data, input_metadata)
  metadata <- .bm_apply_reference(aligned$metadata, reference)
  metadata <- .bm_drop_incomplete(metadata, parsed, say)
  features <- aligned$data[rownames(metadata), , drop = FALSE]

  features <- .bm_normalize(features, normalization, zero_threshold)
  n_before <- ncol(features)
  features <- .bm_filter(features, min_abundance, min_prevalence,
                         max_prevalence, min_variance)
  say(sprintf("%d of %d features kept after filtering.", ncol(features), n_before))
  if (!ncol(features)) stop("No features left after filtering.", call. = FALSE)
  features <- .bm_transform(features, transform)

  if (standardize) metadata <- .bm_standardize(metadata, parsed)
  model <- bm_inla_formula(parsed, metadata, response = "Y", prior = prior)

  # One design per set of parts actually fitted (both, or just one when a
  # feature has no absences or too few non-zero values).
  designs <- list()
  get_design <- function(p) {
    key <- paste(p, collapse = "+")
    if (is.null(designs[[key]])) {
      designs[[key]] <<- .bm_hurdle_design(model, p, prior, link)
    }
    designs[[key]]
  }
  design_full <- get_design(parts)
  n_fixed <- ncol(design_full$X)
  control <- list(num.threads = paste0(cores, ":1"), int_strategy = int_strategy,
                  vb_tol = 0.01, vb_max_iter = 10)

  compare <- c(
    prevalence = median_comparison_prevalence && "prevalence" %in% parts,
    abundance = median_comparison_abundance && "abundance" %in% parts
  )
  thresholds <- c(prevalence = median_comparison_prevalence_threshold,
                  abundance = median_comparison_abundance_threshold)
  draws_needed <- if (shrinkage == "none" && any(compare)) n_draws else 0
  na_draws <- function(r) {
    if (is.null(r) || !draws_needed) return(NULL)
    matrix(NA_real_, nrow(r), draws_needed)
  }
  skip_reasons <- function(y) {
    n_nz <- sum(!is.na(y))
    skipped <- list()
    if ("prevalence" %in% parts && (n_nz == 0 || n_nz == length(y))) {
      skipped$prevalence <- if (n_nz == 0) "Feature absent in all samples" else
        "Feature present in all samples: prevalence model not fitted"
    }
    if ("abundance" %in% parts && n_nz <= n_fixed + 1) {
      skipped$abundance <- "Too few non-zero samples for the abundance model"
    }
    skipped
  }

  fits <- list()
  if (shrinkage == "joint") {
    if (link == "shared") {
      warning("link = \"shared\" is not used with shrinkage = \"joint\".", call. = FALSE)
    }
    skips <- lapply(features, skip_reasons)
    ok_prev <- vapply(skips, function(s) is.null(s$prevalence), logical(1))
    ok_abun <- vapply(skips, function(s) is.null(s$abundance), logical(1))
    pieces <- list()
    for (j in seq_len(ncol(features))) {
      for (part in names(skips[[j]])) {
        pieces[[length(pieces) + 1]] <- .bm_empty_rows(
          design_full, part, colnames(features)[j], features[[j]], skips[[j]][[part]])
      }
    }
    if (sum(ok_prev | ok_abun) >= 2) {
      say(sprintf("Fitting the joint model for %d features (%d prevalence, %d abundance)...",
                  sum(ok_prev | ok_abun), sum(ok_prev), sum(ok_abun)))
      jf <- tryCatch(
        .bm_fit_joint(features, ok_prev, ok_abun, model, design_full, prior,
                      compare, subtract_median, control,
                      n_draws = min(n_draws, 1000)),
        error = function(e) conditionMessage(e))
      if (is.character(jf)) {
        for (j in which(ok_prev | ok_abun)) {
          for (part in parts[c(ok_prev[j], ok_abun[j])]) {
            pieces[[length(pieces) + 1]] <- .bm_empty_rows(
              design_full, part, colnames(features)[j], features[[j]],
              paste("INLA error:", jf))
          }
        }
      } else {
        pieces[[length(pieces) + 1]] <- jf$results
        if (save_models) fits$joint <- jf$fit
        say(sprintf("Prevalence-abundance link (lambda): %s",
                    paste(sprintf("%s = %.2f", names(jf$lambda), jf$lambda), collapse = ", ")))
      }
    }
    pieces <- lapply(pieces, function(p) {
      if (!is.null(p) && !"null_hypothesis" %in% names(p)) p$null_hypothesis <- 0
      p
    })
    res <- do.call(rbind, pieces)
    rownames(res) <- NULL
  } else if (shrinkage == "hierarchical") {
    if (link == "shared") {
      warning("link = \"shared\" is not used with shrinkage = \"hierarchical\".",
              call. = FALSE)
    }
    skips <- lapply(features, skip_reasons)
    pieces <- list()
    for (part in parts) {
      reason <- vapply(skips, function(s) {
        if (is.null(s[[part]])) NA_character_ else s[[part]]
      }, character(1))
      ok <- is.na(reason)
      for (j in which(!ok)) {
        pieces[[length(pieces) + 1]] <- .bm_empty_rows(
          design_full, part, colnames(features)[j], features[[j]], reason[j])
      }
      if (sum(ok) < 2) next
      say(sprintf("Fitting the hierarchical %s model for %d features...",
                  part, sum(ok)))
      hf <- tryCatch(
        .bm_fit_hierarchical_part(part, features, ok, model, design_full,
                                  prior, compare[[part]], subtract_median,
                                  control, n_draws = min(n_draws, 1000),
                                  reference = hierarchical_reference,
                                  residual_variance = hierarchical_variance),
        error = function(e) conditionMessage(e))
      if (is.character(hf)) {
        for (j in which(ok)) {
          pieces[[length(pieces) + 1]] <- .bm_empty_rows(
            design_full, part, colnames(features)[j], features[[j]],
            paste("INLA error:", hf))
        }
      } else {
        pieces[[length(pieces) + 1]] <- hf$results
        if (save_models) fits[[part]] <- hf$fit
      }
    }
    pieces <- lapply(pieces, function(p) {
      if (!is.null(p) && !"null_hypothesis" %in% names(p)) p$null_hypothesis <- 0
      p
    })
    res <- do.call(rbind, pieces)
    rownames(res) <- NULL
  } else {
    results <- list()
    all_draws <- list()
    for (j in seq_len(ncol(features))) {
      feature <- colnames(features)[j]
      y <- features[[j]]
      skipped <- skip_reasons(y)
      fit_parts <- setdiff(parts, names(skipped))

      rows <- lapply(names(skipped), function(p) {
        .bm_empty_rows(design_full, p, feature, y, skipped[[p]])
      })
      row_draws <- lapply(rows, na_draws)
      if (length(fit_parts)) {
        design <- get_design(fit_parts)
        fit <- .bm_fit_feature(y, design, control)
        if (is.character(fit)) {
          failed <- lapply(fit_parts, function(p) {
            .bm_empty_rows(design_full, p, feature, y,
                           paste("INLA error:", fit))
          })
          rows <- c(rows, failed)
          row_draws <- c(row_draws, lapply(failed, na_draws))
        } else {
          summ <- .bm_summarise_fit(fit, design, feature, y, list(),
                                    n_draws = draws_needed)
          rows <- c(rows, list(summ))
          row_draws <- c(row_draws, list(attr(summ, "draws")))
          if (save_models) fits[[feature]] <- fit
        }
      }
      results[[feature]] <- do.call(rbind, rows)
      all_draws[[feature]] <- do.call(rbind, row_draws)
      say(sprintf("Fitted feature %d/%d: %s", j, ncol(features), feature))
    }

    res <- do.call(rbind, unname(results))
    rownames(res) <- NULL
    res$null_hypothesis <- 0
    if (any(compare)) {
      draws <- do.call(rbind, unname(all_draws))
      for (part in names(compare)[compare]) {
        res <- .bm_median_comparison(res, draws, part, thresholds[[part]],
                                     subtract_median)
      }
    }
  }
  res <- .bm_add_qvalues(res, selection, rope_delta, subtract_median)
  res <- res[order(!is.na(res$error), res$qval_individual), c(
    "feature", "metadata", "value", "name", "coef", "stderr", "lower_95",
    "upper_95", "prob_positive", "lfsr", "qval_individual", "lfsr_joint",
    "qval_joint", "null_hypothesis", "error", "model", "N", "N_not_zero")]
  rownames(res) <- NULL

  if (!is.null(output)) .bm_write_results(res, output, max_significance)

  structure(
    list(
      results = res,
      fits = if (save_models) fits else NULL,
      data = features,
      metadata = metadata,
      model = model,
      settings = list(
        output = output, formula = formula,
        min_abundance = min_abundance, min_prevalence = min_prevalence,
        max_prevalence = max_prevalence, zero_threshold = zero_threshold,
        min_variance = min_variance, max_significance = max_significance,
        normalization = normalization, transform = transform,
        standardize = standardize, evaluate_only = evaluate_only,
        prior = prior, link = link, shrinkage = shrinkage,
        hierarchical_reference = hierarchical_reference,
        hierarchical_variance = hierarchical_variance, int_strategy = int_strategy,
        median_comparison_abundance = median_comparison_abundance,
        median_comparison_prevalence = median_comparison_prevalence,
        median_comparison_abundance_threshold = median_comparison_abundance_threshold,
        median_comparison_prevalence_threshold = median_comparison_prevalence_threshold,
        subtract_median = subtract_median, n_draws = n_draws,
        selection = selection, rope_delta = rope_delta,
        cores = cores, verbosity = verbosity
      )
    ),
    class = "bayome"
  )
}

#' @export
print.bayome <- function(x, ...) {
  r <- x$results
  thr <- x$settings$max_significance
  cat("<bayome>\n")
  cat(sprintf("  %d samples, %d features, formula: %s\n", nrow(x$data),
              ncol(x$data), deparse1(.bm_as_rhs_formula(x$settings$formula))))
  for (part in unique(r$model)) {
    i <- r$model == part
    cat(sprintf("  %-10s %d associations with qval_individual <= %g\n", part,
                sum(r$qval_individual[i] <= thr & is.na(r$error[i]), na.rm = TRUE),
                thr))
  }
  invisible(x)
}

# --- internal helpers -------------------------------------------------------

.bm_write_results <- function(res, output, max_significance) {
  if (!dir.exists(output)) dir.create(output, recursive = TRUE)
  utils::write.table(res, file.path(output, "all_results.tsv"), sep = "\t",
                     quote = TRUE, row.names = FALSE)
  sig <- res[is.na(res$error) &
               ((!is.na(res$qval_individual) & res$qval_individual <= max_significance) |
                  (!is.na(res$qval_joint) & res$qval_joint <= max_significance)), ]
  utils::write.table(sig, file.path(output, "significant_results.tsv"),
                     sep = "\t", quote = TRUE, row.names = FALSE)
  invisible(NULL)
}

# Drops samples with missing values in any variable used by the model.
.bm_drop_incomplete <- function(metadata, parsed, say) {
  vars <- unique(c(
    all.vars(stats::as.formula(paste("~", paste(c("1", parsed$fixed),
                                                collapse = " + ")))),
    parsed$group, parsed$ordered, parsed$strata, parsed$random$group,
    stats::na.omit(parsed$random$slope)
  ))
  vars <- intersect(vars, colnames(metadata))
  ok <- stats::complete.cases(metadata[, vars, drop = FALSE])
  if (!all(ok)) {
    say(sprintf("%d sample(s) with missing metadata were dropped.", sum(!ok)))
  }
  metadata[ok, , drop = FALSE]
}

# Puts samples in rows of the feature table and keeps only samples present in
# both tables, in the same order. Like maaslin3, the orientation is detected
# from which dimension's names match the metadata row names.
.bm_align_samples <- function(input_data, input_metadata) {
  data <- as.data.frame(input_data, check.names = FALSE)
  metadata <- as.data.frame(input_metadata)
  meta_ids <- rownames(metadata)

  n_row <- length(intersect(rownames(data), meta_ids))
  n_col <- length(intersect(colnames(data), meta_ids))
  if (n_col > n_row) {
    data <- as.data.frame(t(data), check.names = FALSE)
  }
  common <- intersect(rownames(data), meta_ids)
  if (!length(common)) {
    stop("No sample names are shared between `input_data` and `input_metadata`.",
         call. = FALSE)
  }
  dropped <- length(union(rownames(data), meta_ids)) - length(common)
  if (dropped > 0) {
    message(sprintf("%d sample(s) not present in both tables were dropped.", dropped))
  }
  list(data = data[common, , drop = FALSE],
       metadata = metadata[common, , drop = FALSE])
}

# Applies maaslin3-style reference levels: "var,level" strings, given either
# as a vector or separated by ";".
.bm_apply_reference <- function(metadata, reference) {
  if (is.null(reference) || !length(reference)) return(metadata)
  specs <- trimws(unlist(strsplit(reference, ";", fixed = TRUE)))
  specs <- specs[nzchar(specs)]
  for (s in specs) {
    parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    if (length(parts) != 2) {
      stop(sprintf("Invalid reference '%s': use 'variable,level'.", s), call. = FALSE)
    }
    var <- parts[1]
    lev <- parts[2]
    if (!var %in% colnames(metadata)) {
      stop(sprintf("Reference variable '%s' not found in metadata.", var), call. = FALSE)
    }
    x <- as.factor(metadata[[var]])
    if (!lev %in% levels(x)) {
      stop(sprintf("Level '%s' not found in variable '%s'.", lev, var), call. = FALSE)
    }
    metadata[[var]] <- stats::relevel(x, ref = lev)
  }
  metadata
}
