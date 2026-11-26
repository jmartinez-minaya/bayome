sim_hurdle <- function(n = 80, J = 4, seed = 7) {
  set.seed(seed)
  meta <- data.frame(
    diagnosis = factor(sample(c("nonIBD", "CD"), n, TRUE),
                       levels = c("nonIBD", "CD")),
    age = rnorm(n, 45, 12),
    stage = factor(sample(c("low", "mid", "high"), n, TRUE),
                   levels = c("low", "mid", "high")),
    subject = rep(paste0("s", seq_len(n / 4)), each = 4),
    row.names = paste0("S", seq_len(n))
  )
  dx <- as.numeric(meta$diagnosis == "CD")
  counts <- sapply(seq_len(J), function(j) {
    bp <- if (j == 1) 2 else 0
    pres <- rbinom(n, 1, stats::plogis(0.2 + bp * dx))
    round(pres * exp(4 + rnorm(n, 0, 0.4)))
  })
  colnames(counts) <- paste0("taxon", seq_len(J))
  rownames(counts) <- rownames(meta)
  list(meta = meta, counts = counts)
}

test_that("bm_qvalue is the running mean of sorted lfsr", {
  expect_equal(bm_qvalue(c(0.2, 0, 0.1, NA)), c(0.1, 0, 0.05, NA))
  expect_equal(bm_qvalue(c(NA, NA)), c(NA_real_, NA_real_))
})

test_that("preprocessing mirrors maaslin3 (TSS, NA for zeros, log2)", {
  x <- data.frame(a = c(0, 2, 6), b = c(4, 2, 2), row.names = c("s1", "s2", "s3"))
  nx <- .bm_normalize(x, "TSS", 0)
  expect_true(is.na(nx$a[1]))
  expect_equal(nx$b, c(1, 0.5, 0.25))
  expect_equal(.bm_transform(nx, "LOG")$b, log2(c(1, 0.5, 0.25)))
  expect_equal(ncol(.bm_filter(nx, 0, 0.7, 1.01, 0)), 1)
})

test_that("hurdle model fits both parts and recovers a prevalence effect", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle()
  out <- file.path(tempdir(), "bm_test_out")
  fit <- bayome(s$counts, s$meta, output = out,
                  formula = ~ diagnosis + age + (1 | subject),
                  reference = "diagnosis,nonIBD", shrinkage = "none",
                  verbosity = "quiet")
  expect_s3_class(fit, "bayome")
  r <- fit$results
  expect_setequal(unique(r$model), c("prevalence", "abundance"))
  expect_true(all(c("coef", "stderr", "qval_individual", "qval_joint",
                    "lfsr", "N_not_zero") %in% names(r)))
  hit <- r[r$feature == "taxon1" & r$model == "prevalence" &
             r$name == "diagnosisCD", ]
  expect_equal(nrow(hit), 1)
  expect_gt(hit$coef, 0.5)
  expect_lt(hit$lfsr, 0.05)
  expect_true(all(r$lower_95 <= r$coef & r$coef <= r$upper_95, na.rm = TRUE))
  expect_true(file.exists(file.path(out, "all_results.tsv")))
  expect_true(file.exists(file.path(out, "significant_results.tsv")))
})

test_that("ordered() gives consecutive-level contrasts", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 2)
  fit <- bayome(s$counts, s$meta, formula = ~ diagnosis + ordered(stage),
                  shrinkage = "none",
                  verbosity = "quiet")
  st <- fit$results[fit$results$metadata == "stage", ]
  expect_setequal(unique(st$value), c("mid", "high"))
  expect_equal(nrow(st), 2 * 2 * 2)  # 2 features x 2 parts x 2 contrasts
  expect_true(all(is.finite(st$coef)))
})

test_that("evaluate_only, shared link and skipped parts work", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 2)
  s$counts[, 2] <- 10 + seq_len(nrow(s$counts))  # present everywhere
  fit <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                  link = "shared", shrinkage = "none", verbosity = "quiet")
  r <- fit$results
  skipped <- r[r$feature == "taxon2" & r$model == "prevalence", ]
  expect_match(skipped$error, "present in all samples")
  expect_true(all(is.na(r$error[r$feature == "taxon2" & r$model == "abundance"])))

  fp <- bayome(s$counts, s$meta, formula = ~ diagnosis,
                 evaluate_only = "prevalence", verbosity = "quiet")
  expect_equal(unique(fp$results$model), "prevalence")
})

test_that("unsupported maaslin3 arguments warn", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 2)
  expect_warning(bayome(s$counts, s$meta, formula = ~ diagnosis,
                          augment = TRUE, verbosity = "quiet"),
                 "augment")
})

test_that("median comparison tests against the posterior median", {
  set.seed(3)
  S <- 2000
  res <- data.frame(
    feature = paste0("f", 1:6), name = "diagnosisCD", model = "abundance",
    coef = c(-1, -1, -1, -1, 1, 1), stderr = 0.1, lower_95 = NA,
    upper_95 = NA, prob_positive = NA, lfsr = NA, error = NA_character_,
    stringsAsFactors = FALSE)
  draws <- t(sapply(res$coef, function(m) rnorm(S, m, 0.1)))
  out <- .bm_median_comparison(res, draws, "abundance")
  expect_equal(out$null_hypothesis[1], -1, tolerance = 0.1)
  expect_true(all(out$lfsr[5:6] < 0.01))   # far above the median
  expect_true(all(out$prob_positive[5:6] > 0.99))
  expect_equal(out$coef, res$coef)         # raw coefficients kept

  sub <- .bm_median_comparison(res, draws, "abundance", subtract_median = TRUE)
  expect_equal(sub$coef[5], 2, tolerance = 0.1)

  thr <- .bm_median_comparison(res, draws, "abundance", threshold = 5)
  expect_true(all(thr$lfsr[5:6] > 0.99))   # nothing is 5 units from the median

  few <- .bm_median_comparison(res[1:2, ], draws[1:2, ], "abundance")
  expect_true(all(is.na(few$lfsr)))        # fewer than 3 features: no comparison
})

test_that("bayome reports the median used as null hypothesis", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 4)
  set.seed(1)
  fit <- bayome(s$counts, s$meta, formula = ~ diagnosis, n_draws = 500,
                  verbosity = "quiet")
  r <- fit$results
  expect_true("null_hypothesis" %in% names(r))
  ab <- r[r$model == "abundance" & is.na(r$error), ]
  expect_true(all(ab$null_hypothesis != 0))
  expect_true(all(r$null_hypothesis[r$model == "prevalence"] == 0))

  off <- bayome(s$counts, s$meta, formula = ~ diagnosis,
                  median_comparison_abundance = FALSE, verbosity = "quiet")
  expect_true(all(off$results$null_hypothesis == 0))
})

test_that("shrinkage options run end to end", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 5)
  set.seed(2)
  fa <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                 n_draws = 300, verbosity = "quiet")
  expect_equal(fa$settings$shrinkage, "joint")
  expect_equal(fa$settings$selection, "bh_lfsr")
  set.seed(2)
  fh <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                 shrinkage = "hierarchical", n_draws = 300, verbosity = "quiet")
  rh <- fh$results
  expect_setequal(unique(rh$model), c("prevalence", "abundance"))
  expect_equal(nrow(rh), 2 * 5)
  expect_true(all(rh$null_hypothesis[rh$model == "prevalence"] == 0))
  expect_true(all(rh$null_hypothesis[rh$model == "abundance"] != 0))
  expect_true(all(is.finite(rh$lfsr)), info = paste(unique(rh$error), collapse = " | "))
  expect_equal(fh$settings$hierarchical_reference, "mean")

  set.seed(2)
  fm <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                 shrinkage = "hierarchical", hierarchical_reference = "median",
                 n_draws = 300, save_models = TRUE, verbosity = "quiet")
  rm_ab <- fm$results[fm$results$model == "abundance", ]
  expect_true(all(is.finite(rm_ab$lfsr)))
  expect_length(unique(rm_ab$null_hypothesis), 1)
  expect_false(is.null(fm$fits$abundance))

  set.seed(2)
  fs <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                 shrinkage = "hierarchical", subtract_median = TRUE,
                 n_draws = 300, save_models = TRUE, verbosity = "quiet")
  ab <- fs$results[fs$results$model == "abundance", ]
  u <- fs$fits$abundance$summary.random[[".bm_b_diagnosisCD"]]
  expect_equal(sort(ab$coef), sort(u$mean), tolerance = 1e-8)
  expect_equal(unique(ab$null_hypothesis),
               fs$fits$abundance$summary.fixed["mu_diagnosisCD", "mean"])

  fe <- bayome(s$counts, s$meta, formula = ~ diagnosis + ordered(stage),
                 shrinkage = "hierarchical", verbosity = "quiet")
  expect_true(all(grepl("ordered", fe$results$error)))
})

test_that("hierarchical model can estimate one residual variance per feature", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 5)
  set.seed(2)
  fit <- bayome(s$counts, s$meta, formula = ~ diagnosis,
                  shrinkage = "hierarchical", hierarchical_variance = "estimated",
                  int_strategy = "eb", evaluate_only = "abundance",
                  n_draws = 300, save_models = TRUE, verbosity = "quiet")
  r <- fit$results
  expect_true(all(is.finite(r$lfsr)), info = paste(unique(r$error), collapse = " | "))
  hyp <- rownames(fit$fits$abundance$summary.hyperpar)
  expect_equal(sum(grepl("Gaussian observations", hyp)), 5)
})

test_that("variational variance updates converge and need no lm()", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 6)
  set.seed(2)
  fit <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                  shrinkage = "hierarchical", hierarchical_variance = "vb",
                  evaluate_only = "abundance", n_draws = 300,
                  save_models = TRUE, verbosity = "quiet")
  r <- fit$results
  expect_true(all(is.finite(r$lfsr)), info = paste(unique(r$error), collapse = " | "))
  f <- fit$fits$abundance
  v <- attr(f, "bm_residual_variance")
  expect_length(v, 6)
  expect_true(all(v > 0))
  info <- attr(f, "bm_vb")
  expect_lte(info$iterations, 10)
})

test_that("moment estimator recovers a common variance prior", {
  set.seed(4)
  d <- rep(20, 400)
  sigma2 <- 30 * 0.5 / rchisq(400, df = 30)       # Inv-chi2(d0 = 30, s0^2 = 0.5)
  s2 <- sigma2 * rchisq(400, df = d) / d
  m <- .bm_moderate_variances(s2, d)
  expect_equal(m$s02, 0.5, tolerance = 0.1)
  expect_gt(m$d0, 15)
  expect_lt(m$d0, 60)
  expect_identical(.bm_moderate_variances(rep(1, 50), rep(10, 50))$d0, Inf)
})

test_that("joint model links abundance to prevalence and returns both parts", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 6)
  set.seed(3)
  fit <- bayome(s$counts, s$meta, formula = ~ diagnosis + (1 | subject),
                  n_draws = 300, save_models = TRUE, verbosity = "quiet")
  r <- fit$results
  expect_setequal(unique(r$model), c("prevalence", "abundance"))
  ok <- r[is.na(r$error), ]
  expect_true(all(is.finite(ok$coef) & is.finite(ok$lfsr)))
  expect_true(all(ok$lower_95 <= ok$upper_95))
  lam <- attr(fit$fits$joint, "bm_lambda")
  expect_named(lam, "diagnosisCD")
  expect_true(is.finite(lam))
  expect_true(all(ok$null_hypothesis[ok$model == "prevalence"] == 0))
})

test_that("joint model gives the full structure to the exposure only", {
  skip_if_not_installed("INLA")
  s <- sim_hurdle(J = 6)
  set.seed(3)
  expect_message(
    fit <- bayome(s$counts, s$meta, formula = ~ diagnosis + age + (1 | subject),
                  n_draws = 300, save_models = TRUE, verbosity = "INFO"),
    "exposure = \"diagnosis\"")
  expect_named(fit$lambda, "diagnosisCD")
  expect_identical(fit$settings$exposure, "diagnosis")
  r <- fit$results[is.na(fit$results$error), ]
  # Adjustment covariate: estimated in both parts, finite.
  age <- r[r$name == "age", ]
  expect_setequal(unique(age$model), c("prevalence", "abundance"))
  expect_true(all(is.finite(age$coef) & age$lower_95 <= age$upper_95))
  # Only the exposure adds hyperparameters: lambda and tau for diagnosisCD,
  # plus one random-effect precision per part.
  hp <- rownames(fit$fits$joint$summary.hyperpar)
  expect_length(hp, 4)
  expect_false(any(grepl("age", hp)))

  # Every covariate as exposure: one lambda each.
  fit2 <- bayome(s$counts, s$meta, formula = ~ diagnosis + age,
                 exposure = c("diagnosis", "age"), n_draws = 300,
                 verbosity = "quiet")
  expect_named(fit2$lambda, c("diagnosisCD", "age"))
  expect_error(bayome(s$counts, s$meta, formula = ~ diagnosis, exposure = "nope",
                      verbosity = "quiet"), "not found")
})

test_that("selection rules give the expected q-values", {
  res <- data.frame(
    feature = rep(c("a", "b", "c"), each = 2), name = "x",
    model = rep(c("prevalence", "abundance"), 3),
    coef = c(2, 0.1, 0, 0, 1, 1.5), stderr = 0.2, null_hypothesis = 0,
    lfsr = c(0.001, 0.3, 0.4, 0.45, 0.01, 0.002), stringsAsFactors = FALSE)
  bh <- .bm_add_qvalues(res, "bh_lfsr")
  pr <- res$model == "prevalence"
  expect_equal(bh$qval_individual[pr], p.adjust(pmin(2 * res$lfsr[pr], 1), "BH"))
  pj <- pbeta(pmin(2 * c(0.001, 0.4, 0.002), 1), 1, 2)
  expect_equal(unname(bh$qval_joint[c(1, 3, 5)]), p.adjust(pj, "BH"))
  lf <- .bm_add_qvalues(res, "lfsr")
  expect_equal(lf$qval_individual[pr], bm_qvalue(res$lfsr[pr]))
  rp <- .bm_add_qvalues(res, "rope", rope_delta = 0.5)
  # Large effects have small probability of being practically null.
  expect_lt(rp$qval_individual[1], 0.01)
  expect_gt(rp$qval_individual[3], 0.1)                 # null effect not selected
  expect_gt(rp$qval_individual[3], rp$qval_individual[5])
})
