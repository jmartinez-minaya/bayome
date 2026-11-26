make_plot_fit <- function() {
  set.seed(1)
  n <- 40
  feats <- c("taxon_a", "taxon_b", "taxon_c")
  data <- as.data.frame(matrix(stats::rnorm(n * 3, -5), n, 3,
                               dimnames = list(paste0("s", 1:n), feats)))
  data[sample(n, 10), "taxon_a"] <- NA
  metadata <- data.frame(diagnosis = factor(rep(c("nonIBD", "CD"), n / 2)),
                         age = stats::rnorm(n), reads = stats::rnorm(n),
                         row.names = rownames(data))
  res <- expand.grid(feature = feats, name = c("diagnosisCD", "age", "reads"),
                     model = c("prevalence", "abundance"),
                     stringsAsFactors = FALSE)
  res$metadata <- ifelse(res$name == "diagnosisCD", "diagnosis", res$name)
  res$coef <- stats::rnorm(nrow(res))
  res$lower_95 <- res$coef - 1
  res$upper_95 <- res$coef + 1
  res$qval_individual <- stats::runif(nrow(res))
  res$qval_individual[1] <- 0.001
  res$qval_joint <- res$qval_individual
  res$null_hypothesis <- ifelse(res$model == "abundance", 0.1, 0)
  res$error <- NA_character_
  res$error[nrow(res)] <- "Too few non-zero samples for the abundance model"
  structure(list(results = res, data = data, metadata = metadata,
                 lambda = c(diagnosisCD = 0.4, age = 0, reads = 0),
                 settings = list(max_significance = 0.1)),
            class = "bayome")
}

draws <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  print(p)
  invisible(TRUE)
}

test_that("plot functions return plots that draw", {
  fit <- make_plot_fit()
  plots <- list(
    bm_plot_summary(fit),
    bm_plot_summary(fit$results, coef_names = "diagnosisCD"),
    bm_plot_summary(fit, coef_names = c("diagnosisCD", "age"),
                    heatmap_names = character()),
    bm_plot_summary(fit, n_top = 2),
    bm_plot_effects(fit, "diagnosisCD"),
    bm_plot_link(fit, "diagnosisCD"),
    bm_plot_association(fit, "taxon_a", "diagnosis"),
    bm_plot_association(fit, "taxon_a", "age")
  )
  for (p in plots) {
    expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
    expect_no_error(draws(p))
  }
})

test_that("bm_plot_summary puts the exposure on the left and the rest in the heatmap", {
  fit <- make_plot_fit()
  fit$settings$exposure <- "diagnosis"
  p <- bm_plot_summary(fit)
  expect_s3_class(p, "patchwork")
  # Left: coefficient plot of the exposure; right: heatmap of the others.
  expect_setequal(unique(as.character(p[[1]]$data$name)), "diagnosisCD")
  expect_setequal(unique(as.character(p[[2]]$data$name)), c("age", "reads"))
  q <- bm_plot_summary(fit$results, coef_names = "age")
  expect_setequal(unique(as.character(q[[2]]$data$name)), c("diagnosisCD", "reads"))
})

test_that("plot helpers rank features and drop errors", {
  fit <- make_plot_fit()
  res <- .bm_plot_results(fit)
  expect_false(any(!is.na(res$error)))
  expect_identical(.bm_top_features(fit$results, 1), "taxon_a")
  expect_equal(.bm_plot_threshold(fit, NULL), 0.1)
  expect_equal(.bm_plot_threshold(fit$results, 0.05), 0.05)
})

test_that("bm_plot_link shows lambda and plot inputs are checked", {
  fit <- make_plot_fit()
  p <- bm_plot_link(fit, "diagnosisCD")
  expect_match(ggplot2::get_labs(p)$subtitle, "0.40")
  expect_error(bm_plot_effects(fit, "nope"), "not found")
  expect_error(bm_plot_summary(fit, coef_names = "nope"), "not found")
  expect_error(bm_plot_summary(list(a = 1)), "bayome object")
  expect_error(bm_plot_association(fit, "nope", "diagnosis"), "not found")
  expect_error(bm_plot_association(fit$results, "taxon_a", "diagnosis"),
               "are needed")
})
