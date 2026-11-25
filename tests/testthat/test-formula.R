make_meta <- function(n = 40) {
  set.seed(1)
  data.frame(
    diagnosis = factor(sample(c("nonIBD", "UC", "CD"), n, replace = TRUE)),
    age = round(runif(n, 18, 70)),
    stage = factor(sample(c("low", "mid", "high"), n, replace = TRUE),
                   levels = c("low", "mid", "high")),
    site = sample(c("A", "B", "C"), n, replace = TRUE),
    pair = rep(seq_len(n / 2), each = 2),
    subject = rep(paste0("s", seq_len(n / 4)), each = 4),
    row.names = paste0("sample", seq_len(n))
  )
}

test_that("parser separates fixed, specials and random terms", {
  p <- bm_parse_formula(~ diagnosis + age + group(site) + ordered(stage) +
                          strata(pair) + (1 | subject))
  expect_s3_class(p, "bm_formula")
  expect_equal(p$fixed, c("diagnosis", "age"))
  expect_equal(p$group, "site")
  expect_equal(p$ordered, "stage")
  expect_equal(p$strata, "pair")
  expect_equal(p$random$group, "subject")
  expect_true(is.na(p$random$slope))
  expect_true(p$intercept)
})

test_that("parser accepts strings and ignores a left-hand side", {
  p1 <- bm_parse_formula("~ diagnosis + (1 | subject)")
  p2 <- bm_parse_formula(y ~ diagnosis + (1 | subject))
  expect_equal(p1$fixed, p2$fixed)
  expect_equal(p1$random, p2$random)
})

test_that("double bar expands into intercept and uncorrelated slope", {
  p <- bm_parse_formula(~ diagnosis + (age || subject))
  expect_equal(nrow(p$random), 2)
  expect_equal(sort(p$random$slope, na.last = TRUE), c("age", NA))
})

test_that("correlated random slopes give an informative error", {
  expect_error(bm_parse_formula(~ diagnosis + (age | subject)),
               "not supported yet")
})

test_that("only one strata term is allowed", {
  expect_error(bm_parse_formula(~ diagnosis + strata(a) + strata(b)),
               "Only one strata")
})

test_that("compose_formula mirrors maaslin3 effect arguments", {
  f <- bm_compose_formula(fixed_effects = c("diagnosis", "age"),
                          random_effects = "subject",
                          ordered_effects = "stage")
  p <- bm_parse_formula(f)
  expect_equal(p$fixed, c("diagnosis", "age"))
  expect_equal(p$ordered, "stage")
  expect_equal(p$random$group, "subject")
  expect_error(bm_compose_formula(), "No effects")
})

test_that("INLA formula has the expected f() terms and index columns", {
  meta <- make_meta()
  m <- bm_inla_formula(~ diagnosis + age + group(site) + ordered(stage) +
                         strata(pair) + (age || subject), meta)
  txt <- deparse1(m$formula)
  expect_match(txt, "^y ~ 1 \\+ diagnosis \\+ age \\+ site")
  expect_match(txt, "f\\(.bm_ord_stage, model = \"rw1\"")
  expect_match(txt, "f\\(.bm_strata_pair, model = \"iid\"")
  expect_match(txt, "f\\(.bm_re_subject, model = \"iid\"")
  expect_match(txt, "f\\(.bm_re_subject__age, age, model = \"iid\"")
  expect_equal(m$data$.bm_ord_stage, as.integer(meta$stage))
  expect_equal(max(m$data$.bm_re_subject), 10L)
  expect_true(is.factor(m$data$site))
  expect_true(exists(".bm_hyper", envir = environment(m$formula)))
})

test_that("missing metadata variables are reported", {
  meta <- make_meta()
  expect_error(bm_inla_formula(~ diagnosis + bmi, meta), "bmi")
})


test_that("generated model runs in INLA", {
  skip_if_not_installed("INLA")
  meta <- make_meta(200)
  meta$subject <- rep(paste0("s", 1:20), each = 10)
  u <- rnorm(20, 0, 0.5)
  meta$y <- 1 + 0.03 * meta$age + u[as.integer(factor(meta$subject))] +
    rnorm(200, 0, 0.3)
  m <- bm_inla_formula(~ age + ordered(stage) + (1 | subject), meta)
  fit <- INLA::inla(m$formula, data = m$data, family = "gaussian",
                    control.fixed = list(prec = bm_prior()$fixed$abundance))
  expect_equal(fit$summary.fixed["age", "mean"], 0.03, tolerance = 0.3)
  expect_true(all(c(".bm_re_subject", ".bm_ord_stage") %in% names(fit$summary.random)))
})
