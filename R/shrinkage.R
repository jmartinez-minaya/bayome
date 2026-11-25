# Information sharing across features.
#
# shrinkage = "hierarchical" (one joint model per part, fully Bayesian)
#   All features are stacked and fitted in a single INLA model per part.
#   Without median comparison: beta_jk ~ N(0, tau_k^2), centred at zero.
#   With median comparison: beta_jk = mu_k + u_jk, u_jk ~ N(0, tau_k^2).
#   The compositional reference is either
#     "mean"   (default): the latent common centre mu_k. The corrected effect
#              is u_jk = beta_jk - mu_k, a model parameter with its own
#              posterior marginal;
#     "median": the median over features of beta_jk, computed on joint
#              posterior samples. More robust when the true effects are
#              asymmetric (mostly in one direction), which biases mu_k.
#   Residual variances of the abundance part are one per feature, either
#     "vb" (default): INLA alternating with conjugate variational updates of
#                 each variance (see the "vb" branch below),
#     "plugin":    estimated with lm() per feature and passed as
#                 fixed weights (INLA `scale`); INLA estimates one common
#                 multiplier, or
#     "estimated": one Gaussian likelihood per feature, so INLA estimates
#                 every precision (one hyperparameter per feature).
#   tau_k has a PC prior. Random effects of the formula are replicated per
#   feature with shared variances. The prior is Gaussian (INLA's latent
#   Gaussian models), so it has no point mass at the reference: lfsr
#   controls sign errors, not point nulls.

# --- hierarchical INLA model --------------------------------------------------

# Fits all features of one part in a single INLA model. `features` holds the
# transformed abundances (NA = absent); `ok` is a logical vector of the
# features that can be fitted for this part.
.bm_fit_hierarchical_part <- function(part, features, ok, model, design_full,
                                      prior, compare, subtract_median, control,
                                      n_draws = 1000,
                                      reference = c("mean", "median"),
                                      residual_variance = c("vb", "plugin", "estimated")) {
  reference <- match.arg(reference)
  residual_variance <- match.arg(residual_variance)
  if (is.null(control$int_strategy)) control$int_strategy <- "auto"
  if (is.null(control$vb_tol)) control$vb_tol <- 0.01
  if (is.null(control$vb_max_iter)) control$vb_max_iter <- 10
  if (any(model$terms$type == "ordered")) {
    stop("ordered() is not supported yet with shrinkage = \"hierarchical\".",
         call. = FALSE)
  }
  X <- design_full$X
  n <- nrow(X)
  feats <- colnames(features)[ok]
  J <- length(feats)
  data <- model$data
  terms <- model$terms

  # Long format: one row per sample and feature (abundance: present only).
  Y <- as.matrix(features[, feats, drop = FALSE])
  rows <- list()
  for (j in seq_len(J)) {
    y <- Y[, j]
    keep <- if (part == "abundance") which(!is.na(y)) else seq_len(n)
    rows[[j]] <- data.frame(.i = keep, .j = j,
                            .y = if (part == "abundance") y[keep] else
                              as.numeric(!is.na(y[keep])))
  }
  long <- do.call(rbind, rows)
  N <- nrow(long)

  # Feature-specific intercepts: iid with fixed, very small precision.
  stack <- list(Y = long$.y, .bm_feat = long$.j, .bm_feat_rep = long$.j)
  rhs <- c("1", "f(.bm_feat, model = \"iid\", hyper = .bm_hyper$strata)")
  coef_cols <- setdiff(colnames(X), "(Intercept)")
  safe <- make.unique(gsub("[^A-Za-z0-9_]", "_", coef_cols))
  b_names <- paste0(".bm_b_", safe)
  mu_names <- paste0("mu_", safe)
  for (k in seq_along(coef_cols)) {
    x <- X[long$.i, coef_cols[k]]
    stack[[b_names[k]]] <- long$.j
    stack[[paste0(".bm_w_", safe[k])]] <- x
    rhs <- c(rhs, sprintf(
      "f(%s, %s, model = \"iid\", n = %d, hyper = .bm_hyper$between)",
      b_names[k], paste0(".bm_w_", safe[k]), J))
    if (compare) {
      stack[[mu_names[k]]] <- x
      rhs <- c(rhs, mu_names[k])
    }
  }

  # Formula random effects / strata, replicated per feature.
  ft <- terms[!is.na(terms$index), , drop = FALSE]
  for (i in seq_len(nrow(ft))) {
    idx <- ft$index[i]
    stack[[idx]] <- data[[idx]][long$.i]
    w <- ""
    if (!is.na(ft$weights[i])) {
      stack[[ft$weights[i]]] <- data[[ft$weights[i]]][long$.i]
      w <- paste0(", ", ft$weights[i])
    }
    rhs <- c(rhs, sprintf("f(%s%s, %s, replicate = .bm_feat_rep, hyper = .bm_hyper$%s)",
                          idx, w, ft$f_args[i], ft$hyper[i]))
  }

  hyper <- prior$hyper
  hyper$between <- prior$hyper$random
  formula <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")),
                               env = .bm_formula_env(list(hyper = hyper)))
  ctrl_fixed <- list(mean = 0, prec = prior$fixed[[part]], prec.intercept = 0)
  ctrl_compute <- list(return.marginals = TRUE, config = compare)

  if (part == "abundance" && residual_variance == "estimated") {
    # One Gaussian likelihood per feature, each with its own precision
    # estimated by INLA: the response is an N x J matrix with feature j's
    # values in column j and NA elsewhere.
    Ymat <- matrix(NA_real_, N, J)
    Ymat[cbind(seq_len(N), long$.j)] <- long$.y
    stack$Y <- Ymat
    fit <- INLA::inla(formula, data = stack, family = rep("gaussian", J),
                      control.fixed = ctrl_fixed,
                      control.compute = ctrl_compute,
                      control.inla = list(int.strategy = control$int_strategy),
                      num.threads = control$num.threads, verbose = FALSE)
  } else if (part == "abundance" && residual_variance == "vb") {
    # Mean-field variational updates of the residual variances, alternating
    # with INLA. Given q(sigma2_j), INLA fits everything else with the
    # Gaussian precision fixed at 1 and `scale` = E[1 / sigma2_j]. Given the
    # INLA fit, each sigma2_j has a conjugate scaled-Inv-chi2 update with a
    # common prior Inv-chi2(d0, s0^2) whose parameters are estimated from all
    # features (empirical Bayes, as in limma):
    #   E[1 / sigma2_j] = (d0 + n_j) / (d0 * s0^2 + SS_j),
    #   SS_j = sum_i [(y_ij - E eta_ij)^2 + Var eta_ij].
    n_j <- tabulate(long$.j, nbins = J)
    prec <- rep(1 / stats::var(long$.y), J)      # start: common variance
    ctrl_family <- list(hyper = list(prec = list(initial = 0, fixed = TRUE)))
    vb_iter <- 0
    repeat {
      vb_iter <- vb_iter + 1
      stack$.bm_scale <- prec[long$.j]
      fit <- INLA::inla(formula, data = stack, family = "gaussian",
                        scale = .bm_scale, control.family = ctrl_family,
                        control.fixed = ctrl_fixed,
                        control.compute = ctrl_compute,
                        control.predictor = list(compute = TRUE),
                        num.threads = control$num.threads, verbose = FALSE)
      lp <- fit$summary.linear.predictor[seq_len(N), ]
      ss <- as.vector(tapply((long$.y - lp$mean)^2 + lp$sd^2,
                             factor(long$.j, levels = seq_len(J)), sum))
      mom <- .bm_moderate_variances(ss / n_j, n_j)
      prec_new <- if (is.finite(mom$d0)) {
        (mom$d0 + n_j) / (mom$d0 * mom$s02 + ss)
      } else {
        rep(1 / mom$s02, J)
      }
      change <- max(abs(prec_new - prec) / prec)
      prec <- prec_new
      if (change < control$vb_tol || vb_iter >= control$vb_max_iter) break
    }
    attr(fit, "bm_residual_variance") <- stats::setNames(1 / prec, feats)
    attr(fit, "bm_vb") <- list(iterations = vb_iter, last_change = change,
                               d0 = mom$d0, s02 = mom$s02)
  } else if (part == "abundance") {
    # Per-feature residual variances, plugged in through INLA's `scale`.
    s2 <- vapply(seq_len(J), function(j) {
      r <- long$.j == j
      fit <- stats::lm.fit(X[long$.i[r], , drop = FALSE], long$.y[r])
      max(sum(fit$residuals^2) / max(sum(r) - fit$rank, 1), 1e-6)
    }, numeric(1))
    stack$.bm_scale <- (stats::median(s2) / s2)[long$.j]
    fit <- INLA::inla(formula, data = stack, family = "gaussian",
                      scale = .bm_scale, control.fixed = ctrl_fixed,
                      control.compute = ctrl_compute,
                      num.threads = control$num.threads, verbose = FALSE)
    tau <- fit$summary.hyperpar["Precision for the Gaussian observations", "mean"]
    attr(fit, "bm_residual_variance") <-
      stats::setNames(s2 / (stats::median(s2) * tau), feats)
  } else {
    stack$.bm_ntrials <- rep(1, N)
    fit <- INLA::inla(formula, data = stack, family = "binomial",
                      Ntrials = .bm_ntrials, control.fixed = ctrl_fixed,
                      control.compute = ctrl_compute,
                      num.threads = control$num.threads, verbose = FALSE)
  }

  # Joint posterior samples of beta_jk = mu_k + u_jk: needed for the raw
  # coefficients and, with reference = "median", for the median itself.
  if (compare) {
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop(paste("shrinkage = \"hierarchical\" with median comparison needs the",
                 "sn package: install.packages(\"sn\")"), call. = FALSE)
    }
    sel <- c(stats::setNames(as.list(rep(1, length(mu_names))), mu_names),
             stats::setNames(rep(list(seq_len(J)), length(b_names)), b_names))
    samp <- INLA::inla.posterior.sample(n_draws, fit, selection = sel,
                                        num.threads = control$num.threads)
    lat <- vapply(samp, function(s) s$latent[, 1], numeric(length(samp[[1]]$latent)))
    lat_names <- rownames(samp[[1]]$latent)
  }

  info <- design_full$coef_info
  info <- info[info$column != "(Intercept)", , drop = FALSE]
  out <- list()
  for (k in seq_along(coef_cols)) {
    b_marg <- fit$marginals.random[[b_names[k]]]
    b_sum <- fit$summary.random[[b_names[k]]]
    st <- .bm_marginal_stats(b_marg, b_sum)
    null_h <- 0
    if (compare) {
      mu <- lat[match(paste0(mu_names[k], ":1"), lat_names), ]
      b <- lat[match(paste0(b_names[k], ":", seq_len(J)), lat_names), , drop = FALSE]
      beta <- sweep(b, 2, mu, "+")                 # J x draws
      summ <- function(m) {
        data.frame(
          coef = rowMeans(m), stderr = apply(m, 1, stats::sd),
          lower_95 = apply(m, 1, stats::quantile, 0.025, names = FALSE),
          upper_95 = apply(m, 1, stats::quantile, 0.975, names = FALSE))
      }
      if (reference == "mean") {
        # Corrected effect u_jk = beta_jk - mu_k: exact INLA marginals.
        u_st <- st
        st <- if (subtract_median) u_st else summ(beta)
        st$prob_positive <- u_st$prob_positive
        st$lfsr <- u_st$lfsr
        null_h <- fit$summary.fixed[mu_names[k], "mean"]
      } else {
        med <- apply(beta, 2, stats::median)
        delta <- sweep(beta, 2, med, "-")
        null_h <- mean(med)
        st <- summ(if (subtract_median) delta else beta)
        st$prob_positive <- rowMeans(delta > 0)
        st$lfsr <- pmin(rowMeans(delta <= 0), rowMeans(delta >= 0))
      }
    }
    out[[k]] <- data.frame(
      feature = feats, metadata = info$metadata[k], value = info$value[k],
      name = info$name[k], st, error = NA_character_, model = part,
      N = n, N_not_zero = colSums(!is.na(Y)), null_hypothesis = null_h,
      stringsAsFactors = FALSE, row.names = NULL)
  }
  list(results = do.call(rbind, out), fit = fit)
}

# Method-of-moments estimate of a common scaled-Inv-chi2(d0, s0^2) prior for
# variances s2 with degrees of freedom d (limma's fitFDist, without the
# robust and trend options). d0 = Inf means all variances are equal.
.bm_moderate_variances <- function(s2, d) {
  ok <- is.finite(s2) & s2 > 0 & d > 0
  s2 <- s2[ok]
  d <- d[ok]
  e <- log(s2) - digamma(d / 2) + log(d / 2)
  extra <- stats::var(e) - mean(trigamma(d / 2))   # = trigamma(d0 / 2)
  if (!is.finite(extra) || extra <= 0) {
    return(list(d0 = Inf, s02 = exp(mean(e))))
  }
  # Inverse of trigamma by Newton's method (as in limma::trigammaInverse).
  y <- 0.5 + 1 / extra
  for (i in seq_len(50)) {
    tri <- trigamma(y)
    step <- tri * (1 - tri / extra) / psigamma(y, deriv = 2)
    y <- y + step
    if (-step / y < 1e-8) break
  }
  d0 <- 2 * y
  list(d0 = d0, s02 = exp(mean(e) + digamma(d0 / 2) - log(d0 / 2)))
}
