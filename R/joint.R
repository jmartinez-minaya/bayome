# Joint two-part model (shrinkage = "joint"): prevalence and abundance of all
# features in one INLA fit, with the abundance effects linked to the
# prevalence effects.
#
# For feature j, sample i, exposure columns k (model-matrix columns of the
# exposure variable) and adjustment columns l (all other covariates):
#   prevalence  logit P(present_ij) = a_j^P + sum_k x_ik (mu_k^P + u_jk^P)
#                                     + sum_l z_il g_jl^P + random effects^P
#   abundance   log2 y_ij ~ N(a_j^A + sum_k x_ik (mu_k^A + lambda_k u_jk^P + e_jk)
#                             + sum_l z_il g_jl^A + random effects^A, sigma2_j)
#                                                                (present only)
# with
#   u_jk^P ~ N(0, 1 / prec_prev)   fixed, weakly informative: no pooling
#   e_jk   ~ N(0, 1 / tau_k)       pooled across features (PC prior on tau_k)
#   lambda_k                       one per exposure column, N(0, 1) prior
#                                  (INLA `copy` of the prevalence effect)
#   g_jl   ~ N(0, 1 / prec_adj)    adjustment effects, one per feature and
#                                  part, fixed prior: no pooling, no link and
#                                  no hyperparameters
#   sigma2_j                       "vb" updates with a common Inv-chi2 prior
# Hyperparameters: lambda_k and tau_k per exposure column, plus the
# variances of the formula's random effects.
# Random effects of the formula (e.g. (1 | subject)) get one copy per part,
# replicated per feature with shared variances. The two parts are separated
# with zero weights (INLA crashed with NA indices). Features without both
# presences and absences only enter the abundance part, and features with
# too few non-zero samples only the prevalence part; for the former the
# lambda link is switched off.
#
# Under median comparison, effects are compared with their median over
# features, computed on joint posterior samples.

.bm_fit_joint <- function(features, ok_prev, ok_abun, model, design_full, prior,
                          compare, subtract_median, control, n_draws = 1000,
                          exposure_cols = NULL) {
  if (any(model$terms$type == "ordered")) {
    stop("ordered() is not supported yet with shrinkage = \"joint\".", call. = FALSE)
  }
  keep <- ok_prev | ok_abun
  feats <- colnames(features)[keep]
  ok_prev <- ok_prev[keep]
  ok_abun <- ok_abun[keep]
  Y <- as.matrix(features[, feats, drop = FALSE])
  X <- design_full$X
  n <- nrow(Y)
  J <- length(feats)
  data <- model$data
  coef_cols <- setdiff(colnames(X), "(Intercept)")
  safe <- make.unique(gsub("[^A-Za-z0-9_]", "_", coef_cols))
  K <- length(coef_cols)
  if (!K) stop("The joint model needs at least one covariate.", call. = FALSE)
  if (is.null(exposure_cols)) exposure_cols <- coef_cols
  is_exp <- coef_cols %in% exposure_cols
  if (!any(is_exp)) stop("No model column belongs to the exposure.", call. = FALSE)

  # Long format. Prevalence: every sample of features with presences and
  # absences. Abundance: present samples of features with enough of them.
  pj <- which(ok_prev)
  prev <- data.frame(i = rep(seq_len(n), length(pj)), j = rep(pj, each = n))
  prev$y <- as.numeric(!is.na(Y[cbind(prev$i, prev$j)]))
  obs <- which(!is.na(Y), arr.ind = TRUE)
  obs <- obs[ok_abun[obs[, 2]], , drop = FALSE]
  abun <- data.frame(i = obs[, 1], j = obs[, 2], y = Y[obs])
  NP <- nrow(prev)
  NA_ <- nrow(abun)
  N <- NP + NA_
  rows <- rbind(prev, abun)
  isA <- c(rep(0, NP), rep(1, NA_))
  isP <- 1 - isA
  resp <- matrix(NA_real_, N, 2)
  resp[seq_len(NP), 1] <- prev$y
  resp[NP + seq_len(NA_), 2] <- abun$y

  stack <- list(Y = resp, intP = isP, intA = isA,
                itP = rows$j, witP = isP, itA = rows$j, witA = isA,
                .bm_rep = rows$j)
  hyper <- prior$hyper
  hyper$between <- prior$hyper$random
  hyper$prev <- list(prec = list(initial = log(prior$fixed$prevalence), fixed = TRUE))
  hyper$beta <- list(beta = list(prior = "normal", param = c(0, 1)))
  adj_prec <- if (is.null(prior$fixed$adjustment)) 1 / 2.5^2 else prior$fixed$adjustment
  hyper$adj <- list(prec = list(initial = log(adj_prec), fixed = TRUE))
  rhs <- c("-1", "intP", "intA",
           sprintf("f(itP, witP, model = \"iid\", n = %d, hyper = .bm_hyper$strata)", J),
           sprintf("f(itA, witA, model = \"iid\", n = %d, hyper = .bm_hyper$strata)", J))
  prec <- list(intP = 0, intA = 0)
  link_on <- ok_prev[rows$j]
  for (k in seq_len(K)) {
    x <- X[rows$i, coef_cols[k]]
    s <- safe[k]
    if (!is_exp[k]) {
      # Adjustment covariate: one effect per feature and part, fixed prior.
      stack[[paste0("gP_", s)]] <- rows$j
      stack[[paste0("wgP_", s)]] <- x * isP
      stack[[paste0("gA_", s)]] <- rows$j
      stack[[paste0("wgA_", s)]] <- x * isA
      rhs <- c(rhs,
               sprintf("f(gP_%s, wgP_%s, model = \"iid\", n = %d, hyper = .bm_hyper$adj)",
                       s, s, J),
               sprintf("f(gA_%s, wgA_%s, model = \"iid\", n = %d, hyper = .bm_hyper$adj)",
                       s, s, J))
      next
    }
    stack[[paste0("muP_", s)]] <- x * isP
    stack[[paste0("muA_", s)]] <- x * isA
    stack[[paste0("dP_", s)]] <- rows$j
    stack[[paste0("wdP_", s)]] <- x * isP
    stack[[paste0("dPc_", s)]] <- rows$j
    stack[[paste0("wdPc_", s)]] <- x * isA * link_on
    stack[[paste0("eA_", s)]] <- rows$j
    stack[[paste0("weA_", s)]] <- x * isA
    prec[[paste0("muP_", s)]] <- prior$fixed$prevalence
    prec[[paste0("muA_", s)]] <- prior$fixed$abundance
    rhs <- c(rhs, paste0("muP_", s), paste0("muA_", s),
             sprintf("f(dP_%s, wdP_%s, model = \"iid\", n = %d, hyper = .bm_hyper$prev)",
                     s, s, J),
             sprintf(paste("f(dPc_%s, wdPc_%s, copy = \"dP_%s\", fixed = FALSE,",
                           "hyper = .bm_hyper$beta)"), s, s, s),
             sprintf("f(eA_%s, weA_%s, model = \"iid\", n = %d, hyper = .bm_hyper$between)",
                     s, s, J))
  }

  # Formula random effects and strata: one term per part, replicated per
  # feature, separated by zero weights.
  ft <- model$terms[!is.na(model$terms$index), , drop = FALSE]
  for (t in seq_len(nrow(ft))) {
    idx <- ft$index[t]
    w <- if (is.na(ft$weights[t])) rep(1, N) else data[[ft$weights[t]]][rows$i]
    for (part in c("P", "A")) {
      nm <- paste0(idx, "_", part)
      stack[[nm]] <- data[[idx]][rows$i]
      stack[[paste0("w", nm)]] <- w * if (part == "P") isP else isA
      rhs <- c(rhs, sprintf("f(%s, w%s, %s, replicate = .bm_rep, hyper = .bm_hyper$%s)",
                            nm, nm, ft$f_args[t], ft$hyper[t]))
    }
  }

  formula <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")),
                               env = .bm_formula_env(list(hyper = hyper)))
  stack$.bm_ntrials <- rep(1, N)
  n_j <- tabulate(abun$j, nbins = J)
  ctrl_family <- list(list(), list(hyper = list(prec = list(initial = 0, fixed = TRUE))))
  # One INLA fit. `predictor`: compute the linear predictor (needed only by
  # the "vb" updates); `theta`: hyperparameter mode of the previous fit, used
  # as a warm start.
  run <- function(prec_j, config = FALSE, strategy = control$int_strategy,
                  predictor = TRUE, theta = NULL) {
    stack$.bm_scale <- c(rep(1, NP), prec_j[abun$j])
    INLA::inla(formula, data = stack, family = c("binomial", "gaussian"),
               Ntrials = .bm_ntrials, scale = .bm_scale,
               control.family = ctrl_family,
               control.fixed = list(mean = 0, prec = prec),
               control.predictor = list(compute = predictor),
               control.compute = list(config = config, return.marginals = FALSE,
                                      return.marginals.predictor = FALSE),
               control.inla = list(int.strategy = strategy),
               control.mode = if (is.null(theta)) list() else
                 list(theta = theta, restart = TRUE),
               num.threads = control$num.threads, verbose = FALSE)
  }

  # Residual variances: "vb". Every fit uses the hyperparameter mode only
  # ("eb"), unless the user asks the final fit to integrate over the
  # hyperparameters (int_strategy = "ccd" or "grid"). With many features the
  # hyperparameters are well estimated and integrating over them changes the
  # results little but costs a fit per integration point.
  final_strategy <- if (control$int_strategy == "auto") "eb" else control$int_strategy
  prec_j <- 1 / .bm_initial_variances(abun$y, X[abun$i, , drop = FALSE], abun$j, J)
  vb_iter <- 0
  theta <- NULL
  seconds <- numeric()
  changes <- numeric()
  last_step <- NULL
  say <- if (is.function(control$say)) control$say else function(...) NULL
  t_start <- Sys.time()
  say(sprintf(paste("  %d INLA hyperparameters (%d exposure coefficient(s)",
                    "x lambda and tau, plus random effects); residual variances",
                    "by \"vb\", stopping when they change less than %g%%."),
              2 * sum(is_exp) + 2 * sum(ft$hyper != "strata"), sum(is_exp),
              100 * control$vb_tol))
  repeat {
    vb_iter <- vb_iter + 1
    t0 <- Sys.time()
    fit <- run(prec_j, strategy = "eb", theta = theta)
    seconds[vb_iter] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    theta <- fit$mode$theta
    lp <- fit$summary.linear.predictor[NP + seq_len(NA_), ]
    ss <- as.vector(tapply((abun$y - lp$mean)^2 + lp$sd^2,
                           factor(abun$j, levels = seq_len(J)), sum))
    has <- n_j > 0
    mom <- .bm_moderate_variances(ss[has] / n_j[has], n_j[has])
    new <- prec_j
    new[has] <- if (is.finite(mom$d0)) (mom$d0 + n_j[has]) / (mom$d0 * mom$s02 + ss[has]) else
      1 / mom$s02
    change <- max(abs(new[has] / prec_j[has] - 1))
    changes[vb_iter] <- change
    # The updates converge linearly, each step a near-constant fraction r of
    # the previous one; extrapolating the remaining steps (Aitken), on the log
    # scale, saves iterations. Skipped once converged.
    step <- log(new[has]) - log(prec_j[has])
    if (!is.null(last_step) && change >= control$vb_tol) {
      r <- sqrt(sum(step^2) / sum(last_step^2))
      if (is.finite(r) && r > 0 && r < 0.9) {
        r <- min(r, 0.8)
        new[has] <- exp(log(new[has]) + r / (1 - r) * step)
      }
    }
    last_step <- step
    prec_j <- new
    say(sprintf("  vb iteration %d: largest change in the residual variances %.1f%% (%.0f s)",
                vb_iter, 100 * change, seconds[vb_iter]))
    if (change < control$vb_tol || vb_iter >= control$vb_max_iter) break
  }
  if (change >= control$vb_tol) {
    say(sprintf("  vb stopped after %d iterations without reaching the tolerance.",
                vb_iter))
  }
  say(sprintf("  Final fit (int_strategy = \"%s\")...", final_strategy))
  t0 <- Sys.time()
  fit <- run(prec_j, config = TRUE, strategy = final_strategy, predictor = FALSE,
             theta = theta)
  seconds <- c(seconds, final = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  say(sprintf("  Final fit done (%.0f s). Drawing %d joint posterior samples...",
              seconds[["final"]], n_draws))
  attr(fit, "bm_residual_variance") <- stats::setNames(1 / prec_j, feats)
  attr(fit, "bm_vb") <- list(iterations = vb_iter, last_change = change,
                             changes = changes, d0 = mom$d0, s02 = mom$s02,
                             seconds = seconds)

  sel <- list()
  for (s in safe[is_exp]) {
    sel[[paste0("muP_", s)]] <- 1
    sel[[paste0("muA_", s)]] <- 1
    sel[[paste0("dP_", s)]] <- seq_len(J)
    sel[[paste0("eA_", s)]] <- seq_len(J)
  }
  for (s in safe[!is_exp]) {
    sel[[paste0("gP_", s)]] <- seq_len(J)
    sel[[paste0("gA_", s)]] <- seq_len(J)
  }
  smp <- INLA::inla.posterior.sample(n_draws, fit, selection = sel,
                                     num.threads = control$num.threads)
  nm <- rownames(smp[[1]]$latent)
  lat <- vapply(smp, function(z) z$latent[, 1], numeric(length(nm)))
  pick <- function(name, idx) lat[match(paste0(name, ":", idx), nm), , drop = FALSE]

  summ <- function(m) {
    data.frame(
      coef = rowMeans(m), stderr = apply(m, 1, stats::sd),
      lower_95 = apply(m, 1, stats::quantile, 0.025, names = FALSE),
      upper_95 = apply(m, 1, stats::quantile, 0.975, names = FALSE))
  }
  info <- design_full$coef_info
  info <- info[info$column != "(Intercept)", , drop = FALSE]
  N_nz <- colSums(!is.na(Y))
  out <- list()
  lambda <- stats::setNames(rep(NA_real_, K), coef_cols)
  for (k in seq_len(K)) {
    s <- safe[k]
    if (is_exp[k]) {
      lam <- vapply(smp, function(z) {
        z$hyperpar[grep(paste0("^Beta for dPc_", s, "$"), names(z$hyperpar))]
      }, numeric(1))
      lambda[k] <- mean(lam)
      uP <- pick(paste0("dP_", s), seq_len(J))
      e <- pick(paste0("eA_", s), seq_len(J))
      betaP <- sweep(uP, 2, pick(paste0("muP_", s), 1)[1, ], "+")
      betaA <- sweep(sweep(uP * ok_prev, 2, lam, "*") + e, 2,
                     pick(paste0("muA_", s), 1)[1, ], "+")
    } else {
      betaP <- pick(paste0("gP_", s), seq_len(J))
      betaA <- pick(paste0("gA_", s), seq_len(J))
    }
    for (part in c("prevalence", "abundance")) {
      take <- if (part == "prevalence") ok_prev else ok_abun
      if (!any(take)) next
      beta <- (if (part == "prevalence") betaP else betaA)[take, , drop = FALSE]
      null_h <- 0
      if (compare[[part]] && sum(take) >= 3) {
        med <- apply(beta, 2, stats::median)
        delta <- sweep(beta, 2, med, "-")
        null_h <- mean(med)
      } else {
        delta <- beta
      }
      st <- summ(if (subtract_median) delta else beta)
      st$prob_positive <- rowMeans(delta > 0)
      st$lfsr <- pmin(rowMeans(delta <= 0), rowMeans(delta >= 0))
      out[[length(out) + 1]] <- data.frame(
        feature = feats[take], metadata = info$metadata[k], value = info$value[k],
        name = info$name[k], st, error = NA_character_, model = part, N = n,
        N_not_zero = N_nz[take], null_hypothesis = null_h,
        stringsAsFactors = FALSE, row.names = NULL)
    }
  }
  lambda <- lambda[is_exp]
  attr(fit, "bm_lambda") <- lambda
  say(sprintf("  Joint model done in %.1f minutes.",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  list(results = do.call(rbind, out), fit = fit, lambda = lambda)
}

# Starting values for the residual variances of the "vb" updates: for each
# feature, the residual variance of a least-squares fit of its abundance on
# the fixed effects, moderated towards a common value as in limma (features
# with too few samples get the common value). Random effects are ignored,
# so the values are somewhat too large, but much closer than one global
# variance, which also absorbs the differences between feature means.
.bm_initial_variances <- function(y, X, j, J) {
  s2 <- rep(NA_real_, J)
  d <- rep(0, J)
  for (k in seq_len(J)) {
    i <- which(j == k)
    if (length(i) < 2) next
    q <- qr(X[i, , drop = FALSE])
    d[k] <- length(i) - q$rank
    if (d[k] > 0) s2[k] <- sum(qr.resid(q, y[i])^2) / d[k]
  }
  ok <- is.finite(s2) & s2 > 0 & d > 0
  if (!any(ok)) return(rep(stats::var(y), J))
  mom <- .bm_moderate_variances(s2[ok], d[ok])
  out <- rep(mom$s02, J)
  if (is.finite(mom$d0)) {
    out[ok] <- (mom$d0 * mom$s02 + d[ok] * s2[ok]) / (mom$d0 + d[ok])
  }
  out
}
