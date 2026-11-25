# Joint two-part model (shrinkage = "joint"): prevalence and abundance of all
# features in one INLA fit, with the abundance effects linked to the
# prevalence effects.
#
# For feature j, sample i and covariate column k (model-matrix columns):
#   prevalence  logit P(present_ij) = a_j^P + sum_k x_ik (mu_k^P + u_jk^P)
#                                     + random effects^P
#   abundance   log2 y_ij ~ N(a_j^A + sum_k x_ik (mu_k^A + lambda_k u_jk^P + e_jk)
#                             + random effects^A, sigma2_j)      (present only)
# with
#   u_jk^P ~ N(0, 1 / prec_prev)   fixed, weakly informative: no pooling
#   e_jk   ~ N(0, 1 / tau_k)       pooled across features (PC prior on tau_k)
#   lambda_k                       one per covariate column, N(0, 1) prior
#                                  (INLA `copy` of the prevalence effect)
#   sigma2_j                       "vb" updates with a common Inv-chi2 prior
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
                          max_iter = 30, tol = 0.01) {
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

  # Long format. Prevalence: every sample of features with presences and
  # absences. Abundance: present samples of features with enough of them.
  pj <- which(ok_prev)
  prev <- data.frame(i = rep(seq_len(n), length(pj)), j = rep(pj, each = n))
  prev$y <- as.numeric(!is.na(Y[cbind(prev$i, prev$j)]))
  obs <- which(!is.na(Y), arr.ind = TRUE)
  obs <- obs[ok_abun[obs[, 2]], , drop = FALSE]
  abun <- data.frame(i = obs[, 1], j = obs[, 2], y = Y[obs])
  NP <- nrow(prev); NA_ <- nrow(abun); N <- NP + NA_
  rows <- rbind(prev, abun)
  isA <- c(rep(0, NP), rep(1, NA_)); isP <- 1 - isA
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
  rhs <- c("-1", "intP", "intA",
           sprintf("f(itP, witP, model = \"iid\", n = %d, hyper = .bm_hyper$strata)", J),
           sprintf("f(itA, witA, model = \"iid\", n = %d, hyper = .bm_hyper$strata)", J))
  prec <- list(intP = 0, intA = 0)
  link_on <- ok_prev[rows$j]
  for (k in seq_len(K)) {
    x <- X[rows$i, coef_cols[k]]
    s <- safe[k]
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
             sprintf("f(dP_%s, wdP_%s, model = \"iid\", n = %d, hyper = .bm_hyper$prev)", s, s, J),
             sprintf("f(dPc_%s, wdPc_%s, copy = \"dP_%s\", fixed = FALSE, hyper = .bm_hyper$beta)", s, s, s),
             sprintf("f(eA_%s, weA_%s, model = \"iid\", n = %d, hyper = .bm_hyper$between)", s, s, J))
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
  run <- function(prec_j, config = FALSE, strategy = control$int_strategy) {
    stack$.bm_scale <- c(rep(1, NP), prec_j[abun$j])
    INLA::inla(formula, data = stack, family = c("binomial", "gaussian"),
               Ntrials = .bm_ntrials, scale = .bm_scale,
               control.family = ctrl_family,
               control.fixed = list(mean = 0, prec = prec),
               control.predictor = list(compute = TRUE),
               control.compute = list(config = config, return.marginals = FALSE),
               control.inla = list(int.strategy = strategy),
               num.threads = control$num.threads, verbose = FALSE)
  }

  # Residual variances: "vb" (intermediate iterations use the hyperparameter
  # mode only; the final fit integrates over the hyperparameters).
  prec_j <- rep(1 / stats::var(abun$y), J)
  vb_iter <- 0
  repeat {
    vb_iter <- vb_iter + 1
    fit <- run(prec_j, strategy = "eb")
    lp <- fit$summary.linear.predictor[NP + seq_len(NA_), ]
    ss <- as.vector(tapply((abun$y - lp$mean)^2 + lp$sd^2,
                           factor(abun$j, levels = seq_len(J)), sum))
    has <- n_j > 0
    mom <- .bm_moderate_variances(ss[has] / n_j[has], n_j[has])
    new <- prec_j
    new[has] <- if (is.finite(mom$d0)) (mom$d0 + n_j[has]) / (mom$d0 * mom$s02 + ss[has]) else
      1 / mom$s02
    change <- max(abs(new[has] / prec_j[has] - 1))
    prec_j <- new
    if (change < control$vb_tol || vb_iter >= control$vb_max_iter) break
  }
  fit <- run(prec_j, config = TRUE)
  attr(fit, "bm_residual_variance") <- stats::setNames(1 / prec_j, feats)
  attr(fit, "bm_vb") <- list(iterations = vb_iter, last_change = change,
                             d0 = mom$d0, s02 = mom$s02)

  sel <- list()
  for (s in safe) {
    sel[[paste0("muP_", s)]] <- 1
    sel[[paste0("muA_", s)]] <- 1
    sel[[paste0("dP_", s)]] <- seq_len(J)
    sel[[paste0("eA_", s)]] <- seq_len(J)
  }
  smp <- INLA::inla.posterior.sample(n_draws, fit, selection = sel,
                                     num.threads = control$num.threads)
  nm <- rownames(smp[[1]]$latent)
  lat <- vapply(smp, function(z) z$latent[, 1], numeric(length(nm)))
  pick <- function(name, idx) lat[match(paste0(name, ":", idx), nm), , drop = FALSE]

  summ <- function(m) data.frame(
    coef = rowMeans(m), stderr = apply(m, 1, stats::sd),
    lower_95 = apply(m, 1, stats::quantile, 0.025, names = FALSE),
    upper_95 = apply(m, 1, stats::quantile, 0.975, names = FALSE))
  info <- design_full$coef_info
  info <- info[info$column != "(Intercept)", , drop = FALSE]
  N_nz <- colSums(!is.na(Y))
  out <- list()
  lambda <- stats::setNames(numeric(K), coef_cols)
  for (k in seq_len(K)) {
    s <- safe[k]
    lam <- vapply(smp, function(z) z$hyperpar[grep(paste0("^Beta for dPc_", s, "$"),
                                                   names(z$hyperpar))], numeric(1))
    lambda[k] <- mean(lam)
    uP <- pick(paste0("dP_", s), seq_len(J))
    e <- pick(paste0("eA_", s), seq_len(J))
    betaP <- sweep(uP, 2, pick(paste0("muP_", s), 1)[1, ], "+")
    betaA <- sweep(sweep(uP * ok_prev, 2, lam, "*") + e, 2,
                   pick(paste0("muA_", s), 1)[1, ], "+")
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
  attr(fit, "bm_lambda") <- lambda
  list(results = do.call(rbind, out), fit = fit, lambda = lambda)
}
