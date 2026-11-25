# Two-part (hurdle) model fitted in a single INLA call per feature.
#
# The data are stacked in blocks, one block of n rows per part:
#   prevalence block: presence (0/1) ~ binomial (logit link), all samples
#   abundance block:  log2 abundance ~ gaussian, non-zero samples only
#                     (absent samples have NA response and do not contribute)
# Every fixed effect, ordered(), strata() and random effect gets its own copy
# per part, so the two parts have separate coefficients. With
# link = "shared", random intercepts of the abundance part are a scaled copy
# of the prevalence ones (INLA's `copy` feature), so the two parts share
# subject-level variation.
#
# The design only depends on the metadata, so it is built once per set of
# parts and reused for every feature; only the response changes.

.bm_parts <- c("prevalence", "abundance")
.bm_prefix <- c(prevalence = "prev", abundance = "abun")

# Builds the stacked design for the requested parts.
.bm_hurdle_design <- function(model, parts = .bm_parts, prior = bm_prior(),
                              link = c("none", "shared")) {
  link <- match.arg(link)
  parts <- intersect(.bm_parts, parts)
  data <- model$data
  parsed <- model$parsed
  n <- nrow(data)
  B <- length(parts)

  # Fixed-effect design matrix (fixed terms + group() factors).
  fixed_rhs <- c(if (parsed$intercept) "1" else "-1", parsed$fixed, parsed$group)
  X <- stats::model.matrix(
    stats::as.formula(paste("~", paste(fixed_rhs, collapse = " + "))), data)
  coef_info <- .bm_coef_info(X, data, c(parsed$fixed, parsed$group))

  safe <- make.unique(gsub("[^A-Za-z0-9_]", "_", colnames(X)))
  block_of <- rep(seq_len(B), each = n)

  stack <- list()
  fixed_names <- list()
  prec <- list(default = 0)
  for (b in seq_len(B)) {
    pre <- .bm_prefix[[parts[b]]]
    nm <- paste0(pre, "_", safe)
    fixed_names[[parts[b]]] <- stats::setNames(nm, colnames(X))
    for (k in seq_len(ncol(X))) {
      v <- numeric(n * B)
      v[block_of == b] <- X[, k]
      stack[[nm[k]]] <- v
      is_int <- colnames(X)[k] == "(Intercept)"
      prec[[nm[k]]] <- if (is_int) {
        if (parts[b] == "prevalence") 0.01 else 1e-4
      } else prior$fixed[[parts[b]]]
    }
  }

  # f() terms, one copy per part with NA index outside its block.
  f_terms <- character()
  ft <- model$terms[!is.na(model$terms$index), , drop = FALSE]
  lincombs <- list()
  lincomb_info <- list()
  for (b in seq_len(B)) {
    pre <- .bm_prefix[[parts[b]]]
    for (i in seq_len(nrow(ft))) {
      idx <- ft$index[i]
      name <- paste0(pre, idx)
      v <- rep(NA_integer_, n * B)
      v[block_of == b] <- data[[idx]]
      stack[[name]] <- v
      if (!is.na(ft$weights[i])) stack[[ft$weights[i]]] <- rep(data[[ft$weights[i]]], B)

      prev_name <- paste0(.bm_prefix[["prevalence"]], idx)
      if (link == "shared" && parts[b] == "abundance" &&
          ft$type[i] == "random_intercept" && "prevalence" %in% parts) {
        f_terms <- c(f_terms, sprintf(
          "f(%s, copy = \"%s\", fixed = FALSE, hyper = .bm_hyper$copy)",
          name, prev_name))
      } else {
        f_terms <- c(f_terms, .bm_render_f(ft[i, , drop = FALSE], index = name))
      }

      # ordered(): differences between consecutive levels, as in maaslin3.
      if (ft$type[i] == "ordered") {
        levs <- levels(factor(model$data[[ft$variable[i]]]))
        K <- length(levs)
        for (k in 2:K) {
          w <- numeric(K)
          w[k - 1] <- -1
          w[k] <- 1
          lc_name <- sprintf("%s__%s__%d", pre, ft$variable[i], k)
          lc <- do.call(INLA::inla.make.lincomb, stats::setNames(list(w), name))
          names(lc) <- lc_name
          lincombs <- c(lincombs, lc)
          lincomb_info[[length(lincomb_info) + 1]] <- data.frame(
            lincomb = lc_name, part = parts[b], metadata = ft$variable[i],
            value = levs[k], name = paste0(ft$variable[i], levs[k]),
            stringsAsFactors = FALSE)
        }
      }
    }
  }

  rhs <- c("-1", unlist(fixed_names, use.names = FALSE), f_terms)
  hyper <- prior$hyper
  hyper$copy <- list(beta = list(prior = "normal", param = c(0, 1)))
  env <- .bm_formula_env(list(hyper = hyper))
  formula <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")),
                               env = env)

  list(
    formula = formula, stack = stack, parts = parts, n = n, X = X,
    fixed_names = fixed_names, coef_info = coef_info, prec = prec,
    lincombs = if (length(lincombs)) lincombs else NULL,
    lincomb_info = if (length(lincomb_info)) do.call(rbind, lincomb_info) else NULL,
    link = link
  )
}

# Maps model-matrix columns to maaslin3's metadata / value / name columns.
.bm_coef_info <- function(X, data, term_labels) {
  assign <- attr(X, "assign")
  labs <- c("(Intercept)", term_labels)[assign + 1]
  value <- colnames(X)
  for (k in seq_along(value)) {
    lab <- labs[k]
    if (lab %in% colnames(data) && (is.factor(data[[lab]]) ||
                                    is.character(data[[lab]]) ||
                                    is.logical(data[[lab]]))) {
      value[k] <- substring(colnames(X)[k], nchar(lab) + 1)
    } else if (lab %in% colnames(data)) {
      value[k] <- lab
    }
  }
  data.frame(column = colnames(X), metadata = labs, value = value,
             name = colnames(X), stringsAsFactors = FALSE)
}

# Fits one feature. `abundance` is the transformed abundance with NA meaning
# absent. Returns the INLA fit (or an error string).
.bm_fit_feature <- function(abundance, design, control) {
  n <- design$n
  B <- length(design$parts)
  Y <- matrix(NA_real_, n * B, B)
  for (b in seq_len(B)) {
    rows <- (b - 1) * n + seq_len(n)
    Y[rows, b] <- if (design$parts[b] == "prevalence") {
      as.numeric(!is.na(abundance))
    } else {
      abundance
    }
  }
  stack <- design$stack
  stack$Y <- if (B == 1) as.vector(Y) else Y
  family <- ifelse(design$parts == "prevalence", "binomial", "gaussian")
  # INLA evaluates Ntrials inside `data`, so the vector must live there.
  .bm_ntrials <- rep(1, n * B)
  stack$.bm_ntrials <- .bm_ntrials

  tryCatch(
    INLA::inla(
      design$formula, data = stack, family = family,
      Ntrials = .bm_ntrials,
      control.fixed = list(mean = 0, prec = design$prec),
      lincomb = design$lincombs,
      control.compute = list(return.marginals = TRUE),
      num.threads = control$num.threads,
      verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
}

# Turns one INLA fit into maaslin3-style rows (one per coefficient and part).
# With n_draws > 0 the result carries an attribute "draws": a matrix with
# one row per result row and n_draws posterior samples, used by the median
# comparison.
.bm_summarise_fit <- function(fit, design, feature, abundance, part_errors,
                              n_draws = 0) {
  N <- length(abundance)
  N_nz <- sum(!is.na(abundance))
  rows <- list()
  draws <- list()

  for (part in design$parts) {
    err <- part_errors[[part]]
    nm <- design$fixed_names[[part]]
    info <- design$coef_info
    keep <- info$column != "(Intercept)"
    info <- info[keep, , drop = FALSE]
    nm <- nm[keep]
    pieces <- list()
    piece_draws <- list()
    if (nrow(info)) {
      fo <- data.frame(feature = feature, info[, c("metadata", "value", "name")],
                       stringsAsFactors = FALSE)
      st <- .bm_marginal_stats(fit$marginals.fixed[nm],
                               fit$summary.fixed[nm, , drop = FALSE], n_draws)
      fo <- cbind(fo, st)
      piece_draws$fixed <- attr(st, "draws")
      # Abundance coefficients of a level never seen among non-zero samples
      # are only informed by the prior.
      fo$error <- NA_character_
      if (part == "abundance") {
        Xnz <- design$X[!is.na(abundance), info$column, drop = FALSE]
        unseen <- colSums(abs(Xnz)) == 0
        fo$error[unseen] <-
          "No non-zero samples for this level: estimate is prior-driven"
      }
      pieces$fixed <- fo
    }

    li <- design$lincomb_info
    if (!is.null(li)) li <- li[li$part == part, , drop = FALSE]
    if (!is.null(li) && nrow(li)) {
      lo <- data.frame(feature = feature, li[, c("metadata", "value", "name")],
                       stringsAsFactors = FALSE)
      st <- .bm_marginal_stats(
        fit$marginals.lincomb.derived[li$lincomb],
        fit$summary.lincomb.derived[li$lincomb, , drop = FALSE], n_draws)
      lo <- cbind(lo, st)
      piece_draws$ordered <- attr(st, "draws")
      lo$error <- NA_character_
      pieces$ordered <- lo
    }
    if (!length(pieces)) next
    out <- do.call(rbind, unname(pieces))

    if (!is.null(err)) out$error <- err
    out$model <- part
    out$N <- N
    out$N_not_zero <- N_nz
    rows[[part]] <- out
    draws[[part]] <- do.call(rbind, unname(piece_draws))
  }
  res <- do.call(rbind, unname(rows))
  if (n_draws > 0) attr(res, "draws") <- do.call(rbind, unname(draws))
  res
}

.bm_marginal_stats <- function(marginals, summary, n_draws = 0) {
  p_pos <- vapply(marginals, function(m) {
    1 - INLA::inla.pmarginal(0, m)
  }, numeric(1))
  p_pos <- pmin(pmax(p_pos, 0), 1)
  out <- data.frame(
    coef = summary[, "mean"],
    stderr = summary[, "sd"],
    lower_95 = summary[, "0.025quant"],
    upper_95 = summary[, "0.975quant"],
    prob_positive = unname(p_pos),
    lfsr = unname(pmin(p_pos, 1 - p_pos)),
    row.names = NULL
  )
  if (n_draws > 0) {
    attr(out, "draws") <- t(vapply(marginals, function(m) {
      INLA::inla.rmarginal(n_draws, m)
    }, numeric(n_draws)))
  }
  out
}

# Rows for a part that could not be fitted for a feature.
.bm_empty_rows <- function(design_full, part, feature, abundance, msg) {
  info <- design_full$coef_info
  info <- info[info$column != "(Intercept)", , drop = FALSE]
  li <- design_full$lincomb_info
  if (!is.null(li)) {
    li <- li[li$part == part, c("metadata", "value", "name"), drop = FALSE]
    info <- rbind(info[, c("metadata", "value", "name")], li)
  } else {
    info <- info[, c("metadata", "value", "name")]
  }
  if (!nrow(info)) return(NULL)
  data.frame(feature = feature, info, coef = NA_real_, stderr = NA_real_,
             lower_95 = NA_real_, upper_95 = NA_real_, prob_positive = NA_real_,
             lfsr = NA_real_, error = msg, model = part,
             N = length(abundance), N_not_zero = sum(!is.na(abundance)),
             stringsAsFactors = FALSE)
}

#' Bayesian q-values from local false sign rates
#'
#' For each result, the q-value is the expected proportion of sign errors
#' among all results with a local false sign rate (lfsr) at most as large
#' (the running mean of the sorted lfsr values). Selecting the results with
#' `qval <= 0.1` controls the Bayesian false sign rate at 10%.
#'
#' @param lfsr Numeric vector of local false sign rates; `NA` values get
#'   `NA`.
#' @return Numeric vector of q-values, same length as `lfsr`.
#' @examples
#' bm_qvalue(c(0.001, 0.2, 0.01, NA))
#' @export
bm_qvalue <- function(lfsr) {
  q <- rep(NA_real_, length(lfsr))
  ok <- which(!is.na(lfsr))
  if (!length(ok)) return(q)
  o <- ok[order(lfsr[ok])]
  q[o] <- cumsum(lfsr[o]) / seq_along(o)
  q
}

# Adds individual and joint q-values according to the selection rule.
#   "bh_lfsr" (default): Benjamini-Hochberg on 2 * lfsr, within each model
#             part. Joint ("associated in at least one part"): as MaAsLin3,
#             pbeta(min(p_prevalence, p_abundance), 1, 2), then BH.
#   "rope":   probability that the effect is practically null,
#             P(|delta| <= rope_delta), from a normal approximation of the
#             posterior of the corrected effect; q-value = running mean of
#             the sorted probabilities (Bayesian FDR for relevant effects).
#             Joint: the smaller probability of the two parts.
#   "lfsr":   running mean of the sorted lfsr (Bayesian false sign rate).
#             Joint: the smaller lfsr of the two parts.
# lfsr_joint is always the smaller lfsr of the two parts.
.bm_add_qvalues <- function(res, selection = c("bh_lfsr", "rope", "lfsr"),
                            rope_delta = 0.5, subtract_median = FALSE) {
  selection <- match.arg(selection)
  key <- paste(res$feature, res$name, sep = "\r")
  min_by_key <- function(v) {
    m <- tapply(v, key, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
    m
  }
  lj <- min_by_key(res$lfsr)
  res$lfsr_joint <- unname(lj[key])

  score <- switch(selection,
    bh_lfsr = pmin(2 * res$lfsr, 1),
    lfsr = res$lfsr,
    rope = {
      # Corrected effect: coef already corrected if subtract_median.
      m <- res$coef - if (subtract_median) 0 else res$null_hypothesis
      stats::pnorm(rope_delta, m, res$stderr) - stats::pnorm(-rope_delta, m, res$stderr)
    })
  adjust <- function(p) {
    if (selection == "bh_lfsr") {
      q <- rep(NA_real_, length(p)); ok <- !is.na(p)
      q[ok] <- stats::p.adjust(p[ok], "BH")
      q
    } else {
      bm_qvalue(p)
    }
  }
  res$qval_individual <- NA_real_
  for (part in unique(res$model)) {
    i <- res$model == part
    res$qval_individual[i] <- adjust(score[i])
  }
  sj <- min_by_key(score)
  if (selection == "bh_lfsr") {
    n_parts <- tapply(!is.na(score), key, sum)
    sj <- ifelse(n_parts[names(sj)] >= 2, stats::pbeta(sj, 1, 2), sj)
  }
  qj <- stats::setNames(adjust(unname(sj)), names(sj))
  res$qval_joint <- unname(qj[key])
  res
}
