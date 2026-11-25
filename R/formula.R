# Translation of MaAsLin3 / lme4-style formulas into INLA formulas.
#
# Supported right-hand-side terms:
#   x, x1 * x2, ...        fixed effects (passed to INLA unchanged)
#   group(x)               categorical variable tested jointly (kept as a
#                          fixed factor; joint inference happens afterwards)
#   ordered(x)             ordered categorical -> f(x, model = "rw1")
#   strata(x)              matched design -> iid effect with fixed, very
#                          small precision (conditional-logistic style)
#   (1 | g)                random intercept -> f(g, model = "iid")
#   (0 + x | g), (x || g)  uncorrelated random slope -> f(g_copy, x, model = "iid")
#   (x | g)                correlated intercept + slope: not supported yet

.bm_specials <- c("group", "ordered", "strata")

#' Parse a MaAsLin3-style formula
#'
#' Splits a one-sided formula (as used by `maaslin3()`) into fixed effects,
#' the MaAsLin3 specials `group()`, `ordered()` and `strata()`, and
#' lme4-style random-effect terms.
#'
#' @param formula A formula or a character string such as
#'   `"~ diagnosis + age + (1 | subject)"`. A left-hand side, if present, is
#'   ignored.
#'
#' @return An object of class `bm_formula`: a list with elements `fixed`
#'   (character vector of fixed-effect term labels), `intercept` (logical),
#'   `group`, `ordered`, `strata` (character vectors of variable names) and
#'   `random` (a data frame with one row per random-effect term, with columns
#'   `group` and `slope`; `slope` is `NA` for a random intercept).
#'
#' @examples
#' bm_parse_formula(~ diagnosis + age + ordered(stage) + (1 | subject))
#' @export
bm_parse_formula <- function(formula) {
  formula <- .bm_as_rhs_formula(formula)

  # Random-effect terms; findbars() already expands (x || g) into
  # (1 | g) + (0 + x | g).
  bars <- reformulas::findbars(formula)
  random <- do.call(rbind, c(
    list(data.frame(group = character(), slope = character(),
                    stringsAsFactors = FALSE)),
    lapply(bars, .bm_parse_bar)
  ))
  random <- unique(random)
  rownames(random) <- NULL

  fixed_formula <- reformulas::nobars(formula)
  if (is.null(fixed_formula)) {
    fixed_formula <- ~ 1
  } else if (!inherits(fixed_formula, "formula")) {
    fixed_formula <- stats::as.formula(paste("~", deparse1(fixed_formula)))
  }
  tt <- stats::terms(fixed_formula)
  labels <- attr(tt, "term.labels")

  special <- list(group = character(), ordered = character(),
                  strata = character())
  is_special <- logical(length(labels))
  for (i in seq_along(labels)) {
    m <- regmatches(labels[i],
                    regexec("^(group|ordered|strata)\\((.+)\\)$", labels[i]))[[1]]
    if (length(m) == 3) {
      var <- trimws(m[3])
      if (!grepl("^[.A-Za-z][.A-Za-z0-9_]*$", var)) {
        stop(sprintf("%s() must contain a single variable name, got '%s'.",
                     m[2], var), call. = FALSE)
      }
      special[[m[2]]] <- c(special[[m[2]]], var)
      is_special[i] <- TRUE
    }
  }
  if (length(special$strata) > 1) {
    stop("Only one strata() term is allowed.", call. = FALSE)
  }

  structure(
    list(
      fixed = labels[!is_special],
      intercept = attr(tt, "intercept") == 1,
      group = special$group,
      ordered = special$ordered,
      strata = special$strata,
      random = random
    ),
    class = "bm_formula"
  )
}

#' Build a MaAsLin3-style formula from separate effect arguments
#'
#' Mirrors how `maaslin3()` combines `fixed_effects`, `random_effects`,
#' `group_effects`, `ordered_effects` and `strata_effects` when no
#' `formula` is given.
#'
#' @param fixed_effects,random_effects,group_effects,ordered_effects,strata_effects
#'   Character vectors of metadata column names, or `NULL`.
#' @return A one-sided formula.
#' @examples
#' bm_compose_formula(fixed_effects = c("diagnosis", "age"),
#'                    random_effects = "subject")
#' @export
bm_compose_formula <- function(fixed_effects = NULL, random_effects = NULL,
                               group_effects = NULL, ordered_effects = NULL,
                               strata_effects = NULL) {
  terms <- c(
    fixed_effects,
    if (length(group_effects)) paste0("group(", group_effects, ")"),
    if (length(ordered_effects)) paste0("ordered(", ordered_effects, ")"),
    if (length(strata_effects)) paste0("strata(", strata_effects, ")"),
    if (length(random_effects)) paste0("(1 | ", random_effects, ")")
  )
  if (!length(terms)) {
    stop("No effects given: supply `formula` or at least one effect argument.",
         call. = FALSE)
  }
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

#' Translate a parsed formula into an INLA model
#'
#' Creates the INLA formula and the index columns it needs. Random effects,
#' `ordered()` and `strata()` terms become `f()` terms with integer indices
#' added to `data`.
#'
#' @param formula A formula, a character string, or a `bm_formula` object
#'   from [bm_parse_formula()].
#' @param data Metadata data frame (samples in rows).
#' @param response Name of the response variable in the INLA formula.
#' @param prior A list of priors, as returned by [bm_prior()].
#'
#' @return An object of class `bm_inla_model`: a list with `formula` (the
#'   INLA formula, whose environment holds the hyperprior lists), `data`
#'   (the metadata with the added index columns), `parsed` (the
#'   `bm_formula`) and `terms` (a data frame describing each term).
#'
#' @examples
#' meta <- data.frame(diagnosis = factor(c("CD", "UC", "nonIBD", "CD")),
#'                    age = c(20, 35, 50, 41),
#'                    subject = c("a", "a", "b", "b"))
#' m <- bm_inla_formula(~ diagnosis + age + (1 | subject), meta)
#' m$formula
#' @export
bm_inla_formula <- function(formula, data, response = "y",
                            prior = bm_prior()) {
  parsed <- if (inherits(formula, "bm_formula")) formula else
    bm_parse_formula(formula)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)

  needed <- unique(c(
    all.vars(stats::as.formula(paste("~", paste(c("1", parsed$fixed),
                                                collapse = " + ")))),
    parsed$group, parsed$ordered, parsed$strata,
    parsed$random$group, stats::na.omit(parsed$random$slope)
  ))
  missing <- setdiff(needed, colnames(data))
  if (length(missing)) {
    stop(sprintf("Effect name(s) not found in metadata: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (any(grepl("^\\.bm_", colnames(data)))) {
    stop("Metadata column names must not start with '.bm_'.", call. = FALSE)
  }

  # One row per term. For f() terms we keep the pieces (index column,
  # weights, model arguments, hyperprior key) so other builders, such as the
  # two-part hurdle model, can re-render them with their own index columns.
  term_info <- list()
  add_info <- function(type, variable, index = NA_character_,
                       weights = NA_character_, f_args = NA_character_,
                       hyper = NA_character_) {
    term_info[[length(term_info) + 1]] <<- data.frame( # nolint: assignment_linter.
      type = type, variable = variable, index = index, weights = weights,
      f_args = f_args, hyper = hyper, stringsAsFactors = FALSE)
  }

  # Fixed part: ordinary fixed effects plus group() variables as factors.
  for (v in parsed$group) data[[v]] <- as.factor(data[[v]])
  fixed_terms <- c(parsed$fixed, parsed$group)
  for (v in parsed$fixed) add_info("fixed", v)
  for (v in parsed$group) add_info("group", v)

  for (v in parsed$ordered) {
    x <- data[[v]]
    if (!is.factor(x)) x <- factor(x)
    if (nlevels(x) < 3) {
      stop(sprintf("ordered(%s) needs at least 3 levels; use it as a fixed effect instead.",
                   v), call. = FALSE)
    }
    idx <- paste0(".bm_ord_", v)
    data[[idx]] <- as.integer(x)
    add_info("ordered", v, idx,
             f_args = "model = \"rw1\", scale.model = TRUE, constr = TRUE",
             hyper = "ordered")
  }

  for (v in parsed$strata) {
    idx <- paste0(".bm_strata_", v)
    data[[idx]] <- as.integer(factor(data[[v]]))
    add_info("strata", v, idx, f_args = "model = \"iid\"", hyper = "strata")
  }

  for (i in seq_len(nrow(parsed$random))) {
    g <- parsed$random$group[i]
    s <- parsed$random$slope[i]
    if (is.na(s)) {
      idx <- paste0(".bm_re_", g)
      data[[idx]] <- as.integer(factor(data[[g]]))
      add_info("random_intercept", g, idx, f_args = "model = \"iid\"",
               hyper = "random")
    } else {
      if (!is.numeric(data[[s]])) {
        stop(sprintf("Random slope '%s' must be numeric (factor slopes are not supported yet).",
                     s), call. = FALSE)
      }
      idx <- paste0(".bm_re_", g, "__", s)
      data[[idx]] <- as.integer(factor(data[[g]]))
      add_info("random_slope", paste0(s, " | ", g), idx, weights = s,
               f_args = "model = \"iid\"", hyper = "random")
    }
  }

  terms <- if (length(term_info)) do.call(rbind, term_info) else
    data.frame(type = character(), variable = character(), index = character(),
               weights = character(), f_args = character(),
               hyper = character(), stringsAsFactors = FALSE)
  is_f <- !is.na(terms$index)
  terms$inla_term <- ifelse(is_f, .bm_render_f(terms), terms$variable)

  rhs <- c(if (parsed$intercept) "1" else "-1", fixed_terms,
           terms$inla_term[is_f])
  inla_formula <- stats::as.formula(
    paste(response, "~", paste(rhs, collapse = " + ")),
    env = .bm_formula_env(prior))

  structure(
    list(formula = inla_formula, data = data, parsed = parsed, terms = terms),
    class = "bm_inla_model"
  )
}

# Renders f() calls for the rows of a term table. `index` overrides the index
# column names (e.g. part-specific copies in the hurdle model).
.bm_render_f <- function(terms, index = terms$index) {
  if (!nrow(terms)) return(character())
  w <- ifelse(is.na(terms$weights), "", paste0(", ", terms$weights))
  sprintf("f(%s%s, %s, hyper = .bm_hyper$%s)", index, w, terms$f_args,
          terms$hyper)
}

.bm_formula_env <- function(prior) {
  env <- new.env(parent = asNamespace("stats"))
  env$.bm_hyper <- prior$hyper
  env
}

#' @export
print.bm_formula <- function(x, ...) {
  cat("<bm_formula>\n")
  cat("  intercept:", x$intercept, "\n")
  show <- function(label, v) {
    cat(sprintf("  %-9s %s\n", label,
                if (length(v)) paste(v, collapse = ", ") else "-"))
  }
  show("fixed:", x$fixed)
  show("group:", x$group)
  show("ordered:", x$ordered)
  show("strata:", x$strata)
  re <- if (nrow(x$random)) {
    ifelse(is.na(x$random$slope), paste0("(1 | ", x$random$group, ")"),
           paste0("(0 + ", x$random$slope, " | ", x$random$group, ")"))
  } else {
    character()
  }
  show("random:", re)
  invisible(x)
}

#' @export
print.bm_inla_model <- function(x, ...) {
  cat("<bm_inla_model>\n")
  cat("INLA formula:\n  ", deparse1(x$formula), "\n", sep = "")
  invisible(x)
}

# --- internal helpers -------------------------------------------------------

.bm_as_rhs_formula <- function(formula) {
  if (is.character(formula)) {
    formula <- tryCatch(stats::as.formula(paste(formula, collapse = " ")),
                        error = function(e) {
                          stop(sprintf("Invalid formula: %s", formula), call. = FALSE)
                        })
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula or a character string.", call. = FALSE)
  }
  # Drop any left-hand side: bayome adds the response per feature.
  stats::as.formula(paste("~", deparse1(formula[[length(formula)]])),
                    env = environment(formula))
}

.bm_parse_bar <- function(bar) {
  grp <- bar[[3]]
  if (!is.name(grp)) {
    stop(sprintf("Random-effect grouping must be a single variable, got '%s'.",
                 deparse1(grp)), call. = FALSE)
  }
  grp <- as.character(grp)
  lhs_terms <- stats::terms(stats::as.formula(paste("~", deparse1(bar[[2]]))))
  slopes <- attr(lhs_terms, "term.labels")
  has_int <- attr(lhs_terms, "intercept") == 1

  if (!has_int && !length(slopes)) {
    stop(sprintf("Empty random-effect term '(%s | %s)'.", deparse1(bar[[2]]), grp),
         call. = FALSE)
  }
  if (has_int && length(slopes)) {
    stop(sprintf(paste0("Correlated random effects '(%s | %s)' are not supported yet. ",
                        "Use '(%s || %s)' for uncorrelated intercept and slopes."),
                 deparse1(bar[[2]]), grp, deparse1(bar[[2]]), grp),
         call. = FALSE)
  }
  if (any(!grepl("^[.A-Za-z][.A-Za-z0-9_]*$", slopes))) {
    stop(sprintf("Random slopes must be single variable names in '(%s | %s)'.",
                 deparse1(bar[[2]]), grp), call. = FALSE)
  }
  data.frame(
    group = grp,
    slope = if (has_int) NA_character_ else slopes,
    stringsAsFactors = FALSE
  )
}
