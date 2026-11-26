#' Plot bayome results
#'
#' Plots built with ggplot2 from a [bayome()] fit or from its results table
#' (for example read back from `all_results.tsv`). The layout and colours
#' follow the MaAsLin3 plots, so that results of both look familiar; the
#' intervals are 95% posterior credible intervals.
#'
#' * `bm_plot_summary()`: summary of the top features. Left, a coefficient
#'   plot for the main coefficients (`coef_names`): one point per feature and
#'   model part (circles: prevalence; triangles: abundance), filled when the
#'   association is selected, with its 95% credible interval and the
#'   reference it is compared with (dashed: prevalence; solid: abundance).
#'   Right, a heatmap of the other coefficients (`heatmap_names`), marked
#'   `*` when `qval_individual <= threshold` and `**` when
#'   `qval_individual <= threshold / 10`.
#' * `bm_plot_effects()`: coefficient plot of one coefficient, prevalence
#'   and abundance in separate panels.
#' * `bm_plot_link()`: prevalence effect against abundance effect of one
#'   coefficient, one point per feature, with the estimated link `lambda`
#'   in the subtitle when available.
#' * `bm_plot_association()`: the data behind one feature and one metadata
#'   variable. Prevalence: share of samples where the feature is present in
#'   each level (or, for a continuous variable, the variable in samples
#'   where the feature is present or absent). Abundance: the feature's
#'   abundance where present, by level or against the variable.
#'
#' Features are ranked by `qval_joint` when available, otherwise by their
#' smallest `qval_individual`. Rows with an `error` are left out.
#'
#' @param x A `bayome` object, or a data frame with the bayome results
#'   columns.
#' @param n_top Number of features to show.
#' @param coef_names Coefficients (values of the `name` column, such as
#'   `"diagnosisCD"`) for the coefficient plot. Default: the coefficients
#'   of the `exposure` for a joint fit; otherwise the two with the smallest
#'   mean log q-value among the top features.
#' @param heatmap_names Coefficients for the heatmap. Default: all the other
#'   coefficients.
#' @param name One coefficient (value of the `name` column).
#' @param threshold Q-value threshold for selected associations. Default:
#'   `max_significance` of the fit, or 0.1 for a results table.
#' @param feature Feature to plot.
#' @param variable Metadata variable to plot (a column of `metadata`).
#' @param data Abundance table used by the model: samples in rows, features
#'   in columns, absences as `NA`. Default: `x$data`, the processed table of
#'   the fit (log2 scale with the default `transform = "LOG"`).
#' @param metadata Sample metadata, rows matching `data`. Default:
#'   `x$metadata` (continuous variables are standardised there when
#'   `standardize = TRUE`).
#' @return A `ggplot` object; `bm_plot_summary()` and
#'   `bm_plot_association()` return two plots combined with patchwork.
#' @examples
#' res <- data.frame(
#'   feature = rep(c("taxon_a", "taxon_b", "taxon_c"), each = 4),
#'   name = rep(c("diagnosisCD", "age"), each = 2, times = 3),
#'   model = rep(c("prevalence", "abundance"), 6),
#'   coef = c(1.2, -0.8, 0.1, 0.2, 0.1, 0.05, -0.3, 0.1, -2, -1.5, 0.4, 0.6),
#'   qval_individual = c(0.01, 0.04, 0.5, 0.3, 0.8, 0.9, 0.6, 0.7, 0.001,
#'                       0.002, 0.2, 0.05),
#'   error = NA_character_)
#' res$lower_95 <- res$coef - 0.7
#' res$upper_95 <- res$coef + 0.7
#' bm_plot_summary(res)
#' bm_plot_effects(res, "diagnosisCD")
#' bm_plot_link(res, "diagnosisCD")
#' @name bm_plot
NULL

# Colours of the MaAsLin3 plots.
.bm_col <- c(prevalence = "#008B8B", abundance = "#8B008B")
.bm_shape <- c(prevalence = 21, abundance = 24)
.bm_heat <- c("#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B")

#' @rdname bm_plot
#' @export
bm_plot_summary <- function(x, n_top = 30, coef_names = NULL,
                            heatmap_names = NULL, threshold = NULL) {
  res <- .bm_plot_results(x)
  threshold <- .bm_plot_threshold(x, threshold)
  if (!is.null(coef_names) || !is.null(heatmap_names)) {
    wanted <- c(coef_names, heatmap_names)
    miss <- setdiff(wanted, res$name)
    if (length(miss)) {
      stop(sprintf("Coefficient(s) not found in the results: %s.",
                   paste(miss, collapse = ", ")), call. = FALSE)
    }
    # With only coef_names, the other coefficients still go to the heatmap.
    if (!is.null(heatmap_names)) res <- res[res$name %in% wanted, , drop = FALSE]
  }
  top <- .bm_top_features(res, n_top)
  res <- res[res$feature %in% top, , drop = FALSE]

  # Default: the exposure of a joint fit goes to the coefficient plot.
  exposure <- if (inherits(x, "bayome")) x$settings$exposure else NULL
  if (is.null(coef_names) && length(exposure) && "metadata" %in% colnames(res)) {
    coef_names <- setdiff(unique(res$name[res$metadata %in% exposure]), heatmap_names)
    if (!length(coef_names)) coef_names <- NULL
  }
  if (is.null(coef_names)) {
    q <- pmax(.bm_rank_q(res), 1e-300)
    mlq <- sort(tapply(log(q), res$name, mean, na.rm = TRUE))
    coef_names <- utils::head(setdiff(names(mlq), heatmap_names), 2)
  }
  if (is.null(heatmap_names)) heatmap_names <- setdiff(unique(res$name), coef_names)

  # Features ordered by their effect on the first coefficient plotted.
  first <- res[res$name == c(coef_names, heatmap_names)[1], , drop = FALSE]
  ord <- tapply(first$coef, first$feature, mean)
  ord <- c(names(sort(ord)), setdiff(top, names(ord)))
  res$feature <- factor(res$feature, levels = ord)

  p1 <- if (length(coef_names)) {
    .bm_coef_plot(res[res$name %in% coef_names, , drop = FALSE], coef_names,
                  threshold)
  }
  p2 <- if (length(heatmap_names)) {
    .bm_heatmap_plot(res[res$name %in% heatmap_names, , drop = FALSE],
                     heatmap_names, threshold)
  }
  if (is.null(p1)) return(p2)
  if (is.null(p2)) return(p1)
  p2 <- p2 + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_blank(),
                            axis.ticks.y = ggplot2::element_blank())
  # The heatmap needs room for its two panel titles ("Prevalence",
  # "Abundance") even with few coefficients.
  patchwork::wrap_plots(p1, p2, widths = c(1.2 * length(coef_names),
                                           max(1.6, length(heatmap_names) / 2.5)),
                        guides = "collect")
}

#' @rdname bm_plot
#' @export
bm_plot_effects <- function(x, name, n_top = 25, threshold = NULL) {
  res <- .bm_plot_results(x, name)
  threshold <- .bm_plot_threshold(x, threshold)
  top <- .bm_top_features(res, n_top)
  res <- res[res$feature %in% top, , drop = FALSE]
  ord <- tapply(res$coef, res$feature, mean)
  res$feature <- factor(res$feature, levels = names(sort(ord)))
  res$panel <- factor(res$model, levels = c("prevalence", "abundance"),
                      labels = c("Prevalence (log-odds)", "Abundance (log2)"))
  .bm_coef_plot(res, name, threshold, facet = "panel") +
    ggplot2::labs(x = sprintf("%s coefficient (95%% credible interval)", name))
}

#' @rdname bm_plot
#' @export
bm_plot_link <- function(x, name, threshold = NULL) {
  res <- .bm_plot_results(x, name)
  threshold <- .bm_plot_threshold(x, threshold)
  keep <- c("feature", "coef", "qval_individual")
  prev <- res[res$model == "prevalence", keep]
  abun <- res[res$model == "abundance", keep]
  d <- merge(prev, abun, by = "feature", suffixes = c("_prev", "_abun"))
  if (!nrow(d)) {
    stop("No feature has both a prevalence and an abundance estimate.", call. = FALSE)
  }
  sel_p <- !is.na(d$qval_individual_prev) & d$qval_individual_prev <= threshold
  sel_a <- !is.na(d$qval_individual_abun) & d$qval_individual_abun <= threshold
  d$selected <- factor(ifelse(sel_p & sel_a, "Both", ifelse(sel_p, "Prevalence",
                              ifelse(sel_a, "Abundance", "None"))),
                       levels = c("Both", "Prevalence", "Abundance", "None"))
  lambda <- if (inherits(x, "bayome") && !is.null(x$lambda)) x$lambda[name] else NA
  subtitle <- NULL
  if (!is.na(lambda)) subtitle <- sprintf("Estimated link lambda = %.2f", lambda)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$coef_prev, y = .data$coef_abun,
                                  fill = .data$selected)) +
    ggplot2::geom_hline(yintercept = 0, colour = "darkgray", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "darkgray", linetype = "dashed") +
    ggplot2::geom_point(shape = 21, size = 2.5, colour = "black", stroke = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(Both = "#B2182B", Prevalence = .bm_col[["prevalence"]],
                 Abundance = .bm_col[["abundance"]], None = "white"),
      name = sprintf("Selected (q <= %g)", threshold), drop = FALSE) +
    ggplot2::labs(x = sprintf("Prevalence coefficient of %s (log-odds)", name),
                  y = sprintf("Abundance coefficient of %s (log2)", name),
                  subtitle = subtitle) +
    .bm_theme()
}

#' @rdname bm_plot
#' @export
bm_plot_association <- function(x, feature, variable,
                                data = x$data, metadata = x$metadata) {
  if (is.null(data) || is.null(metadata)) {
    stop("`data` and `metadata` are needed: pass a bayome object or give them.",
         call. = FALSE)
  }
  if (!feature %in% colnames(data)) {
    stop(sprintf("Feature '%s' not found in `data`.", feature), call. = FALSE)
  }
  if (!variable %in% colnames(metadata)) {
    stop(sprintf("Variable '%s' not found in `metadata`.", variable), call. = FALSE)
  }
  common <- intersect(rownames(data), rownames(metadata))
  y <- data[common, feature]
  v <- metadata[common, variable]
  categorical <- is.factor(v) || is.character(v) || is.logical(v)
  if (categorical) v <- factor(v)
  status <- factor(ifelse(is.na(y), "Absent", "Present"),
                   levels = c("Present", "Absent"))

  res <- tryCatch(.bm_plot_results(x), error = function(e) NULL)
  label <- function(part) {
    if (is.null(res) || !"metadata" %in% colnames(res)) return(NULL)
    r <- res[res$feature == feature & res$metadata == variable &
               res$model == part, , drop = FALSE]
    if (!nrow(r)) return(NULL)
    paste(sprintf("%s: coef = %.2f, q = %s", r$name, r$coef,
                  formatC(r$qval_individual, format = "e", digits = 1)),
          collapse = "\n")
  }
  note <- function(p, part) {
    lab <- label(part)
    if (is.null(lab)) return(p)
    p + ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.1,
                          label = lab, size = 3, fontface = "italic")
  }

  # Prevalence
  if (categorical) {
    d <- as.data.frame(table(level = v, status = status))
    d$share <- d$Freq / stats::ave(d$Freq, d$level, FUN = sum)
    # Present at the bottom of each bar, labelled in its middle.
    d$status <- factor(d$status, levels = c("Absent", "Present"))
    pres <- d[d$status == "Present", ]
    p_prev <- ggplot2::ggplot(d, ggplot2::aes(x = .data$level, y = .data$share,
                                              fill = .data$status)) +
      ggplot2::geom_col(colour = "black", linewidth = 0.2, width = 0.7) +
      ggplot2::geom_text(data = pres[pres$share > 0.06, ],
                         ggplot2::aes(y = .data$share / 2,
                                      label = sprintf("%.0f%%", 100 * .data$share)),
                         size = 3, colour = "white") +
      ggplot2::scale_fill_manual(values = c(Present = .bm_col[["prevalence"]],
                                            Absent = "grey85"), name = NULL) +
      ggplot2::scale_y_continuous(labels = function(b) paste0(100 * b, "%"),
                                  breaks = seq(0, 1, 0.25),
                                  expand = ggplot2::expansion(mult = c(0, 0.25))) +
      ggplot2::labs(x = variable, y = "Samples")
  } else {
    d <- data.frame(status = status, value = v)
    p_prev <- ggplot2::ggplot(d, ggplot2::aes(x = .data$status, y = .data$value,
                                              fill = .data$status)) +
      ggplot2::geom_boxplot(outlier.alpha = 0, alpha = 0.5, show.legend = FALSE) +
      ggplot2::geom_jitter(shape = 21, size = 1, stroke = 0.15, alpha = 0.75,
                           width = 0.2, height = 0, show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c(Present = .bm_col[["prevalence"]],
                                            Absent = "grey85")) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.6, 1.5))) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = variable)
  }
  p_prev <- note(p_prev + ggplot2::ggtitle("Prevalence") + .bm_theme(), "prevalence")

  # Abundance where present
  d <- data.frame(value = v, y = y)[!is.na(y), ]
  if (categorical) {
    p_abun <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$y,
                                              fill = .data$value)) +
      ggplot2::geom_boxplot(outlier.alpha = 0, alpha = 0.5, show.legend = FALSE) +
      ggplot2::geom_jitter(shape = 21, size = 1, stroke = 0.15, alpha = 0.75,
                           width = 0.2, height = 0, show.legend = FALSE) +
      ggplot2::scale_fill_brewer(palette = "Spectral")
  } else {
    p_abun <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$y)) +
      ggplot2::geom_point(fill = "darkolivegreen4", colour = "black", alpha = 0.5,
                          shape = 21, size = 1, stroke = 0.15) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.5,
                           colour = "blue")
  }
  p_abun <- p_abun +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25))) +
    ggplot2::labs(x = variable, y = "Abundance where present") +
    ggplot2::ggtitle("Abundance") + .bm_theme()
  p_abun <- note(p_abun, "abundance")

  patchwork::wrap_plots(p_prev, p_abun, ncol = 2) +
    patchwork::plot_annotation(title = feature)
}

# --- internal helpers -------------------------------------------------------

.bm_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "transparent"),
                   legend.background = ggplot2::element_rect(fill = "transparent"))
}

# Coefficient plot: points by model part (shape and colour), filled when
# selected, with 95% credible intervals and the reference line.
.bm_coef_plot <- function(res, names, threshold, facet = "name") {
  res$model <- factor(res$model, levels = c("prevalence", "abundance"))
  res$selected <- !is.na(res$qval_individual) & res$qval_individual <= threshold
  res$fill <- ifelse(res$selected, as.character(res$model), "none")
  res[[facet]] <- if (facet == "name") factor(res$name, levels = names) else res[[facet]]
  ref <- if ("null_hypothesis" %in% colnames(res)) res$null_hypothesis else 0
  ref[is.na(ref)] <- 0
  refs <- unique(data.frame(facet = res[[facet]], model = res$model, ref = ref))
  names(refs)[1] <- facet
  has_ci <- all(c("lower_95", "upper_95") %in% colnames(res))
  dodge <- ggplot2::position_dodge(width = 0.6)

  p <- ggplot2::ggplot(res, ggplot2::aes(x = .data$coef, y = .data$feature,
                                         group = .data$model)) +
    ggplot2::geom_vline(data = refs, ggplot2::aes(xintercept = .data$ref,
                                                  linetype = .data$model),
                        colour = "darkgray")
  if (has_ci) {
    p <- p + ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lower_95,
                                                  xmax = .data$upper_95,
                                                  colour = .data$model),
                                     position = dodge, show.legend = FALSE)
  }
  p +
    ggplot2::geom_point(ggplot2::aes(shape = .data$model, fill = .data$fill,
                                     colour = .data$model),
                        size = 3, stroke = 0.6, position = dodge) +
    ggplot2::scale_shape_manual(values = .bm_shape, name = "Association",
                                labels = c("Prevalence", "Abundance"), drop = FALSE) +
    ggplot2::scale_colour_manual(values = .bm_col, guide = "none") +
    ggplot2::scale_fill_manual(values = c(.bm_col, none = "white"),
                               breaks = c("prevalence", "abundance"),
                               labels = c("Prevalence", "Abundance"),
                               name = sprintf("Selected (q <= %g)", threshold),
                               drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c(prevalence = "dashed",
                                              abundance = "solid"),
                                   labels = c("Prevalence", "Abundance"),
                                   name = "Reference", drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(
      shape = unname(.bm_shape), colour = unname(.bm_col)))) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), scales = "free_x",
                        nrow = 1) +
    ggplot2::labs(x = "Coefficient (95% credible interval)", y = "Feature") +
    .bm_theme()
}

# Heatmap with binned coefficients and q-value stars, as in MaAsLin3.
.bm_heatmap_plot <- function(res, names, threshold) {
  q <- res$qval_individual
  res$star <- ifelse(!is.na(q) & q <= threshold / 10, "**",
                     ifelse(!is.na(q) & q <= threshold, "*", ""))
  thr <- round(max(abs(stats::quantile(res$coef, c(0.1, 0.9), na.rm = TRUE))) / 10,
               1) * 5
  if (!is.finite(thr) || thr == 0) thr <- 0.5
  br <- c(-Inf, -thr, -thr / 2, 0, thr / 2, thr, Inf)
  labs <- c(sprintf("(-Inf, %g]", -thr), sprintf("(%g, %g]", -thr, -thr / 2),
            sprintf("(%g, 0]", -thr / 2), sprintf("(0, %g]", thr / 2),
            sprintf("(%g, %g]", thr / 2, thr), sprintf("(%g, Inf)", thr))
  res$bin <- cut(res$coef, br, labels = labs)
  grid <- expand.grid(feature = levels(res$feature), name = names,
                      model = c("prevalence", "abundance"),
                      stringsAsFactors = FALSE)
  d <- merge(grid, res[, c("feature", "name", "model", "bin", "star")],
             by = c("feature", "name", "model"), all.x = TRUE)
  d$feature <- factor(d$feature, levels = levels(res$feature))
  d$name <- factor(d$name, levels = names)
  d$model <- factor(d$model, levels = c("prevalence", "abundance"),
                    labels = c("Prevalence", "Abundance"))
  d$star[is.na(d$star)] <- ""

  ggplot2::ggplot(d, ggplot2::aes(x = .data$name, y = .data$feature)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$bin), colour = "white",
                       linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$star), size = 4, vjust = 0.75) +
    ggplot2::scale_fill_manual(values = stats::setNames(.bm_heat, labs),
                               na.value = "#EEEEEE", drop = FALSE,
                               name = "Coefficient") +
    ggplot2::facet_grid(~ model) +
    ggplot2::labs(x = NULL, y = "Feature",
                  caption = sprintf("* q <= %g, ** q <= %g", threshold,
                                    threshold / 10)) +
    .bm_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1),
                   panel.spacing = ggplot2::unit(0, "lines"))
}

.bm_plot_results <- function(x, name = NULL) {
  res <- if (inherits(x, "bayome")) x$results else x
  if (!is.data.frame(res)) {
    stop("`x` must be a bayome object or a results data frame.", call. = FALSE)
  }
  need <- c("feature", "name", "model", "coef", "qval_individual")
  miss <- setdiff(need, colnames(res))
  if (length(miss)) {
    stop(sprintf("Results are missing column(s): %s.", paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  if ("error" %in% colnames(res)) {
    res <- res[is.na(res$error) | res$error == "", , drop = FALSE]
  }
  res <- res[is.finite(res$coef), , drop = FALSE]
  if (!is.null(name)) {
    res <- res[res$name == name, , drop = FALSE]
    if (!nrow(res)) {
      stop(sprintf("Coefficient '%s' not found in the results.", name), call. = FALSE)
    }
  }
  res
}

.bm_plot_threshold <- function(x, threshold) {
  if (!is.null(threshold)) return(threshold)
  if (inherits(x, "bayome") && !is.null(x$settings$max_significance)) {
    return(x$settings$max_significance)
  }
  0.1
}

# q-value used to rank: qval_joint when present, otherwise qval_individual.
.bm_rank_q <- function(res) {
  q <- if ("qval_joint" %in% colnames(res)) res$qval_joint else res$qval_individual
  q[is.na(q)] <- res$qval_individual[is.na(q)]
  q[is.na(q)] <- Inf
  q
}

.bm_top_features <- function(res, n_top) {
  q <- .bm_rank_q(res)
  best <- tapply(q, res$feature, min)
  big <- tapply(abs(res$coef), res$feature, max)
  o <- order(best, -big[names(best)])
  utils::head(names(best)[o], n_top)
}
