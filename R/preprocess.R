# Pre-processing that mirrors maaslin3: normalise -> filter -> transform.
# As in maaslin3, values at or below `zero_threshold` become NA after
# normalisation. NA therefore means "absent": the prevalence model uses
# !is.na(x) and the abundance model uses the non-NA values only.

.bm_normalize <- function(features, normalization, zero_threshold) {
  x <- as.matrix(features)
  if (any(x < 0, na.rm = TRUE) && normalization != "NONE") {
    stop("Input data contains negative values; use normalization = \"NONE\".",
         call. = FALSE)
  }
  present <- !is.na(x) & x > zero_threshold
  out <- switch(normalization,
    TSS = {
      tot <- pmax(rowSums(x, na.rm = TRUE), .Machine$double.eps)
      x / tot
    },
    CLR = {
      lgm <- apply(x, 1, function(r) mean(log(r[!is.na(r) & r > zero_threshold])))
      log(x) - lgm
    },
    NONE = x
  )
  out[!present] <- NA
  as.data.frame(out, check.names = FALSE)
}

.bm_filter <- function(features, min_abundance, min_prevalence, max_prevalence,
                       min_variance) {
  x <- as.matrix(features)
  n <- nrow(x)
  n_above <- colSums(!is.na(x) & x > min_abundance)
  keep <- n_above > n * min_prevalence & n_above <= n * max_prevalence
  x <- x[, keep, drop = FALSE]
  v <- apply(x, 2, stats::var, na.rm = TRUE)
  v[is.na(v)] <- 0
  x <- x[, v > min_variance, drop = FALSE]
  as.data.frame(x, check.names = FALSE)
}

.bm_transform <- function(features, transform) {
  x <- as.matrix(features)
  out <- switch(transform,
    LOG = {
      if (any(x <= 0, na.rm = TRUE)) {
        stop("LOG transform needs positive values (CLR output is not; use transform = \"NONE\").",
             call. = FALSE)
      }
      log2(x)
    },
    PLOG = {
      if (any(x < 0, na.rm = TRUE)) {
        stop("PLOG transform needs non-negative values.", call. = FALSE)
      }
      log2(x)
    },
    NONE = x
  )
  as.data.frame(out, check.names = FALSE)
}

# z-scores continuous metadata used as fixed effects (as maaslin3's
# standardize = TRUE). Random-effect groups, ordered() and strata() variables
# are left untouched; random slopes are standardised too.
.bm_standardize <- function(metadata, parsed) {
  vars <- unique(c(
    all.vars(stats::as.formula(paste("~", paste(c("1", parsed$fixed),
                                                collapse = " + ")))),
    stats::na.omit(parsed$random$slope)
  ))
  for (v in vars) {
    x <- metadata[[v]]
    if (is.numeric(x) && length(unique(x[!is.na(x)])) > 2) {
      metadata[[v]] <- as.numeric(scale(x))
    }
  }
  metadata
}
