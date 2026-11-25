#' bayome: Bayesian Occurrence-Abundance Microbiome Estimation
#'
#' Bayesian models for microbiome association studies, fitted with 'R-INLA'.
#' The occurrence (presence/absence) and the abundance of each feature are
#' modelled together: abundance effects are linked to occurrence effects and
#' information is shared across features. Users write the same formulas as
#' in `maaslin3()`; bayome translates them into INLA models. See [bayome()],
#' [bm_prior()], [bm_parse_formula()] and [bm_inla_formula()].
#'
#' bayome is an independent package, not part of 'MaAsLin3' or 'bioBakery'.
#' Its interface, pre-processing, two-part structure, median comparison and
#' output format follow 'MaAsLin3' (MIT license); the Bayesian models, the
#' joint occurrence-abundance link, the information sharing across features
#' and the selection rules are new.
#'
#' @references
#' Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA,
#' Thompson KN, Nearing JT, Huttenhower C (2026). MaAsLin 3: refining and
#' extending generalized multivariable linear models for meta-omic
#' association discovery. Nature Methods, 23, 554-564.
#' \doi{10.1038/s41592-025-02923-9}
#'
#' @keywords internal
"_PACKAGE"

# Variables that INLA evaluates inside its `data` argument.
utils::globalVariables(c(".bm_scale", ".bm_ntrials"))
