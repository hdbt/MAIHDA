#' @keywords internal
"_PACKAGE"

#' MAIHDA: Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy
#'
#' The MAIHDA package provides tools for conducting Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy. This approach is useful for examining
#' intersectional inequalities in health and other outcomes.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{make_strata}}: Create intersectional strata from multiple variables
#'   \item \code{\link{fit_maihda}}: Fit multilevel models using lme4 or brms
#'   \item \code{\link{summary_maihda}}: Summarize models with variance partition and stratum estimates
#'   \item \code{\link{predict_maihda}}: Make predictions at individual or stratum level
#'   \item \code{\link{plot_maihda}}: Visualize results with various plot types
#'   \item \code{\link{compare_maihda}}: Compare models with bootstrap confidence intervals
#' }
#'
#' @docType package
#' @name MAIHDA-package
#' @aliases MAIHDA
NULL
