#' Reconstruct the adjusted-model main-effect terms for a MAIHDA decomposition
#'
#' For each stratum-defining variable, returns the model term to add as an additive
#' fixed main effect in the adjusted model, plus the data augmented with any
#' reconstructed binned factors. A categorical dimension is used directly; a numeric
#' dimension that \code{\link{make_strata}} tertile-binned is reconstructed as the
#' \emph{same} binned factor (using the stored breaks/labels), because
#' \code{make_strata()} bins only a temporary copy and leaves the original numeric
#' column intact -- adding the raw numeric column would wrongly enter a linear term
#' instead of the factor that defines the strata.
#'
#' @param strata_vars Character vector of stratum-defining variable names.
#' @param autobin_info Named list of \code{list(breaks, labels)} per auto-binned
#'   variable (the \code{strata_autobin_info} stored on a \code{maihda_model}).
#' @param data Data frame containing the original stratum-defining columns.
#' @return A list with \code{terms} (character vector of RHS term names) and
#'   \code{data} (the input augmented with any \code{.maihda_dim_*} binned columns).
#' @keywords internal
maihda_adjusted_terms <- function(strata_vars, autobin_info, data) {
  terms <- character(0)
  for (v in strata_vars) {
    if (!is.null(autobin_info) && v %in% names(autobin_info)) {
      # The dimension was auto-binned for the strata; the additive main effect must
      # be the SAME tertile factor (make_strata left the original column numeric).
      info <- autobin_info[[v]]
      new_col <- paste0(".maihda_dim_", v)
      data[[new_col]] <- cut(data[[v]], breaks = info$breaks,
                             include.lowest = TRUE, labels = info$labels)
      terms <- c(terms, new_col)
    } else {
      terms <- c(terms, v)
    }
  }
  list(terms = terms, data = data)
}

#' Build the adjusted-model formula and data for a MAIHDA decomposition
#'
#' Given a fitted null model's formula (in \code{(1 | stratum)} form) and its stored
#' strata metadata, returns the adjusted formula -- the null formula plus the additive
#' main effects of the stratum dimensions -- and the data carrying any reconstructed
#' binned factors. Returns \code{NULL} when fewer than two dimensions are available,
#' because there is no intersection to decompose and the single dimension's main effect
#' would render the stratum random intercept redundant (singular).
#'
#' @param null_formula The null model formula using \code{(1 | stratum)}.
#' @param strata_vars Character vector of stratum-defining variables.
#' @param autobin_info Auto-binning recipe (\code{strata_autobin_info}).
#' @param data The null model's data (\code{original_data}) with the dimension columns.
#' @return A list with \code{formula} and \code{data}, or \code{NULL} if fewer than two
#'   dimensions are available.
#' @keywords internal
#' @importFrom stats update as.formula
maihda_adjusted_formula <- function(null_formula, strata_vars, autobin_info, data) {
  if (is.null(strata_vars) || length(strata_vars) < 2) {
    return(NULL)
  }
  adj <- maihda_adjusted_terms(strata_vars, autobin_info, data)
  rhs <- paste(sprintf("`%s`", adj$terms), collapse = " + ")
  adjusted_formula <- stats::update(null_formula,
                                    stats::as.formula(paste(". ~ . +", rhs)))
  list(formula = adjusted_formula, data = adj$data)
}
