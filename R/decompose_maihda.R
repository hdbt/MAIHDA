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
  # Quote via maihda_quote_name() (deparse(as.name())) rather than a manual
  # sprintf("`%s`"): the helper escapes the rare legal column name that itself
  # contains a backtick, which naive backtick-wrapping would turn into a broken
  # formula.
  rhs <- paste(vapply(adj$terms, maihda_quote_name, character(1)), collapse = " + ")
  adjusted_formula <- stats::update(null_formula,
                                    stats::as.formula(paste(". ~ . +", rhs)))
  list(formula = adjusted_formula, data = adj$data)
}

#' Build the crossed-dimensions-model formula and data for a MAIHDA decomposition
#'
#' The crossed-dimensions alternative to the two-model (fixed-effects PCV)
#' decomposition (the function name keeps the historical "cross_classified"
#' spelling). Given a null model's formula (in \code{(1 | stratum)} form, carrying
#' only the covariates) and the stratum metadata, returns the single crossed formula --
#' the covariates plus an \emph{additive random intercept for each stratum dimension}
#' plus the intersection (\code{stratum}) random intercept -- together with the data
#' carrying any reconstructed binned factors. In the fitted model each dimension's RE
#' variance is that dimension's additive main-effect variance and the \code{stratum} RE
#' variance is the interaction beyond additive; see \code{\link{maihda}}.
#'
#' Returns \code{NULL} when fewer than two dimensions are available (there is no
#' intersection to decompose). The dimension grouping factor reuses the dimension's own
#' column for a categorical dimension and the reconstructed \code{.maihda_dim_*} tertile
#' factor for an auto-binned numeric dimension (via \code{\link{maihda_adjusted_terms}}),
#' so the additive REs are crossed on exactly the levels that define the strata.
#'
#' @param null_formula The null model formula using \code{(1 | stratum)} (covariates
#'   only -- any dimension main effects written as fixed terms should be removed first,
#'   because they enter as random effects here).
#' @param strata_vars Character vector of stratum-defining variables.
#' @param autobin_info Auto-binning recipe (\code{strata_autobin_info}).
#' @param data The null model's data (\code{original_data}) with the \code{stratum}
#'   column and the dimension columns.
#' @param interaction_group Name of the intersection grouping factor (the column whose
#'   random intercept captures the interaction). Default \code{"stratum"}.
#' @return A list with \code{formula}, \code{data}, \code{dim_groups} (a named character
#'   vector mapping each \code{strata_var} to its random-effect grouping-factor name) and
#'   \code{interaction_group} (\code{"stratum"}); or \code{NULL} if fewer than two
#'   dimensions are available.
#' @keywords internal
#' @importFrom stats update as.formula
maihda_cross_classified_formula <- function(null_formula, strata_vars, autobin_info,
                                            data, interaction_group = "stratum") {
  if (is.null(strata_vars) || length(strata_vars) < 2) {
    return(NULL)
  }
  adj <- maihda_adjusted_terms(strata_vars, autobin_info, data)
  # One additive random intercept per dimension (on the dimension's own grouping
  # factor) plus the intersection random intercept. nobars() on null_formula keeps the
  # covariates; we re-add the stratum RE so the builder is idempotent w.r.t. it.
  fixed_formula <- reformulas::nobars(null_formula)
  re_terms <- c(
    sprintf("(1 | %s)", vapply(adj$terms, maihda_quote_name, character(1))),
    sprintf("(1 | %s)", maihda_quote_name(interaction_group))
  )
  cc_formula <- stats::update(
    fixed_formula,
    stats::as.formula(paste(". ~ . +", paste(re_terms, collapse = " + ")))
  )
  dim_groups <- stats::setNames(adj$terms, strata_vars)
  list(formula = cc_formula, data = adj$data, dim_groups = dim_groups,
       interaction_group = interaction_group)
}
