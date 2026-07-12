# Reserved column-name prefix for the internal binned-dimension factors that the
# decomposition machinery writes when a numeric stratum dimension was auto-binned:
# maihda_adjusted_terms() adds them at fit time (the adjusted / crossed formulae
# reference them) and maihda_add_binned_dim_columns() rebuilds them at predict time.
# Centralised so the two writers cannot drift on the naming, and so make_strata()
# can reserve it up front. Mirrors the reserved weight-column constants in
# design_weights.R.
.maihda_dim_prefix <- ".maihda_dim_"

# Internal column name carrying the reconstructed tertile factor for auto-binned
# dimension `v`.
maihda_dim_col <- function(v) paste0(.maihda_dim_prefix, v)

# Guard the reserved '.maihda_dim_<v>' column an auto-binned dimension will write
# against silently overwriting a user column of the same name. Called from
# make_strata() -- the SOLE producer of the auto-bin recipe -- on the raw user data,
# so it fires before any internal augmentation and never false-positives on the
# '.maihda_dim_*' columns the package itself adds to a fitted model's `original_data`
# downstream (those are re-augmented idempotently, not user data). Errors with a
# rename hint, mirroring maihda_guard_reserved_weight_col().
maihda_guard_reserved_dim_col <- function(v, data) {
  col <- maihda_dim_col(v)
  if (col %in% names(data)) {
    stop("Auto-binning the numeric stratum variable '", v, "' would create the ",
         "reserved internal column '", col, "', but 'data' already contains a column ",
         "of that name; auto-binning would overwrite it. Rename or remove '", col,
         "', or pass autobin = FALSE and bin '", v, "' yourself before creating strata.",
         call. = FALSE)
  }
  invisible(TRUE)
}

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
      new_col <- maihda_dim_col(v)
      data[[new_col]] <- cut(data[[v]], breaks = info$breaks,
                             include.lowest = TRUE, labels = info$labels)
      terms <- c(terms, new_col)
    } else {
      terms <- c(terms, v)
    }
  }
  list(terms = terms, data = data)
}

#' Warn when a numeric stratum dimension enters the adjusted model as a linear term
#'
#' In MAIHDA the stratum-defining dimensions are categorical, but the adjusted /
#' decomposition model adds an un-binned numeric dimension by its raw column name --
#' i.e. as a single linear slope, not categorical main effects (see the \code{else}
#' branch of \code{\link{maihda_adjusted_terms}}). That silently changes the
#' PCV/decomposition interpretation. This helper flags the offending dimensions and
#' warns; it does \emph{not} alter the model (the term is still entered as-is). The
#' flagged case is a stratum variable that is \emph{numeric} in the data, \emph{not} a
#' key in \code{autobin_info}, and has \strong{three or more} distinct values -- either
#' a category code with few unique values, or a many-valued numeric fitted with
#' \code{autobin = FALSE}. Categorical (factor/character) dimensions and auto-binned
#' numerics -- reconstructed as the \code{.maihda_dim_*} tertile factors -- are excluded,
#' and so is a binary (two-level) numeric dimension: with only two levels a linear term
#' and the two-level factor span the same design space, so the fit and the PCV coincide
#' and there is nothing to warn about.
#'
#' Called once per decomposition at each public entry point (not inside
#' \code{maihda_adjusted_terms()}, which runs several times per decomposition), so the
#' warning fires exactly once.
#'
#' @param strata_vars Character vector of stratum-defining variable names.
#' @param autobin_info Auto-binning recipe (\code{strata_autobin_info}); its names are
#'   the auto-binned numeric dimensions, which are excluded from the warning.
#' @param data Data frame carrying the original stratum-defining columns.
#' @param fn Name of the calling function, used in the warning prefix.
#' @return Invisibly, the character vector of flagged dimension names (empty if none).
#' @keywords internal
maihda_warn_linear_strata_dims <- function(strata_vars, autobin_info, data,
                                           fn = "maihda") {
  # names(NULL) is NULL, so an absent/empty autobin_info flags every numeric dimension.
  # Require >= 3 distinct values: a binary (two-level) numeric enters identically as a
  # linear term or a two-level factor, so there is no linear-vs-categorical discrepancy.
  binned <- names(autobin_info)
  offenders <- strata_vars[vapply(strata_vars, function(v) {
    !(v %in% binned) && v %in% names(data) && is.numeric(data[[v]]) &&
      length(unique(stats::na.omit(data[[v]]))) >= 3L
  }, logical(1))]
  if (length(offenders) == 0L) {
    return(invisible(character(0)))
  }
  warning(fn, "(): the numeric stratum dimension(s) ",
          paste(offenders, collapse = ", "),
          " enter the adjusted/decomposition model as linear fixed effect(s), not ",
          "categorical main effects. In MAIHDA the stratum dimensions are categorical, ",
          "so a linear term changes the PCV/decomposition interpretation. If these hold ",
          "category codes, convert them with factor() before fitting (e.g. data$",
          offenders[1], " <- factor(data$", offenders[1], ")); if a variable is ",
          "genuinely continuous, use it as a covariate rather than a stratum dimension. ",
          "(Numeric dimensions with >10 values are auto-binned into tertiles when ",
          "autobin = TRUE.)", call. = FALSE)
  invisible(offenders)
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
#' @param context Optional character vector of higher-level grouping variables that the
#'   caller re-appends as contextual random intercepts (via \code{context =} on the fit).
#'   Named here only so the extra-random-effect guard treats them as legitimate rather
#'   than flagging them; the builder itself does not add them.
#' @return A list with \code{formula}, \code{data}, \code{dim_groups} (a named character
#'   vector mapping each \code{strata_var} to its random-effect grouping-factor name) and
#'   \code{interaction_group} (\code{"stratum"}); or \code{NULL} if fewer than two
#'   dimensions are available.
#' @keywords internal
#' @importFrom stats update as.formula
maihda_cross_classified_formula <- function(null_formula, strata_vars, autobin_info,
                                            data, interaction_group = "stratum",
                                            context = NULL) {
  if (is.null(strata_vars) || length(strata_vars) < 2) {
    return(NULL)
  }
  adj <- maihda_adjusted_terms(strata_vars, autobin_info, data)

  # Stripping the bars keeps the covariates and lets the builder re-add exactly the
  # dimension + intersection intercepts (the caller re-appends any `context` random
  # intercept). A random effect written directly in the formula that is NOT the
  # intersection group, a stratum dimension, or a supplied `context` variable would
  # be dropped here without a trace, silently misallocating its variance to the
  # strata or the residual. Reject it with a directed error rather than drop it.
  allowed_re_vars <- unique(c(interaction_group, strata_vars, adj$terms, context))
  extra_re <- character(0)
  for (b in reformulas::findbars(null_formula)) {
    if (!all(all.vars(b[[3]]) %in% allowed_re_vars)) {
      extra_re <- c(extra_re, deparse(b))
    }
  }
  if (length(extra_re) > 0) {
    stop("decomposition = \"crossed-dimensions\" cannot carry the extra random ",
         "effect(s) ", paste(sprintf("(%s)", extra_re), collapse = ", "),
         ": the crossed model replaces the random part with one intercept per ",
         "stratum dimension plus the intersection (stratum) intercept, so these ",
         "terms would be silently dropped. Supply a higher-level grouping through ",
         "context = instead (it composes with the crossed model), or use ",
         "decomposition = \"two-model\".", call. = FALSE)
  }

  # One additive random intercept per dimension (on the dimension's own grouping
  # factor) plus the intersection random intercept. Stripping the bars from
  # null_formula keeps the covariates; we re-add the stratum RE so the builder is
  # idempotent w.r.t. it.
  fixed_formula <- maihda_nobars(null_formula)
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
