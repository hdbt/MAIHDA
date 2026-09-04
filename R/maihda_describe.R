# Pre-model MAIHDA sample description.
#
# The natural counterpart to the post-model maihda_table(): where maihda_table()
# reports the fitted null/adjusted results and ranked strata, maihda_describe()
# reports the sample the model is about to see -- the "Table 1" of an applied
# MAIHDA write-up (total N, the complete-case analytic sample, observed vs.
# expected intersectional strata, dimension distributions, per-stratum outcome
# summaries, missingness accounting, and context units). It introduces no new
# machinery: the strata come from make_strata() via the same shorthand
# resolution as fit_maihda() (maihda_resolve_strata_formula), the analytic
# sample from the same row mask the engines apply (maihda_analytic_keep_mask),
# the family from the same auto-detection as the fitters, and the outcome
# values from the same family-aware extraction as the plots
# (maihda_observed_outcome_for_plot) -- so the description is guaranteed to
# match a subsequent make_strata() / fit_maihda() / maihda() exactly.

# Documented thresholds for the data-quality flags in $warnings.
.maihda_describe_high_missing <- 10     # % outcome missingness that trips the overall flag
.maihda_describe_conc_floor <- 10       # % floor for the per-stratum concentration flag
.maihda_describe_context_levels <- 10   # fewer context levels weakly identify the variance
.maihda_describe_small_unit <- 5        # a context unit with fewer obs is "small"
.maihda_describe_id_levels <- 50        # a dimension with >= this many levels looks like an ID
.maihda_describe_id_ratio <- 0.5        # ... or with more unique levels than half the rows
.maihda_describe_max_enum <- 100000     # cap on Cartesian enumeration of empty strata

#' Describe the MAIHDA sample before fitting
#'
#' @description
#' Builds the standard \dQuote{Table 1} descriptives of a MAIHDA analytic sample
#' \emph{before} any model is fitted: the total and complete-case sample sizes,
#' the observed vs. expected intersectional strata (and which are empty or
#' small), the distribution of each stratum-defining dimension, family-aware
#' per-stratum outcome summaries, missing-data accounting, and -- for a
#' contextual cross-classified design -- the context units. It is the pre-model
#' counterpart of \code{\link{maihda_table}}.
#'
#' The description comes from the \emph{same machinery} as the model: the strata
#' are built by \code{\link{make_strata}} via the same formula-shorthand
#' resolution as \code{\link{fit_maihda}}, the analytic sample reproduces the
#' engines' own row selection (missing outcome, covariates, stratum dimensions,
#' context, and -- for a weighted design -- missing/non-positive sampling
#' weights), and the outcome family uses the same auto-detection as the fitters.
#' Stratum IDs, labels, and counts therefore match a subsequent
#' \code{make_strata()} / \code{fit_maihda()} / \code{\link{maihda}} on the same
#' formula and data exactly.
#'
#' @details
#' \strong{Family-aware outcome summaries.} The per-stratum (and overall)
#' outcome summary depends on the resolved family, and the applied summary type
#' is recorded in \code{$outcome_summary} -- a 0/1 outcome is never described
#' with a Gaussian mean/SD as if continuous:
#' \itemize{
#'   \item gaussian (and other continuous families): mean, SD, median, min, max;
#'   \item binomial (including an aggregated outcome written either as
#'     \code{cbind(success, failure)} or, for \code{brms}, as
#'     \code{success | trials(n)}): event count, trials, and the observed
#'     proportion. Trials are the binomial denominator, so a row with a missing
#'     or non-positive trial count has no observed outcome;
#'   \item poisson / negative binomial: the same numeric summaries, read as the
#'     observed count mean (rate);
#'   \item cumulative (ordinal): the mean and (lower) median category score
#'     (categories scored 1..K in level order, the same convention as the
#'     package's plots), with the full category distribution in
#'     \code{$outcome_levels}.
#' }
#'
#' \strong{Intersectional sparsity.} \code{$overview} contrasts
#' \code{n_strata_observed} with \code{n_strata_expected} (the full Cartesian
#' product of the observed dimension levels) and \code{n_empty_strata}; with
#' \code{include_empty_strata = TRUE} the empty combinations are enumerated as
#' rows of \code{$strata} (with \code{stratum = NA}, so real stratum IDs stay
#' identical to a subsequent fit). This sparsity is exactly what motivates the
#' Bayesian engine for small-cell designs (\code{engine = "brms"}).
#'
#' \strong{Missing data.} \code{$overview} separates the rows lost to a missing
#' stratum dimension (\code{make_strata()} routes them to the NA stratum) from
#' rows with a missing outcome, and reports the resulting complete-case
#' \code{n_analytic} -- computed with the same row mask the engines apply, so it
#' equals \code{nrow(fit_maihda(...)$data)}.
#'
#' \strong{Weights.} With \code{sampling_weights}, weighted counts and weighted
#' outcome means/proportions are reported alongside the unweighted ones, and
#' rows with a missing or non-positive weight are excluded from the analytic
#' sample (as the weighted engines do). As elsewhere in the package, these are
#' population-weighted \emph{point} summaries -- not a complex survey design; no
#' design-based variances are computed.
#'
#' \strong{Units.} \code{maihda_describe()} treats rows as the unit of
#' description. For longitudinal (long-format) data a row is a measurement
#' occasion, not a person; describing a fitted longitudinal model
#' (\code{fit_maihda(id =, time =)}) additionally reports \code{n_individuals}
#' and the median occasions per individual in \code{$overview}.
#'
#' \strong{Fitted-model input.} A \code{maihda_model} (or the
#' \code{maihda_analysis} bundle from \code{\link{maihda}}) can be described
#' post hoc; the formula, data, family, context, and sampling weights are taken
#' from the fit, so the description covers the exact analytic sample the model
#' used. For the \code{wemix}/\code{ordinal} engines the stored data already
#' \emph{is} the pre-filtered analytic sample, so total and analytic counts
#' coincide there.
#'
#' @param x A model formula (with the intersectional shorthand
#'   \code{outcome ~ covars + (1 | var1:var2)} or \code{... + (1 | stratum)}
#'   after \code{\link{make_strata}}), or a fitted \code{maihda_model} /
#'   \code{maihda_analysis} to describe post hoc.
#' @param data A data frame with the model variables; required when \code{x} is
#'   a formula, and must be omitted for a fitted-model input.
#' @param context Optional character vector naming higher-level context
#'   column(s), exactly as in \code{\link{fit_maihda}}; adds the per-unit
#'   \code{$context} table and the contextual identification checks. Must be
#'   omitted for a fitted-model input (taken from the fit).
#' @param family \code{NULL} (default) auto-detects the family the same way the
#'   fitters do -- a binary outcome becomes \code{"binomial"}, an ordered factor
#'   with 3+ levels becomes the cumulative (ordinal) model, anything else is
#'   described as \code{"gaussian"}. Otherwise any family specification
#'   \code{\link{fit_maihda}} accepts (e.g. \code{"poisson"},
#'   \code{"negbinomial"}, a family object). Must be omitted for a fitted-model
#'   input.
#' @param sampling_weights Optional single character string naming a numeric
#'   column of \code{data} with individual sampling (design) weights, as in
#'   \code{\link{fit_maihda}}; adds weighted counts and outcome summaries. Must
#'   be omitted for a fitted-model input.
#' @param flag_stratum_n Strata with \code{n} at or below this threshold are
#'   \emph{flagged} as small (column \code{small} of \code{$strata} and a
#'   \code{$warnings} entry) -- never dropped. This deliberately differs from
#'   \code{make_strata()}'s \code{min_n}, which drops. Default 20.
#' @param include_empty_strata Logical; when \code{TRUE} (default) the
#'   zero-count combinations of the observed dimension levels are appended to
#'   \code{$strata} with \code{stratum = NA} and \code{n = 0}.
#' @param autobin Passed to \code{\link{make_strata}} when the formula shorthand
#'   builds the strata, so the description matches a fit with the same
#'   \code{autobin} setting. Ignored for a fitted-model input (the strata are
#'   already built). Default \code{TRUE}.
#' @param digits Decimal places used by the \code{print()} method. Default 3.
#' @param weights Optional precision weights, given as a bare column name of
#'   \code{data} exactly as in \code{\link{fit_maihda}}, so the same call
#'   describes and fits the same sample. Rows whose weight is missing,
#'   zero, negative or non-finite are excluded from the analytic sample, as the
#'   engines exclude them. For a \strong{binomial} model they also supply the
#'   denominator of R's second aggregated-binomial idiom (see \code{?glm}): a
#'   proportion response with the trial counts passed as \code{weights}, which
#'   is then summarised as events out of trials rather than as one trial per row.
#'   That reading is taken only when the weights are non-unit whole numbers, the
#'   same rule \code{\link{maihda_discriminatory_accuracy}} applies, so the two
#'   never report different sample sizes for one model; it covers the frequency-cell
#'   spelling (a 0/1 response whose counts ride in \code{weights}) as well as a
#'   proportion response. Non-integral weights are not counts and are left alone.
#'   There is no \code{binomial_weights} override here as there is on the AUC:
#'   those select an estimand for the concordance, not a way to count a sample.
#'   Mutually
#'   exclusive with \code{sampling_weights}, which are design weights and mean
#'   something different. Must be omitted for a fitted-model input (the fit
#'   already carries its weights). Placed last in the argument list for backward
#'   compatibility; supply it by name.
#'
#' @return An object of class \code{maihda_describe}: a list of export-ready
#'   data frames (pass to \code{write.csv()} or \code{knitr::kable()}) plus
#'   metadata:
#'   \item{overview}{one row: \code{n_total}, \code{n_missing_outcome},
#'     \code{pct_missing_outcome}, \code{n_rows_missing_dimensions} (rows
#'     outside every stratum), \code{n_analytic}, \code{n_strata_observed},
#'     \code{n_strata_expected}, \code{n_empty_strata}; plus
#'     \code{n_invalid_weights}/\code{sum_weights} for a weighted description
#'     and \code{n_individuals}/\code{median_occasions_per_individual} for a
#'     fitted longitudinal model}
#'   \item{dimensions}{counts and percentages for each level of each MAIHDA
#'     dimension (\code{NULL} when the dimensions are unknown, e.g. a hand-built
#'     \code{stratum} column)}
#'   \item{strata}{one row per stratum: \code{stratum}, \code{label}, the
#'     dimension values, \code{n}, \code{n_analytic}, \code{n_missing_outcome},
#'     \code{pct_missing_outcome}, the family-aware outcome summary columns, and
#'     the \code{small}/\code{empty} flags}
#'   \item{outcome_overall}{the same family-aware outcome summary over the whole
#'     sample}
#'   \item{outcome_levels}{the outcome's category distribution (binomial /
#'     ordinal outcomes; \code{NULL} otherwise)}
#'   \item{context}{one row per context unit: \code{context}, \code{level},
#'     \code{n}, \code{n_analytic}, \code{pct} (\code{NULL} without a context)}
#'   \item{missingness}{per variable (outcome, each dimension, each context,
#'     the weight column): \code{n_missing} and \code{pct_missing}}
#'   \item{warnings}{data-quality flags as a data frame of \code{check} /
#'     \code{message} rows (zero rows when clean): auto-binned or ID-like or
#'     linear-numeric dimensions, empty/small strata, rows lost to missing
#'     dimensions, high or concentrated outcome missingness, and weakly
#'     identified contexts}
#'   \item{observations}{a slim per-row frame (stratum and the outcome
#'     numerator/denominator on the summary scale) kept so \code{plot()} is
#'     self-contained}
#'   \item{outcome, family, family_detected, outcome_summary, event_level,
#'     strata_vars, context_vars, sampling_weights, source, engine,
#'     longitudinal, response_recoding, call, ...}{metadata used by
#'     \code{print()}}
#'
#' @seealso \code{\link{maihda_table}} for the post-model tables,
#'   \code{\link{make_strata}}, \code{\link{fit_maihda}},
#'   \code{\link[=plot.maihda_describe]{plot}} for the descriptive plots.
#'
#' @examples
#' data(maihda_health_data)
#' desc <- maihda_describe(BMI ~ Age + (1 | Gender:Race:Education),
#'                         data = maihda_health_data, flag_stratum_n = 20)
#' desc
#' desc$overview
#' head(desc$strata)
#'
#' # Binary outcome: auto-detected binomial, proportions instead of means
#' desc_bin <- maihda_describe(Obese ~ (1 | Gender:Race:Education),
#'                             data = maihda_health_data)
#' desc_bin$outcome_levels
#'
#' # Contextual cross-classified design (country as context)
#' data(maihda_country_data)
#' desc_ctx <- maihda_describe(math ~ escs + (1 | gender:ses),
#'                             data = maihda_country_data, context = "country")
#' desc_ctx$context
#'
#' # Export-ready tables
#' # write.csv(desc$strata, "table1_strata.csv", row.names = FALSE)
#'
#' @export
#' @importFrom stats median quantile sd complete.cases setNames
maihda_describe <- function(x, data = NULL, context = NULL, family = NULL,
                            sampling_weights = NULL, flag_stratum_n = 20,
                            include_empty_strata = TRUE, autobin = TRUE,
                            digits = 3, weights = NULL) {
  # `weights` is deliberately the LAST formal rather than sitting beside
  # `sampling_weights` where it reads better: inserting it mid-signature would shift
  # every positional argument after it in code written against 0.2.1.
  # Non-standard evaluation, as in fit_maihda(weights = ), so the same call works for
  # both; captured before any use because on the fitted-model branch there is no
  # `data` to evaluate it against, only a request to reject it.
  weights_quo <- rlang::enquo(weights)
  weights_supplied <- !rlang::quo_is_null(weights_quo)
  if (!is.numeric(flag_stratum_n) || length(flag_stratum_n) != 1 ||
      !is.finite(flag_stratum_n) || flag_stratum_n < 0) {
    stop("'flag_stratum_n' must be a single non-negative number.", call. = FALSE)
  }
  if (!isTRUE(include_empty_strata) && !isFALSE(include_empty_strata)) {
    stop("'include_empty_strata' must be TRUE or FALSE.", call. = FALSE)
  }
  call <- match.call()

  # --- Fitted-model input: describe the exact sample the fit used -------------
  if (inherits(x, "maihda_analysis") || inherits(x, "maihda_model")) {
    if (!is.null(data) || !is.null(context) || !is.null(family) ||
        !is.null(sampling_weights) || weights_supplied) {
      stop("'x' is a fitted MAIHDA object and already carries its data, ",
           "context, family, and weights; do not also supply 'data', ",
           "'context', 'family', 'sampling_weights', or 'weights'.", call. = FALSE)
    }
    source <- if (inherits(x, "maihda_analysis")) "maihda_analysis" else "maihda_model"
    model <- if (inherits(x, "maihda_analysis")) x$model else x
    if (!inherits(model, "maihda_model")) {
      stop("The maihda_analysis object carries no fitted model to describe.",
           call. = FALSE)
    }
    return(maihda_describe_from_model(model, source = source,
                                      flag_stratum_n = flag_stratum_n,
                                      include_empty_strata = include_empty_strata,
                                      digits = digits, call = call))
  }

  # --- Formula + data input ----------------------------------------------------
  if (!inherits(x, "formula")) {
    stop("'x' must be a model formula, a maihda_model (from fit_maihda()), or ",
         "a maihda_analysis (from maihda()).", call. = FALSE)
  }
  if (is.null(data)) {
    stop("'data' is required when 'x' is a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  context <- maihda_validate_context(context, data)
  if (!is.null(sampling_weights)) {
    sampling_weights <- maihda_validate_sampling_weights(sampling_weights, data)
  }
  # Precision weights, resolved against `data` exactly as fit_maihda() resolves its
  # own: the same column name gives the same rows and the same binomial denominator
  # in both calls. They are mutually exclusive with design weights there, so they are
  # here -- the two mean different things and only one can be `weights =`.
  weights_value <- NULL
  if (weights_supplied) {
    if (!is.null(sampling_weights)) {
      stop("Supply either 'sampling_weights' (design weights) or 'weights' ",
           "(precision weights), not both.", call. = FALSE)
    }
    weights_value <- tryCatch(
      rlang::eval_tidy(weights_quo, data = data),
      error = function(e) stop("maihda_describe(): could not evaluate 'weights' ",
                               "in the data: ", conditionMessage(e), call. = FALSE))
    if (!is.numeric(weights_value) || length(weights_value) != nrow(data)) {
      stop("'weights' must be a numeric vector with one value per row of 'data'.",
           call. = FALSE)
    }
    # Zero / negative / non-finite -> NA, the rows lme4 actually drops. Shared with
    # the fitters, so the analytic sample described here is the one that would be fit.
    weights_value <- maihda_normalize_precision_weights(
      weights_value, nrow(data), "lme4", "maihda_describe()")
  }

  # Resolve the strata exactly as fit_maihda() would: the shared helper parses
  # the shorthand, calls make_strata(), and rewrites the grouping to
  # (1 | stratum) -- identical IDs, labels, counts, and validation errors.
  strata_res <- maihda_resolve_strata_formula(x, data, autobin)
  formula <- strata_res$formula
  data <- strata_res$data
  if (!"stratum" %in% names(data)) {
    stop("The formula has no intersectional stratum term. Use the shorthand ",
         "(1 | var1:var2), or build the strata with make_strata() and include ",
         "(1 | stratum).", call. = FALSE)
  }

  # Family: auto-detect with the fitters' own detection (on the same analytic
  # sample, incl. the weighted engines' weight-based row drops), or resolve an
  # explicit specification with the fitters' own resolver.
  detect_weights <- if (!is.null(sampling_weights)) {
    maihda_sampling_weight_mask(data[[sampling_weights]])
  } else {
    # NULL when no `weights` was given; otherwise the NA-normalised precision
    # weights, whose NA rows maihda_row_mask() drops as the engines drop them.
    weights_value
  }
  family_detected <- FALSE
  if (is.null(family)) {
    if (isTRUE(tryCatch(
      maihda_response_is_binary(formula, data, weights = detect_weights),
      error = function(e) FALSE))) {
      fam_obj <- stats::binomial()
      family_detected <- TRUE
    } else if (isTRUE(tryCatch(
      maihda_response_is_ordinal(formula, data, weights = detect_weights),
      error = function(e) FALSE))) {
      fam_obj <- maihda_cumulative("logit")
      family_detected <- TRUE
    } else {
      fam_obj <- stats::gaussian()
    }
  } else {
    fam_obj <- maihda_resolve_family_spec(family)
  }

  # Contextual design: same validation + formula append as fit_maihda(), so the
  # analytic row mask below also drops rows with a missing context.
  if (!is.null(context)) {
    formula <- maihda_apply_context_formula(formula, context,
                                            strata_res$strata_vars)
  }

  # The analytic (complete-case) sample: exactly the rows an engine would fit --
  # transformed terms honoured, NA stratum (missing dims), missing outcome /
  # covariates / context dropped, and missing/non-positive sampling weights
  # excluded the way the weighted engines exclude them.
  keep <- maihda_analytic_keep_mask(formula, data, weights = detect_weights)
  if (is.null(keep)) {
    # Fallback when the model frame cannot be built (mirrors the detection
    # helpers' raw fallback): complete cases over the formula's columns.
    vars_needed <- intersect(all.vars(formula), names(data))
    keep <- stats::complete.cases(data[, vars_needed, drop = FALSE])
    if (!is.null(detect_weights)) keep <- keep & !is.na(detect_weights)
  }

  maihda_describe_build(
    formula = formula, data = data, analytic_data = data[keep, , drop = FALSE],
    model_agg = maihda_describe_formula_agg(formula, data, keep, fam_obj,
                                            weights_value),
    fam_obj = fam_obj, family_detected = family_detected,
    strata_info = strata_res$strata_info, strata_vars = strata_res$strata_vars,
    sep = strata_res$strata_sep, autobin_info = strata_res$autobin_info,
    context_vars = context, sampling_weights = sampling_weights,
    flag_stratum_n = flag_stratum_n,
    include_empty_strata = include_empty_strata,
    source = "formula", engine = NULL, longitudinal_info = NULL,
    response_recoding = NULL, digits = digits, call = call
  )
}

# Describe a fitted maihda_model post hoc: formula/data/family/context/weights
# come from the fit, per-stratum analytic counts from the stored analytic frame
# (model$data, the rows the engine actually used), and totals from the pre-fit
# data (model$original_data). For the wemix/ordinal engines original_data
# already IS the pre-filtered analytic sample, so totals == analytic there.
maihda_describe_from_model <- function(model, source, flag_stratum_n,
                                       include_empty_strata, digits, call) {
  data <- model$original_data
  if (!is.data.frame(data)) {
    stop("The fitted model carries no original data to describe.", call. = FALSE)
  }
  maihda_describe_build(
    formula = model$formula, data = data, analytic_data = model$data,
    model_agg = maihda_describe_model_agg(model, data),
    fam_obj = model$family, family_detected = FALSE,
    strata_info = model$strata_info, strata_vars = model$strata_vars,
    sep = model$strata_sep, autobin_info = model$strata_autobin_info,
    context_vars = model$context_vars,
    sampling_weights = model$sampling_weights,
    flag_stratum_n = flag_stratum_n,
    include_empty_strata = include_empty_strata,
    source = source, engine = model$engine,
    longitudinal_info = model$longitudinal_info,
    response_recoding = model$response_recoding,
    digits = digits, call = call
  )
}

# Per-row aggregated-binomial success / trial counts implied by a PRE-FIT
# `weights =` argument, aligned to all of `data`'s rows, or NULL when the response
# and weights are not an aggregated binomial. The formula-path counterpart of
# maihda_describe_model_agg(); both defer to maihda_agg_counts_from_weights(), so
# maihda_describe(f, d, weights = n) and maihda_describe(fit_maihda(f, d, weights = n))
# report the same events and trials.
#
# The rule is applied to the ANALYTIC rows (`keep`) rather than to all of `data`, for
# the same reason the fitted-model path reads it off the analytic frame: those are the
# rows an engine would see, and a single NA-weight or NA-covariate row would otherwise
# decide for the whole sample whether the weights look like trial counts at all. Rows
# outside the mask get NA trials and are counted as having no observed outcome, as on
# the fitted-model path.
maihda_describe_formula_agg <- function(formula, data, keep, fam_obj, weights_value) {
  n_total <- nrow(data)
  if (is.null(weights_value) || length(weights_value) != n_total ||
      !is.logical(keep) || length(keep) != n_total || !any(keep)) {
    return(NULL)
  }
  # Gate on the family before the detector runs, so a Gaussian proportion outcome
  # with integral precision weights is never warned about as "successes".
  fam_name <- tryCatch(maihda_normalize_family_name(fam_obj$family),
                       error = function(e) NA_character_)
  if (!isTRUE(fam_name %in% c("binomial", "quasibinomial"))) {
    return(NULL)
  }
  y <- tryCatch(eval(maihda_describe_response_expr(formula), envir = data,
                     enclos = environment(formula)),
                error = function(e) NULL)
  # A cbind() matrix response is numeric too, and carries its own denominator; it is
  # summarised row-wise upstream and must not be routed through the weights rule.
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) != n_total) {
    return(NULL)
  }
  agg <- maihda_agg_counts_from_weights(y[keep], weights_value[keep],
                                        context = "Aggregated-binomial outcome")
  if (is.null(agg)) {
    return(NULL)
  }
  successes <- rep(NA_real_, n_total)
  trials <- rep(NA_real_, n_total)
  successes[keep] <- agg$successes
  trials[keep] <- agg$trials
  list(successes = successes, trials = trials)
}

# Per-row aggregated-binomial success / trial counts of a FITTED model, aligned to
# its pre-fit `data` (model$original_data) rows, or NULL when the fit is not an
# aggregated binomial that only the fitted object can reveal.
#
# R has two spellings of an aggregated binomial, and only one of them is visible in
# the formula: cbind(successes, failures) (a matrix response, which
# maihda_describe_build() already summarises row-wise) and -- per ?glm -- a
# PROPORTION response whose trial counts are supplied as `weights =`. The second
# leaves no trace on the formula, so maihda_trials_from_formula() cannot see it and
# describe took the denominator to be 1 per ROW: a 12-row, 340-of-617 sample was
# reported as "6.42 events / 12 trials (53.5%)" -- a fractional event count over a
# row count -- while maihda_discriminatory_accuracy() read the same fit as 340 cases
# / 277 controls. Two numbers for one model.
#
# The rule is not re-implemented here: maihda_da_aggregated_counts() IS the rule, so
# describe and the AUC cannot drift apart -- including through changes to the rule.
# It currently reads any non-unit integral prior weights as trial counts whatever the
# response, so the frequency-cell spelling (a 0/1 response whose counts ride in
# `weights =`) reports the same events and trials here as cbind(successes, failures)
# does. describe always takes that "auto" reading; the AUC's binomial_weights =
# "trials" / "analytic" overrides are estimand choices for the concordance and have no
# counterpart in a sample description. Non-lme4 engines return NULL: brms trials come
# off the formula's trials() term, and the wemix/ordinal engines carry no aggregated
# response.
#
# ROW ALIGNMENT: the counts live on the analytic frame (model$data, the rows the
# engine kept), while describe evaluates the outcome on the pre-fit frame, so the two
# differ in length whenever a row was dropped for a missing covariate. model$data is
# a model frame and therefore carries the pre-fit row NAMES; map through them, and
# return NULL rather than guess if that mapping is not exact -- reporting counts
# against the wrong rows would be worse than the row-denominator this fixes.
#
# A pre-fit row that is NOT in the analytic frame gets NA trials, which the caller
# then counts as a missing outcome. That is a real cost of this spelling and is
# stated rather than hidden: the trial counts ARE the response's denominator here and
# are only recoverable from the fit, so an excluded row has no summarisable outcome
# left -- unlike cbind(successes, failures), whose denominator sits in the pre-fit
# data and is summarised for every row. The alternative, totals that silently cover
# fewer rows than the "non-missing" count beside them, would not reconcile.
maihda_describe_model_agg <- function(model, data) {
  if (!identical(model$engine, "lme4") || !is.data.frame(model$data)) {
    return(NULL)
  }
  # Gate on the family BEFORE running the detector, not just on its result. The
  # detector's rule -- a response inside [0, 1] with non-unit integral weights -- is
  # a perfectly ordinary Gaussian fit too (a proportion outcome with precision
  # weights), and while describe would discard the counts a few lines later under
  # `outcome_kind != "binomial"`, the detector would already have warned about
  # "successes" at a model that has no successes in it.
  if (!isTRUE(maihda_model_family_name(model) %in%
              c("binomial", "quasibinomial"))) {
    return(NULL)
  }
  agg <- tryCatch(
    maihda_da_aggregated_counts(model, context = "Aggregated-binomial outcome"),
    error = function(e) NULL)
  if (is.null(agg) || length(agg$trials) != nrow(model$data)) {
    return(NULL)
  }
  n_total <- nrow(data)
  idx <- if (nrow(model$data) == n_total) {
    # No row was dropped, so the analytic frame IS the pre-fit frame, in order.
    seq_len(n_total)
  } else {
    rn <- rownames(model$data)
    dn <- rownames(data)
    m <- if (is.null(rn) || is.null(dn)) NULL else match(rn, dn)
    if (is.null(m) || anyNA(m) || anyDuplicated(m)) return(NULL)
    m
  }
  successes <- rep(NA_real_, n_total)
  trials <- rep(NA_real_, n_total)
  successes[idx] <- as.numeric(agg$successes)
  trials[idx] <- as.numeric(agg$trials)
  list(successes = successes, trials = trials)
}

# The outcome expression of a resolved MAIHDA formula. brms addition terms
# (y | weights(w), y | trials(n)) wrap the response in `|` calls; the outcome is
# the leftmost leaf.
maihda_describe_response_expr <- function(formula) {
  if (length(formula) != 3L) {
    stop("The formula must have an outcome on its left-hand side.", call. = FALSE)
  }
  resp <- formula[[2]]
  while (is.call(resp) && identical(resp[[1]], as.name("|"))) {
    resp <- resp[[2]]
  }
  resp
}

# Assemble the maihda_describe object from resolved pieces. `data` is the full
# pre-model data (with the stratum column; NA stratum = outside every stratum),
# `analytic_data` the complete-case rows an engine fits (formula path: the
# keep-mask rows; model path: the stored analytic frame).
maihda_describe_build <- function(formula, data, analytic_data, fam_obj,
                                  family_detected, strata_info, strata_vars,
                                  sep, autobin_info, context_vars,
                                  sampling_weights, flag_stratum_n,
                                  include_empty_strata, source, engine,
                                  longitudinal_info, response_recoding,
                                  digits, call, model_agg = NULL) {
  if (is.null(sep)) sep <- " \u00d7 "
  n_total <- nrow(data)

  # --- Outcome values and missingness ----------------------------------------
  resp_expr <- maihda_describe_response_expr(formula)
  outcome_name <- paste(deparse(resp_expr), collapse = " ")
  outcome_vals <- tryCatch(
    eval(resp_expr, envir = data, enclos = environment(formula)),
    error = function(e) stop("maihda_describe(): could not evaluate the outcome '",
                             outcome_name, "' in the data: ",
                             conditionMessage(e), call. = FALSE))
  # An aggregated-binomial brms response keeps its denominator in a `trials()`
  # addition term (`y | trials(n)`), which the leftmost-leaf extraction above
  # strips off with the rest of the `|` expression. Recover the trial counts from
  # the same formula so the denominator is the trials, not 1 -- otherwise the
  # "proportion" is a mean success COUNT (26 events / 4 trials = 6.5 for a
  # four-row 26-of-62 sample) and $observations feeds the outcome histogram raw
  # counts. NULL for every other response, which then takes the branches below
  # unchanged.
  outcome_trials <- maihda_trials_from_formula(formula, data, n = n_total)
  # The other aggregated-binomial spelling, which the formula cannot show: a
  # proportion response with the trial counts passed as `weights =` (?glm). Only a
  # fitted model can reveal it (maihda_describe_model_agg(), which defers to the same
  # detector the AUC uses), and only when the formula carries no trials() term and the
  # response is not already a cbind() matrix -- both of those are exact and take
  # precedence. `outcome_for_extract` then holds SUCCESS COUNTS, which is what the
  # extractor's trials branch expects; `outcome_vals` stays the raw response so the
  # missingness and category-level reporting below are unchanged.
  outcome_for_extract <- outcome_vals
  if (is.null(outcome_trials) && !is.null(model_agg) &&
      !is.matrix(outcome_vals) && !is.data.frame(outcome_vals) &&
      length(model_agg$trials) == n_total) {
    outcome_trials <- model_agg$trials
    outcome_for_extract <- model_agg$successes
  }
  if (is.matrix(outcome_vals) || is.data.frame(outcome_vals)) {
    if (nrow(outcome_vals) != n_total) {
      stop("The outcome '", outcome_name, "' does not have one row per data row.",
           call. = FALSE)
    }
    totals <- rowSums(as.matrix(outcome_vals))
    # The same per-row validity rule the plots apply to a matrix binomial
    # outcome: no finite trials (or zero trials) means no observed outcome.
    outcome_missing <- !is.finite(totals) | totals <= 0
  } else {
    if (length(outcome_vals) != n_total) {
      stop("The outcome '", outcome_name, "' does not have one value per data row.",
           call. = FALSE)
    }
    outcome_missing <- is.na(outcome_vals)
  }

  # --- Family kind and the family-aware numerator/denominator -----------------
  fam_name <- maihda_normalize_family_name(fam_obj$family)
  outcome_kind <- if (fam_name %in% c("binomial", "quasibinomial", "bernoulli")) {
    "binomial"
  } else if (identical(fam_name, "cumulative")) {
    "ordinal"
  } else if (fam_name %in% c("poisson", "negbinomial")) {
    "count"
  } else {
    "continuous"
  }
  # brms spells Bernoulli "bernoulli"; hand the extractor the canonical binomial
  # family so factor/character/matrix outcomes take the binomial branches.
  fam_for_extract <- if (outcome_kind == "binomial") {
    list(family = "binomial", link = fam_obj$link)
  } else {
    fam_obj
  }
  # The trial counts are the binomial denominator and nothing else; a trials()
  # term under some other family is not a denominator, so it is dropped rather
  # than allowed to redefine the outcome (or its missingness) there.
  if (outcome_kind != "binomial") {
    outcome_trials <- NULL
    outcome_for_extract <- outcome_vals
  }
  # A row with no finite, positive trial count has no observed outcome -- the same
  # rule the matrix branch applies to a zero/NA row total.
  if (!is.null(outcome_trials)) {
    outcome_missing <- outcome_missing |
      !is.finite(outcome_trials) | outcome_trials <= 0
  }
  od <- tryCatch(
    maihda_observed_outcome_for_plot(outcome_for_extract, fam_for_extract,
                                     trials = outcome_trials),
    error = function(e) stop("maihda_describe(): cannot summarise the outcome '",
                             outcome_name, "' under family '", fam_name, "': ",
                             conditionMessage(e), call. = FALSE))
  num <- od$numerator
  den <- od$denominator
  complete <- maihda_observed_complete(num, den)
  # On the trials path the extractor can also reject a row the checks above
  # accepted -- successes outside [0, trials], which it warns about and drops.
  # Fold that back into the missingness accounting, or the report would not
  # reconcile: "23 events / 52 trials" over "4 non-missing" rows and 0 missing
  # outcomes, when only 3 rows contributed. Restricted to the trials path so no
  # other family's missingness rule moves (the matrix branch already derives
  # outcome_missing from the same row totals the extractor uses).
  if (!is.null(outcome_trials)) {
    outcome_missing <- outcome_missing | !complete
  }

  # Category levels (binomial / ordinal vector outcomes) for $outcome_levels and
  # the event-level report; the second binary level is the modeled event, the
  # same convention as the fitters' 0/1 recoding.
  lv_out <- NULL
  if (outcome_kind == "ordinal" && is.factor(outcome_vals)) {
    lv_out <- levels(droplevels(outcome_vals))
  } else if (outcome_kind == "binomial" && is.null(dim(outcome_vals)) &&
             maihda_is_binary_vector(outcome_vals)) {
    lv_out <- as.character(maihda_binary_levels(outcome_vals))
  }
  event_level <- if (outcome_kind == "binomial" && length(lv_out) == 2) {
    lv_out[2]
  } else {
    NULL
  }
  outcome_levels <- NULL
  if (!is.null(lv_out)) {
    ov_chr <- as.character(outcome_vals)
    cnt <- table(factor(ov_chr[!outcome_missing], levels = lv_out))
    outcome_levels <- data.frame(
      level = lv_out,
      n = as.integer(cnt),
      pct = if (sum(cnt) > 0) 100 * as.integer(cnt) / sum(cnt) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  outcome_summary_label <- switch(outcome_kind,
    binomial = "events/proportion (binomial)",
    ordinal = sprintf("category scores 1..%d (cumulative)",
                      if (!is.null(lv_out)) length(lv_out) else NA_integer_),
    count = sprintf("count mean/rate (%s)", fam_name),
    sprintf("mean/SD (%s)", fam_name))

  # --- Sampling weights --------------------------------------------------------
  w <- NULL
  if (!is.null(sampling_weights) && sampling_weights %in% names(data)) {
    # Missing/non-positive weights -> NA, the weighted engines' exclusion rule.
    w <- maihda_sampling_weight_mask(data[[sampling_weights]])
  }

  # --- Strata table ------------------------------------------------------------
  grp <- as.character(data$stratum)
  valid <- !is.na(data$stratum)
  have_dims <- !is.null(strata_vars) && length(strata_vars) > 0 &&
    all(strata_vars %in% names(data)) && !is.null(strata_info)

  if (have_dims) {
    ids <- as.character(strata_info$stratum)
    labels <- as.character(strata_info$label)
    dim_cols <- lapply(strata_vars, function(v) as.character(strata_info[[v]]))
    names(dim_cols) <- strata_vars
  } else {
    # Hand-built stratum column with no recorded dimensions: describe the strata
    # by their IDs; the dimension tables and the expected-strata Cartesian are
    # not defined.
    ids <- as.character(sort(unique(data$stratum[valid])))
    labels <- ids
    dim_cols <- list()
  }
  n_obs_strata <- length(ids)

  count_by <- function(g) {
    tab <- table(g)
    out <- as.integer(tab[match(ids, names(tab))])
    out[is.na(out)] <- 0L
    out
  }
  st_n <- count_by(grp[valid])
  st_analytic <- count_by(as.character(analytic_data$stratum))
  miss_by <- tapply(outcome_missing[valid], grp[valid], sum)
  st_missing <- as.integer(miss_by[match(ids, names(miss_by))])
  st_missing[is.na(st_missing)] <- 0L

  # Family-aware per-stratum outcome summaries over the complete-outcome rows.
  rows_by <- split(which(valid), factor(grp[valid], levels = ids))
  summ <- maihda_describe_outcome_summaries(rows_by, num, den, complete,
                                            outcome_kind, lv_out, w)

  strata <- data.frame(
    stratum = if (have_dims) strata_info$stratum else ids,
    label = labels, stringsAsFactors = FALSE)
  for (v in names(dim_cols)) strata[[v]] <- dim_cols[[v]]
  strata$n <- st_n
  strata$n_analytic <- st_analytic
  strata$n_missing_outcome <- st_missing
  strata$pct_missing_outcome <- ifelse(st_n > 0, 100 * st_missing / st_n, NA_real_)
  if (!is.null(w)) {
    wn <- vapply(rows_by, function(r) sum(w[r], na.rm = TRUE), numeric(1))
    strata$n_weighted <- as.numeric(wn)
  }
  for (cl in names(summ)) strata[[cl]] <- summ[[cl]]
  strata$small <- strata$n > 0 & strata$n <= flag_stratum_n
  strata$empty <- strata$n == 0L

  # --- Overall outcome summary --------------------------------------------------
  all_rows <- seq_len(n_total)
  overall <- maihda_describe_outcome_summaries(list(all = all_rows), num, den,
                                               complete, outcome_kind, lv_out, w)
  outcome_overall <- data.frame(n_nonmissing = sum(!outcome_missing))
  for (cl in names(overall)) outcome_overall[[cl]] <- overall[[cl]]

  # --- Dimensions, expected strata, empty enumeration ---------------------------
  dimensions <- NULL
  n_expected <- NA_real_
  n_empty <- NA_real_
  binned_dims <- NULL
  if (have_dims) {
    # The strata-building columns: numeric dimensions that make_strata()
    # auto-binned are rebuilt from the stored breaks so the levels match the
    # strata exactly.
    binned_dims <- maihda_apply_autobin_info(
      data[, strata_vars, drop = FALSE], autobin_info)
    levels_by_var <- lapply(strata_vars, function(v)
      as.character(maihda_dim_levels(binned_dims[[v]])))
    names(levels_by_var) <- strata_vars

    dim_rows <- lapply(strata_vars, function(v) {
      col <- binned_dims[[v]]
      lv <- levels_by_var[[v]]
      cnt <- table(factor(as.character(col), levels = lv))
      out <- data.frame(dimension = v, level = lv, n = as.integer(cnt),
                        pct = if (sum(cnt) > 0) 100 * as.integer(cnt) / sum(cnt)
                              else NA_real_,
                        stringsAsFactors = FALSE)
      if (!is.null(w)) {
        wsum <- tapply(w[!is.na(col)], factor(as.character(col[!is.na(col)]),
                                              levels = lv), sum, na.rm = TRUE)
        wsum <- as.numeric(wsum)
        wsum[is.na(wsum)] <- 0
        out$n_weighted <- wsum
        tot <- sum(wsum)
        out$pct_weighted <- if (tot > 0) 100 * wsum / tot else NA_real_
      }
      out
    })
    dimensions <- do.call(rbind, dim_rows)
    rownames(dimensions) <- NULL

    n_expected <- prod(vapply(levels_by_var, length, numeric(1)))
    n_empty <- n_expected - n_obs_strata
    if (include_empty_strata && is.finite(n_expected) && n_empty > 0 &&
        n_expected <= .maihda_describe_max_enum) {
      grid <- expand.grid(levels_by_var, stringsAsFactors = FALSE,
                          KEEP.OUT.ATTRS = FALSE)
      # Anti-join against the observed strata with the same matcher the
      # prediction machinery uses, so "observed" means exactly make_strata()'s
      # combinations.
      matched <- maihda_match_strata_rows(grid, strata_info, strata_vars)
      empty_grid <- grid[is.na(matched), , drop = FALSE]
      if (nrow(empty_grid) > 0) {
        empty_rows <- strata[0, , drop = FALSE][seq_len(nrow(empty_grid)), , drop = FALSE]
        empty_rows$stratum <- NA
        empty_rows$label <- maihda_paste_label_rows(empty_grid, sep)
        for (v in strata_vars) empty_rows[[v]] <- as.character(empty_grid[[v]])
        empty_rows$n <- 0L
        empty_rows$n_analytic <- 0L
        empty_rows$n_missing_outcome <- 0L
        empty_rows$pct_missing_outcome <- NA_real_
        if (!is.null(w)) empty_rows$n_weighted <- 0
        empty_rows$small <- FALSE
        empty_rows$empty <- TRUE
        strata <- rbind(strata, empty_rows)
        rownames(strata) <- NULL
      }
    }
  }

  # --- Context units -------------------------------------------------------------
  context_df <- NULL
  if (!is.null(context_vars) && length(context_vars) > 0) {
    ctx_rows <- lapply(context_vars, function(cv) {
      col <- data[[cv]]
      lv <- as.character(maihda_dim_levels(col))
      cnt <- table(factor(as.character(col), levels = lv))
      an_col <- analytic_data[[cv]]
      an_cnt <- table(factor(as.character(an_col), levels = lv))
      data.frame(context = cv, level = lv, n = as.integer(cnt),
                 n_analytic = as.integer(an_cnt),
                 pct = if (sum(cnt) > 0) 100 * as.integer(cnt) / sum(cnt)
                       else NA_real_,
                 stringsAsFactors = FALSE)
    })
    context_df <- do.call(rbind, ctx_rows)
    rownames(context_df) <- NULL
  }

  # --- Missingness by variable -----------------------------------------------------
  miss_rows <- list(data.frame(variable = outcome_name, role = "outcome",
                               n_missing = sum(outcome_missing),
                               stringsAsFactors = FALSE))
  if (have_dims) {
    miss_rows <- c(miss_rows, lapply(strata_vars, function(v) {
      data.frame(variable = v, role = "dimension",
                 n_missing = sum(is.na(data[[v]])), stringsAsFactors = FALSE)
    }))
  }
  if (!is.null(context_vars)) {
    miss_rows <- c(miss_rows, lapply(context_vars, function(cv) {
      data.frame(variable = cv, role = "context",
                 n_missing = sum(is.na(data[[cv]])), stringsAsFactors = FALSE)
    }))
  }
  if (!is.null(w)) {
    miss_rows <- c(miss_rows, list(data.frame(
      variable = sampling_weights, role = "sampling_weights",
      n_missing = sum(is.na(w)), stringsAsFactors = FALSE)))
  }
  missingness <- do.call(rbind, miss_rows)
  missingness$pct_missing <- if (n_total > 0) {
    100 * missingness$n_missing / n_total
  } else {
    NA_real_
  }
  rownames(missingness) <- NULL

  # --- Overview ------------------------------------------------------------------
  n_analytic <- nrow(analytic_data)
  overview <- data.frame(
    n_total = n_total,
    n_missing_outcome = sum(outcome_missing),
    pct_missing_outcome = if (n_total > 0) 100 * sum(outcome_missing) / n_total
                          else NA_real_,
    n_rows_missing_dimensions = sum(!valid),
    n_analytic = n_analytic,
    n_strata_observed = n_obs_strata,
    n_strata_expected = n_expected,
    n_empty_strata = n_empty
  )
  if (!is.null(w)) {
    overview$n_invalid_weights <- sum(is.na(w))
    overview$sum_weights <- sum(w, na.rm = TRUE)
  }
  longitudinal_meta <- NULL
  if (!is.null(longitudinal_info) && !is.null(longitudinal_info$id) &&
      longitudinal_info$id %in% names(analytic_data)) {
    occ <- table(as.character(analytic_data[[longitudinal_info$id]]))
    longitudinal_meta <- list(id = longitudinal_info$id,
                              time = longitudinal_info$time,
                              n_individuals = length(occ),
                              median_occasions = stats::median(as.integer(occ)))
    overview$n_individuals <- longitudinal_meta$n_individuals
    overview$median_occasions_per_individual <- longitudinal_meta$median_occasions
  }

  # --- Data-quality flags -----------------------------------------------------------
  warnings_df <- maihda_describe_warnings(
    data = data, strata = strata, strata_vars = strata_vars,
    have_dims = have_dims, autobin_info = autobin_info,
    outcome_name = outcome_name, outcome_missing = outcome_missing,
    n_total = n_total, n_expected = n_expected, n_empty = n_empty,
    flag_stratum_n = flag_stratum_n, context_df = context_df,
    n_rows_missing_dimensions = sum(!valid))

  # Slim per-row frame (complete-outcome rows) so plot() needs no re-derivation.
  observations <- data.frame(stratum = grp[complete],
                             numerator = num[complete],
                             denominator = den[complete],
                             stringsAsFactors = FALSE)

  structure(
    list(
      overview = overview,
      dimensions = dimensions,
      strata = strata,
      outcome_overall = outcome_overall,
      outcome_levels = outcome_levels,
      context = context_df,
      missingness = missingness,
      warnings = warnings_df,
      observations = observations,
      outcome = outcome_name,
      family = fam_name,
      family_link = fam_obj$link,
      family_detected = family_detected,
      outcome_kind = outcome_kind,
      outcome_summary = outcome_summary_label,
      event_level = event_level,
      strata_vars = if (have_dims) strata_vars else NULL,
      sep = sep,
      context_vars = context_vars,
      sampling_weights = sampling_weights,
      flag_stratum_n = flag_stratum_n,
      include_empty_strata = include_empty_strata,
      source = source,
      engine = engine,
      longitudinal = longitudinal_meta,
      response_recoding = response_recoding,
      digits = digits,
      call = call
    ),
    class = "maihda_describe"
  )
}

# Family-aware outcome summaries for a list of row-index groups. Returns a named
# list of equal-length column vectors (one value per group, in the groups'
# order): numeric summaries for continuous/count outcomes, event counts and
# proportions for binomial ones, and category scores for ordinal ones -- plus
# the weighted mean/proportion when weights are supplied.
maihda_describe_outcome_summaries <- function(rows_by, num, den, complete,
                                              outcome_kind, lv_out, w = NULL) {
  cc_by <- lapply(rows_by, function(r) r[complete[r]])
  out <- list()
  if (outcome_kind == "binomial") {
    ev <- vapply(cc_by, function(r) sum(num[r]), numeric(1))
    tr <- vapply(cc_by, function(r) sum(den[r]), numeric(1))
    out$outcome_events <- ev
    out$outcome_trials <- tr
    out$outcome_proportion <- ifelse(tr > 0, ev / tr, NA_real_)
    if (!is.null(w)) {
      out$outcome_proportion_weighted <- vapply(cc_by, function(r) {
        r <- r[!is.na(w[r])]
        d <- sum(w[r] * den[r])
        if (d > 0) sum(w[r] * num[r]) / d else NA_real_
      }, numeric(1))
    }
  } else if (outcome_kind == "ordinal") {
    out$outcome_mean_score <- vapply(cc_by, function(r)
      if (length(r)) mean(num[r]) else NA_real_, numeric(1))
    med_idx <- vapply(cc_by, function(r)
      if (length(r)) as.integer(stats::quantile(num[r], 0.5, type = 1))
      else NA_integer_, integer(1))
    out$outcome_median_score <- as.numeric(med_idx)
    # The lower-median category label (quantile type 1 always returns an
    # observed score, so the label is a real category). An ordinal outcome
    # supplied as raw integer scores (possible on the brms cumulative path)
    # carries no level labels; the label column is then NA.
    out$outcome_median_category <- if (is.null(lv_out)) {
      rep(NA_character_, length(med_idx))
    } else {
      ifelse(!is.na(med_idx) & med_idx >= 1 & med_idx <= length(lv_out),
             lv_out[pmax(med_idx, 1L)], NA_character_)
    }
    if (!is.null(w)) {
      out$outcome_mean_score_weighted <- vapply(cc_by, function(r) {
        r <- r[!is.na(w[r])]
        d <- sum(w[r])
        if (d > 0) sum(w[r] * num[r]) / d else NA_real_
      }, numeric(1))
    }
  } else {
    # continuous and count outcomes share the numeric summary; for counts the
    # mean is the observed rate.
    out$outcome_mean <- vapply(cc_by, function(r)
      if (length(r)) mean(num[r]) else NA_real_, numeric(1))
    out$outcome_sd <- vapply(cc_by, function(r)
      if (length(r) > 1) stats::sd(num[r]) else NA_real_, numeric(1))
    out$outcome_median <- vapply(cc_by, function(r)
      if (length(r)) stats::median(num[r]) else NA_real_, numeric(1))
    out$outcome_min <- vapply(cc_by, function(r)
      if (length(r)) min(num[r]) else NA_real_, numeric(1))
    out$outcome_max <- vapply(cc_by, function(r)
      if (length(r)) max(num[r]) else NA_real_, numeric(1))
    if (!is.null(w)) {
      out$outcome_mean_weighted <- vapply(cc_by, function(r) {
        r <- r[!is.na(w[r])]
        d <- sum(w[r])
        if (d > 0) sum(w[r] * num[r]) / d else NA_real_
      }, numeric(1))
    }
  }
  lapply(out, unname)
}

# The data-quality flags of $warnings. Each check is a concrete,
# package-grounded rule (thresholds are the .maihda_describe_* constants at the
# top of this file); a clean sample returns a zero-row frame.
maihda_describe_warnings <- function(data, strata, strata_vars, have_dims,
                                     autobin_info, outcome_name,
                                     outcome_missing, n_total, n_expected,
                                     n_empty, flag_stratum_n, context_df,
                                     n_rows_missing_dimensions) {
  rows <- list()
  add <- function(check, message) {
    rows[[length(rows) + 1]] <<- data.frame(check = check, message = message,
                                            stringsAsFactors = FALSE)
  }

  if (have_dims) {
    # Auto-binned numeric dimensions: the strata are data-dependent.
    for (v in intersect(names(autobin_info), strata_vars)) {
      add("autobinned_dimension", sprintf(
        "Numeric dimension '%s' was auto-binned into 3 groups by make_strata(); the cut-points are data-dependent. Convert it to a factor (or bin it yourself) for explicit, reproducible strata.",
        v))
    }
    for (v in strata_vars) {
      col <- data[[v]]
      if (!is.numeric(col)) next
      k <- length(unique(col[!is.na(col)]))
      if (v %in% names(autobin_info) || k < 3) next
      # Raw numeric dimension: strata by distinct value, but an adjusted model
      # would take it as one linear fixed effect (see ?make_strata).
      add("linear_numeric_dimension", sprintf(
        "Numeric dimension '%s' has %d distinct values; distinct values define distinct strata, but in an adjusted model (maihda()) it would enter as a single linear fixed effect. Wrap it in factor() to treat it categorically.",
        v, k))
    }
    # ID-like dimension: far too many levels to be an intersectional dimension.
    for (v in strata_vars) {
      col <- data[[v]]
      nn <- sum(!is.na(col))
      k <- length(unique(col[!is.na(col)]))
      if (nn >= 20 && (k >= .maihda_describe_id_levels ||
                       k > .maihda_describe_id_ratio * nn)) {
        add("id_like_dimension", sprintf(
          "Dimension '%s' has %d levels (%.0f%% of its non-missing rows are unique); it looks like an identifier rather than an intersectional dimension.",
          v, k, 100 * k / nn))
      }
    }
  }

  if (is.finite(n_empty) && n_empty > 0) {
    add("empty_strata", sprintf(
      "%d of %s expected strata are empty (no observations): the strata space is sparse. Sparse designs are where the Bayesian engine (engine = \"brms\") is most useful.",
      as.integer(n_empty), format(n_expected, big.mark = ",", trim = TRUE)))
  }
  small <- strata[strata$small, , drop = FALSE]
  if (nrow(small) > 0) {
    small <- small[order(small$n), , drop = FALSE]
    shown <- utils::head(small, 3)
    add("small_strata", sprintf(
      "%d strata have n <= %g; smallest: %s. Partial pooling shrinks small-stratum estimates toward the mean; expect wide intervals for these cells.",
      nrow(small), flag_stratum_n,
      paste(sprintf("'%s' (n=%d)", shown$label, shown$n), collapse = ", ")))
  }

  if (n_rows_missing_dimensions > 0) {
    per_dim <- ""
    if (have_dims) {
      cnts <- vapply(strata_vars, function(v) sum(is.na(data[[v]])), integer(1))
      cnts <- cnts[cnts > 0]
      if (length(cnts) > 0) {
        per_dim <- sprintf(" (per dimension: %s)",
                           paste(sprintf("%s %d", names(cnts), cnts),
                                 collapse = ", "))
      }
    }
    add("missing_dimensions", sprintf(
      "%d rows (%.1f%%) have a missing stratum dimension and fall outside every stratum%s. They are excluded from the model's analytic sample.",
      n_rows_missing_dimensions, 100 * n_rows_missing_dimensions / n_total,
      per_dim))
  }

  overall_pct <- if (n_total > 0) 100 * sum(outcome_missing) / n_total else 0
  if (overall_pct > .maihda_describe_high_missing) {
    add("outcome_missingness", sprintf(
      "The outcome '%s' is missing for %d rows (%.1f%%).",
      outcome_name, sum(outcome_missing), overall_pct))
  }
  # Missingness concentrated in particular strata: at least twice the overall
  # rate (with an absolute floor), in cells large enough to mean something.
  conc <- strata[!strata$empty & strata$n >= 10 & strata$n_missing_outcome >= 5 &
                   !is.na(strata$pct_missing_outcome) &
                   strata$pct_missing_outcome >= pmax(2 * overall_pct,
                                                      .maihda_describe_conc_floor), ,
                 drop = FALSE]
  if (nrow(conc) > 0) {
    conc <- conc[order(-conc$pct_missing_outcome), , drop = FALSE]
    shown <- utils::head(conc, 3)
    add("outcome_missingness_concentrated", sprintf(
      "Outcome missingness is concentrated in particular strata (overall %.1f%%): %s.",
      overall_pct,
      paste(sprintf("'%s' %.1f%% (n=%d)", shown$label,
                    shown$pct_missing_outcome, shown$n), collapse = ", ")))
  }

  if (!is.null(context_df)) {
    for (cv in unique(context_df$context)) {
      units <- context_df[context_df$context == cv, , drop = FALSE]
      k <- nrow(units)
      if (k < .maihda_describe_context_levels) {
        add("context_few_levels", sprintf(
          "Context '%s' has %d levels; fewer than ~%d levels weakly identify the context variance and often yield a singular lme4 fit. Consider engine = \"brms\".",
          cv, k, .maihda_describe_context_levels))
      }
      n_small <- sum(units$n < .maihda_describe_small_unit)
      if (n_small > 0) {
        add("context_small_units", sprintf(
          "%d of %d '%s' context units have fewer than %d observations (median %s obs/unit).",
          n_small, k, cv, .maihda_describe_small_unit,
          format(stats::median(units$n))))
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(check = character(), message = character(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Print a MAIHDA sample description
#'
#' @param x A \code{maihda_describe} object from \code{\link{maihda_describe}}.
#' @param digits Decimal places (defaults to the value stored on \code{x}).
#' @param ... Additional arguments (not used).
#' @return Invisibly, \code{x}.
#' @export
print.maihda_describe <- function(x, digits = x$digits, ...) {
  pal <- maihda_palette()
  fnum <- function(v, d = digits) formatC(v, format = "f", digits = d)
  fpct <- function(v) sprintf("%.1f%%", v)
  fint <- function(v) format(v, big.mark = ",", trim = TRUE)

  cat(pal$bold("MAIHDA Sample Description"), "\n", sep = "")
  cat("=========================\n\n")

  cat(sprintf("Outcome: %s | Family: %s%s | Summary: %s\n",
              x$outcome, x$family,
              if (isTRUE(x$family_detected)) " (auto-detected)" else "",
              x$outcome_summary))
  if (!is.null(x$strata_vars)) {
    cat("Dimensions: ", paste(x$strata_vars, collapse = ", "), "\n", sep = "")
  } else {
    cat("Dimensions: (not recorded; strata described by their IDs)\n")
  }
  if (!is.null(x$context_vars)) {
    cat("Context: ", paste(x$context_vars, collapse = ", "), "\n", sep = "")
  }
  if (!identical(x$source, "formula")) {
    cat(sprintf("Described from a fitted %s (engine %s): counts refer to the data supplied to that fit.\n",
                x$source, if (!is.null(x$engine)) x$engine else "?"))
  }
  if (!is.null(x$response_recoding)) {
    rr <- x$response_recoding
    cat(sprintf("Binary outcome recoded at fit time: '%s' = 0 (reference), '%s' = 1 (event).\n",
                rr$level[rr$value == 0][1], rr$level[rr$value == 1][1]))
  }
  if (!is.null(x$longitudinal)) {
    cat(sprintf("Longitudinal fit: rows are measurement occasions; %s individuals, median %s occasions each.\n",
                fint(x$longitudinal$n_individuals),
                format(x$longitudinal$median_occasions)))
  }

  ov <- x$overview
  cat("\n", pal$bold("Sample:"), "\n", sep = "")
  cat(sprintf("  Rows (total):                    %s\n", fint(ov$n_total)))
  cat(sprintf("  Missing outcome:                 %s (%s)\n",
              fint(ov$n_missing_outcome), fpct(ov$pct_missing_outcome)))
  cat(sprintf("  Outside strata (missing dims):   %s\n",
              fint(ov$n_rows_missing_dimensions)))
  cat(sprintf("  Analytic sample (complete case): %s\n", fint(ov$n_analytic)))
  if (!is.null(x$sampling_weights)) {
    cat(sprintf("  Sampling weights: '%s' (sum %s; %s rows with missing/non-positive weights excluded).\n",
                x$sampling_weights, format(ov$sum_weights, big.mark = ","),
                fint(ov$n_invalid_weights)))
    cat("  Weighted summaries are population-weighted point estimates; no design-based variances.\n")
  }

  oo <- x$outcome_overall
  cat("\n", pal$bold("Observed outcome"),
      sprintf(" (%s non-missing): ", fint(oo$n_nonmissing)), sep = "")
  if (x$outcome_kind == "binomial") {
    cat(sprintf("%s events / %s trials (%s%s)\n",
                fint(oo$outcome_events), fint(oo$outcome_trials),
                fpct(100 * oo$outcome_proportion),
                if (!is.null(x$event_level)) paste0(" '", x$event_level, "'")
                else ""))
  } else if (x$outcome_kind == "ordinal") {
    lv_txt <- if (!is.null(x$outcome_levels)) {
      paste(sprintf("%s %s", x$outcome_levels$level,
                    fpct(x$outcome_levels$pct)), collapse = ", ")
    } else {
      ""
    }
    cat(sprintf("mean score %s (%s)\n", fnum(oo$outcome_mean_score), lv_txt))
  } else {
    cat(sprintf("mean %s (SD %s), median %s, range %s to %s\n",
                fnum(oo$outcome_mean), fnum(oo$outcome_sd),
                fnum(oo$outcome_median), fnum(oo$outcome_min),
                fnum(oo$outcome_max)))
  }

  st <- x$strata
  obs <- st[!st$empty, , drop = FALSE]
  cat("\n", pal$bold("Strata"), sep = "")
  if (is.finite(ov$n_strata_expected)) {
    cat(sprintf(" (%s observed / %s expected, %s empty):\n",
                fint(ov$n_strata_observed), fint(ov$n_strata_expected),
                fint(ov$n_empty_strata)))
  } else {
    cat(sprintf(" (%s observed):\n", fint(ov$n_strata_observed)))
  }
  if (nrow(obs) > 0) {
    cat(sprintf("  Cell sizes: min %s, median %s, max %s; %s strata at/below n = %g%s\n",
                fint(min(obs$n)), format(stats::median(obs$n)), fint(max(obs$n)),
                fint(sum(obs$small)), x$flag_stratum_n,
                if (any(obs$small)) " (flagged)" else ""))
    smallest <- utils::head(obs[order(obs$n), , drop = FALSE], 5)
    if (any(obs$small)) {
      cat("  Smallest strata:\n")
      for (i in seq_len(nrow(smallest))) {
        cat(sprintf("    %-40s n = %s\n", smallest$label[i], fint(smallest$n[i])))
      }
    }
    cat("  (full table in $strata)\n")
  }

  if (!is.null(x$dimensions)) {
    cat("\n", pal$bold("Dimensions:"), "\n", sep = "")
    for (v in unique(x$dimensions$dimension)) {
      d <- x$dimensions[x$dimensions$dimension == v, , drop = FALSE]
      if (nrow(d) > 8) {
        shown <- utils::head(d, 6)
        lv_txt <- paste0(paste(sprintf("%s %s (%s)", shown$level,
                                       fint(shown$n), fpct(shown$pct)),
                               collapse = ", "),
                         sprintf(", ... (%d levels)", nrow(d)))
      } else {
        lv_txt <- paste(sprintf("%s %s (%s)", d$level, fint(d$n), fpct(d$pct)),
                        collapse = ", ")
      }
      cat(strwrap(sprintf("%s: %s", v, lv_txt), width = 78, indent = 2,
                  exdent = 4), sep = "\n")
    }
  }

  if (!is.null(x$context)) {
    cat("\n", pal$bold("Context:"), "\n", sep = "")
    for (cv in unique(x$context$context)) {
      u <- x$context[x$context$context == cv, , drop = FALSE]
      cat(sprintf("  %s: %s units; obs/unit min %s, median %s, max %s (per-unit table in $context)\n",
                  cv, fint(nrow(u)), fint(min(u$n)), format(stats::median(u$n)),
                  fint(max(u$n))))
    }
  }

  if (nrow(x$warnings) > 0) {
    cat("\n", pal$bold("Warnings:"), "\n", sep = "")
    for (i in seq_len(nrow(x$warnings))) {
      msg <- strwrap(x$warnings$message[i], width = 76, exdent = 4)
      cat(pal$warn(paste0("  ! ", msg[1])), "\n", sep = "")
      if (length(msg) > 1) {
        for (m in msg[-1]) cat(pal$warn(m), "\n", sep = "")
      }
    }
  } else {
    cat("\n", pal$muted("No data-quality warnings."), "\n", sep = "")
  }
  invisible(x)
}

#' Plot a MAIHDA sample description
#'
#' Descriptive plots for a \code{\link{maihda_describe}} object: the
#' stratum-size distribution (the sparsity of the strata space), the
#' family-aware outcome distribution, and per-variable missingness.
#'
#' @param x A \code{maihda_describe} object.
#' @param type One of \code{"stratum_size"} (default; histogram of stratum
#'   sizes with the small-stratum threshold marked, empty strata included as
#'   zeros when they were enumerated), \code{"outcome"} (histogram of a
#'   continuous/count outcome, or category bars for a binomial/ordinal one), or
#'   \code{"missingness"} (percent missing per variable).
#' @param ... Additional arguments (not used).
#' @return A \pkg{ggplot2} object, extendable with the usual \code{+} grammar.
#' @examples
#' desc <- maihda_describe(BMI ~ Age + (1 | Gender:Race:Education),
#'                         data = maihda_health_data)
#' plot(desc, type = "stratum_size")
#' plot(desc, type = "outcome")
#' plot(desc, type = "missingness")
#' @export
plot.maihda_describe <- function(x, type = c("stratum_size", "outcome",
                                             "missingness"), ...) {
  type <- match.arg(type)
  switch(type,
         stratum_size = maihda_describe_plot_sizes(x),
         outcome = maihda_describe_plot_outcome(x),
         missingness = maihda_describe_plot_missingness(x))
}

maihda_describe_plot_sizes <- function(x) {
  st <- x$strata
  ov <- x$overview
  subtitle <- if (is.finite(ov$n_strata_expected)) {
    sprintf("%d observed of %s expected strata; %s empty; %d at/below n = %g",
            ov$n_strata_observed, format(ov$n_strata_expected, trim = TRUE),
            format(ov$n_empty_strata, trim = TRUE),
            sum(st$small), x$flag_stratum_n)
  } else {
    sprintf("%d observed strata; %d at/below n = %g",
            ov$n_strata_observed, sum(st$small), x$flag_stratum_n)
  }
  ggplot2::ggplot(st, ggplot2::aes(x = .data$n)) +
    ggplot2::geom_histogram(bins = min(30L, max(5L, length(unique(st$n)))),
                            fill = "#0072B2", colour = "white") +
    ggplot2::geom_vline(xintercept = x$flag_stratum_n, colour = "#D55E00",
                        linetype = "dashed") +
    ggplot2::labs(title = "Stratum sizes", subtitle = subtitle,
                  x = "Individuals per stratum", y = "Strata") +
    theme_maihda()
}

maihda_describe_plot_outcome <- function(x) {
  if (!is.null(x$outcome_levels)) {
    d <- x$outcome_levels
    d$level <- factor(d$level, levels = d$level)
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data$level, y = .data$pct)) +
        ggplot2::geom_col(fill = "#0072B2") +
        ggplot2::labs(title = sprintf("Observed outcome: %s", x$outcome),
                      subtitle = sprintf("Family: %s", x$family),
                      x = x$outcome, y = "% of non-missing") +
        theme_maihda()
    )
  }
  obs <- x$observations
  vals <- data.frame(value = obs$numerator / obs$denominator)
  ggplot2::ggplot(vals, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = 30, fill = "#0072B2", colour = "white") +
    ggplot2::labs(title = sprintf("Observed outcome: %s", x$outcome),
                  subtitle = sprintf("Family: %s", x$family),
                  x = x$outcome, y = "Observations") +
    theme_maihda()
}

maihda_describe_plot_missingness <- function(x) {
  m <- x$missingness
  m$variable <- stats::reorder(m$variable, m$pct_missing)
  all_zero <- all(m$n_missing == 0)
  ggplot2::ggplot(m, ggplot2::aes(x = .data$pct_missing, y = .data$variable)) +
    ggplot2::geom_col(fill = "#0072B2") +
    # Anchor the axis at 0 (and give it some span) so a fully complete sample
    # shows flat zero bars instead of an empty panel around 0.
    ggplot2::expand_limits(x = c(0, 1)) +
    ggplot2::labs(title = "Missingness by variable",
                  subtitle = if (all_zero) {
                    "No missing values"
                  } else {
                    sprintf("Of %s rows",
                            format(x$overview$n_total, big.mark = ","))
                  },
                  x = "% of rows missing", y = NULL) +
    theme_maihda()
}
