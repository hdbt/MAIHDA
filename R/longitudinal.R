# Longitudinal / growth-curve MAIHDA (3-level: occasions within individuals
# within intersectional strata).
#
# The cross-sectional MAIHDA models assume one intercept-only stratum random
# effect and read a single between-stratum variance. A longitudinal MAIHDA
# instead fits a growth curve with random INTERCEPTS AND SLOPES on time at both
# the individual (level 2) and stratum (level 3) levels, so the between-stratum
# variance -- and hence the VPC -- becomes a function of time:
#
#   y_tij = beta0 + beta1 * time + ... + (u0j + u1j * time)        [stratum  L3]
#                                      + (v0ij + v1ij * time)       [person   L2]
#                                      + e_tij                       [occasion L1]
#
#   VarS(t) = a(t)' Sigma_s a(t),  a(t) = (1, t, t^2, ...)'
#   VPC_S(t) = VarS(t) / (VarS(t) + VarI(t) + level1var)
#
# This file implements the validation, the 3-level growth formula builder, the
# random-effect covariance-block extractors (the longitudinal generalisation of
# the single-cell maihda_stratum_variance_lme4()), the time-varying VPC summary
# (lme4 + brms), the additive-vs-multiplicative PCV, and a small print helper.
# The intercept-only guards elsewhere are left untouched: a longitudinal fit is
# tagged with $longitudinal_info and routed to these helpers, so every other
# model still rejects random slopes.
#
# Method follows Bell, Evans, Holman & Leckie (2024, Soc Sci Med 351:116955,
# <doi:10.1016/j.socscimed.2024.116955>).
#
# INTERNAL TIME CENTERING. The growth terms are fit on internally centered time
# (t - min(t)) whenever the time axis does not start at 0 (age, calendar year,
# waves coded 10, 11, ...). The raw polynomial basis (1, t, t^2, ...) over an
# offset range is near-collinear and its random-effect covariance scales are
# wildly heterogeneous, so lme4 can converge to a false optimum -- observed on a
# quadratic fit anchored at wave 10: ~128 log-likelihood units below the true
# optimum with NO convergence warning, and a baseline between-stratum variance
# three orders of magnitude off. Centering makes the offset-axis fit the SAME
# optimization problem as the equivalent 0-anchored coding, and anchors the
# covariance blocks at the observed baseline (Sigma[1,1] is the baseline
# intercept variance). Every user-facing time -- ref_time, reporting grids,
# plots, prediction newdata -- stays on the original scale; only the model terms
# use the derived column below, and a(t)' Sigma a(t) evaluations subtract the
# stored centering offset first.

# ---- internal time centering -------------------------------------------------

# Reserved column name for the internally centered time variable, mirroring the
# reserved-column pattern of .maihda_dim_* (decompose_maihda.R) and .maihda_sw
# (design_weights.R). Written into the analytic data by fit_maihda() when
# centering applies; prediction newdata rebuilds it from the original time column
# (maihda_prepare_prediction_data), so callers always work in original units.
.maihda_ctime_col <- ".maihda_ctime"

# The model variable the growth terms were built on: the original time column
# when no centering was needed, else the derived centered column. NULL-safe for
# maihda_model objects stored by package versions that predate internal
# centering (they fall back to the original column, matching their fits).
maihda_lng_time_term <- function(lng) {
  if (!is.null(lng$time_term)) lng$time_term else lng$time
}

# The centering offset (0 when the fit used raw time). An original-scale time t
# maps to the model's coefficient coordinates as t - center; every
# a(t)' Sigma a(t) evaluation of a (possibly centered) covariance block must
# subtract this first. NULL-safe like maihda_lng_time_term().
maihda_lng_time_center <- function(lng) {
  if (!is.null(lng$time_center) && is.numeric(lng$time_center) &&
      length(lng$time_center) == 1 && is.finite(lng$time_center)) {
    lng$time_center
  } else {
    0
  }
}

# Centering offset for a time vector: its minimum finite value, anchoring the
# growth terms at the observed baseline. 0 (no centering; the historical raw-time
# path, byte-identical results) when the axis already starts at 0.
maihda_longitudinal_center <- function(time_values) {
  tv <- time_values[is.finite(time_values)]
  if (length(tv) == 0) {
    return(0)
  }
  m <- min(tv)
  if (m == 0) 0 else m
}

# ---- validation -------------------------------------------------------------

#' Validate a longitudinal (id / time) specification
#'
#' @param id Single column name: the person/unit identifier (level 2).
#' @param time Single column name: a numeric measurement-time variable (level 1).
#' @param time_degree Integer >= 1: polynomial degree of the growth curve (1 =
#'   linear). brms supports degree 1 only.
#' @param data The model data.
#' @param engine The fitting engine; only "lme4"/"brms" support the 3-level
#'   growth structure.
#' @param sampling_weights,context Must be NULL -- design-weighted and contextual
#'   longitudinal models are out of scope.
#' @param formula Optional model formula, used only to tell a package-derived
#'   refit (whose formula already references the reserved \code{.maihda_ctime}
#'   column) from a fresh user call when guarding that reserved name.
#' @return A list \code{list(id, time, time_degree)}.
#' @keywords internal
maihda_validate_longitudinal <- function(id, time, time_degree, data,
                                         engine = "lme4",
                                         sampling_weights = NULL,
                                         context = NULL,
                                         formula = NULL) {
  if (!is.character(time) || length(time) != 1 || is.na(time) || !nzchar(time)) {
    stop("'time' must be a single column name (a character string) naming the ",
         "measurement-time variable for a longitudinal MAIHDA.", call. = FALSE)
  }
  if (is.null(id) || !is.character(id) || length(id) != 1 || is.na(id) || !nzchar(id)) {
    stop("A longitudinal MAIHDA needs 'id', a single column name identifying the ",
         "person/unit measured repeatedly (level 2). Supply it alongside 'time'.",
         call. = FALSE)
  }
  missing_cols <- setdiff(c(id, time), names(data))
  if (length(missing_cols) > 0) {
    stop("Longitudinal column(s) not found in data: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  if (identical(id, time)) {
    stop("'id' and 'time' must name different columns.", call. = FALSE)
  }
  if (identical(id, "stratum") || identical(time, "stratum")) {
    stop("'id'/'time' may not be named 'stratum' (reserved for the intersectional ",
         "grouping).", call. = FALSE)
  }
  if (.maihda_ctime_col %in% c(id, time)) {
    stop("'id'/'time' may not use the reserved internal column name '",
         .maihda_ctime_col, "' (it carries the internally centered time).",
         call. = FALSE)
  }
  if (!is.numeric(data[[time]])) {
    stop("The 'time' column '", time, "' must be numeric (the growth curve is a ",
         "polynomial in time). Code occasions/waves as 0, 1, 2, ... or use age.",
         call. = FALSE)
  }
  if (!is.numeric(time_degree) || length(time_degree) != 1 || is.na(time_degree) ||
      time_degree < 1 || time_degree != floor(time_degree)) {
    stop("'time_degree' must be a single whole number >= 1 (1 = linear growth).",
         call. = FALSE)
  }
  time_degree <- as.integer(time_degree)

  # Genuinely repeated measures: at least one id must appear more than once, else
  # the level-2 (person) random effects are unidentified and this is not longitudinal.
  ids <- data[[id]]
  if (!any(duplicated(ids[!is.na(ids)]))) {
    stop("The data do not look longitudinal: every '", id, "' value is unique, so ",
         "there are no repeated measurements to model. Supply long-format data ",
         "(one row per measurement occasion).", call. = FALSE)
  }

  if (!engine %in% c("lme4", "brms")) {
    stop("A longitudinal MAIHDA (id/time) is supported only by engine = \"lme4\" ",
         "or \"brms\"; the 3-level random-slope growth structure has no ",
         "wemix/ordinal representation. This model uses engine = \"", engine, "\".",
         call. = FALSE)
  }
  if (!is.null(sampling_weights)) {
    stop("A design-weighted longitudinal MAIHDA is out of scope: 'sampling_weights' ",
         "is not supported with 'id'/'time'. Fit unweighted (lme4/brms).",
         call. = FALSE)
  }
  if (!is.null(context)) {
    stop("A contextual-and-longitudinal model (context x stratum x time) is out of ",
         "scope: 'context' is not supported with 'id'/'time'.", call. = FALSE)
  }
  if (identical(engine, "brms") && time_degree > 1L) {
    stop("The brms longitudinal engine currently supports linear growth only ",
         "(time_degree = 1). Use engine = \"lme4\" for higher-degree growth.",
         call. = FALSE)
  }

  # Guard the reserved centered-time column against silently overwriting a user
  # variable of the same name. Only fires when centering will actually occur
  # (min(time) != 0 writes the column) AND the formula does not already reference
  # it -- a formula that does is a package-derived refit (maihda()'s null/adjusted
  # pair re-entering fit_maihda with the first fit's original_data), whose column
  # is the package's own and is recomputed idempotently.
  if (.maihda_ctime_col %in% names(data) &&
      maihda_longitudinal_center(data[[time]]) != 0 &&
      (is.null(formula) || !.maihda_ctime_col %in% all.vars(formula))) {
    stop("'", .maihda_ctime_col, "' is a reserved internal column name for the ",
         "longitudinal MAIHDA fit (it carries the internally centered time, ",
         "needed because '", time, "' does not start at 0), but 'data' already ",
         "contains a column of that name; fitting would overwrite it. Rename ",
         "your '", .maihda_ctime_col, "' column before fitting.", call. = FALSE)
  }

  list(id = id, time = time, time_degree = time_degree)
}

#' Guard against longitudinal ids reused across strata
#'
#' \code{(time | id)} treats every row sharing an id value as the SAME person,
#' so ids that are only unique within a site or group (person "1" in every
#' site) silently merge different people's trajectories into one level-2 unit.
#' An id appearing in more than one stratum is the observable symptom of that
#' (a person's intersectional stratum is a person-level attribute, so a
#' genuinely unique id maps to exactly one stratum); reject it with guidance
#' rather than guessing a nesting. Called after strata resolution, when the
#' \code{stratum} column exists. Rows missing id or stratum are ignored (the
#' engines drop them).
#'
#' @param data The model data carrying the resolved \code{stratum} column.
#' @param id Name of the person/unit identifier column.
#' @return \code{NULL}, invisibly; stops on ambiguous ids.
#' @keywords internal
maihda_check_longitudinal_ids <- function(data, id) {
  ids <- data[[id]]
  strat <- data[["stratum"]]
  ok <- !is.na(ids) & !is.na(strat)
  pairs <- unique(data.frame(id = as.character(ids[ok]),
                             stratum = as.character(strat[ok]),
                             stringsAsFactors = FALSE))
  dup <- unique(pairs$id[duplicated(pairs$id)])
  if (length(dup) > 0) {
    shown <- paste(utils::head(sort(dup), 5), collapse = ", ")
    if (length(dup) > 5) shown <- paste0(shown, ", ...")
    stop("Longitudinal '", id, "' values appear in more than one stratum (",
         length(dup), " id value(s): ", shown, "). The growth model treats ",
         "every row sharing an id as the same person, so ids that are only ",
         "unique within a site or group would merge different people's ",
         "trajectories. Give every person a globally unique id (e.g. ",
         "paste(site, ", id, ")) before fitting. If instead a person's ",
         "stratum variables changed between occasions, fix the stratum at a ",
         "single reference occasion (e.g. baseline) first -- the strata are ",
         "person-level groupings and must be constant within a person.",
         call. = FALSE)
  }
  invisible(NULL)
}

# Polynomial-in-time term labels, e.g. c("wave", "I(wave^2)") for degree 2. The
# first random/fixed term is the linear time; higher degrees use I(time^k) so the
# design vector at time t is a(t) = (1, t, t^2, ..., t^degree).
maihda_time_terms <- function(time, time_degree) {
  qtime <- maihda_quote_name(time)
  terms <- qtime
  if (time_degree >= 2) {
    # Use the already-quoted name inside I() too, so a non-syntactic time column
    # (e.g. "time point") yields valid formula text I(`time point`^2) rather than
    # the unparseable I(time point^2). For a syntactic name maihda_quote_name() is
    # a no-op, so this also matches the (Intercept), time, I(time^2), ... names
    # lme4 records in VarCorr (see maihda_re_block_lme4()).
    terms <- c(terms, sprintf("I(%s^%d)", qtime, 2:time_degree))
  }
  terms
}

#' Build the 3-level growth formula for a longitudinal MAIHDA
#'
#' Given a base formula already carrying the covariates and the resolved stratum
#' grouping (\code{y ~ covars + (1 | stratum)}), returns the growth formula
#' \code{y ~ covars + time(+ I(time^2)...) + (time... | id) + (time... | stratum)}:
#' the time polynomial enters the fixed part (if absent) and a random
#' intercept+slope block is placed at both the individual and stratum levels. Any
#' random effects in the base formula are replaced by this canonical structure.
#'
#' When the growth terms are built on internally centered time (\code{time} is
#' the derived \code{.maihda_ctime} column; see the file header), any bare
#' raw-time polynomial the user wrote in the fixed part (\code{orig_time},
#' \code{I(orig_time^2)}, ...) is \emph{replaced} by the centered terms rather
#' than kept alongside them, which would be perfectly collinear.
#'
#' @param base_formula The resolved formula (fixed part + stratum grouping).
#' @param id,time,time_degree The longitudinal specification; \code{time} is the
#'   model variable the growth terms are built on (the centered column when
#'   centering applies).
#' @param orig_time The user's original time column name; differs from
#'   \code{time} only when centering applies.
#' @return The growth formula (same environment as \code{base_formula}).
#' @keywords internal
#' @importFrom stats update as.formula terms
maihda_longitudinal_formula <- function(base_formula, id, time, time_degree,
                                        orig_time = time) {
  poly_terms <- maihda_time_terms(time, time_degree)
  ptime <- paste(poly_terms, collapse = " + ")

  fixed <- maihda_nobars(base_formula)
  fixed_labels <- attr(stats::terms(fixed), "term.labels")
  if (!identical(orig_time, time)) {
    # Centering active: drop any bare raw-time polynomial from the fixed part --
    # the centered terms below replace it (keeping both would be collinear).
    raw_terms <- intersect(maihda_time_terms(orig_time, time_degree), fixed_labels)
    if (length(raw_terms) > 0) {
      fixed <- stats::update(fixed, stats::as.formula(
        paste(". ~ . -", paste(raw_terms, collapse = " - "))))
      fixed_labels <- attr(stats::terms(fixed), "term.labels")
    }
  }
  add_fixed <- setdiff(poly_terms, fixed_labels)
  if (length(add_fixed) > 0) {
    fixed <- stats::update(fixed, stats::as.formula(
      paste(". ~ . +", paste(add_fixed, collapse = " + "))))
  }

  re <- sprintf("(%s | %s) + (%s | %s)",
                ptime, maihda_quote_name(id),
                ptime, maihda_quote_name("stratum"))
  stats::update(fixed, stats::as.formula(paste(". ~ . +", re)))
}

#' Build the adjusted-model formula for a longitudinal MAIHDA decomposition
#'
#' The longitudinal analogue of \code{\link{maihda_adjusted_formula}}: the null
#' growth model plus the stratum dimensions' additive main effects AND their
#' interactions with the time polynomial (\code{dim:time}), so the remaining
#' stratum-level intercept/slope variance is the intersectional interaction beyond
#' additive. Auto-binned numeric dimensions reuse their reconstructed tertile
#' factor (\code{.maihda_dim_*}, via \code{\link{maihda_adjusted_terms}}).
#'
#' @param null_formula The fitted null growth formula.
#' @param strata_vars,autobin_info,data Stratum metadata (as for
#'   \code{maihda_adjusted_formula}).
#' @param time,time_degree The longitudinal specification; \code{time} must be
#'   the variable the growth terms were built on (the internally centered column
#'   for a centered fit, \code{maihda_lng_time_term()}), so the \code{dim:time}
#'   interactions reference the same terms as the null formula.
#' @return A list with \code{formula} and \code{data}, or \code{NULL} if fewer
#'   than two dimensions are available.
#' @keywords internal
#' @importFrom stats update as.formula
maihda_longitudinal_adjusted_formula <- function(null_formula, strata_vars,
                                                 autobin_info, data, time,
                                                 time_degree) {
  if (is.null(strata_vars) || length(strata_vars) < 2) {
    return(NULL)
  }
  adj <- maihda_adjusted_terms(strata_vars, autobin_info, data)
  main <- vapply(adj$terms, maihda_quote_name, character(1))
  poly_terms <- maihda_time_terms(time, time_degree)
  # dim main effects + every dim x (time polynomial) interaction.
  inter <- as.vector(outer(main, poly_terms, function(d, p) paste0(d, ":", p)))
  rhs <- paste(c(main, inter), collapse = " + ")
  adjusted_formula <- stats::update(null_formula,
                                    stats::as.formula(paste(". ~ . +", rhs)))
  list(formula = adjusted_formula, data = adj$data)
}

# ---- random-effect covariance blocks ----------------------------------------

# The ordered k x k random-effect covariance matrix of a grouping factor (lme4),
# rows/cols in the order (Intercept), time, I(time^2), ..., so the time design
# vector a(t) = (1, t, t^2, ...) lines up with it. This is the longitudinal
# generalisation of maihda_stratum_variance_lme4()'s single (Intercept,Intercept)
# cell.
maihda_re_block_lme4 <- function(model, group, time, time_degree) {
  vc <- lme4::VarCorr(model)
  if (!group %in% names(vc)) {
    stop("No '", group, "' random effect found in the model.", call. = FALSE)
  }
  m <- as.matrix(vc[[group]])
  want <- c("(Intercept)", maihda_time_terms(time, time_degree))
  idx <- match(want, rownames(m))
  if (anyNA(idx)) {
    stop("The '", group, "' random effect is missing the growth term(s): ",
         paste(want[is.na(idx)], collapse = ", "),
         ". A longitudinal MAIHDA needs a random intercept and slope on '", time,
         "' at this level.", call. = FALSE)
  }
  m[idx, idx, drop = FALSE]
}

# Posterior-mean ordered covariance block of a grouping factor (brms), for the
# point-estimate PCV / components table. Built from the SD and correlation draws
# (maihda_re_cov_draws_brms), avoiding any dependence on the dimension order of
# brms::VarCorr()'s $cov array. The brms longitudinal engine is restricted to
# linear growth (time_degree 1), so the block is the 2x2 (intercept, slope).
maihda_re_block_brms <- function(model, group, time, time_degree) {
  if (time_degree != 1L) {
    stop("The brms longitudinal engine supports linear growth only ",
         "(time_degree = 1).", call. = FALSE)
  }
  draws <- maihda_posterior_draws_brms(model)
  blk <- maihda_re_cov_draws_brms(draws, group, time)
  v0 <- mean(blk$v0); v1 <- mean(blk$v1); cv <- mean(blk$cov)
  matrix(c(v0, cv, cv, v1), nrow = 2,
         dimnames = list(c("(Intercept)", time), c("(Intercept)", time)))
}

# Engine-agnostic ordered covariance block (point estimate) for a maihda_model.
# The block is in the model's coefficient coordinates -- CENTERED time when the
# fit used the internal centering -- so evaluations at an original-scale time t
# must subtract maihda_lng_time_center() first.
maihda_re_block <- function(object, group) {
  lng <- object$longitudinal_info
  time_term <- maihda_lng_time_term(lng)
  if (identical(object$engine, "lme4")) {
    maihda_re_block_lme4(object$model, group, time_term, lng$time_degree)
  } else if (identical(object$engine, "brms")) {
    maihda_re_block_brms(object$model, group, time_term, lng$time_degree)
  } else {
    stop("Longitudinal MAIHDA is supported only for lme4/brms.", call. = FALSE)
  }
}

# Between-level variance implied by a covariance block at time(s) t:
# a(t)' Sigma a(t) with a(t) = (1, t, t^2, ...). Vectorised over t. t must be on
# the block's own coefficient scale: for a fit on internally centered time,
# convert an original-scale time first (t - maihda_lng_time_center(lng)).
maihda_var_at_time <- function(Sigma, t) {
  degree <- nrow(Sigma) - 1L
  vapply(t, function(ti) {
    a <- ti^(0:degree)
    as.numeric(crossprod(a, Sigma %*% a))
  }, numeric(1))
}

# Between-level variance of the INSTANTANEOUS SLOPE implied by a covariance
# block at time(s) t: b(t)' Sigma b(t) with b(t) = d a(t)/dt = (0, 1, 2t, 3t^2,
# ...). For a linear (2x2) block this is Sigma[2,2] at every t; for higher
# degrees the slope variance is itself time-varying, so -- exactly as for
# maihda_var_at_time() -- t must be on the block's (possibly centered)
# coefficient scale. Vectorised over t.
maihda_slope_var_at_time <- function(Sigma, t) {
  degree <- nrow(Sigma) - 1L
  vapply(t, function(ti) {
    b <- c(0, seq_len(degree) * ti^(seq_len(degree) - 1L))
    as.numeric(crossprod(b, Sigma %*% b))
  }, numeric(1))
}

# Stop with a consistent message when a cross-sectional, single-value-per-stratum
# summary (a scalar BLUP ranking, predicted value, or BLUP-based plot) is
# requested for a longitudinal MAIHDA. In a growth model each stratum's estimand
# is a TRAJECTORY (random intercept + slope(s)), so collapsing it to one number
# produces a cross-sectional-looking result that is not the right quantity. Point
# the user to the trajectory tools instead.
maihda_stop_longitudinal_scalar <- function(what) {
  stop(what, " is not defined for a longitudinal MAIHDA: each stratum is a ",
       "trajectory (random intercept + slope), not a single value. Use ",
       "predict(type = \"strata\") for the per-stratum intercept and slope, ",
       "plot(type = \"trajectories\") for the stratum mean trajectories, or ",
       "plot(type = \"vpc_trajectory\") for the time-varying VPC.", call. = FALSE)
}

# A time grid for reporting VPC(t): the observed unique times when few, else a
# 25-point grid spanning their range.
maihda_longitudinal_time_grid <- function(time_values) {
  u <- sort(unique(time_values[is.finite(time_values)]))
  if (length(u) <= 12L) {
    return(u)
  }
  seq(min(u), max(u), length.out = 25L)
}

# Move a reporting frame to a single time, setting EVERY column the model reads as
# time -- the centered term the growth block is fit on AND the user's original time
# column.
#
# Under internal centering the two are different columns, and the original one does
# not simply disappear from the fixed part: maihda_longitudinal_formula() drops
# only the BARE raw-time polynomial that the centered terms replace, so anything
# else the user wrote on the raw axis survives into the fitted model -- a covariate
# interaction (x:wave), a transformation (poly(wave, 2), log(wave)), or a formula
# offset. Moving only the centered column leaves those terms evaluated at each
# row's OWN observed time while the result is reported as "at time t", which biases
# any quantity built from the linear predictor on the grid (the count VPC's
# marginal lambda). The columns must move together, so this is the one place that
# does it.
maihda_longitudinal_set_time <- function(newdata, time_term, t_c,
                                         orig_time = time_term, center = 0) {
  newdata[[time_term]] <- t_c
  if (!identical(orig_time, time_term)) {
    # Assigned unconditionally, not only when the column is already present: a
    # formula offset such as offset(0.5 * wave) is re-evaluated against this frame
    # and needs the raw time even when the model frame carries only the derived
    # 'offset(...)' column.
    newdata[[orig_time]] <- t_c + center
  }
  newdata
}

# Level-1 (residual) variance of a longitudinal VPC evaluated on a time grid
# (lme4). For Gaussian and binomial families the level-1 variance is a constant
# (attr(sc)^2, pi^2/3, or 1), so the single maihda_residual_variance_lme4() value
# is recycled to length(t_c) and the trajectory is unaffected. For Poisson /
# negative-binomial the latent-scale level-1 variance log1p(1/lambda[ + 1/theta])
# depends on the MARGINAL expected count lambda, which changes with time in a
# growth model, so a single sample-wide average reused at every time point biases
# VPC(t). Evaluate it at each grid time instead: set the model's time term to that
# value (averaging the per-row marginal counts over the observed covariate rows
# before the transform, as the cross-sectional path averages over the sample --
# see maihda_count_level1_variance()) and apply the
# log-normal mean correction v(t)/2 from the total growth-block variance v_t at
# that time. `t_c` is on the model's (possibly centered) coefficient scale; `v_t`
# is VarStratum(t) + VarId(t) at the same times, aligned to `t_c`.
maihda_longitudinal_resid_grid_lme4 <- function(model, time_term, t_c, v_t,
                                                orig_time = time_term, center = 0) {
  fam <- maihda_family(model)
  is_count <- !is.null(fam) && identical(fam$link, "log") &&
    fam$family %in% c("poisson", "negbinomial")
  if (!is_count) {
    return(rep(maihda_residual_variance_lme4(model), length(t_c)))
  }
  theta <- if (identical(fam$family, "negbinomial")) {
    maihda_negbin_theta_lme4(model)
  } else {
    Inf
  }
  w <- maihda_fit_prior_weights(model)
  frame <- maihda_model_frame(model)
  # Offset on the modified-time frame. predict.merMod here would DROP an external offset=
  # (biasing the marginal count lambda, hence VPC(t)) and ERROR on a formula offset()
  # term, so build the fixed-effects linear predictor directly. A FORMULA offset such as
  # offset(0.5 * time) is RE-EVALUATED at each grid time -- it shifts the marginal count
  # and must track the modified time -- whereas an EXTERNAL offset= vector has no
  # expression to re-evaluate and is kept per fitted row (a fixed per-row exposure). Under
  # internal time centering both the offset AND any surviving raw-time fixed term
  # reference the ORIGINAL time column, so maihda_longitudinal_set_time() moves it in
  # lockstep with the centered term. off_ext is NULL for a fit with no external
  # offset; a no-offset fit reproduces predict(re.form = NA).
  off_ext <- maihda_fitted_offset_external(model)
  has_formula_off <- !is.null(
    attr(stats::terms(maihda_nobars(stats::formula(model))), "offset"))
  vapply(seq_along(t_c), function(j) {
    nd <- maihda_longitudinal_set_time(frame, time_term, t_c[j],
                                       orig_time = orig_time, center = center)
    off_j <- off_ext
    if (has_formula_off) {
      fo <- maihda_lme4_formula_offset_at(model, nd)
      off_j <- if (is.null(off_j)) fo else off_j + fo
    }
    eta <- maihda_lme4_fixed_link(model, nd, offset = off_j)
    mu <- pmax(exp(eta + v_t[j] / 2), .Machine$double.eps)
    maihda_count_level1_variance(mu, theta = theta, w = w)
  }, numeric(1))
}

# Level-1 (residual) variance of a longitudinal count VPC at a SINGLE time `t_c`
# (brms). The brms analogue of the count branch of maihda_longitudinal_resid_grid_lme4:
# the marginal expected count lambda is the posterior-mean plug-in exp(eta(t) +
# v(t)/2) with eta the posterior-mean fixed-part linear predictor at time `t_c`
# (marginalized over the observed covariate rows) -- matching the "lambda held at
# its posterior mean" treatment the per-draw residual path documents
# (maihda_residual_variance_draws_brms). Returns a scalar for Poisson; for the
# negative binomial the 'shape' (theta) draws are propagated per draw exactly as in
# that path, so the return is a per-draw vector. `v_t` is VarStratum(t) + VarId(t)
# (posterior-mean) at `t_c`. Only called for count families.
maihda_longitudinal_resid_at_brms <- function(model, draws, time_term, t_c, v_t,
                                              orig_time = time_term, center = 0) {
  nd <- maihda_longitudinal_set_time(model$data, time_term, t_c,
                                     orig_time = orig_time, center = center)
  eta <- maihda_brms_linpred_mean(model, newdata = nd, re_formula = NA)
  mu <- pmax(exp(eta + v_t / 2), .Machine$double.eps)
  w <- maihda_fit_prior_weights(model)
  fam <- maihda_family(model)
  if (identical(fam$family, "negbinomial")) {
    if (!"shape" %in% names(draws)) {
      stop("Could not extract the negative-binomial 'shape' (theta) draws from ",
           "the brms posterior.", call. = FALSE)
    }
    shape_d <- as.numeric(draws[["shape"]])
    return(vapply(shape_d,
                  function(s) maihda_count_level1_variance(mu, theta = s, w = w),
                  numeric(1)))
  }
  maihda_count_level1_variance(mu, w = w)
}

# ---- time-varying VPC summary -----------------------------------------------

# Percentile band for one grid time's bootstrap VPC column. Applies the same
# absolute floor (10) AND majority-of-eligible rule as maihda_bootstrap_ci() so a
# time whose refits mostly produced non-finite VPCs yields NA rather than a band from
# a biased handful of survivors. No draws are legitimately excluded per time, so the
# eligible count is n_boot; a failing time returns NA (soft) rather than erroring the
# whole trajectory the way maihda_bootstrap_ci() does for the headline interval.
maihda_longitudinal_vpc_band <- function(col, n_boot, conf_level) {
  col <- col[is.finite(col)]
  if (length(col) < 10L || length(col) < ceiling(0.5 * n_boot)) {
    return(c(NA_real_, NA_real_))
  }
  a <- 1 - conf_level
  stats::quantile(col, c(a / 2, 1 - a / 2), names = FALSE)
}

#' Time-varying VPC summary for a longitudinal MAIHDA (lme4)
#'
#' @param object A longitudinal \code{maihda_model} (lme4 engine).
#' @param bootstrap,n_boot,conf_level Parametric-bootstrap controls for the VPC(t)
#'   band.
#' @return A list with \code{vpc_result} (the reference-time VPC, for the headline
#'   print), \code{variance_components}, and \code{longitudinal} (the trajectory).
#' @keywords internal
#' @importFrom lme4 refit
maihda_longitudinal_summary_lme4 <- function(object, bootstrap = FALSE,
                                             n_boot = 1000, conf_level = 0.95) {
  lng <- object$longitudinal_info
  model <- object$model
  time_term <- maihda_lng_time_term(lng)
  center <- maihda_lng_time_center(lng)
  Sigma_s <- maihda_re_block_lme4(model, "stratum", time_term, lng$time_degree)
  Sigma_i <- maihda_re_block_lme4(model, lng$id, time_term, lng$time_degree)

  grid <- maihda_longitudinal_time_grid(object$data[[lng$time]])
  ref_time <- lng$ref_time
  # Times are user-facing (original scale); the covariance blocks are in the
  # (possibly centered) coefficient coordinates, so every a(t)' Sigma a(t)
  # evaluation subtracts the centering offset first.
  grid_c <- grid - center
  ref_c <- ref_time - center

  # Between-level variances on the grid, and the level-1 residual evaluated AT
  # each grid time. For Gaussian/binomial the residual is constant; for count
  # families it tracks the time-varying marginal mean, so it is a vector, not one
  # sample-wide average reused at every time (see maihda_longitudinal_resid_grid_lme4).
  var_s_grid <- maihda_var_at_time(Sigma_s, grid_c)
  var_i_grid <- maihda_var_at_time(Sigma_i, grid_c)
  var_s_ref <- maihda_var_at_time(Sigma_s, ref_c)
  var_i_ref <- maihda_var_at_time(Sigma_i, ref_c)
  resid_grid <- maihda_longitudinal_resid_grid_lme4(model, time_term, grid_c,
                                                    var_s_grid + var_i_grid,
                                                    orig_time = lng$time, center = center)
  resid_ref <- maihda_longitudinal_resid_grid_lme4(model, time_term, ref_c,
                                                   var_s_ref + var_i_ref,
                                                   orig_time = lng$time, center = center)
  # Scalar residual for the components table / headline: the reference-time value
  # (a single latent-scale number that the trajectory generalises over time).
  var_resid <- as.numeric(resid_ref)

  vpc_t_est <- var_s_grid / (var_s_grid + var_i_grid + resid_grid)
  ref_vpc <- as.numeric(var_s_ref / (var_s_ref + var_i_ref + resid_ref))

  vpc_lower <- rep(NA_real_, length(grid))
  vpc_upper <- rep(NA_real_, length(grid))
  ref_ci <- NULL
  # Count contributing bootstrap draws whose refit optimiser did not converge, so the
  # reported n_boot_ok does not silently imply convergence (see bootstrap_vpc()).
  n_nonconv <- 0L
  if (bootstrap) {
    boot <- matrix(NA_real_, nrow = n_boot, ncol = length(grid))
    ref_boot <- rep(NA_real_, n_boot)
    sim <- stats::simulate(model, nsim = n_boot)
    for (i in seq_len(n_boot)) {
      tryCatch({
        bm <- lme4::refit(model, newresp = sim[[i]])
        Ss <- maihda_re_block_lme4(bm, "stratum", time_term, lng$time_degree)
        Si <- maihda_re_block_lme4(bm, lng$id, time_term, lng$time_degree)
        vs <- maihda_var_at_time(Ss, grid_c); vi <- maihda_var_at_time(Si, grid_c)
        vr_grid <- maihda_longitudinal_resid_grid_lme4(bm, time_term, grid_c, vs + vi,
                                                       orig_time = lng$time, center = center)
        boot[i, ] <- vs / (vs + vi + vr_grid)
        rs <- maihda_var_at_time(Ss, ref_c); ri <- maihda_var_at_time(Si, ref_c)
        vr_ref <- maihda_longitudinal_resid_grid_lme4(bm, time_term, ref_c, rs + ri,
                                                      orig_time = lng$time, center = center)
        ref_boot[i] <- rs / (rs + ri + vr_ref)
        if (maihda_lme4_optimizer_failed(bm)) n_nonconv <- n_nonconv + 1L
      }, error = function(e) NULL)
    }
    # Non-converged draws are retained but reported (see bootstrap_pcv()).
    if (n_nonconv > 0) {
      warning(sprintf(paste0(
        "%d of %d contributing longitudinal VPC bootstrap draw(s) had an lme4 ",
        "optimiser that did not converge; they are retained in the intervals. ",
        "n_boot_ok counts converged and non-converged refits alike -- interpret the ",
        "intervals accordingly."), n_nonconv, sum(is.finite(ref_boot))), call. = FALSE)
    }
    for (j in seq_along(grid)) {
      band <- maihda_longitudinal_vpc_band(boot[, j], n_boot, conf_level)
      vpc_lower[j] <- band[1]
      vpc_upper[j] <- band[2]
    }
    ref_ci <- maihda_bootstrap_ci(ref_boot, n_boot, conf_level, "VPC")
  }

  vpc_result <- if (bootstrap && !is.null(ref_ci)) {
    list(estimate = ref_vpc, ci_lower = ref_ci[1], ci_upper = ref_ci[2],
         conf_level = conf_level, bootstrap = TRUE, method = "bootstrap",
         ref_time = ref_time, n_boot_ok = attr(ref_ci, "n_ok"),
         n_boot_nonconverged = n_nonconv,
         mc_se = attr(ref_ci, "mc_se"))
  } else {
    list(estimate = ref_vpc, bootstrap = FALSE, ref_time = ref_time)
  }

  longitudinal <- list(
    vpc_t = data.frame(time = grid, estimate = vpc_t_est,
                       lower = vpc_lower, upper = vpc_upper),
    var_stratum_t = var_s_grid,
    var_id_t = var_i_grid,
    var_resid = var_resid,
    var_resid_t = as.numeric(resid_grid),
    Sigma_stratum = Sigma_s,
    Sigma_id = Sigma_i,
    time_grid = grid,
    ref_time = ref_time,
    time = lng$time,
    time_term = time_term,
    time_center = center,
    time_degree = lng$time_degree,
    id = lng$id,
    bootstrap = isTRUE(bootstrap),
    conf_level = conf_level
  )

  list(
    vpc_result = vpc_result,
    variance_components = maihda_longitudinal_components_table(
      Sigma_s, Sigma_i, var_resid, lng$time, lng$id, center = center),
    longitudinal = longitudinal
  )
}

#' Time-varying VPC summary for a longitudinal MAIHDA (brms, linear growth)
#'
#' @param object A longitudinal \code{maihda_model} (brms engine, time_degree 1).
#' @param conf_level Credible-interval level.
#' @return As \code{maihda_longitudinal_summary_lme4}, with posterior bands.
#' @keywords internal
maihda_longitudinal_summary_brms <- function(object, conf_level = 0.95) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to summarize brms models.", call. = FALSE)
  }
  lng <- object$longitudinal_info
  model <- object$model
  time_term <- maihda_lng_time_term(lng)
  center <- maihda_lng_time_center(lng)
  draws <- maihda_posterior_draws_brms(model)
  sig_s <- maihda_re_cov_draws_brms(draws, "stratum", time_term)
  sig_i <- maihda_re_cov_draws_brms(draws, lng$id, time_term)

  # Posterior-mean covariance blocks (components table, and the log-normal mean
  # correction v(t)/2 in the count residual below).
  Sigma_s <- maihda_re_block_brms(model, "stratum", time_term, lng$time_degree)
  Sigma_i <- maihda_re_block_brms(model, lng$id, time_term, lng$time_degree)

  # Per-draw VarS(t) / VarI(t) for the linear growth block:
  # a(t)' Sigma a(t) = v0 + 2 t cov + t^2 v1, t on the coefficient scale. Defined
  # ahead of the residual so the count residual reuses the SAME per-draw
  # random-effect variance for its log-normal mean correction.
  var_at <- function(blk, t_c) blk$v0 + 2 * t_c * blk$cov + t_c^2 * blk$v1

  # Level-1 residual as a function of time. For count families it tracks the
  # time-varying marginal mean instead of one own-time average reused at every time
  # (which would bias VPC(t)); for other families it is the per-draw
  # var_resid_draws, constant over time.
  fam <- maihda_family(model)
  is_count <- !is.null(fam) && identical(fam$link, "log") &&
    fam$family %in% c("poisson", "negbinomial")
  resid_at <- if (is_count) {
    # Count residual computed DRAW-BY-DRAW rather than from posterior-mean
    # predictors: the per-draw fixed-effect eta at each time and the per-draw total
    # random-effect variance both feed the observation-level variance
    # E_i[log1p(1 / lambda_i(t))], so its posterior uncertainty propagates into the
    # VPC band. This matches the cross-sectional count VPC
    # (maihda_count_resid_var_from_linpred(), reused here) and no longer understates
    # the band for low or strongly time-varying counts. Falls back to the documented
    # posterior-mean plug-in only if the posterior_linpred draw axis cannot be
    # aligned with the SD draws.
    nb_extra <- NULL
    if (identical(fam$family, "negbinomial")) {
      if (!"shape" %in% names(draws)) {
        stop("Could not extract the negative-binomial 'shape' (theta) draws from ",
             "the brms posterior.", call. = FALSE)
      }
      nb_extra <- 1 / as.numeric(draws[["shape"]])
    }
    resid_w <- maihda_fit_prior_weights(model)
    n_draws <- length(sig_s$v0)
    function(t_c) {
      v_t_draws <- var_at(sig_s, t_c) + var_at(sig_i, t_c)
      nd <- maihda_longitudinal_set_time(model$data, time_term, t_c,
                                         orig_time = lng$time, center = center)
      eta_link <- tryCatch(
        as.matrix(brms::posterior_linpred(model, newdata = nd, re_formula = NA)),
        error = function(e) NULL)
      if (!is.null(eta_link) && length(dim(eta_link)) == 2L &&
          nrow(eta_link) == n_draws && length(v_t_draws) == n_draws) {
        maihda_count_resid_var_from_linpred(eta_link, v_t_draws, w = resid_w,
                                            extra = nb_extra)
      } else {
        v_t_mean <- maihda_var_at_time(Sigma_s, t_c) +
          maihda_var_at_time(Sigma_i, t_c)
        maihda_longitudinal_resid_at_brms(model, draws, time_term, t_c, v_t_mean,
                                          orig_time = lng$time, center = center)
      }
    }
  } else {
    var_resid_draws <- maihda_residual_variance_draws_brms(model, draws)
    if (length(var_resid_draws) == 1L) {
      var_resid_draws <- rep(var_resid_draws, length(sig_s$v0))
    }
    function(t_c) var_resid_draws
  }

  grid <- maihda_longitudinal_time_grid(object$data[[lng$time]])
  ref_time <- lng$ref_time
  # As in the lme4 sibling: reporting times are on the original scale, the
  # posterior blocks are in (possibly centered) coefficient coordinates.
  grid_c <- grid - center
  ref_c <- ref_time - center
  a <- 1 - conf_level

  # Level-1 residual evaluated once per reporting time (a scalar or a per-draw
  # vector), reused for both the VPC bands and the components-table scalar.
  resid_grid_draws <- lapply(grid_c, resid_at)
  resid_ref_draws <- resid_at(ref_c)

  vpc_draws_at <- function(t_c, resid) {
    vs <- var_at(sig_s, t_c); vi <- var_at(sig_i, t_c)
    vs / (vs + vi + resid)
  }

  summ <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA_real_, NA_real_, NA_real_))
    c(stats::median(v), stats::quantile(v, a / 2, names = FALSE),
      stats::quantile(v, 1 - a / 2, names = FALSE))
  }
  mat <- vapply(seq_along(grid_c),
                function(j) summ(vpc_draws_at(grid_c[j], resid_grid_draws[[j]])),
                numeric(3))
  ref <- summ(vpc_draws_at(ref_c, resid_ref_draws))

  vpc_result <- list(estimate = ref[1], ci_lower = ref[2], ci_upper = ref[3],
                     conf_level = conf_level, bootstrap = FALSE,
                     method = "posterior", ref_time = ref_time)

  # Scalar residual for the components table: the reference-time value (mean over
  # draws), matching mean(var_resid_draws) for the non-count families.
  var_resid <- mean(resid_ref_draws)
  resid_grid <- vapply(resid_grid_draws, mean, numeric(1))

  longitudinal <- list(
    vpc_t = data.frame(time = grid, estimate = mat[1, ],
                       lower = mat[2, ], upper = mat[3, ]),
    var_stratum_t = maihda_var_at_time(Sigma_s, grid_c),
    var_id_t = maihda_var_at_time(Sigma_i, grid_c),
    var_resid = var_resid,
    var_resid_t = resid_grid,
    Sigma_stratum = Sigma_s,
    Sigma_id = Sigma_i,
    time_grid = grid,
    ref_time = ref_time,
    time = lng$time,
    time_term = time_term,
    time_center = center,
    time_degree = lng$time_degree,
    id = lng$id,
    bootstrap = FALSE,
    conf_level = conf_level
  )

  list(
    vpc_result = vpc_result,
    variance_components = maihda_longitudinal_components_table(
      Sigma_s, Sigma_i, var_resid, lng$time, lng$id, center = center),
    longitudinal = longitudinal
  )
}

# Per-draw 2x2 covariance pieces (v0 = intercept var, v1 = slope var,
# cov = intercept-slope covariance) of a group's linear growth block (brms).
maihda_re_cov_draws_brms <- function(draws, group, time) {
  sd0 <- draws[[paste0("sd_", group, "__Intercept")]]
  sd1 <- draws[[paste0("sd_", group, "__", time)]]
  cor01 <- draws[[paste0("cor_", group, "__Intercept__", time)]]
  if (is.null(sd0) || is.null(sd1) || is.null(cor01)) {
    stop("Could not find the intercept/slope SD and correlation draws for the '",
         group, "' random effect in the brms posterior (expected sd_", group,
         "__Intercept, sd_", group, "__", time, ", cor_", group, "__Intercept__",
         time, ").", call. = FALSE)
  }
  sd0 <- as.numeric(sd0); sd1 <- as.numeric(sd1); cor01 <- as.numeric(cor01)
  list(v0 = sd0^2, v1 = sd1^2, cov = cor01 * sd0 * sd1)
}

# Variance-components table for a longitudinal summary. The variances are
# time-varying, so this lists the covariance-block pieces (intercept var, slope
# var, intercept-slope covariance) for the stratum and individual levels plus the
# residual -- it is NOT a single proportion stack (use the VPC trajectory for the
# share over time). Tagged kind = "longitudinal" so plot_vpc()/print route around
# the proportion-stack logic.
maihda_longitudinal_components_table <- function(Sigma_s, Sigma_i, var_resid,
                                                 time, id, center = 0) {
  block_rows <- function(Sigma, level) {
    deg <- nrow(Sigma) - 1L
    # Diagonal (variance) rows. The first is the random-INTERCEPT variance -- the
    # between-level variance at the model's coefficient origin, labelled with
    # that origin: raw time 0 for a zero-anchored fit, or the centering offset
    # (the observed baseline) when the growth terms were fit on internally
    # centered time. The baseline variance the VPC summary and the longitudinal
    # PCV report is a(t)'Sigma a(t) evaluated at ref_time (see
    # maihda_var_at_time()); under centering the two coincide whenever
    # ref_time equals the centering offset (the usual case).
    diag_names <- c(sprintf("intercept (time = %g)", center),
                    if (deg >= 1) paste0("slope (", time, ")"),
                    if (deg >= 2) paste0("slope^", 2:deg, " (", time, ")"))
    vars <- diag(Sigma)
    out <- data.frame(component = sprintf("%s: %s", level, diag_names),
                      variance = as.numeric(vars),
                      sd = sqrt(pmax(as.numeric(vars), 0)),
                      stringsAsFactors = FALSE)
    # Off-diagonal covariance rows: EVERY unique pair (i < j) of the block, not
    # just intercept-slope. The time-varying variance a(t)'Sigma a(t) behind the
    # VPC uses the whole matrix, so a quadratic (3x3) block also carries the
    # intercept-quadratic and slope-quadratic covariances. (For the linear 2x2
    # block this reduces to the single intercept-slope covariance as before.)
    if (deg >= 1) {
      short <- c("intercept", "slope", if (deg >= 2) paste0("slope^", 2:deg))
      pairs <- utils::combn(deg + 1L, 2L)
      out <- rbind(out, data.frame(
        component = sprintf("%s: %s-%s covariance", level,
                            short[pairs[1, ]], short[pairs[2, ]]),
        variance = as.numeric(Sigma[t(pairs)]), sd = NA_real_,
        stringsAsFactors = FALSE))
    }
    out
  }
  tab <- rbind(
    block_rows(Sigma_s, "Between-stratum"),
    block_rows(Sigma_i, sprintf("Between-individual (%s)", id)),
    data.frame(component = "Within (residual)", variance = var_resid,
               sd = sqrt(max(var_resid, 0)), stringsAsFactors = FALSE)
  )
  attr(tab, "kind") <- "longitudinal"
  tab
}

# ---- proportional change in variance (additive vs multiplicative) -----------

# REML lmer between-stratum (co)variance estimates are not comparable across
# models with different fixed effects -- exactly the longitudinal null-vs-adjusted
# pair, where the adjusted growth model adds the dimensions' main effects and
# their dim:time interactions. This applies maihda_pcv_refit_ml()'s rule (see
# calculate_pcv.R) to the growth models that helper deliberately skips: refit a
# REML lmer growth fit with ML (lme4::refitML) before maihda_longitudinal_pcv()
# reads its stratum covariance block. Left on REML, the comparison biases both
# PCVs downward (overstating the multiplicative/interaction share) and can flip
# their sign outright with few strata. glmer (GLMM) fits and the brms engine are
# already on the ML / posterior scale and are returned unchanged, as is a fit
# whose stratum growth block sits on the boundary (see below). The single-model
# summaries (the time-varying VPC, the components table) deliberately keep their
# REML fit, exactly like the cross-sectional VPC/ICC summaries.
maihda_longitudinal_refit_ml <- function(model) {
  if (!inherits(model, "maihda_model") || !identical(model$engine, "lme4") ||
      is.null(model$longitudinal_info)) {
    return(model)
  }
  is_reml <- tryCatch(isTRUE(lme4::isREML(model$model)), error = function(e) FALSE)
  if (!is_reml) return(model)
  # Skip the ML refit when the stratum growth block is itself on the boundary
  # (every variance effectively zero): a(t)' Sigma a(t) is ~0 under either
  # criterion, the PCV degrades to NA (null at boundary) or ~1 (adjusted at
  # boundary) either way, and re-optimising would only nudge an exact-zero block
  # off the boundary. Mirrors maihda_pcv_refit_ml()'s boundary guard. If
  # refitML() itself fails, the tryCatch below keeps the original (REML) fit.
  if (maihda_stratum_growth_at_boundary_lme4(model$model)) return(model)
  refit <- tryCatch(lme4::refitML(model$model), error = function(e) NULL)
  if (!is.null(refit)) {
    model$model <- refit
  } else {
    # The ML refit was warranted (a non-boundary REML growth fit) but failed: keep the
    # REML fit and SAY SO, rather than silently leaving a REML block the comparison would
    # otherwise present as ML. maihda_longitudinal_pcv() then sees the fit is still REML
    # and sets ml_refit = FALSE, so its print method makes no ML claim -- but without this
    # warning the user is never told the requested ML basis could not be applied.
    warning("estimation = \"ML\" was requested, but lme4::refitML() failed for a ",
            "longitudinal growth model; its REML fit is used instead, so this is not a ",
            "pure ML-basis comparison. Compare against estimation = \"fitted\".",
            call. = FALSE)
  }
  model
}

# TRUE when EVERY variance of the stratum growth block sits on the lower boundary
# (effectively zero), reproducing lme4::isSingular()'s relative tolerance per
# component: each scaled SD sqrt(var)/sigma_resid is below `tol`. Zero diagonals
# bound the covariances (|cov| <= sqrt(v_i v_j)), so the whole block -- and hence
# a(t)' Sigma a(t) at every t -- is ~0. A block with ANY non-negligible variance
# is worth the ML refit. The growth analogue of
# maihda_stratum_at_boundary_lme4() (calculate_pcv.R), which reads only the
# single intercept cell. Returns TRUE (skip the refit) if the block cannot be
# read at all; maihda_longitudinal_pcv()'s own NA guards then surface the problem.
maihda_stratum_growth_at_boundary_lme4 <- function(model, tol = 1e-4) {
  d <- tryCatch({
    vc <- lme4::VarCorr(model)
    if (!"stratum" %in% names(vc)) {
      NA_real_
    } else {
      as.numeric(diag(as.matrix(vc[["stratum"]])))
    }
  }, error = function(e) NA_real_)
  if (length(d) == 0 || anyNA(d)) return(TRUE)
  sigma_resid <- tryCatch(stats::sigma(model), error = function(e) NA_real_)
  if (!is.finite(sigma_resid) || sigma_resid <= 0) {
    # No usable residual scale: fall back to an absolute near-zero test.
    return(all(d <= .Machine$double.eps))
  }
  all(sqrt(pmax(d, 0)) / sigma_resid < tol)
}

#' Longitudinal MAIHDA proportional change in variance (PCV)
#'
#' Compares the stratum-level random-effect covariance block of the null growth
#' model with that of the adjusted model (null + dimension main effects + their
#' \code{dim:time} interactions). Reports the PCV in the baseline variance and in
#' the instantaneous-slope variance at baseline -- the additive-vs-multiplicative
#' split of the intersectional trajectory inequality (Bell, Evans, Holman &
#' Leckie 2024) -- and the time-specific PCV over the supplied times. Both are
#' evaluated at the observed baseline time (\code{ref_time}), so they are
#' invariant to how the time axis is coded (for linear growth the slope variance
#' is the same at every time, so this reduces to the slope-variance cell).
#'
#' As in \code{\link{calculate_pcv}}, the \code{estimation} argument selects the
#' variance-estimation basis. With \code{estimation = "ML"}, REML \code{lmer} growth
#' fits are refitted with maximum likelihood (\code{\link[lme4]{refitML}}) before the
#' comparison -- the null and adjusted models differ in fixed effects (the dimensions'
#' main effects and their \code{dim:time} interactions), across which REML applies a
#' model-specific correction -- for a correction-free comparison. With
#' \code{estimation = "fitted"} (the default) each fit's own REML covariance block is
#' used, matching the single-model summaries and avoiding ML's finite-sample bias. The
#' stored models (and the single-model summaries computed from them, e.g. the
#' time-varying VPC) always keep their REML fit; \code{ml_refit} on the result records
#' whether an ML refit fully applied, and \code{estimation_used} records the basis
#' ACTUALLY used -- \code{"fitted"}, \code{"ML"}, or \code{"mixed"} when an ML refit
#' was requested but a model kept its REML fit (a boundary skip or a failed refit),
#' so a mixed REML/ML comparison is not mistaken for a clean one. glmer (GLMM) and
#' brms fits are already on the ML / posterior scale, so the choice does not affect
#' them.
#'
#' @param null_model,adjusted_model Longitudinal \code{maihda_model}s from a
#'   \code{maihda(decomposition = "longitudinal")} pair.
#' @param times Optional numeric times for the time-specific PCV; defaults to the
#'   null model's reporting grid.
#' @param estimation Variance-estimation basis, \code{"fitted"} (default) or
#'   \code{"ML"}; see \code{\link{calculate_pcv}}.
#' @return An object of class \code{maihda_long_pcv}.
#' @keywords internal
maihda_longitudinal_pcv <- function(null_model, adjusted_model, times = NULL,
                                    estimation = c("fitted", "ML")) {
  estimation <- match.arg(estimation)
  lng <- null_model$longitudinal_info
  center <- maihda_lng_time_center(lng)
  adj_center <- maihda_lng_time_center(adjusted_model$longitudinal_info)
  if (!isTRUE(all.equal(center, adj_center))) {
    stop("The null and adjusted longitudinal models were fit with different ",
         "internal time centerings (", center, " vs ", adj_center, "), so their ",
         "covariance blocks are in different coordinates and cannot be compared.",
         call. = FALSE)
  }

  # REML vs ML basis (see calculate_pcv()'s `estimation` argument). With
  # estimation = "ML" the two growth models are refit with ML before their
  # between-stratum covariance blocks are read (the fits differ in fixed effects, so
  # REML applies a model-specific correction); with "fitted" (default) each fit's own
  # REML block is used. Either way the caller's stored models are untouched (copy
  # semantics), so summary()'s time-varying VPC keeps each fit's own REML estimate.
  # ml_refit records whether the comparison is on the ML scale -- FALSE for "fitted",
  # and also when a boundary skip or a failed refitML() left a REML fit in place (the
  # print method then makes no ML claim).
  is_reml_lng <- function(m) {
    identical(m$engine, "lme4") &&
      tryCatch(isTRUE(lme4::isREML(m$model)), error = function(e) FALSE)
  }
  do_ml <- identical(estimation, "ML")
  null_was_reml <- is_reml_lng(null_model)
  adj_was_reml  <- is_reml_lng(adjusted_model)
  refit_needed <- do_ml && (null_was_reml || adj_was_reml)
  if (do_ml) {
    null_model <- maihda_longitudinal_refit_ml(null_model)
    adjusted_model <- maihda_longitudinal_refit_ml(adjusted_model)
  }
  # A model we WANTED on ML that is STILL REML after the refit -- a boundary skip
  # (maihda_longitudinal_refit_ml()'s guard) or a failed refitML() -- makes the
  # comparison MIXED: one covariance block on ML, the other on REML, which is
  # exactly the cross-model REML incomparability estimation = "ML" was meant to
  # remove. `ml_refit` alone cannot record this (it is FALSE for a plain fitted/REML
  # request too), so it silently read as a clean fitted comparison. Track the basis
  # ACTUALLY used ("fitted"/"ML"/"mixed") on the result, mirroring calculate_pcv().
  still_reml <- is_reml_lng(null_model) || is_reml_lng(adjusted_model)
  ml_incomplete <- refit_needed && still_reml
  # ml_refit stays TRUE only for a fully-applied ML refit (kept for backward
  # compatibility with objects/tests that read it); estimation_used carries the
  # three-way basis and is what print() keys on.
  ml_refit <- refit_needed && !still_reml
  estimation_used <- maihda_pcv_estimation_used(estimation, ml_incomplete,
                                                engine = null_model$engine)

  Sn <- maihda_re_block(null_model, "stratum")
  Sa <- maihda_re_block(adjusted_model, "stratum")
  deg <- nrow(Sn) - 1L

  # A between-stratum variance a(t)' Sn a(t) that is only boundary-level positive --
  # e.g. 2.5e-13 from a singular null growth fit with no stratum trajectory signal --
  # passes a plain > 0 test yet makes the PCV ratio explode (PCV(t) in the hundreds
  # or thousands, with the wrong sign). Gate every denominator with the SAME relative
  # singularity tolerance calculate_pcv() / stepwise_pcv() apply
  # (maihda_variance_at_boundary(): sqrt(var) / sqrt(resid) < 1e-4), scaling the null
  # between-stratum variance against the null model's own residual/latent variance; a
  # boundary-level denominator yields NA rather than a spurious ratio. For lme4 this is
  # sigma^2 (Gaussian) or the latent-scale constant; where it is unavailable (e.g. brms)
  # the helper falls back to an absolute near-zero test.
  resid_var_null <- tryCatch(
    if (identical(null_model$engine, "lme4"))
      maihda_residual_variance_lme4(null_model$model) else NA_real_,
    error = function(e) NA_real_)
  denom_at_boundary <- function(v) maihda_variance_at_boundary(v, resid_var_null)
  # TRUE when the WHOLE null stratum growth block is degenerate (every variance on the
  # boundary): a(t)' Sn a(t) is then ~0 at every t, so no PCV cell is defined and all
  # are NA below. Recorded on the result (and surfaced by print()) so the all-NA output
  # reads as an honest "no between-stratum trajectory variance to explain", not a silent
  # gap -- the longitudinal analogue of calculate_pcv()'s degenerate-null guard.
  null_at_boundary <- identical(null_model$engine, "lme4") &&
    isTRUE(tryCatch(maihda_stratum_growth_at_boundary_lme4(null_model$model),
                    error = function(e) FALSE))

  pcv_cell <- function(vn, va) if (!denom_at_boundary(vn)) (vn - va) / vn else NA_real_

  # The "baseline" PCV is the proportional change in the between-stratum variance
  # at the OBSERVED baseline time (lng$ref_time = min(time)), not the raw
  # intercept-variance cell Sn[1, 1]. The covariance blocks are in the model's
  # (possibly centered) coefficient coordinates, so a(t)'Sigma a(t) is evaluated
  # at ref_time - center; under the internal centering (center = min(time)) this
  # is usually 0, making the baseline the intercept cell itself. Matches how the
  # VPC summary reports its baseline.
  ref_time <- lng$ref_time
  ref_c <- ref_time - center
  var_baseline_null <- maihda_var_at_time(Sn, ref_c)
  var_baseline_adjusted <- maihda_var_at_time(Sa, ref_c)
  pcv_intercept <- pcv_cell(var_baseline_null, var_baseline_adjusted)
  # The slope PCV compares the INSTANTANEOUS-slope variance at the baseline,
  # b(t)'Sigma b(t) with b = da/dt -- not the raw Sn[2, 2] cell, which is the
  # slope variance at the coefficient origin and is NOT invariant to where time
  # is zeroed once time_degree >= 2 (for linear growth the two coincide at every
  # t, so this reduces to Sn[2, 2] there).
  var_slope_null <- if (deg >= 1) maihda_slope_var_at_time(Sn, ref_c) else NA_real_
  var_slope_adjusted <- if (deg >= 1) maihda_slope_var_at_time(Sa, ref_c) else NA_real_
  pcv_slope <- if (deg >= 1) pcv_cell(var_slope_null, var_slope_adjusted) else NA_real_

  if (is.null(times)) {
    times <- maihda_longitudinal_time_grid(null_model$data[[lng$time]])
  }
  vn_t <- maihda_var_at_time(Sn, times - center)
  va_t <- maihda_var_at_time(Sa, times - center)
  pcv_t <- ifelse(vapply(vn_t, denom_at_boundary, logical(1)), NA_real_,
                  (vn_t - va_t) / vn_t)

  structure(
    list(
      pcv_intercept = pcv_intercept,
      pcv_slope = pcv_slope,
      var_baseline_null = var_baseline_null,
      var_baseline_adjusted = var_baseline_adjusted,
      var_slope_null = var_slope_null,
      var_slope_adjusted = var_slope_adjusted,
      ref_time = ref_time,
      time_center = center,
      pcv_t = data.frame(time = times, var_null = vn_t, var_adjusted = va_t,
                         pcv = pcv_t),
      Sigma_stratum_null = Sn,
      Sigma_stratum_adjusted = Sa,
      time = lng$time,
      time_degree = lng$time_degree,
      ml_refit = ml_refit,
      estimation = estimation,
      estimation_used = estimation_used,
      null_at_boundary = null_at_boundary
    ),
    class = "maihda_long_pcv"
  )
}

#' Per-stratum trajectory parameters for a longitudinal MAIHDA
#'
#' The stratum-level random-effect estimates as a wide table, one row per stratum:
#' the stratum's deviation at the baseline time (\code{baseline}, the longitudinal
#' analogue of a cross-sectional stratum BLUP), the random intercept
#' (\code{intercept}; the deviation at the model's coefficient origin -- the
#' internal centering offset, i.e. the observed baseline, when centering applied,
#' or raw time 0 otherwise) and the random slope(s) on time (\code{slope}, ...).
#' This is the longitudinal shape of \code{predict_maihda(type = "strata")} -- a
#' stratum is now a \emph{trajectory}, not a single value.
#'
#' \code{baseline} is \eqn{a(t_0 - c)' coef} with \eqn{a(t) = (1, t, t^2, ...)},
#' \eqn{t_0 = } the reference (baseline) time \code{ref_time = min(time)} and
#' \eqn{c} the internal centering offset; it equals \code{intercept} whenever
#' \code{ref_time} coincides with the centering origin (the usual case for a
#' centered fit, and the zero-anchored case for a raw fit).
#'
#' @param object A longitudinal \code{maihda_model}.
#' @return A data frame: \code{stratum}, \code{stratum_id}, optional \code{label},
#'   \code{baseline}, \code{intercept}, \code{slope}(, \code{slope2}, ...).
#' @keywords internal
maihda_longitudinal_strata_predictions <- function(object) {
  re <- maihda_longitudinal_stratum_re(object)
  lng <- object$longitudinal_info
  deg <- lng$time_degree
  # The coefficients are in the model's (possibly centered) coordinates, so the
  # baseline deviation is evaluated at ref_time - center.
  ref_c <- lng$ref_time - maihda_lng_time_center(lng)
  mat <- do.call(rbind, lapply(re$coef, function(co) {
    out <- rep(NA_real_, deg + 1L)
    out[seq_along(co)] <- co
    out
  }))
  colnames(mat) <- c("intercept",
                     if (deg >= 1) "slope",
                     if (deg >= 2) paste0("slope", 2:deg))
  # Deviation at the baseline time, a(ref_time - center)' coef.
  baseline <- as.numeric(mat %*% ref_c^(0:deg))
  df <- data.frame(stratum = re$stratum, stratum_id = re$stratum_id,
                   stringsAsFactors = FALSE)
  if (!is.null(re$label)) df$label <- re$label
  df$baseline <- baseline
  cbind(df, as.data.frame(mat, stringsAsFactors = FALSE))
}

#' Print a longitudinal MAIHDA PCV
#'
#' @param x A \code{maihda_long_pcv} object.
#' @param ... Unused.
#' @return The object, invisibly.
#' @export
print.maihda_long_pcv <- function(x, ...) {
  pal <- maihda_palette()
  fmt <- function(v) if (isTRUE(is.finite(v))) pal$accent(sprintf("%.1f%%", 100 * v)) else "NA"
  cat(pal$bold("Longitudinal PCV (additive vs. multiplicative intersectionality)"), "\n", sep = "")
  cat("================================================================\n\n")
  # Baseline = the between-stratum variance at the observed baseline time
  # (ref_time), evaluated in the model's (possibly centered) coefficient
  # coordinates (see maihda_longitudinal_pcv).
  cat(sprintf("Baseline (%s = %g) variance: %.4f (null) -> %.4f (adjusted)\n",
              x$time, x$ref_time, x$var_baseline_null, x$var_baseline_adjusted))
  cat(sprintf("  PCV_intercept: %s of the baseline between-stratum inequality is additive.\n",
              fmt(x$pcv_intercept)))
  if (nrow(x$Sigma_stratum_null) >= 2) {
    # Instantaneous-slope variance at the baseline (see maihda_slope_var_at_time);
    # objects stored by older package versions carry no var_slope_* fields, so
    # fall back to the slope-variance cell they reported.
    vsn <- if (!is.null(x$var_slope_null)) x$var_slope_null else x$Sigma_stratum_null[2, 2]
    vsa <- if (!is.null(x$var_slope_adjusted)) x$var_slope_adjusted else x$Sigma_stratum_adjusted[2, 2]
    cat(sprintf("Slope (%s) variance at baseline: %.4f (null) -> %.4f (adjusted)\n",
                x$time, vsn, vsa))
    cat(sprintf("  PCV_slope:     %s of the *trajectory* between-stratum inequality is additive\n",
                fmt(x$pcv_slope)))
    cat("                 (the remainder is the multiplicative/interaction part).\n")
  }
  cat(pal$muted(paste0(
      "\nThe PCV is the share of the null model's between-stratum (trajectory) variance\n",
      "explained by the dimensions' additive main effects and their time interactions;\n",
      "a high PCV_slope means trajectory inequalities are 'mostly additive'.\n")))
  if (isTRUE(x$null_at_boundary)) {
    # The null growth model's between-stratum block sits on the singularity boundary
    # (effectively no stratum trajectory variance), so every PCV is a 0/0 ratio and is
    # reported as NA rather than an exploded value from a near-zero denominator. Say so,
    # mirroring calculate_pcv()'s boundary caveat.
    cat(pal$warn(paste0(
        "\nNote: the null model's between-stratum (trajectory) variance is at the singularity\n",
        "boundary (a singular fit -- effectively no between-stratum trajectory variation to\n",
        "explain), so the PCV is undefined (0/0) and reported as NA. This is consistent with\n",
        "genuinely additive strata as well as a degenerate fit; the two cannot be told apart.\n")))
  }
  if (identical(x$estimation_used, "mixed")) {
    # estimation = "ML" was requested, but a growth model kept its REML fit (a
    # between-stratum trajectory variance at the singularity boundary, or a failed
    # refitML()), so the comparison mixes an ML block with a REML one -- NOT the
    # pure, correction-free ML basis requested. Surface it, mirroring
    # calculate_pcv()'s "mixed" basis note, so the output is not mistaken for a
    # clean fitted (REML) comparison.
    cat(pal$warn(paste0(
        "\nNote: estimation = \"ML\" was requested, but a growth model kept its REML fit\n",
        "(a between-stratum trajectory variance at the singularity boundary, or a failed\n",
        "refitML()), so this compares an ML variance block against a REML one -- not the\n",
        "pure, correction-free ML basis requested. Compare against estimation = \"fitted\".\n",
        "See ?calculate_pcv.\n")))
  } else if (isTRUE(x$ml_refit)) {
    # REML growth fits were ML-refitted for this comparison (see
    # maihda_longitudinal_pcv), so the variances above are on the ML scale while
    # summary()'s time-varying VPC keeps each fit's own REML estimate -- say so,
    # mirroring maihda_table()'s basis note for the cross-sectional PCV.
    cat(pal$muted(paste0(
        "(REML growth fits were refitted with maximum likelihood for this comparison --\n",
        "REML variances are not comparable across the null vs. adjusted fixed effects --\n",
        "so these variances are on the ML scale; summary()'s time-varying VPC keeps\n",
        "each fit's own REML estimate. See ?calculate_pcv.)\n")))
  }
  invisible(x)
}
