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

  fixed <- reformulas::nobars(base_formula)
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

# ---- time-varying VPC summary -----------------------------------------------

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
  var_resid <- maihda_residual_variance_lme4(model)

  grid <- maihda_longitudinal_time_grid(object$data[[lng$time]])
  ref_time <- lng$ref_time
  # Times are user-facing (original scale); the covariance blocks are in the
  # (possibly centered) coefficient coordinates, so every a(t)' Sigma a(t)
  # evaluation subtracts the centering offset first.
  grid_c <- grid - center
  ref_c <- ref_time - center

  vpc_fun <- function(Ss, Si, t_c) {
    vs <- maihda_var_at_time(Ss, t_c)
    vi <- maihda_var_at_time(Si, t_c)
    vs / (vs + vi + var_resid)
  }
  vpc_t_est <- vpc_fun(Sigma_s, Sigma_i, grid_c)
  ref_vpc <- vpc_fun(Sigma_s, Sigma_i, ref_c)

  vpc_lower <- rep(NA_real_, length(grid))
  vpc_upper <- rep(NA_real_, length(grid))
  ref_ci <- NULL
  if (bootstrap) {
    boot <- matrix(NA_real_, nrow = n_boot, ncol = length(grid))
    ref_boot <- rep(NA_real_, n_boot)
    sim <- stats::simulate(model, nsim = n_boot)
    for (i in seq_len(n_boot)) {
      tryCatch({
        bm <- lme4::refit(model, newresp = sim[[i]])
        Ss <- maihda_re_block_lme4(bm, "stratum", time_term, lng$time_degree)
        Si <- maihda_re_block_lme4(bm, lng$id, time_term, lng$time_degree)
        vr <- maihda_residual_variance_lme4(bm)
        vs <- maihda_var_at_time(Ss, grid_c); vi <- maihda_var_at_time(Si, grid_c)
        boot[i, ] <- vs / (vs + vi + vr)
        rs <- maihda_var_at_time(Ss, ref_c); ri <- maihda_var_at_time(Si, ref_c)
        ref_boot[i] <- rs / (rs + ri + vr)
      }, error = function(e) NULL)
    }
    a <- 1 - conf_level
    for (j in seq_along(grid)) {
      col <- boot[, j][is.finite(boot[, j])]
      if (length(col) >= 10L) {
        vpc_lower[j] <- stats::quantile(col, a / 2, names = FALSE)
        vpc_upper[j] <- stats::quantile(col, 1 - a / 2, names = FALSE)
      }
    }
    ref_ci <- maihda_bootstrap_ci(ref_boot, n_boot, conf_level, "VPC")
  }

  vpc_result <- if (bootstrap && !is.null(ref_ci)) {
    list(estimate = ref_vpc, ci_lower = ref_ci[1], ci_upper = ref_ci[2],
         conf_level = conf_level, bootstrap = TRUE, method = "bootstrap",
         ref_time = ref_time, n_boot_ok = attr(ref_ci, "n_ok"),
         mc_se = attr(ref_ci, "mc_se"))
  } else {
    list(estimate = ref_vpc, bootstrap = FALSE, ref_time = ref_time)
  }

  longitudinal <- list(
    vpc_t = data.frame(time = grid, estimate = vpc_t_est,
                       lower = vpc_lower, upper = vpc_upper),
    var_stratum_t = maihda_var_at_time(Sigma_s, grid_c),
    var_id_t = maihda_var_at_time(Sigma_i, grid_c),
    var_resid = var_resid,
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
  var_resid_draws <- maihda_residual_variance_draws_brms(model, draws)
  if (length(var_resid_draws) == 1L) {
    var_resid_draws <- rep(var_resid_draws, length(sig_s$v0))
  }

  grid <- maihda_longitudinal_time_grid(object$data[[lng$time]])
  ref_time <- lng$ref_time
  # As in the lme4 sibling: reporting times are on the original scale, the
  # posterior blocks are in (possibly centered) coefficient coordinates.
  grid_c <- grid - center
  ref_c <- ref_time - center
  a <- 1 - conf_level

  # Per-draw VarS(t) / VarI(t) for the linear block:
  # a(t)'Sigma a(t) = v0 + 2 t cov + t^2 v1, t on the coefficient scale.
  var_at <- function(blk, t_c) blk$v0 + 2 * t_c * blk$cov + t_c^2 * blk$v1
  vpc_draws_at <- function(t_c) {
    vs <- var_at(sig_s, t_c); vi <- var_at(sig_i, t_c)
    vs / (vs + vi + var_resid_draws)
  }

  summ <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA_real_, NA_real_, NA_real_))
    c(stats::median(v), stats::quantile(v, a / 2, names = FALSE),
      stats::quantile(v, 1 - a / 2, names = FALSE))
  }
  mat <- vapply(grid_c, function(t_c) summ(vpc_draws_at(t_c)), numeric(3))
  ref <- summ(vpc_draws_at(ref_c))

  vpc_result <- list(estimate = ref[1], ci_lower = ref[2], ci_upper = ref[3],
                     conf_level = conf_level, bootstrap = FALSE,
                     method = "posterior", ref_time = ref_time)

  Sigma_s <- maihda_re_block_brms(model, "stratum", time_term, lng$time_degree)
  Sigma_i <- maihda_re_block_brms(model, lng$id, time_term, lng$time_degree)
  var_resid <- mean(var_resid_draws)

  longitudinal <- list(
    vpc_t = data.frame(time = grid, estimate = mat[1, ],
                       lower = mat[2, ], upper = mat[3, ]),
    var_stratum_t = maihda_var_at_time(Sigma_s, grid_c),
    var_id_t = maihda_var_at_time(Sigma_i, grid_c),
    var_resid = var_resid,
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
#' @param null_model,adjusted_model Longitudinal \code{maihda_model}s from a
#'   \code{maihda(decomposition = "longitudinal")} pair.
#' @param times Optional numeric times for the time-specific PCV; defaults to the
#'   null model's reporting grid.
#' @return An object of class \code{maihda_long_pcv}.
#' @keywords internal
maihda_longitudinal_pcv <- function(null_model, adjusted_model, times = NULL) {
  lng <- null_model$longitudinal_info
  center <- maihda_lng_time_center(lng)
  adj_center <- maihda_lng_time_center(adjusted_model$longitudinal_info)
  if (!isTRUE(all.equal(center, adj_center))) {
    stop("The null and adjusted longitudinal models were fit with different ",
         "internal time centerings (", center, " vs ", adj_center, "), so their ",
         "covariance blocks are in different coordinates and cannot be compared.",
         call. = FALSE)
  }
  Sn <- maihda_re_block(null_model, "stratum")
  Sa <- maihda_re_block(adjusted_model, "stratum")
  deg <- nrow(Sn) - 1L

  pcv_cell <- function(vn, va) if (is.finite(vn) && vn > 0) (vn - va) / vn else NA_real_

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
  pcv_t <- ifelse(vn_t > 0, (vn_t - va_t) / vn_t, NA_real_)

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
      time_degree = lng$time_degree
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
  invisible(x)
}
