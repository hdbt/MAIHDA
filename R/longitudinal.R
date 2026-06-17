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
#' @return A list \code{list(id, time, time_degree)}.
#' @keywords internal
maihda_validate_longitudinal <- function(id, time, time_degree, data,
                                         engine = "lme4",
                                         sampling_weights = NULL,
                                         context = NULL) {
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

  list(id = id, time = time, time_degree = time_degree)
}

# Polynomial-in-time term labels, e.g. c("wave", "I(wave^2)") for degree 2. The
# first random/fixed term is the linear time; higher degrees use I(time^k) so the
# design vector at time t is a(t) = (1, t, t^2, ..., t^degree).
maihda_time_terms <- function(time, time_degree) {
  terms <- maihda_quote_name(time)
  if (time_degree >= 2) {
    terms <- c(terms, sprintf("I(%s^%d)", time, 2:time_degree))
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
#' @param base_formula The resolved formula (fixed part + stratum grouping).
#' @param id,time,time_degree The longitudinal specification.
#' @return The growth formula (same environment as \code{base_formula}).
#' @keywords internal
#' @importFrom stats update as.formula terms
maihda_longitudinal_formula <- function(base_formula, id, time, time_degree) {
  poly_terms <- maihda_time_terms(time, time_degree)
  ptime <- paste(poly_terms, collapse = " + ")

  fixed <- reformulas::nobars(base_formula)
  fixed_labels <- attr(stats::terms(fixed), "term.labels")
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
#' @param time,time_degree The longitudinal specification.
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
maihda_re_block <- function(object, group) {
  lng <- object$longitudinal_info
  if (identical(object$engine, "lme4")) {
    maihda_re_block_lme4(object$model, group, lng$time, lng$time_degree)
  } else if (identical(object$engine, "brms")) {
    maihda_re_block_brms(object$model, group, lng$time, lng$time_degree)
  } else {
    stop("Longitudinal MAIHDA is supported only for lme4/brms.", call. = FALSE)
  }
}

# Between-level variance implied by a covariance block at time(s) t:
# a(t)' Sigma a(t) with a(t) = (1, t, t^2, ...). Vectorised over t.
maihda_var_at_time <- function(Sigma, t) {
  degree <- nrow(Sigma) - 1L
  vapply(t, function(ti) {
    a <- ti^(0:degree)
    as.numeric(crossprod(a, Sigma %*% a))
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
  Sigma_s <- maihda_re_block_lme4(model, "stratum", lng$time, lng$time_degree)
  Sigma_i <- maihda_re_block_lme4(model, lng$id, lng$time, lng$time_degree)
  var_resid <- maihda_residual_variance_lme4(model)

  grid <- maihda_longitudinal_time_grid(object$data[[lng$time]])
  ref_time <- lng$ref_time

  vpc_fun <- function(Ss, Si, t) {
    vs <- maihda_var_at_time(Ss, t)
    vi <- maihda_var_at_time(Si, t)
    vs / (vs + vi + var_resid)
  }
  vpc_t_est <- vpc_fun(Sigma_s, Sigma_i, grid)
  ref_vpc <- vpc_fun(Sigma_s, Sigma_i, ref_time)

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
        Ss <- maihda_re_block_lme4(bm, "stratum", lng$time, lng$time_degree)
        Si <- maihda_re_block_lme4(bm, lng$id, lng$time, lng$time_degree)
        vr <- maihda_residual_variance_lme4(bm)
        vs <- maihda_var_at_time(Ss, grid); vi <- maihda_var_at_time(Si, grid)
        boot[i, ] <- vs / (vs + vi + vr)
        rs <- maihda_var_at_time(Ss, ref_time); ri <- maihda_var_at_time(Si, ref_time)
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
    var_stratum_t = maihda_var_at_time(Sigma_s, grid),
    var_id_t = maihda_var_at_time(Sigma_i, grid),
    var_resid = var_resid,
    Sigma_stratum = Sigma_s,
    Sigma_id = Sigma_i,
    time_grid = grid,
    ref_time = ref_time,
    time = lng$time,
    time_degree = lng$time_degree,
    id = lng$id,
    bootstrap = isTRUE(bootstrap),
    conf_level = conf_level
  )

  list(
    vpc_result = vpc_result,
    variance_components = maihda_longitudinal_components_table(
      Sigma_s, Sigma_i, var_resid, lng$time, lng$id),
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
  draws <- maihda_posterior_draws_brms(model)
  sig_s <- maihda_re_cov_draws_brms(draws, "stratum", lng$time)
  sig_i <- maihda_re_cov_draws_brms(draws, lng$id, lng$time)
  var_resid_draws <- maihda_residual_variance_draws_brms(model, draws)
  if (length(var_resid_draws) == 1L) {
    var_resid_draws <- rep(var_resid_draws, length(sig_s$v0))
  }

  grid <- maihda_longitudinal_time_grid(object$data[[lng$time]])
  ref_time <- lng$ref_time
  a <- 1 - conf_level

  # Per-draw VarS(t) / VarI(t) for the linear block:
  # a(t)'Sigma a(t) = v0 + 2 t cov + t^2 v1.
  var_at <- function(blk, t) blk$v0 + 2 * t * blk$cov + t^2 * blk$v1
  vpc_draws_at <- function(t) {
    vs <- var_at(sig_s, t); vi <- var_at(sig_i, t)
    vs / (vs + vi + var_resid_draws)
  }

  summ <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(c(NA_real_, NA_real_, NA_real_))
    c(stats::median(v), stats::quantile(v, a / 2, names = FALSE),
      stats::quantile(v, 1 - a / 2, names = FALSE))
  }
  mat <- vapply(grid, function(t) summ(vpc_draws_at(t)), numeric(3))
  ref <- summ(vpc_draws_at(ref_time))

  vpc_result <- list(estimate = ref[1], ci_lower = ref[2], ci_upper = ref[3],
                     conf_level = conf_level, bootstrap = FALSE,
                     method = "posterior", ref_time = ref_time)

  Sigma_s <- maihda_re_block_brms(model, "stratum", lng$time, lng$time_degree)
  Sigma_i <- maihda_re_block_brms(model, lng$id, lng$time, lng$time_degree)
  var_resid <- mean(var_resid_draws)

  longitudinal <- list(
    vpc_t = data.frame(time = grid, estimate = mat[1, ],
                       lower = mat[2, ], upper = mat[3, ]),
    var_stratum_t = maihda_var_at_time(Sigma_s, grid),
    var_id_t = maihda_var_at_time(Sigma_i, grid),
    var_resid = var_resid,
    Sigma_stratum = Sigma_s,
    Sigma_id = Sigma_i,
    time_grid = grid,
    ref_time = ref_time,
    time = lng$time,
    time_degree = lng$time_degree,
    id = lng$id,
    bootstrap = FALSE,
    conf_level = conf_level
  )

  list(
    vpc_result = vpc_result,
    variance_components = maihda_longitudinal_components_table(
      Sigma_s, Sigma_i, var_resid, lng$time, lng$id),
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
                                                 time, id) {
  block_rows <- function(Sigma, level) {
    deg <- nrow(Sigma) - 1L
    # Diagonal (variance) rows. The first is the random-INTERCEPT variance, i.e.
    # the between-level variance at time = 0 -- NOT the variance at the baseline
    # (ref_time = min(time)); the two coincide only when time is zero-referenced.
    # The baseline variance is reported by the VPC summary and the longitudinal
    # PCV as a(t)'Sigma a(t) evaluated at ref_time (see maihda_var_at_time()).
    diag_names <- c("intercept (time = 0)",
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
#' \code{dim:time} interactions). Reports the PCV in the baseline (intercept)
#' variance and in the slope variance -- the additive-vs-multiplicative split of
#' the intersectional trajectory inequality (Bell, Evans, Holman & Leckie 2024) --
#' and the time-specific PCV over the supplied times.
#'
#' @param null_model,adjusted_model Longitudinal \code{maihda_model}s from a
#'   \code{maihda(decomposition = "longitudinal")} pair.
#' @param times Optional numeric times for the time-specific PCV; defaults to the
#'   null model's reporting grid.
#' @return An object of class \code{maihda_long_pcv}.
#' @keywords internal
maihda_longitudinal_pcv <- function(null_model, adjusted_model, times = NULL) {
  lng <- null_model$longitudinal_info
  Sn <- maihda_re_block(null_model, "stratum")
  Sa <- maihda_re_block(adjusted_model, "stratum")
  deg <- nrow(Sn) - 1L

  pcv_cell <- function(vn, va) if (is.finite(vn) && vn > 0) (vn - va) / vn else NA_real_

  # The "baseline" PCV is the proportional change in the between-stratum variance
  # at the OBSERVED baseline time (lng$ref_time = min(time)), not the raw-time-0
  # intercept variance Sn[1, 1]. The two coincide only when time is centred so the
  # baseline is 0; for waves like 10:12 the time-0 intercept variance is an
  # extrapolation and its PCV is meaningless (it can even go negative). Evaluating
  # a(t)'Sigma a(t) at ref_time matches how the VPC summary reports its baseline.
  ref_time <- lng$ref_time
  var_baseline_null <- maihda_var_at_time(Sn, ref_time)
  var_baseline_adjusted <- maihda_var_at_time(Sa, ref_time)
  pcv_intercept <- pcv_cell(var_baseline_null, var_baseline_adjusted)
  # The slope variance is invariant to where time is zeroed, so Sn[2, 2]/Sa[2, 2]
  # are already the right cells.
  pcv_slope <- if (deg >= 1) pcv_cell(Sn[2, 2], Sa[2, 2]) else NA_real_

  if (is.null(times)) {
    times <- maihda_longitudinal_time_grid(null_model$data[[lng$time]])
  }
  vn_t <- maihda_var_at_time(Sn, times)
  va_t <- maihda_var_at_time(Sa, times)
  pcv_t <- ifelse(vn_t > 0, (vn_t - va_t) / vn_t, NA_real_)

  structure(
    list(
      pcv_intercept = pcv_intercept,
      pcv_slope = pcv_slope,
      var_baseline_null = var_baseline_null,
      var_baseline_adjusted = var_baseline_adjusted,
      ref_time = ref_time,
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
#' analogue of a cross-sectional stratum BLUP), the raw random intercept at
#' time 0 (\code{intercept}) and the random slope(s) on time (\code{slope}, ...).
#' This is the longitudinal shape of \code{predict_maihda(type = "strata")} -- a
#' stratum is now a \emph{trajectory}, not a single value.
#'
#' \code{baseline} is \eqn{a(t_0)' coef} with \eqn{a(t) = (1, t, t^2, ...)} and
#' \eqn{t_0 = } the reference (baseline) time \code{ref_time = min(time)}; the
#' package defines the baseline at \code{ref_time}, so it equals the raw
#' \code{intercept} (deviation at time 0) only when time is zero-referenced.
#'
#' @param object A longitudinal \code{maihda_model}.
#' @return A data frame: \code{stratum}, \code{stratum_id}, optional \code{label},
#'   \code{baseline}, \code{intercept}, \code{slope}(, \code{slope2}, ...).
#' @keywords internal
maihda_longitudinal_strata_predictions <- function(object) {
  re <- maihda_longitudinal_stratum_re(object)
  deg <- object$longitudinal_info$time_degree
  ref_time <- object$longitudinal_info$ref_time
  mat <- do.call(rbind, lapply(re$coef, function(co) {
    out <- rep(NA_real_, deg + 1L)
    out[seq_along(co)] <- co
    out
  }))
  colnames(mat) <- c("intercept",
                     if (deg >= 1) "slope",
                     if (deg >= 2) paste0("slope", 2:deg))
  # Deviation at the baseline time, a(ref_time)' coef, NOT the raw time-0
  # intercept -- these differ whenever time is not zero-referenced.
  baseline <- as.numeric(mat %*% ref_time^(0:deg))
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
  fmt <- function(v) if (isTRUE(is.finite(v))) sprintf("%.1f%%", 100 * v) else "NA"
  cat("Longitudinal PCV (additive vs. multiplicative intersectionality)\n")
  cat("================================================================\n\n")
  # Baseline = the between-stratum variance at the observed baseline time
  # (ref_time), not the raw time-0 intercept variance Sigma[1, 1] (see
  # maihda_longitudinal_pcv); the two coincide only when time is centred.
  cat(sprintf("Baseline (%s = %g) variance: %.4f (null) -> %.4f (adjusted)\n",
              x$time, x$ref_time, x$var_baseline_null, x$var_baseline_adjusted))
  cat(sprintf("  PCV_intercept: %s of the baseline between-stratum inequality is additive.\n",
              fmt(x$pcv_intercept)))
  if (nrow(x$Sigma_stratum_null) >= 2) {
    cat(sprintf("Slope (%s) variance:          %.4f (null) -> %.4f (adjusted)\n",
                x$time, x$Sigma_stratum_null[2, 2], x$Sigma_stratum_adjusted[2, 2]))
    cat(sprintf("  PCV_slope:     %s of the *trajectory* between-stratum inequality is additive\n",
                fmt(x$pcv_slope)))
    cat("                 (the remainder is the multiplicative/interaction part).\n")
  }
  cat("\nThe PCV is the share of the null model's between-stratum (trajectory) variance\n")
  cat("explained by the dimensions' additive main effects and their time interactions;\n")
  cat("a high PCV_slope means trajectory inequalities are 'mostly additive'.\n")
  invisible(x)
}
