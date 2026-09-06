# Design-weighted MAIHDA (sampling / survey weights).
#
# lme4's `weights=` are PRECISION weights (they scale the residual variance), not
# sampling weights, so feeding survey weights to lmer/glmer gives the wrong
# objective and invalid standard errors for population-representative estimates.
# The design-weighted MAIHDA of the literature (Evans et al.; following
# Rabe-Hesketh & Skrondal 2006) instead maximises a weighted pseudo-likelihood.
# This file implements that path via the 'wemix' engine (WeMix::mix(), the
# pseudo-maximum-likelihood mixed-model fitter built for NAEP/PISA analysis) and a
# pseudo-posterior path for the brms engine (sampling weights as likelihood
# weights). The intersectional strata are exhaustive population cells -- every
# stratum is "sampled" with probability 1 -- so the level-2 weights are 1 and the
# individual sampling weights enter at level 1 unchanged (conditional and
# unconditional level-1 weights coincide).

# Reserved column names added to the analytic data by the weighted engines.
# These are package-internal: '.maihda_l2wt' is the constant level-2 weight WeMix
# is handed, and '.maihda_sw' is the normalized likelihood-weight column injected
# into the brms formula. They must not collide with a user variable the model
# uses (see maihda_guard_reserved_weight_col()).
.maihda_wemix_l2_col <- ".maihda_l2wt"
.maihda_brms_weights_col <- ".maihda_sw"

# Gradient-criterion threshold below which a WeMix fit counts as having reached a
# maximum; see maihda_wemix_convergence() for how it is calibrated.
.maihda_wemix_grad_tol <- 1e-2

# Guard a reserved internal weight column against silently overwriting a user
# variable. The weighted engines write '.maihda_sw' / '.maihda_l2wt' into the
# analytic data; if the user's data already carries a column of that name AND the
# model formula references it, fitting would clobber their variable with the
# internal weight (a brms covariate '.maihda_sw' or a WeMix covariate
# '.maihda_l2wt' would change value mid-fit). Reject that case with a rename hint.
# A reserved column merely present but NOT referenced by the formula -- e.g. one
# carried along in a prior fit's '$original_data' when maihda() refits the null /
# adjusted models -- is overwritten harmlessly and allowed.
maihda_guard_reserved_weight_col <- function(col, data, formula, engine) {
  if (col %in% names(data) && col %in% all.vars(formula)) {
    stop("'", col, "' is a reserved internal column name for the design-weighted ",
         engine, " path, but your 'data' already contains it and the model formula ",
         "references it; fitting would overwrite that variable with the internal ",
         "weight column. Rename your '", col, "' column before fitting.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a sampling-weights specification
#'
#' @param sampling_weights A single character string naming a numeric column of
#'   \code{data} holding the individual sampling (design) weights.
#' @param data The data frame the weights must live in.
#' @return The validated column name.
#' @keywords internal
maihda_validate_sampling_weights <- function(sampling_weights, data) {
  if (!is.character(sampling_weights) || length(sampling_weights) != 1 ||
      is.na(sampling_weights) || !nzchar(sampling_weights)) {
    stop("'sampling_weights' must be a single column name (a character string) ",
         "identifying the sampling-weight variable in 'data'. To use an external ",
         "vector, add it to 'data' as a column first.", call. = FALSE)
  }
  if (sampling_weights %in% c(.maihda_wemix_l2_col, .maihda_brms_weights_col)) {
    stop("'sampling_weights' may not use the reserved column name '",
         sampling_weights, "'.", call. = FALSE)
  }
  if (!sampling_weights %in% names(data)) {
    stop("Sampling-weight column not found in data: ", sampling_weights,
         call. = FALSE)
  }
  w <- data[[sampling_weights]]
  if (!is.numeric(w)) {
    stop("Sampling-weight column '", sampling_weights, "' must be numeric.",
         call. = FALSE)
  }
  if (!any(is.finite(w) & w > 0)) {
    stop("Sampling-weight column '", sampling_weights,
         "' has no positive finite values.", call. = FALSE)
  }
  sampling_weights
}

# Fingerprint of a maihda_model's SAMPLING weights (design-weighted fits), so the
# PCV and VPC comparisons do not silently mix fits with different design weights
# -- or one weighted and one unweighted fit -- whose variance estimates are not
# comparable. Unweighted fits map to "none"; a weighted fit is keyed by the weight
# column name and its values on the analytic rows. The companion to
# maihda_weight_fingerprint(), which covers lme4 PRECISION weights (and degrades
# to "unit" for engines whose prior weights are not recoverable, wemix included).
maihda_sampling_weight_fingerprint <- function(model) {
  sw <- model$sampling_weights
  if (is.null(sw)) {
    return("none")
  }
  w <- if (is.data.frame(model$data) && sw %in% names(model$data)) {
    model$data[[sw]]
  } else if (is.data.frame(model$data) && .maihda_brms_weights_col %in% names(model$data)) {
    # A brms fit's analytic frame carries the normalized weight column instead.
    model$data[[.maihda_brms_weights_col]]
  } else {
    NULL
  }
  if (is.null(w)) {
    return(paste0("col:", sw))
  }
  # Align to the analytic frame's row names so a reordered-but-identical fit gets
  # the same fingerprint (the values live in model$data, keyed by its row names).
  w <- maihda_order_by_ids(as.numeric(w), rownames(model$data))
  paste0(sw, ":",
         paste(formatC(w, format = "g", digits = 12), collapse = "\r"))
}

# Stop early with an installation hint when WeMix is unavailable.
maihda_require_wemix <- function() {
  if (!requireNamespace("WeMix", quietly = TRUE)) {
    stop("Package 'WeMix' is required for the design-weighted (engine = \"wemix\") ",
         "fit. Please install it with: install.packages('WeMix') -- or use ",
         "engine = \"brms\" for the pseudo-posterior alternative.", call. = FALSE)
  }
  invisible(TRUE)
}

# The wemix engine fits the canonical MAIHDA structure only: one intercept-only
# (1 | stratum) random effect. WeMix has no support for crossed random effects, so
# the crossed-dimensions decomposition and contextual cross-classified models must
# use lme4/brms.
maihda_wemix_check_formula <- function(formula) {
  re_terms <- reformulas::findbars(formula)
  ok <- length(re_terms) == 1 &&
    identical(paste(deparse(re_terms[[1]][[2]]), collapse = " "), "1") &&
    identical(all.vars(re_terms[[1]][[3]]), "stratum")
  if (!ok) {
    stop("engine = \"wemix\" supports the canonical MAIHDA structure only: a ",
         "single intercept-only random effect (1 | stratum) (or the (1 | var1:var2) ",
         "shorthand that resolves to it). For crossed or additional random effects ",
         "(context =, decomposition = \"crossed-dimensions\", extra (1 | g) terms), ",
         "use engine = \"lme4\" or \"brms\".", call. = FALSE)
  }
  invisible(TRUE)
}

# WeMix::mix() supports linear and binomial-logit models; the MAIHDA variance
# summaries additionally need a defined level-1 variance, so restrict to exactly
# those two families up front rather than failing inside WeMix.
maihda_wemix_check_family <- function(family) {
  ok <- (family$family == "gaussian" && family$link == "identity") ||
    (family$family == "binomial" && family$link == "logit")
  if (!ok) {
    stop("engine = \"wemix\" supports gaussian(identity) and binomial(logit) ",
         "models; this model uses ", family$family, "(", family$link, "). ",
         "Use engine = \"brms\" with sampling weights for other families.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Fit a design-weighted MAIHDA model via WeMix
#'
#' Internal engine call for \code{fit_maihda(engine = "wemix")}. Builds the
#' analytic sample (complete cases on the model variables and the weight column,
#' positive weights only) so the stored \code{data} matches the rows WeMix fits,
#' attaches the constant level-2 weight column (strata are exhaustive population
#' cells, sampled with certainty), and calls \code{WeMix::mix()} with the
#' unconditional weights \code{c(level1, level2)}.
#'
#' @param formula The resolved model formula (with \code{(1 | stratum)}).
#' @param data The data (after strata creation / response recoding).
#' @param family The resolved family object (gaussian-identity or binomial-logit).
#' @param sampling_weights Name of the level-1 sampling-weight column.
#' @param dot_vals Named list of evaluated \code{...} arguments forwarded to
#'   \code{WeMix::mix()} (e.g. \code{nQuad}, \code{verbose}, \code{fast}).
#' @return A list with \code{model} (the \code{WeMixResults}) and \code{data}
#'   (the analytic data frame actually fitted, including the weight columns).
#' @keywords internal
maihda_fit_wemix <- function(formula, data, family, sampling_weights, dot_vals) {
  maihda_guard_reserved_weight_col(.maihda_wemix_l2_col, data, formula, "wemix")
  w <- as.numeric(data[[sampling_weights]])
  # Keep exactly the rows WeMix::mix() will fit: the evaluated analytic frame (fixed-
  # effect transformations applied, rows missing AFTER them dropped) intersected with
  # finite, positive sampling weights. Passing the weight mask makes the analytic-frame
  # helper drop non-positive/NA-weight rows too (maihda_sampling_weight_mask maps them
  # to NA). A raw complete.cases() over the formula's columns misses NAs introduced by
  # a transformed term (e.g. log(x) of x <= 0), which left the weight vector longer than
  # the model matrix and triggered a weight/X row-count mismatch in mix(). Fall back to
  # the raw check only when the analytic frame cannot be built.
  keep <- maihda_analytic_keep_mask(formula, data,
                                    weights = maihda_sampling_weight_mask(w))
  if (is.null(keep)) {
    model_vars <- intersect(unique(c(all.vars(formula), sampling_weights)),
                            names(data))
    keep <- stats::complete.cases(data[, model_vars, drop = FALSE]) &
      is.finite(w) & w > 0
  }
  if (!any(keep)) {
    stop("No usable rows remain for the wemix fit after dropping rows with ",
         "missing model variables or non-positive sampling weights.", call. = FALSE)
  }
  if (sum(!keep) > 0) {
    warning(sprintf(paste0("fit_maihda(): dropped %d row(s) with missing model ",
                           "variables or non-positive sampling weights before the ",
                           "wemix fit."), sum(!keep)), call. = FALSE)
    data <- data[keep, , drop = FALSE]
  }
  data[[.maihda_wemix_l2_col]] <- 1

  # WeMix silently returns whatever it has when the Newton loop hits
  # 'max_iteration' (its own source carries an EMPTY `if (iteration >=
  # max_iteration) {}` block where a non-convergence warning belongs), and its
  # loop guard `iteration < max_iteration` accepts anything comparable: 0 and
  # negative values skip the optimisation entirely and return the starting
  # values, while a string is compared ALPHABETICALLY. Since the returned object
  # records no iteration count, a nonsense limit is unrecoverable after the fact
  # -- so reject it here, where the user can still act on it.
  dot_vals <- maihda_validate_wemix_max_iteration(dot_vals)

  args <- list(
    formula = formula,
    data = data,
    # Unconditional weights, level 1 first. The level-2 (stratum) weight is 1:
    # intersectional strata are population cells included with certainty, so the
    # level-1 conditional and unconditional weights coincide.
    weights = c(sampling_weights, .maihda_wemix_l2_col)
  )
  # mix() fits a linear mixed model unless a (binomial) family is supplied.
  if (family$family == "binomial") {
    args$family <- stats::binomial(link = "logit")
  }
  model <- do.call(WeMix::mix, c(args, dot_vals))

  list(model = model, data = data)
}

# Reject a 'max_iteration' forwarded to WeMix::mix() that would abandon (or never
# enter) the optimisation. WeMix's loop condition is `iteration < max_iteration`,
# so any value <= 0 returns the starting values untouched, and a non-numeric one
# is coerced by `<` into a silent string comparison. Returns dot_vals unchanged
# when the argument is absent or valid.
maihda_validate_wemix_max_iteration <- function(dot_vals) {
  if (!"max_iteration" %in% names(dot_vals)) {
    return(dot_vals)
  }
  mi <- dot_vals$max_iteration
  bad <- !is.numeric(mi) || length(mi) != 1L || !is.finite(mi) ||
    mi < 1 || mi != round(mi)
  if (bad) {
    stop("'max_iteration' passed to the wemix engine must be a single whole ",
         "number >= 1; got ", paste(deparse(mi), collapse = ""),
         ". WeMix stops at that many iterations WITHOUT reporting that it did, ",
         "so a value below 1 would silently return the unoptimized starting ",
         "values as if they were the fit.", call. = FALSE)
  }
  dot_vals$max_iteration <- as.integer(mi)
  dot_vals
}

# Convergence criterion of a fitted WeMixResults, or NA when it cannot be judged.
#
# WeMix reports NOTHING about how its optimisation terminated: the returned object
# carries no iteration count, no convergence code and no gradient, and the fitter
# neither warns nor errors when it exhausts 'max_iteration' (linear path) or its
# bobyqa return code is non-zero (whose `opt` object it discards). The one thing it
# does return is the log-likelihood FUNCTION, so re-derive the evidence: evaluate
# the gradient at the reported estimates and apply WeMix's own exit test,
# max(pmin(|est * g|, |g|)).
#
# The two engine paths parameterise that function differently -- the adaptive
# (non-Gaussian) path takes c(beta, variances) directly, the linear path takes
# (v = named theta, sigma, beta) and returns a list -- and neither signature is
# documented API. So the objective is SELF-CHECKED first: unless it reproduces the
# fit's own reported log-likelihood, this returns NA and the caller reports
# "unknown" rather than risk a verdict built on a misread parameterisation.
maihda_wemix_gradient_criterion <- function(model) {
  tryCatch({
    lnlf <- model$lnlf
    lnl0 <- suppressWarnings(as.numeric(model$lnl))
    if (!is.function(lnlf) || length(lnl0) != 1L || !is.finite(lnl0)) {
      return(NA_real_)
    }
    if (isTRUE(model$is_adaptive)) {
      pars <- c(as.numeric(model$coef), as.numeric(model$vars))
      f <- function(p) sum(as.numeric(lnlf(p)))
    } else {
      theta <- model$theta
      k <- length(model$coef)
      q <- length(theta)
      pars <- c(as.numeric(model$coef), as.numeric(theta), as.numeric(model$sigma))
      if (q < 1L || !is.finite(pars[k + q + 1L])) {
        return(NA_real_)
      }
      f <- function(p) {
        v <- p[k + seq_len(q)]
        # The linear-path objective looks its variance terms up BY NAME.
        names(v) <- names(theta)
        as.numeric(lnlf(v = v, sigma = p[k + q + 1L],
                        beta = p[seq_len(k)])$lnl)
      }
    }
    if (length(pars) < 1L || !all(is.finite(pars))) {
      return(NA_real_)
    }
    base <- f(pars)
    # Self-check: refuse to judge unless the rebuilt objective agrees with the
    # likelihood the fit itself reports.
    if (!is.finite(base) ||
        abs(base - lnl0) > 1e-4 * max(1, abs(lnl0))) {
      return(NA_real_)
    }
    grad <- vapply(seq_along(pars), function(j) {
      h <- max(1e-6, 1e-6 * abs(pars[j]))
      up <- down <- pars
      up[j] <- up[j] + h
      down[j] <- down[j] - h
      (f(up) - f(down)) / (2 * h)
    }, numeric(1))
    if (!all(is.finite(grad))) {
      return(NA_real_)
    }
    max(pmin(abs(pars * grad), abs(grad)))
  }, error = function(e) NA_real_)
}

# Convergence verdict + messages for a WeMixResults, mirroring the structure the
# lme4 / ordinal / brms branches of maihda_fit_diagnostics() produce.
#
# `singular` is the caller's boundary flag. A variance pinned at zero sits on the
# EDGE of the parameter space, where the gradient legitimately does not vanish
# (the optimality condition is one-sided), so the stationarity test cannot
# distinguish a converged boundary fit from an abandoned one -- report "unknown"
# there and let the separately-surfaced singular flag carry the warning.
maihda_wemix_convergence <- function(model, singular = NA) {
  fail <- function(msg) list(converged = FALSE, messages = msg)

  # A field that is PRESENT but non-finite is failure evidence; a field that is
  # ABSENT is no evidence at all, and must not be read as one.
  lnl <- suppressWarnings(as.numeric(model$lnl))
  if (length(lnl) == 1L && !is.finite(lnl)) {
    return(fail("The weighted pseudo-log-likelihood is not finite: the fit did not converge."))
  }
  coefs <- suppressWarnings(as.numeric(model$coef))
  if (length(coefs) > 0 && !all(is.finite(coefs))) {
    return(fail("One or more fixed-effect estimates are not finite: the fit did not converge."))
  }
  ses <- suppressWarnings(as.numeric(model$SE))
  if (length(ses) > 0 && !all(is.finite(ses))) {
    return(fail(paste0("One or more fixed-effect standard errors are not finite: the ",
                       "fit did not converge, or its information matrix is singular.")))
  }

  crit <- maihda_wemix_gradient_criterion(model)
  # Deliberately far looser than WeMix's own 1e-5 exit test. This gradient is a
  # central DIFFERENCE taken in the variance parameterisation, not the analytic one
  # the fitter optimises in, so the two do not agree to WeMix's precision: across
  # healthy gaussian and binomial fits the criterion here ran from 0 to 8e-04,
  # while fits whose optimisation was abandoned sat between 0.1 and 5. The
  # threshold is placed in that gap. The purpose is to catch a grossly unfinished
  # optimisation, NOT to second-guess where WeMix chose to stop -- a false alarm on
  # a sound fit would cost more than missing a marginal one.
  if (!is.na(crit) && crit < .maihda_wemix_grad_tol) {
    return(list(converged = TRUE, messages = character(0)))
  }
  if (!is.na(crit) && !isTRUE(singular)) {
    return(fail(sprintf(
      paste0("The optimisation stopped well away from a maximum (gradient ",
             "criterion %.3g, against a %g threshold): WeMix returns the last ",
             "iterate without reporting this, so treat the estimates as ",
             "provisional. Consider raising 'max_iteration' or 'nQuad'."),
      crit, .maihda_wemix_grad_tol)))
  }
  # No usable evidence either way (unreadable likelihood function, an unfamiliar
  # WeMix version, or a boundary fit): say so rather than assume success.
  list(converged = NA, messages = character(0))
}

#' Variance components of a wemix MAIHDA fit
#'
#' Reads the between-stratum variance (and, for a linear model, the residual
#' variance) from the \code{WeMixResults} variance table. For a binomial-logit
#' model the level-1 variance is the usual latent-scale \eqn{\pi^2/3}, matching
#' the lme4/brms summaries.
#'
#' @param object A \code{maihda_model} with engine \code{"wemix"}.
#' @return A list with \code{stratum} and \code{residual} variances.
#' @keywords internal
maihda_wemix_variances <- function(object) {
  vd <- object$model$varDF
  if (is.null(vd) || !all(c("grp", "vcov") %in% names(vd))) {
    stop("Could not read the variance components from the WeMix fit.", call. = FALSE)
  }
  s_rows <- vd$grp == "stratum" &
    (is.na(vd$var1) | vd$var1 %in% c("(Intercept)", "Intercept"))
  if (!any(s_rows)) {
    stop("No 'stratum' random-effect variance found in the WeMix fit.", call. = FALSE)
  }
  var_stratum <- as.numeric(vd$vcov[s_rows][1])

  if (object$family$family == "gaussian") {
    r_rows <- vd$grp == "Residual"
    if (!any(r_rows)) {
      stop("No residual variance found in the WeMix fit.", call. = FALSE)
    }
    var_residual <- as.numeric(vd$vcov[r_rows][1])
  } else {
    # binomial-logit: latent-scale level-1 variance, as in the other engines.
    var_residual <- (pi^2) / 3
  }

  list(stratum = var_stratum, residual = var_residual)
}

#' Fixed-part (and optionally full) linear predictor of a wemix fit
#'
#' WeMix's own \code{predict()} method needs the grouping structure re-resolved
#' and offers no fixed-only form, so predictions are built directly from the
#' coefficient vector and the stored stratum effects: the fixed design matrix is
#' constructed with the training data's factor levels AND transformation basis (so a
#' data-dependent term such as \code{scale(x)} uses the fit's centre and scale rather
#' than recomputing them from \code{newdata}) and multiplied by
#' \code{coef}, any formula offset term is evaluated on \code{newdata} and added,
#' and \code{include_re} adds each row's stratum effect (conditional
#' mode; an unseen stratum contributes 0 -- the population-average fallback that
#' \code{\link{predict_maihda}} only reaches when \code{allow_new_levels = TRUE},
#' having otherwise rejected unseen strata upstream). Everything is on the link
#' scale.
#'
#' @param object A \code{maihda_model} with engine \code{"wemix"}.
#' @param newdata Data to predict for; defaults to the analytic data.
#' @param include_re Add the stratum random effect (conditional mode)?
#' @return A numeric vector of link-scale predictions.
#' @keywords internal
maihda_wemix_linpred <- function(object, newdata = NULL, include_re = TRUE) {
  if (is.null(newdata)) {
    newdata <- object$data
  }
  # Terms rebuilt from the FITTED data, so a data-dependent transformation such as
  # scale(x) / poly(x, 2) / ns(x, 3) evaluates on the fit's basis instead of being
  # recomputed from the prediction batch (see maihda_fitted_predict_terms()).
  basis <- maihda_fitted_predict_terms(object$formula, object$data)
  tt <- basis$terms
  mf <- stats::model.frame(tt, newdata, xlev = basis$xlev,
                           na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf)
  beta <- object$model$coef
  missing_cols <- setdiff(names(beta), colnames(X))
  if (length(missing_cols) > 0) {
    stop("Could not rebuild the wemix design matrix; missing column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  eta <- drop(X[, names(beta), drop = FALSE] %*% beta)

  # A formula offset term (offset(.) in the model formula) is part of the linear
  # predictor WeMix::mix() fits but is NOT a column of X, so model.matrix() never
  # rebuilds it -- add it back explicitly from the model frame, evaluated on
  # newdata, or response-scale predictions would be off by the offset.
  off <- stats::model.offset(mf)
  if (!is.null(off)) {
    eta <- eta + off
  }

  if (include_re) {
    re <- maihda_wemix_ranef_vector(object)
    u <- re[as.character(newdata$stratum)]
    u[is.na(u)] <- 0
    eta <- eta + unname(u)
  }
  eta
}

# Named vector of stratum random-effect estimates (conditional modes) from the
# WeMix fit, keyed by stratum label.
maihda_wemix_ranef_vector <- function(object) {
  rm <- object$model$ranefMat
  if (is.null(rm) || !"stratum" %in% names(rm)) {
    stop("No 'stratum' random effects found in the WeMix fit.", call. = FALSE)
  }
  tab <- rm[["stratum"]]
  cols <- intersect(c("(Intercept)", "Intercept"), colnames(tab))
  if (length(cols) == 0) {
    stop("The 'stratum' random effect must include an intercept for MAIHDA ",
         "stratum estimates.", call. = FALSE)
  }
  stats::setNames(as.numeric(tab[[cols[1]]]), rownames(tab))
}

#' Stratum random-effect table for a wemix fit
#'
#' Mirrors \code{maihda_stratum_ranef_lme4()}: one row per stratum with the
#' random-effect estimate (conditional mode), a conditional standard error, and a
#' 95\% interval. WeMix reports no conditional variances, so the SE is computed
#' analytically from the weighted pseudo-likelihood: for a Gaussian model the
#' conditional precision of \eqn{u_j} is \eqn{1/\tau^2 + \sum_j w_{ij}/\sigma^2}
#' (the design-weighted analogue of lme4's \code{condVar}, to which it reduces at
#' unit weights), and for a binomial-logit model the Laplace curvature at the
#' conditional mode, \eqn{1/\tau^2 + \sum_j w_{ij}\,\hat p_{ij}(1-\hat p_{ij})}.
#' These are model-based approximations, not design-based (replicate-weight)
#' uncertainty.
#'
#' @param object A \code{maihda_model} with engine \code{"wemix"}.
#' @return A data frame with \code{stratum}, \code{stratum_id},
#'   \code{random_effect}, \code{se}, \code{lower_95}, \code{upper_95}.
#' @keywords internal
maihda_wemix_stratum_ranef <- function(object) {
  re <- maihda_wemix_ranef_vector(object)
  vars <- maihda_wemix_variances(object)
  tau2 <- vars$stratum

  data <- object$data
  w <- maihda_prior_weights(object)
  strata <- as.character(data$stratum)

  if (object$family$family == "gaussian") {
    sigma2 <- vars$residual
    info <- vapply(names(re), function(s) {
      sum(w[strata == s], na.rm = TRUE) / sigma2
    }, numeric(1))
  } else {
    # Curvature of the weighted Bernoulli log-likelihood at the conditional mode.
    p <- stats::plogis(maihda_wemix_linpred(object, include_re = TRUE))
    info <- vapply(names(re), function(s) {
      sel <- strata == s
      sum(w[sel] * p[sel] * (1 - p[sel]), na.rm = TRUE)
    }, numeric(1))
  }

  se <- if (is.finite(tau2) && tau2 > 0) {
    sqrt(1 / (1 / tau2 + info))
  } else {
    # Boundary fit (zero between-stratum variance): the conditional distribution
    # collapses on 0, so the SE is 0 rather than undefined.
    rep(0, length(re))
  }

  data.frame(
    stratum = names(re),
    stratum_id = suppressWarnings(as.integer(names(re))),
    random_effect = unname(re),
    se = unname(se),
    lower_95 = unname(re - 1.96 * se),
    upper_95 = unname(re + 1.96 * se),
    stringsAsFactors = FALSE
  )
}

#' Per-stratum predictions for a wemix fit
#'
#' wemix counterpart of \code{maihda_stratum_predictions_lme4()}: per-stratum
#' means of the fixed-part prediction plus the stratum effect, aggregated with
#' the SAMPLING weights so the stratum-level summaries are design-weighted
#' (population-representative under the weights), unlike the lme4 prior-weight
#' aggregation.
#'
#' @param object A \code{maihda_model} with engine \code{"wemix"}.
#' @param summary_obj Its \code{maihda_summary} (for the stratum estimates).
#' @param scale "response" or "link".
#' @return A data frame as from \code{maihda_weighted_stratum_aggregate()}.
#' @keywords internal
maihda_stratum_predictions_wemix <- function(object, summary_obj,
                                             scale = c("response", "link")) {
  scale <- match.arg(scale)
  data <- object$data
  if (!"stratum" %in% names(data)) {
    stop("'stratum' variable not found in fitted model data.")
  }

  linkinv <- maihda_linkinv(object$family)
  prior_w <- maihda_prediction_weights(object)
  eta_fixed <- maihda_wemix_linpred(object, include_re = FALSE)

  stratum_est <- summary_obj$stratum_estimates
  if (is.null(stratum_est) || nrow(stratum_est) == 0) {
    stop("No stratum estimates available.")
  }

  key <- as.character(data$stratum)
  idx <- match(key, as.character(stratum_est$stratum))
  transform_eta <- function(eta) {
    if (scale == "response") linkinv(eta) else eta
  }

  pred_df <- data.frame(
    stratum = key,
    predicted_row = transform_eta(eta_fixed + stratum_est$random_effect[idx]),
    lower_row = transform_eta(eta_fixed + stratum_est$lower_95[idx]),
    upper_row = transform_eta(eta_fixed + stratum_est$upper_95[idx]),
    fixed_row = transform_eta(eta_fixed),
    weight = prior_w,
    stringsAsFactors = FALSE
  )

  maihda_weighted_stratum_aggregate(
    pred_df, c("predicted_row", "lower_row", "upper_row", "fixed_row")
  )
}

#' Inject sampling weights into a brms formula
#'
#' Rewrites \code{y ~ ...} as \code{y | weights(w) ~ ...}. An existing addition
#' term (e.g. an aggregated-binomial \code{y | trials(n)}) is extended with
#' \code{+ weights(w)}; a formula that already carries a \code{weights()} addition
#' term is rejected (the two weight specifications would conflict).
#'
#' @param formula The model formula.
#' @param wcol Name of the (normalized) weight column.
#' @return The rewritten formula (same environment).
#' @keywords internal
maihda_brms_weights_formula <- function(formula, wcol) {
  lhs <- formula[[2]]
  weights_call <- call("weights", as.name(wcol))
  if (is.call(lhs) && identical(lhs[[1]], as.name("|"))) {
    if (grepl("weights\\s*\\(", paste(deparse(lhs[[3]]), collapse = " "))) {
      stop("The formula already carries a weights() addition term; supply the ",
           "sampling weights either there or via 'sampling_weights', not both.",
           call. = FALSE)
    }
    lhs[[3]] <- call("+", lhs[[3]], weights_call)
  } else {
    lhs <- call("|", lhs, weights_call)
  }
  formula[[2]] <- lhs
  formula
}

#' Remove the package-internal weights(.maihda_sw) addition term from a formula
#'
#' The inverse of the injection \code{maihda_brms_weights_formula()} performs for
#' the reserved \code{.maihda_sw} column. \code{maihda()} and
#' \code{compare_maihda_groups()} derive the null / adjusted models from a
#' sampling-weighted fit's \emph{stored} formula, which already carries the
#' injected \code{weights(.maihda_sw)} term; that term must be stripped before the
#' derived model is re-prepared, or the reserved-column guard (\code{.maihda_sw}
#' both in the formula and in the carried-over data) and the weights-formula
#' rewrite (which rejects an existing \code{weights()} term) abort the package's
#' own refit. Only the internal \code{weights(.maihda_sw)} call is removed -- a
#' user's own \code{weights()} term (a genuine conflict) is left in place so it is
#' still rejected -- and any other addition term (e.g. \code{trials(n)}) is
#' preserved.
#'
#' @param formula A model formula.
#' @return The formula with the internal \code{weights(.maihda_sw)} term removed,
#'   or the input unchanged when it carries no such term.
#' @keywords internal
maihda_strip_brms_weights_term <- function(formula) {
  if (length(formula) < 3L) {
    return(formula)                    # one-sided formula: no LHS to strip
  }
  lhs <- formula[[2]]
  if (!is.call(lhs) || !identical(lhs[[1]], as.name("|"))) {
    return(formula)                    # response carries no addition term
  }
  stripped <- maihda_drop_weights_call(lhs[[3]], .maihda_brms_weights_col)
  if (is.null(stripped)) {
    # The addition was solely weights(.maihda_sw): drop the `|`, leaving the bare
    # response.
    formula[[2]] <- lhs[[2]]
  } else if (!identical(stripped, lhs[[3]])) {
    lhs[[3]] <- stripped
    formula[[2]] <- lhs
  }
  formula
}

# Recursively drop a weights(<wcol>) call from a brms addition-term expression (a
# `+`-tree of addition terms, e.g. `trials(n) + weights(.maihda_sw)`). Returns the
# pruned expression, or NULL when nothing remains.
maihda_drop_weights_call <- function(expr, wcol) {
  is_target <- is.call(expr) && identical(expr[[1]], as.name("weights")) &&
    length(expr) == 2L && identical(expr[[2]], as.name(wcol))
  if (is_target) {
    return(NULL)
  }
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    left <- maihda_drop_weights_call(expr[[2]], wcol)
    right <- maihda_drop_weights_call(expr[[3]], wcol)
    if (is.null(left) && is.null(right)) {
      return(NULL)
    }
    if (is.null(left)) {
      return(right)
    }
    if (is.null(right)) {
      return(left)
    }
    return(call("+", left, right))
  }
  expr
}

#' Rows complete on every variable a brms model will use
#'
#' Mirrors the rows brms retains after its own NA exclusion: complete on the
#' response and its addition-term variables (\code{y | trials(n)}), on the
#' fixed-effect terms -- evaluated through the model frame, so a transformed
#' predictor whose transformation yields \code{NA}/\code{NaN} counts as
#' incomplete -- and on every variable in the random-effect terms.
#'
#' @param formula The (pre-weights-injection) model formula.
#' @param data The model data.
#' @return A logical vector over the rows of \code{data}.
#' @keywords internal
maihda_brms_complete_rows <- function(formula, data) {
  ok <- rep(TRUE, nrow(data))
  fixed <- maihda_nobars(formula)
  tt <- tryCatch(stats::delete.response(stats::terms(fixed, data = data)),
                 error = function(e) NULL)
  if (!is.null(tt)) {
    mf <- tryCatch(stats::model.frame(tt, data, na.action = stats::na.pass),
                   error = function(e) NULL)
    if (!is.null(mf) && nrow(mf) == nrow(data)) {
      ok <- ok & stats::complete.cases(mf)
    }
  }
  raw_vars <- unique(c(
    if (length(formula) == 3) all.vars(formula[[2]]),
    unlist(lapply(reformulas::findbars(formula), all.vars), use.names = FALSE)
  ))
  raw_vars <- intersect(raw_vars, names(data))
  if (length(raw_vars) > 0) {
    ok <- ok & stats::complete.cases(data[raw_vars])
  }
  ok
}

#' Prepare data and formula for a sampling-weighted brms fit
#'
#' Drops rows with missing or non-positive sampling weights and rows incomplete
#' on the model variables (each with a warning), normalizes the remaining
#' weights to mean 1 -- likelihood weights scale the effective sample size, so
#' unnormalized expansion weights (summing to the population) would massively
#' overstate the information in the data -- and rewrites the formula with a
#' \code{weights()} addition term.
#'
#' Incomplete rows must go BEFORE the normalization: brms silently excludes
#' them after receiving the data, so normalizing over all weight-valid rows
#' would hand the sampler surviving weights with an arbitrary mean, scaling the
#' pseudo-posterior's effective sample size by that mean.
#'
#' @param data The model data.
#' @param formula The model formula.
#' @param sampling_weights Name of the sampling-weight column.
#' @return A list with \code{data} (weights column \code{.maihda_sw} added),
#'   \code{formula} (rewritten), and \code{keep} (a logical mask over the input
#'   rows marking those retained, so the caller can re-slice any row-aligned
#'   forwarded arguments to the same rows).
#' @keywords internal
maihda_prepare_brms_sampling_weights <- function(data, formula, sampling_weights) {
  # A maihda()- / compare_maihda_groups()-derived refit re-enters here with the
  # internal weights term already on the formula and the .maihda_sw column already
  # in `data` -- both copied from the prior sampling-weighted fit the null /
  # adjusted model was derived from. Strip that term so this re-prepares cleanly
  # from `sampling_weights`: otherwise the reserved-column guard (below) sees
  # .maihda_sw referenced in the formula AND present in data and aborts, and the
  # formula rewrite (further down) rejects the already-present weights() term. The
  # weight column is re-normalized and re-injected from the ORIGINAL weight column
  # (still carried in the derived data), which is the intended behaviour.
  formula <- maihda_strip_brms_weights_term(formula)
  maihda_guard_reserved_weight_col(.maihda_brms_weights_col, data, formula, "brms")
  w <- as.numeric(data[[sampling_weights]])
  keep <- is.finite(w) & w > 0
  if (!any(keep)) {
    stop("No usable rows remain after dropping missing or non-positive sampling ",
         "weights.", call. = FALSE)
  }
  if (sum(!keep) > 0) {
    warning(sprintf(paste0("fit_maihda(): dropped %d row(s) with missing or ",
                           "non-positive sampling weights before the brms fit."),
                    sum(!keep)), call. = FALSE)
  }
  incomplete <- keep & !maihda_brms_complete_rows(formula, data)
  if (sum(incomplete) > 0) {
    keep <- keep & !incomplete
    if (!any(keep)) {
      stop("No usable rows remain after dropping rows incomplete on the model ",
           "variables.", call. = FALSE)
    }
    warning(sprintf(paste0("fit_maihda(): dropped %d row(s) incomplete on the ",
                           "model variables before normalizing the sampling ",
                           "weights (brms would drop them only after ",
                           "normalization, leaving the fitted weights with a ",
                           "mean different from 1)."),
            sum(incomplete)), call. = FALSE)
  }
  if (sum(!keep) > 0) {
    data <- data[keep, , drop = FALSE]
    w <- w[keep]
  }
  data[[.maihda_brms_weights_col]] <- w * length(w) / sum(w)
  list(
    data = data,
    formula = maihda_brms_weights_formula(formula, .maihda_brms_weights_col),
    keep = keep
  )
}
