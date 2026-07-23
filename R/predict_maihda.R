#' Predict from MAIHDA Model
#'
#' Makes predictions from a fitted MAIHDA model, either at the stratum level
#' or individual level.
#'
#' @param object A maihda_model object from \code{fit_maihda()}.
#' @param newdata Optional data frame for making predictions. If NULL, uses the
#'   original data from model fitting.
#' @param type Character string specifying prediction type:
#'   \itemize{
#'     \item "individual": Individual-level predictions including random effects
#'     \item "strata": Stratum-level predictions (random effects only). For a
#'       longitudinal (growth-curve) fit a stratum is a \emph{trajectory}, so this
#'       returns the per-stratum trajectory parameters (baseline deviation, random
#'       intercept and random slope(s)) rather than a single random effect. A
#'       \emph{non-longitudinal} fit whose stratum random effects include random
#'       slopes (a hand-written \code{(1 + x | stratum)}) is an error here, as in
#'       \code{summary()}: the single value this would return is the intercept
#'       alone -- the stratum effect only where the slope variables are zero.
#'   }
#'   For backward compatibility, "link" or "response" may also be passed here
#'   and will be interpreted as individual-level predictions on that scale.
#' @param scale Character string specifying the prediction scale for
#'   individual-level predictions: "response" (default) or "link". For a
#'   cumulative (ordinal) model the "link" scale is the latent location
#'   \eqn{\eta} and the "response" scale is the \emph{expected category score}
#'   \eqn{\sum_k k P(Y = k)} (categories scored 1..K in their declared order).
#'   For an aggregated-binomial fit (an lme4 \code{cbind(success, failure)} or a
#'   brms \code{success | trials(n)}) the "response" scale is the per-trial
#'   \emph{probability} on both engines (not the expected success count).
#' @param allow_new_levels Logical. By default (\code{FALSE}) a stratum in
#'   \code{newdata} that the model never saw -- whether supplied directly as a
#'   \code{stratum} column or rebuilt from the grouping variables -- is an error,
#'   for every engine, matching \pkg{lme4}'s default. Set \code{TRUE} to instead
#'   predict unseen strata with the stratum random effect dropped (treated as
#'   zero), while keeping any \emph{other} random effect the row participates in
#'   (e.g. a contextual \code{(1 | school)} intercept from
#'   \code{fit_maihda(context = )}, or a longitudinal growth term) -- the same
#'   behaviour as \pkg{lme4}'s \code{allow.new.levels}, which zeroes only the unseen
#'   level's effect and keeps seen ones. For the usual single-stratum model the
#'   stratum is the only random effect, so this is the \emph{population-average}
#'   (fixed-effects-only) prediction. This affects \code{type = "individual"} only:
#'   a stratum-level prediction (\code{type = "strata"}) has no random effect to
#'   report for an unseen stratum, so unseen strata remain an error there
#'   regardless.
#' @param ... Additional arguments passed to predict method of underlying model.
#'
#' @return Depending on type:
#'   \itemize{
#'     \item For "individual": A numeric vector of predicted values on the
#'       requested scale
#'     \item For "strata": A data frame with stratum ID and predicted random
#'       effect. When \code{newdata} is supplied, the result is restricted to the
#'       strata present in \code{newdata} (and a stratum the model never saw is an
#'       error, as for "individual"); when \code{newdata} is \code{NULL}, every
#'       training stratum is returned.
#'   }
#'
#' @examples
#' \donttest{
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' model <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata_result$data)
#' 
#' # Individual predictions
#' pred_ind <- predict_maihda(model, type = "individual")
#' 
#' # Stratum predictions
#' pred_strata <- predict_maihda(model, type = "strata")
#' }
#'
#' @export
#' @importFrom stats predict
#' @importFrom lme4 ranef
predict_maihda <- function(object, newdata = NULL,
                           type = c("individual", "strata", "response", "link"),
                           scale = c("response", "link"),
                           allow_new_levels = FALSE, ...) {
  if (!inherits(object, "maihda_model")) {
    stop("'object' must be a maihda_model object from fit_maihda()")
  }
  if (!is.logical(allow_new_levels) || length(allow_new_levels) != 1 ||
      is.na(allow_new_levels)) {
    stop("'allow_new_levels' must be TRUE or FALSE.", call. = FALSE)
  }

  type <- match.arg(type)
  if (type %in% c("response", "link")) {
    scale <- type
    type <- "individual"
  } else {
    scale <- match.arg(scale)
  }
  engine <- object$engine
  model <- object$model

  newdata_supplied <- !is.null(newdata)
  if (!newdata_supplied) {
    newdata <- object$data
  } else {
    newdata <- maihda_prepare_prediction_data(object, newdata, type = type,
                                              allow_new_levels = allow_new_levels)
  }

  # Longitudinal (growth-curve) model: a stratum is a TRAJECTORY, so the
  # stratum-level prediction is the random intercept and slope(s), not a single
  # value. Individual-level predictions flow through the engine branches below
  # (predict()/posterior_linpred handle the random slopes).
  if (type == "strata" && !is.null(object$longitudinal_info)) {
    res <- maihda_longitudinal_strata_predictions(object)
    return(maihda_filter_strata_predictions(res, newdata))
  }

  if (engine == "lme4") {
    if (type == "individual") {
      # Individual-level predictions including random effects. Unseen strata are
      # already rejected upstream unless allow_new_levels = TRUE, in which case
      # lme4 must be told to permit them (it then sets their random effect to 0,
      # the same population-average fallback the other engines use).
      dots <- maihda_dots_default(list(...), "allow.new.levels",
                                  isTRUE(allow_new_levels))
      # Training-data predictions (no newdata supplied) call predict() WITHOUT
      # newdata so lme4 reuses its stored linear predictor, which includes any
      # offset. An external offset (offset = ... passed to fit_maihda()) survives
      # only in the fitted (offset) column, which predict.merMod ignores on the
      # newdata path; substituting object$data would therefore silently drop it
      # from individual predictions -- and from the AUC/MOR, tables, and plots
      # built on them.
      if (!newdata_supplied) {
        return(do.call(stats::predict, c(list(model, type = scale), dots)))
      }
      # For genuine newdata, predict.merMod re-evaluates a formula offset() term
      # but cannot recover an external offset (only its fitted values were stored,
      # not the generating expression), so those predictions would be silently
      # wrong. Reject them with a directed error rather than return them.
      if (maihda_lme4_has_external_offset(object)) {
        stop("This model was fit with an external offset (offset = ... passed to ",
             "fit_maihda()), which cannot be reconstructed for new data. Refit with ",
             "the offset written into the formula (e.g. ... + offset(log(exposure))) ",
             "to predict on newdata.", call. = FALSE)
      }
      predictions <- do.call(stats::predict,
                             c(list(model, newdata = newdata, type = scale), dots))
      return(predictions)

    } else if (type == "strata") {
      # Stratum-level predictions (random effects)
      result <- maihda_stratum_ranef_lme4(model)
      result$predicted <- result$random_effect
      result <- result[, c("stratum", "predicted", "se", "lower_95", "upper_95")]
      return(maihda_filter_strata_predictions(result, newdata))
    }
    
  } else if (engine == "wemix") {
    if (type == "individual") {
      # Built from coef + the stored stratum effects (WeMix's own predict() needs
      # the grouping re-resolved and offers no scale argument).
      eta <- maihda_wemix_linpred(object, newdata = newdata, include_re = TRUE)
      if (scale == "response") {
        return(maihda_linkinv(object$family)(eta))
      }
      return(eta)

    } else if (type == "strata") {
      result <- maihda_wemix_stratum_ranef(object)
      result$predicted <- result$random_effect
      result <- result[, c("stratum", "predicted", "se", "lower_95", "upper_95")]
      return(maihda_filter_strata_predictions(result, newdata))
    }

  } else if (engine == "ordinal") {
    if (type == "individual") {
      # predict.clmm does not exist; the latent location eta = x'beta + u is
      # built from the stored coefficients and stratum conditional modes. The
      # "link" scale is that latent location; the "response" scale is the
      # expected category score sum_k k * P(Y = k) (categories scored 1..K in
      # order), the package's response-scale summary of a cumulative model.
      eta <- maihda_clmm_linpred(object, newdata = newdata, include_re = TRUE)
      if (scale == "response") {
        return(maihda_ordinal_eta_to_score(eta, object$model$alpha,
                                           object$family$link))
      }
      return(eta)

    } else if (type == "strata") {
      result <- maihda_clmm_stratum_ranef(object)
      result$predicted <- result$random_effect
      result <- result[, c("stratum", "predicted", "se", "lower_95", "upper_95")]
      return(maihda_filter_strata_predictions(result, newdata))
    }

  } else if (engine == "brms") {
    # Verify brms is available
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required to predict from brms models. Please install it with: install.packages('brms')")
    }

    # A sampling-weighted brms fit carries a weights() addition term, and brms
    # requires its column in newdata even though predictions do not depend on it;
    # supply a unit weight when the caller's newdata lacks it.
    if (!is.null(object$sampling_weights) && ".maihda_sw" %in% all.vars(object$formula) &&
        !".maihda_sw" %in% names(newdata)) {
      newdata$.maihda_sw <- 1
    }

    if (type == "individual") {
      # Individual-level predictions. As for lme4, unseen levels are rejected
      # upstream unless allow_new_levels = TRUE. Forwarding allow_new_levels to brms
      # is NOT enough to honour the documented zero-effect fallback -- brms's default
      # sample_new_levels = "uncertainty" DRAWS a new effect from the estimated
      # random-effects distribution rather than treating it as zero -- so
      # maihda_brms_individual_prediction() zeroes each unseen grouping level (via a
      # per-row re_formula dropping only the terms that row has not seen), stratum
      # AND context/longitudinal alike, matching lme4's allow.new.levels.
      dots <- maihda_dots_default(list(...), "allow_new_levels",
                                  isTRUE(allow_new_levels))
      predictions <- maihda_brms_individual_prediction(object, newdata, scale,
                                                       allow_new_levels, dots)
      return(predictions)

    } else if (type == "strata") {
      # Stratum-level predictions
      result <- maihda_stratum_ranef_brms(model)
      result$predicted <- result$random_effect
      result <- result[, c("stratum", "predicted", "se", "lower_95", "upper_95")]
      return(maihda_filter_strata_predictions(result, newdata))
    }
  }
}

#' Predict method for maihda_model objects
#'
#' S3 \code{\link[stats]{predict}} method for MAIHDA models: a thin alias of
#' \code{\link{predict_maihda}}, so \code{predict(model, type = "strata")} and
#' \code{predict_maihda(model, type = "strata")} are interchangeable. See
#' \code{\link{predict_maihda}} for the full documentation of the arguments
#' and return value.
#'
#' @inheritParams predict_maihda
#' @inherit predict_maihda return
#' @seealso \code{\link{predict_maihda}}
#' @export
predict.maihda_model <- function(object, newdata = NULL,
                                 type = c("individual", "strata", "response",
                                          "link"),
                                 scale = c("response", "link"),
                                 allow_new_levels = FALSE, ...) {
  predict_maihda(object, newdata = newdata, type = type, scale = scale,
                 allow_new_levels = allow_new_levels, ...)
}

# Inject a named default into a list of forwarded `...` arguments unless the
# caller already supplied it, so an engine's new-levels switch can be set from
# `allow_new_levels` without clashing with a user-supplied value of the same name.
maihda_dots_default <- function(dots, name, value) {
  if (!name %in% names(dots)) {
    dots[[name]] <- value
  }
  dots
}

# TRUE when an lme4 fit carries an EXTERNAL offset (offset = ... passed to
# fit_maihda()) rather than a formula offset() term. The distinction matters for
# newdata predictions: predict.merMod re-evaluates a formula offset() from newdata
# but silently ignores an external one (it lives only as the fitted (offset)
# column). The model frame names that column "(offset)" for either kind of offset,
# so a formula offset() term -- identifiable from the fixed-part terms -- is what
# separates the two.
maihda_lme4_has_external_offset <- function(object) {
  if (!identical(object$engine, "lme4")) {
    return(FALSE)
  }
  mf <- object$data
  if (is.null(mf) || !"(offset)" %in% names(mf)) {
    return(FALSE)
  }
  fixed <- tryCatch(maihda_nobars(object$formula),
                    error = function(e) object$formula)
  offset_terms <- tryCatch(attr(stats::terms(fixed), "offset"),
                           error = function(e) NULL)
  is.null(offset_terms) || length(offset_terms) == 0
}

# Individual-level brms predictions, honouring the documented unseen-stratum
# fallback. allow_new_levels = TRUE promises a prediction that drops the stratum
# random effect (treating it as zero) for a stratum the model never saw. brms does
# NOT do this by simply receiving allow_new_levels = TRUE: its default
# sample_new_levels = "uncertainty" SAMPLES a new random effect from the random-
# effects distribution. lme4's allow.new.levels instead ZEROES the effect of any
# unseen grouping level -- stratum OR a context/longitudinal grouping -- while
# keeping the effects of the levels a row HAS seen. To reproduce that, each row is
# predicted with an re_formula that drops exactly the grouping terms whose level
# that row never saw, and keeps the rest; rows sharing a kept-term signature are
# predicted together. A row that keeps every term is the ordinary full-random-
# effects prediction; one that keeps none is the fixed-effects-only population
# average (re_formula = NA). Unseen levels are only possible when
# allow_new_levels = TRUE (otherwise upstream validation / brms rejects them), so
# without it every row takes the full path unchanged. Previously only the STRATUM
# term was zeroed for unseen strata; an unseen context/longitudinal level was left
# in the re_formula and thus SAMPLED, diverging from lme4 (a different, marginal
# estimand, visible on the response scale) -- this generalisation closes that gap.
maihda_brms_individual_prediction <- function(object, newdata, scale,
                                              allow_new_levels, dots) {
  bars <- maihda_brms_re_bars(object)

  # No random effects, no new-levels request, or no rows: every row uses the full
  # model. brms rejects a genuinely new level upstream, so there is nothing to zero.
  if (!isTRUE(allow_new_levels) || length(bars) == 0 || nrow(newdata) == 0) {
    return(maihda_brms_predict_rows(object, newdata, scale, dots))
  }

  # For each bar, the grouping levels seen in training and the level each newdata
  # row takes; a bar is KEPT for a row only where that row's level is known. The
  # check is conservative: a bar whose training levels or newdata column cannot be
  # resolved is kept (never silently dropped), so at worst the previous behaviour is
  # retained. An NA grouping value is likewise left to the normal path.
  known   <- lapply(bars, function(b) maihda_brms_bar_known_levels(object, b))
  row_lab <- lapply(bars, function(b) maihda_brms_bar_row_levels(newdata, b))
  keep <- matrix(TRUE, nrow = nrow(newdata), ncol = length(bars))
  for (j in seq_along(bars)) {
    kn <- known[[j]]
    rl <- row_lab[[j]]
    if (length(kn) > 0 && !is.null(rl)) {
      keep[, j] <- is.na(rl) | rl %in% kn
    }
  }

  # The caller's own random-effect scope, if they set one. Zeroing an unseen level
  # may only NARROW that scope -- never widen or replace it. Overwriting the scope
  # outright reintroduced grouping terms the caller had explicitly excluded (most
  # visibly re_formula = NA, the fixed-effects-only population average, which came
  # back as "~ (1 | <other group>)"), and for a partial re_formula it substituted a
  # different term for the requested one.
  scope <- maihda_brms_requested_re(dots)
  if (!scope$understood) {
    # An re_formula/re.form this helper cannot interpret is left strictly alone --
    # deferring to the caller rather than silently overriding them. brms still
    # receives allow_new_levels via dots and applies its own handling.
    return(maihda_brms_predict_rows(object, newdata, scale, dots))
  }

  # Predict each distinct kept-bar signature once, under the scope that keeps exactly
  # the requested bars minus the ones this row never saw (unchanged dots when the
  # row keeps every bar, or when none of the dropped bars was requested anyway).
  key <- apply(keep, 1L, function(z) paste0(which(z), collapse = ","))
  pred <- rep(NA_real_, nrow(newdata))
  for (k in unique(key)) {
    idx <- which(key == k)
    kept_j <- keep[idx[1L], ]
    grp_dots <- dots
    if (!all(kept_j)) {
      dropped <- vapply(bars[!kept_j], maihda_brms_bar_key, character(1))
      if (is.null(scope$bars)) {
        # No caller restriction: keep every model bar this row has seen.
        grp_dots[[scope$name]] <- maihda_brms_re_formula_from_bars(bars[kept_j])
      } else {
        effective <- Filter(
          function(b) !(maihda_brms_bar_key(b) %in% dropped), scope$bars)
        if (length(effective) < length(scope$bars)) {
          grp_dots[[scope$name]] <- maihda_brms_re_formula_from_bars(effective)
        }
      }
    }
    pred[idx] <- maihda_brms_predict_rows(
      object, newdata[idx, , drop = FALSE], scale, grp_dots)
  }
  pred
}

# A random-effect bar's deparsed text, the key used to tell two bars apart (and the
# same form maihda_brms_re_formula_from_bars() writes out).
maihda_brms_bar_key <- function(bar) {
  paste(deparse(bar, width.cutoff = 500L), collapse = " ")
}

# The random-effect scope a caller requested through predict_maihda()'s `...`.
# Returns the argument NAME to write any narrowed scope back under (so a caller who
# used one spelling never ends up with both), plus the requested `bars`: NULL for
# brms's unrestricted default (re_formula = NULL or absent), an empty list for
# re_formula = NA (fixed effects only), otherwise the bars of the supplied formula.
# brms honours the lme4 spelling re.form as an alias on some prediction paths but not
# others, and lets re.form win when both are given; this resolves the same way.
# `understood = FALSE` marks a value outside brms's documented {NULL, NA, formula}
# contract, which the caller then owns untouched.
maihda_brms_requested_re <- function(dots) {
  nm <- intersect(c("re.form", "re_formula"), names(dots))
  if (length(nm) == 0) {
    return(list(understood = TRUE, name = "re_formula", bars = NULL))
  }
  nm <- nm[1L]
  v <- dots[[nm]]
  if (is.null(v)) {
    return(list(understood = TRUE, name = nm, bars = NULL))
  }
  if (inherits(v, "formula")) {
    b <- tryCatch(reformulas::findbars(v), error = function(e) NULL)
    return(list(understood = TRUE, name = nm,
                bars = if (is.null(b)) list() else b))
  }
  if (is.logical(v) && length(v) == 1L && is.na(v)) {
    return(list(understood = TRUE, name = nm, bars = list()))
  }
  list(understood = FALSE, name = nm, bars = NULL)
}

# The random-effect bars of a brms MAIHDA fit's stored formula (brmsformula-aware),
# or an empty list when there are none or the formula is unusable. Shared by the
# unseen-level prediction logic so bar extraction cannot drift.
maihda_brms_re_bars <- function(object) {
  f <- object$formula
  if (inherits(f, "brmsformula") && inherits(f$formula, "formula")) {
    f <- f$formula
  }
  if (!inherits(f, "formula")) {
    return(list())
  }
  bars <- tryCatch(reformulas::findbars(f), error = function(e) NULL)
  if (is.null(bars)) list() else bars
}

# Build a brms re_formula "~ (bar1) + (bar2) + ..." from a set of random-effect
# bars. Returns NA -- brms's re_formula that drops ALL group terms (the fixed-
# effects-only population average) -- when the set is empty.
maihda_brms_re_formula_from_bars <- function(bars) {
  if (length(bars) == 0) {
    return(NA)
  }
  terms_chr <- vapply(bars, function(b)
    paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"),
    character(1))
  stats::as.formula(paste("~", paste(terms_chr, collapse = " + ")),
                    env = baseenv())
}

# A brms re_formula that keeps every grouping term EXCEPT the stratum one (kept for
# backward compatibility and the Stan-free unit test). Now a thin wrapper over the
# shared bar helpers.
maihda_brms_unseen_re_formula <- function(object) {
  bars <- maihda_brms_re_bars(object)
  keep <- Filter(function(b) !("stratum" %in% all.vars(b[[3]])), bars)
  maihda_brms_re_formula_from_bars(keep)
}

# The grouping levels a random-effect bar's grouping factor took in the TRAINING
# data (character labels), for telling a seen level from a new one. An empty result
# means "cannot determine" -- the caller then keeps the bar rather than risk
# dropping a real effect.
maihda_brms_bar_known_levels <- function(object, bar) {
  data <- object$data
  grp <- bar[[3]]
  if (!is.data.frame(data) || !all(all.vars(grp) %in% names(data))) {
    return(character(0))
  }
  lev <- tryCatch(as.character(eval(grp, data)), error = function(e) NULL)
  if (is.null(lev)) character(0) else unique(lev[!is.na(lev)])
}

# The grouping level each newdata row takes for a random-effect bar's grouping
# factor, or NULL when the grouping column(s) are absent from newdata (the caller
# then keeps the bar so brms raises its usual missing-column error rather than the
# effect being silently dropped).
maihda_brms_bar_row_levels <- function(newdata, bar) {
  grp <- bar[[3]]
  if (!all(all.vars(grp) %in% names(newdata))) {
    return(NULL)
  }
  lab <- tryCatch(as.character(eval(grp, newdata)), error = function(e) NULL)
  if (is.null(lab) || length(lab) != nrow(newdata)) NULL else lab
}

# Response- or link-scale brms predictions for a block of rows. On the response
# scale an aggregated-binomial `y | trials(n)` fit returns the expected COUNT
# (trials * p), so it is normalised to a probability (matching lme4's cbind() fit);
# a categorical-likelihood fit returns an nobs x summary x category array, collapsed
# to the expected category score (categories scored 1..K in order). The link scale
# is the posterior-mean linear predictor (maihda_brms_linpred_mean; brms's
# posterior_linpred returns draws and ignores a summary= argument).
maihda_brms_predict_rows <- function(object, nd, scale, dots) {
  model <- object$model
  if (scale == "response") {
    f <- do.call(stats::fitted,
                 c(list(model, newdata = nd, summary = TRUE), dots))
    if (length(dim(f)) == 3) {
      return(maihda_brms_fitted_array_scores(f))
    }
    return(maihda_brms_response_to_prob(object, f[, "Estimate"], nd))
  }
  do.call(maihda_brms_linpred_mean, c(list(model, newdata = nd), dots))
}

# Expected category scores sum_k k * P(Y = k) from a brms fitted() SUMMARY array
# (nobs x summary-statistics x category), as returned for a cumulative /
# categorical likelihood. The Estimate slice is rebuilt as an explicit
# nobs x ncat matrix: plain f[, "Estimate", ] drops BOTH unit margins for a
# single prediction row, leaving a bare category vector whose ncol() is NULL
# (seq_len() then errors), so one-row newdata used to fail here.
maihda_brms_fitted_array_scores <- function(f) {
  est <- matrix(f[, "Estimate", ], nrow = dim(f)[1], ncol = dim(f)[3])
  drop(est %*% seq_len(ncol(est)))
}

# Restrict a per-stratum prediction table to the strata present in `newdata` so
# type = "strata" respects newdata the way type = "individual" does (instead of
# always returning every training stratum). A stratum in newdata that the model
# never saw is an error, matching the individual path. When newdata carries no
# stratum column the table is returned unchanged. When it carries a stratum column
# whose values are all missing it names no training stratum to keep, so the result
# is empty -- not silently every training stratum.
maihda_filter_strata_predictions <- function(result, newdata) {
  if (is.null(newdata) || !"stratum" %in% names(newdata)) {
    return(result)
  }
  wanted <- unique(as.character(newdata$stratum))
  wanted <- wanted[!is.na(wanted)]
  if (length(wanted) == 0) {
    return(result[0, , drop = FALSE])
  }
  known <- as.character(result$stratum)
  unknown <- setdiff(wanted, known)
  if (length(unknown) > 0) {
    stop("newdata contains strata not present in the fitted model: ",
         paste(utils::head(unknown, 5), collapse = ", "),
         if (length(unknown) > 5) ", ..." else "",
         call. = FALSE)
  }
  result[known %in% wanted, , drop = FALSE]
}
