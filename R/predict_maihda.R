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
#'       intercept and random slope(s)) rather than a single random effect.
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

  if (is.null(newdata)) {
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
      # Individual-level predictions. As for lme4, unseen strata are rejected
      # upstream unless allow_new_levels = TRUE. Forwarding allow_new_levels to brms
      # is NOT enough to honour the documented zero-effect fallback for unseen
      # strata -- brms's default sample_new_levels = "uncertainty" DRAWS a new
      # stratum effect from the estimated random-effects distribution rather than
      # treating it as zero -- so the unseen rows are split off and predicted with
      # re_formula = NA (see maihda_brms_individual_prediction()). allow_new_levels
      # is still forwarded for the seen rows so an unseen *context* level (a
      # different kind of new level) keeps working as before.
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

# Inject a named default into a list of forwarded `...` arguments unless the
# caller already supplied it, so an engine's new-levels switch can be set from
# `allow_new_levels` without clashing with a user-supplied value of the same name.
maihda_dots_default <- function(dots, name, value) {
  if (!name %in% names(dots)) {
    dots[[name]] <- value
  }
  dots
}

# Individual-level brms predictions, honouring the documented unseen-stratum
# fallback. allow_new_levels = TRUE promises a prediction that drops the stratum
# random effect (treating it as zero) for a stratum the model never saw. brms does
# NOT do this by simply receiving allow_new_levels = TRUE: its default
# sample_new_levels = "uncertainty" SAMPLES a new stratum effect from the random-
# effects distribution. So the unseen-stratum rows are predicted separately with an
# re_formula that excludes the stratum term (see maihda_brms_unseen_re_formula());
# the seen-stratum rows keep their estimated stratum effect. A blanket re_formula on
# all rows would wrongly drop the seen strata's effects too, which is why the rows
# are split. Any OTHER random effect the unseen row participates in -- a contextual
# (1 | school) intercept, a longitudinal (time | id) growth term -- is retained,
# matching lme4's allow.new.levels (which zeroes only the unseen level's effect); for
# the usual single-stratum model the excluding re_formula is NA (fixed effects only),
# the population average. Unseen rows are only possible when allow_new_levels = TRUE
# (otherwise upstream validation rejected them), so without it every row takes the
# ordinary full-random-effects path.
maihda_brms_individual_prediction <- function(object, newdata, scale,
                                              allow_new_levels, dots) {
  known <- maihda_known_strata(object)
  unseen <- if (isTRUE(allow_new_levels) && !is.null(known) &&
                "stratum" %in% names(newdata)) {
    !as.character(newdata$stratum) %in% known
  } else {
    rep(FALSE, nrow(newdata))
  }

  pred <- rep(NA_real_, nrow(newdata))
  if (any(!unseen)) {
    pred[!unseen] <- maihda_brms_predict_rows(
      object, newdata[!unseen, , drop = FALSE], scale, dots)
  }
  if (any(unseen)) {
    unseen_dots <- dots
    # Drop the stratum effect (treat as zero) but keep any non-stratum random effect.
    unseen_dots$re_formula <- maihda_brms_unseen_re_formula(object)
    pred[unseen] <- maihda_brms_predict_rows(
      object, newdata[unseen, , drop = FALSE], scale, unseen_dots)
  }
  pred
}

# A brms re_formula that keeps every grouping (random-effect) term EXCEPT the
# stratum one, for predicting an unseen stratum: the stratum random effect is
# dropped (treated as zero -- the population-average fallback) while any other
# random effect the row participates in is retained, matching lme4's
# allow.new.levels (which zeroes only the unseen level's effect, keeping seen ones).
# Returns NA -- the brms re_formula that drops ALL group terms -- when the stratum is
# the only random effect, i.e. the ordinary fixed-effects-only population average.
# A bar's grouping factor is read with all.vars() so a context (1 | school) or a
# longitudinal (poly | id) term is kept while (... | stratum) is removed.
maihda_brms_unseen_re_formula <- function(object) {
  f <- object$formula
  if (inherits(f, "brmsformula") && inherits(f$formula, "formula")) {
    f <- f$formula
  }
  if (!inherits(f, "formula")) {
    return(NA)
  }
  bars <- tryCatch(reformulas::findbars(f), error = function(e) NULL)
  if (is.null(bars) || length(bars) == 0) {
    return(NA)
  }
  keep <- Filter(function(b) !("stratum" %in% all.vars(b[[3]])), bars)
  if (length(keep) == 0) {
    return(NA)
  }
  terms_chr <- vapply(keep, function(b)
    paste0("(", paste(deparse(b, width.cutoff = 500L), collapse = " "), ")"),
    character(1))
  stats::as.formula(paste("~", paste(terms_chr, collapse = " + ")),
                    env = baseenv())
}

# Response- or link-scale brms predictions for a block of rows. On the response
# scale an aggregated-binomial `y | trials(n)` fit returns the expected COUNT
# (trials * p), so it is normalised to a probability (matching lme4's cbind() fit);
# a categorical-likelihood fit returns an nobs x summary x category array, collapsed
# to the expected category score (categories scored 1..K in order).
maihda_brms_predict_rows <- function(object, nd, scale, dots) {
  model <- object$model
  if (scale == "response") {
    f <- do.call(stats::fitted,
                 c(list(model, newdata = nd, summary = TRUE), dots))
    if (length(dim(f)) == 3) {
      est <- f[, "Estimate", ]
      return(drop(est %*% seq_len(ncol(est))))
    }
    return(maihda_brms_response_to_prob(object, f[, "Estimate"], nd))
  }
  do.call(brms::posterior_linpred,
          c(list(model, newdata = nd, summary = TRUE), dots))[, "Estimate"]
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
