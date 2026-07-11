# Discriminatory accuracy for binary MAIHDA models.
#
# The "DA" in MAIHDA is discriminatory accuracy: how well the intersectional
# strata separate individuals who do and do not have a binary outcome. The VPC
# summarises *variation* between strata; discriminatory accuracy summarises
# *prediction* at the individual level, and the two can diverge sharply (a high
# between-stratum VPC can still translate into only modest individual-level AUC).
# This file promotes the helpers previously sketched in the "binary_outcomes"
# vignette to first-class, tested, exported functions.

#' Area under the ROC curve (C-statistic), rank-based
#'
#' @description
#' Computes the AUC / C-statistic as the Mann-Whitney U statistic: the
#' probability that a randomly chosen case (\code{y == 1}) is assigned a higher
#' predicted value than a randomly chosen non-case (\code{y == 0}), with ties
#' counting as one half. This needs no external package. An AUC of 0.5 is chance;
#' 1 is perfect separation.
#'
#' @param prob Numeric vector of predicted probabilities (or any score where
#'   larger means more case-like).
#' @param y Observed binary outcome as 0/1 numeric or logical, the same length as
#'   \code{prob}.
#'
#' @return A single number in \code{[0, 1]}, or \code{NA_real_} if either class is
#'   absent.
#'
#' @references
#' Merlo, J., Wagner, P., Ghith, N., & Leckie, G. (2016). An original stepwise
#' multilevel logistic regression analysis of discriminatory accuracy: the case of
#' neighbourhoods and health. \emph{PLOS ONE}, 11(4), e0153778.
#'
#' @examples
#' maihda_auc(c(0.1, 0.4, 0.35, 0.8), c(0, 0, 1, 1))
#'
#' @export
maihda_auc <- function(prob, y) {
  if (!is.numeric(prob)) {
    stop("'prob' must be a numeric vector of predicted probabilities/scores.",
         call. = FALSE)
  }
  if (is.logical(y)) {
    y <- as.integer(y)
  }
  y <- suppressWarnings(as.numeric(y))
  if (length(prob) != length(y)) {
    stop("'prob' and 'y' must have the same length (", length(prob), " vs ",
         length(y), ").", call. = FALSE)
  }

  keep <- !(is.na(prob) | is.na(y))
  prob <- prob[keep]
  y <- y[keep]

  if (!all(y %in% c(0, 1))) {
    stop("'y' must be a binary 0/1 (or logical) outcome. For a factor outcome, ",
         "convert it to 0/1 first (e.g. as.numeric(factor) - 1).", call. = FALSE)
  }

  # Coerce to double: sum(<logical>) is integer, so n1 * n0 (and n1 * (n1 + 1))
  # overflow the 32-bit integer limit for large samples -- n1 * n0 exceeds
  # .Machine$integer.max once both classes are sizeable (around n ~ 130k at a 13%
  # event rate) -- silently producing NA. Double arithmetic avoids the overflow.
  n1 <- as.double(sum(y == 1))
  n0 <- as.double(sum(y == 0))
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  r <- rank(prob)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Median Odds Ratio (MOR) for a logistic MAIHDA model
#'
#' @description
#' The Median Odds Ratio translates the between-stratum variance of a logistic
#' MAIHDA model onto the odds-ratio scale: the median relative change in the odds
#' of the outcome when comparing two individuals from randomly chosen strata
#' (higher- vs lower-risk). \code{MOR = exp(sqrt(2 * V_A) * qnorm(0.75))}, where
#' \code{V_A} is the between-stratum (latent, logit-scale) variance. An MOR of 1
#' indicates no between-stratum heterogeneity. The MOR is defined only for the
#' \strong{logit} link (it is the median \emph{odds} ratio); a non-logit binomial
#' fit such as \code{probit} is rejected, because its latent variance is on a
#' different scale and the \code{exp(...)} above would not be an odds ratio.
#'
#' For a \strong{cumulative-logit} (ordinal) MAIHDA model the same formula
#' applies to the latent logit-scale between-stratum variance and is the
#' \emph{median cumulative odds ratio}: the median relative change in the odds
#' of being at or below any given outcome category between two randomly chosen
#' strata (under the model's proportional-odds assumption it is the same for
#' every category split).
#'
#' \strong{Scope.} \code{V_A} is the \emph{between-stratum} variance. For a
#' crossed-dimensions fit (from \code{maihda(decomposition =
#' "crossed-dimensions")}) that is the \emph{total} intersectional variance --
#' the additive dimension variances plus the interaction component, since the
#' independent crossed random effects sum at the intersection level. Contextual
#' (\code{context = }) and other non-stratum random effects are never included:
#' the MOR compares two individuals from randomly chosen strata \emph{within
#' the same context}.
#'
#' @param model A \code{maihda_model} from \code{\link{fit_maihda}} fitted with a
#'   \code{binomial} (lme4), \code{bernoulli} (brms), or \code{cumulative}
#'   (ordinal) family and a \strong{logit} link.
#'
#' @return A single number (the MOR, \eqn{\ge 1}), or \code{NA_real_} if the
#'   between-stratum variance is unavailable.
#'
#' @references
#' Larsen, K., & Merlo, J. (2005). Appropriate assessment of neighborhood effects
#' on individual health: integrating random and fixed effects in multilevel
#' logistic regression. \emph{American Journal of Epidemiology}, 161(1), 81-88.
#'
#' @seealso \code{\link{maihda_discriminatory_accuracy}}
#'
#' @export
maihda_mor <- function(model) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  fam <- maihda_model_family_name(model)
  # Accept the lme4 "binomial" family, the brms "bernoulli" family a binary 0/1
  # outcome is fit with (fit_maihda(engine = "brms") routes Bernoulli data to
  # bernoulli()), and the "cumulative" (ordinal) family, whose logit-scale
  # latent variance yields the median cumulative odds ratio.
  if (!isTRUE(fam %in% c("binomial", "bernoulli", "cumulative"))) {
    stop("The Median Odds Ratio is only defined for binomial/Bernoulli (logistic) ",
         "and cumulative-logit (ordinal) MAIHDA models; this model uses family = '",
         fam, "'.", call. = FALSE)
  }
  link <- maihda_model_link_name(model)
  if (!identical(link, "logit")) {
    stop("The Median Odds Ratio is defined only for the logit link -- it is the ",
         "median *odds* ratio, derived from the logistic latent variance -- but this ",
         "model uses a '", fam, "' model with the '", link, "' link, whose latent ",
         "variance is on a different scale.", call. = FALSE)
  }

  v_a <- tryCatch(maihda_mor_between_variance(model), error = function(e) NA_real_)
  if (!is.numeric(v_a) || length(v_a) != 1 || !is.finite(v_a) || v_a < 0) {
    return(NA_real_)
  }

  exp(sqrt(2 * v_a) * stats::qnorm(0.75))
}

# Between-stratum (intersectional) variance for the MOR. For the canonical fit
# this is the stratum intercept variance. For a crossed-dimensions fit
# ($cc_info) the between-stratum effect is the SUM of the independent crossed
# REs -- each additive dimension effect plus the interaction -- so its variance
# is their total; reading only the "stratum" (interaction) component would
# understate the between-stratum heterogeneity the MOR translates. Contextual
# and other non-intersectional REs are never included: the MOR compares two
# individuals from different strata within the same context.
maihda_mor_between_variance <- function(model) {
  if (is.null(model$cc_info)) {
    return(extract_between_variance(model))
  }
  groups <- unique(c(model$cc_info$interaction_group,
                     unname(model$cc_info$dim_groups)))
  var_named <- if (identical(model$engine, "lme4")) {
    maihda_random_variances_lme4(model$model)
  } else if (identical(model$engine, "brms")) {
    maihda_random_variances_brms(model$model)
  } else {
    stop("Crossed-dimensions MOR is supported for the lme4 and brms engines ",
         "only.", call. = FALSE)
  }
  missing_re <- setdiff(groups, names(var_named))
  if (length(missing_re) > 0) {
    stop("Crossed-dimensions MOR is missing the random effect(s): ",
         paste(missing_re, collapse = ", "), ".", call. = FALSE)
  }
  sum(var_named[groups])
}

#' Discriminatory accuracy of a binary MAIHDA model
#'
#' @description
#' Bundles the individual-level discriminatory-accuracy summaries for a binomial
#' MAIHDA model: the AUC / C-statistic (how well the model's predicted
#' probabilities separate cases from non-cases) and the Median Odds Ratio. Applied
#' to a strata-only (null) model, the AUC is the discriminatory accuracy of the
#' intersectional strata themselves -- Merlo's central quantity; comparing it with
#' an adjusted model shows whether individual covariates beyond stratum membership
#' sharpen classification. The AUC is computed for any binomial link; the Median
#' Odds Ratio is reported only for the logit link and is \code{NA} otherwise (e.g.
#' for a probit fit), since the MOR is an odds-ratio-scale quantity.
#'
#' \strong{The AUC is apparent (in-sample).} It scores the same observations used
#' to estimate the model -- the fixed effects, the variance components, and the
#' shrunken stratum BLUPs -- so it is an \emph{apparent} (resubstitution) AUC and
#' is optimistically biased, the more so with small or sparse strata. It is
#' reported as the conventional descriptive MAIHDA discriminatory accuracy
#' (Merlo 2018), \strong{not} as a cross-validated estimate of out-of-sample
#' predictive discrimination; the returned object carries \code{apparent = TRUE}
#' and \code{print()} says so. For genuine predictive accuracy, validate with an
#' out-of-fold (group-aware) scheme or an optimism correction.
#'
#' Aggregated-binomial fits are supported on both engines that fit them -- an lme4
#' \code{cbind(success, failure)} response and a brms \code{y | trials(n)} response:
#' the AUC is the count-weighted C-statistic over the implied individual-level 0/1
#' data, and \code{n_case} / \code{n_control} are the total successes / failures.
#'
#' \strong{Scope.} When the model carries random effects \emph{beyond} the
#' intersectional partition -- a contextual \code{(1 | school)} from
#' \code{fit_maihda(context = )} or an explicit extra grouping such as
#' \code{(1 | site)} -- the headline \code{auc} is the concordance of the
#' \emph{intersectional strata} (fixed effects plus the stratum random effect
#' and, for a crossed-dimensions fit, the additive dimension effects), matching
#' the scope of the MOR; the concordance of the full model including the other
#' random effects is reported separately as \code{auc_full}
#' (\code{auc_scope = "strata"}). For the canonical single-\code{(1 | stratum)}
#' model the two coincide and only \code{auc} is reported
#' (\code{auc_scope = "model"}, \code{auc_full} absent).
#'
#' @param model A \code{maihda_model} from \code{\link{fit_maihda}} fitted with a
#'   \code{binomial} family -- including an aggregated response (an lme4
#'   \code{cbind(success, failure)} or a brms \code{y | trials(n)}) -- or the
#'   \code{bernoulli} family a binary 0/1 outcome is fit with under
#'   \code{engine = "brms"}.
#'
#' @return An object of class \code{maihda_da}: a list with \code{auc},
#'   \code{auc_scope}, \code{auc_full}, \code{mor},
#'   \code{n_case}, \code{n_control}, \code{family}, \code{link}, \code{engine},
#'   \code{weighted}, \code{weight_type} and \code{apparent} (always \code{TRUE} --
#'   the AUC is in-sample; see Description). \code{mor} is \code{NA} for a
#'   non-logit binomial link, where the AUC is still reported. For an
#'   aggregated-binomial fit \code{n_case} / \code{n_control} are the total
#'   successes / failures. \code{weighted} is \code{TRUE} when the AUC is a
#'   weighted concordance -- \code{weight_type} \code{"sampling"} for a
#'   design-weighted fit (each observation contributes its sampling weight as
#'   case/control mass, estimating the population discriminatory accuracy) or
#'   \code{"precision"} for an lme4 fit with non-integer precision weights;
#'   \code{n_case} / \code{n_control} stay unweighted observation counts either
#'   way. \code{weight_type} is \code{NULL} for an unweighted AUC.
#'
#' @references
#' Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
#' discriminatory accuracy (MAIHDA) within an intersectional framework.
#' \emph{Social Science & Medicine}, 203, 74-80.
#'
#' @seealso \code{\link{maihda_auc}}, \code{\link{maihda_mor}}
#'
#' @examples
#' \donttest{
#' # Obese (Yes/No) by intersectional strata of Gender x Race
#' strata <- make_strata(maihda_health_data, vars = c("Gender", "Race"))
#' d <- maihda_health_data
#' d$stratum <- strata$data$stratum
#' m <- fit_maihda(Obese ~ (1 | stratum), data = d, family = "binomial")
#' maihda_discriminatory_accuracy(m)
#' }
#'
#' @export
maihda_discriminatory_accuracy <- function(model) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  fam <- maihda_model_family_name(model)
  # Accept both the lme4 "binomial" family and the brms "bernoulli" family a binary
  # 0/1 outcome is fit with (see fit_maihda(engine = "brms")); both are logistic
  # MAIHDA models for which AUC / MOR are defined.
  if (!isTRUE(fam %in% c("binomial", "bernoulli"))) {
    stop("Discriminatory accuracy (AUC / MOR) is only defined for binomial/Bernoulli ",
         "(logistic) MAIHDA models; this model uses family = '", fam, "'.",
         call. = FALSE)
  }

  prob <- predict_maihda(model, type = "individual", scale = "response")
  resp <- maihda_da_observed_response(model)

  # fit_maihda() supports aggregated-binomial responses on both engines that fit them
  # -- an lme4 cbind(success, failure) and a brms `y | trials(n)`. For lme4, detect
  # aggregation STRUCTURALLY from the fitted model frame -- glmer stores a cbind()
  # response as a two-column [successes, failures] matrix -- and read the trial counts
  # straight from it. For brms there is no weights.brmsfit, so the trials come from the
  # same path the prediction-weighting uses, maihda_brms_trial_counts(), which parses
  # the trials() addition term off the stored formula. Either way we then compute a
  # count-weighted AUC over the implied individual-level 0/1 data instead of passing a
  # non-0/1 response to maihda_auc() (which errors). The earlier lme4 heuristic inferred
  # aggregation from the response carrying values outside {0, 1}, which silently failed
  # when every aggregated row was all-success or all-failure (proportions exactly 0/1):
  # such a fit fell through to an unweighted row-level AUC with one pseudo-observation
  # per stratum and wrong case/control totals. agg_counts is NULL for a true Bernoulli
  # fit (which takes the ordinary rank-based path) and is aligned to prob's rows.
  agg_counts <- if (identical(model$engine, "lme4")) {
    maihda_da_aggregated_counts(model)
  } else if (identical(model$engine, "brms")) {
    maihda_da_brms_aggregated_counts(model)
  } else {
    NULL
  }
  aggregated <- !is.null(agg_counts) && length(agg_counts$successes) == length(prob)
  # Design-weighted fit (sampling_weights supplied): compute the design-weighted
  # AUC, where each observation contributes its sampling weight as case (y = 1) or
  # control (y = 0) mass -- the weighted Mann-Whitney concordance, estimating the
  # POPULATION discriminatory accuracy rather than the sample's. The reported
  # case/control totals stay unweighted observation counts.
  sw <- if (!is.null(model$sampling_weights)) maihda_prior_weights(model) else NULL
  design_weighted <- !is.null(sw) && length(sw) == length(prob) &&
    any(is.finite(sw)) && !isTRUE(all(abs(sw - 1) < sqrt(.Machine$double.eps)))
  # lme4 precision-weighted Bernoulli fit (non-unit, non-integer prior weights,
  # e.g. 1.5): NOT an aggregated binomial -- integer-trial fits were caught by
  # agg_counts above. Each observation contributes its weight as case/control
  # mass (the weighted Mann-Whitney concordance) instead of rounding y * weight
  # into pseudo-counts against fractional trials, which produced an AUC above 1
  # and inflated case totals.
  pw <- NULL
  if (!aggregated && !design_weighted && identical(model$engine, "lme4")) {
    pw_try <- tryCatch(as.numeric(stats::weights(model$model, type = "prior")),
                       error = function(e) NULL)
    if (!is.null(pw_try) && length(pw_try) == length(prob) &&
        all(is.finite(pw_try)) &&
        any(abs(pw_try - 1) > sqrt(.Machine$double.eps))) {
      pw <- pw_try
    }
  }
  # One AUC evaluator over an arbitrary ranking score, so the intersectional-
  # scope AUC below reuses exactly the same case/control weighting as the
  # full-model AUC (the AUC is rank-based, so any monotone score works).
  if (aggregated) {
    successes <- agg_counts$successes
    trials <- agg_counts$trials
    # The count-weighted AUC ranks each row by its per-trial predicted probability.
    # predict_maihda(scale = "response") returns that probability on both engines now
    # (lme4's cbind() fit and a brms `y | trials(n)` fit, which normalises its expected
    # success COUNT by the trial counts internally), so the rows rank directly by it.
    # Reported case/control totals stay unweighted observation counts (matching the
    # design_weighted Bernoulli branch below), so read them off the raw counts before
    # any weighting.
    n_case <- sum(successes, na.rm = TRUE)
    n_control <- sum(trials - successes, na.rm = TRUE)
    # When the same fit ALSO carries sampling weights (a brms `y | trials(n)` fit with
    # sampling_weights -- lme4 never reaches here, as sampling_weights routes to wemix),
    # fold them into the per-row case/control mass so the AUC is the design-weighted
    # (population) concordance. Without this, `weighted = design_weighted` is reported
    # (line below) and print.maihda_da() claims a design-weighted AUC while the number
    # is actually the unweighted one. The ranking score stays unweighted -- it is a
    # per-trial probability, not a mass.
    if (design_weighted) {
      successes <- sw * successes
      trials <- sw * trials
    }
    auc_for <- function(score) maihda_auc_weighted(score, successes, trials)
  } else if (design_weighted) {
    auc_for <- function(score) {
      maihda_auc_weighted(score, successes = sw * resp, trials = sw)
    }
    n_case <- sum(resp == 1, na.rm = TRUE)
    n_control <- sum(resp == 0, na.rm = TRUE)
  } else if (!is.null(pw)) {
    auc_for <- function(score) {
      maihda_auc_weighted(score, successes = pw * resp, trials = pw)
    }
    n_case <- sum(resp == 1, na.rm = TRUE)
    n_control <- sum(resp == 0, na.rm = TRUE)
  } else {
    auc_for <- function(score) maihda_auc(score, resp)
    n_case <- sum(resp == 1, na.rm = TRUE)
    n_control <- sum(resp == 0, na.rm = TRUE)
  }

  # Scope. With non-intersectional random effects in the model (a contextual
  # (1 | school) or an explicit (1 | site)), the full-model prediction folds
  # their effects into the score while the MOR is a between-STRATUM quantity --
  # the two would summarize different models (a strong site effect can carry a
  # high full-model AUC over a negligible stratum effect). The headline AUC is
  # therefore the intersectional-scope concordance -- fixed effects plus the
  # stratum RE (and, for a crossed-dimensions fit, the additive dimension REs)
  # -- matching the MOR's scope, with the full-model AUC reported alongside as
  # auc_full. For the canonical single-(1 | stratum) model the two coincide and
  # only auc is reported.
  scopes <- maihda_da_re_scopes(model)
  auc_full <- NULL
  auc_scope <- "model"
  if (length(scopes$other) > 0) {
    score <- tryCatch(maihda_da_scope_scores(model, scopes$intersectional),
                      error = function(e) NULL)
    if (!is.null(score) && length(score) == length(prob)) {
      auc_full <- auc_for(prob)
      auc <- auc_for(score)
      auc_scope <- "strata"
    } else {
      # The scoped score could not be built; report the full-model AUC rather
      # than nothing, labelled as such.
      auc <- auc_for(prob)
    }
  } else {
    auc <- auc_for(prob)
  }

  # The AUC is link-agnostic (rank-based on predicted probabilities), but the MOR is
  # defined only for the logit link. For other binomial links (e.g. probit) report
  # the AUC with mor = NA rather than an odds ratio that is off the model's scale.
  link <- maihda_model_link_name(model)
  mor <- if (identical(link, "logit")) maihda_mor(model) else NA_real_

  structure(
    list(
      auc = auc,
      auc_scope = auc_scope,
      auc_full = auc_full,
      mor = mor,
      n_case = n_case,
      n_control = n_control,
      family = fam,
      link = link,
      engine = model$engine,
      weighted = design_weighted || !is.null(pw),
      weight_type = if (design_weighted) "sampling" else if (!is.null(pw)) "precision",
      # The AUC is APPARENT (in-sample / resubstitution): it scores the same rows
      # used to estimate the fixed effects, variance components, and stratum BLUPs,
      # so it is optimistically biased -- more so with small/sparse strata, where
      # the shrunken BLUPs still track their own rows. This is the conventional
      # MAIHDA discriminatory accuracy (Merlo 2018), reported as a descriptive
      # in-sample measure; it is NOT a cross-validated estimate of out-of-sample
      # predictive discrimination. Flagged so print()/tidiers can label it and a
      # future out-of-fold variant can set apparent = FALSE.
      apparent = TRUE
    ),
    class = "maihda_da"
  )
}

#' @export
print.maihda_da <- function(x, ...) {
  pal <- maihda_palette()
  cat(pal$bold("Discriminatory accuracy (binomial MAIHDA)"), "\n", sep = "")
  cat(sprintf("  AUC (C-statistic): %s%s\n",
              if (is.finite(x$auc)) pal$accent(sprintf("%.3f", x$auc)) else "NA",
              if (identical(x$auc_scope, "strata")) " [intersectional strata]" else ""))
  if (!is.null(x$auc_full) && is.finite(x$auc_full)) {
    cat(sprintf("  AUC, full model:   %s\n",
                pal$accent(sprintf("%.3f", x$auc_full))))
    cat(pal$muted(paste0(
      "  (The model carries non-stratum random effects; the headline AUC is\n",
      "  the concordance of the intersectional strata -- matching the MOR's\n",
      "  scope -- while the full-model AUC includes the other random effects.)\n")))
  }
  mor_str <- if (is.finite(x$mor)) {
    pal$accent(sprintf("%.3f", x$mor))
  } else if (!is.null(x$link) && !identical(x$link, "logit")) {
    sprintf("NA (requires the logit link; model uses '%s')", x$link)
  } else {
    "NA"
  }
  cat(sprintf("  Median Odds Ratio: %s\n", mor_str))
  cat(sprintf("  Cases / controls:  %d / %d\n", x$n_case, x$n_control))
  if (!isFALSE(x$apparent)) {
    cat(pal$muted(paste0(
      "  (AUC is apparent / in-sample: scored on the same rows used to fit the\n",
      "  model, so it is optimistically biased -- more so with sparse strata. It\n",
      "  is a descriptive measure, not cross-validated out-of-sample discrimination.)\n")))
  }
  if (isTRUE(x$weighted)) {
    msg <- if (identical(x$weight_type, "precision")) {
      paste0(
        "  (AUC is precision-weighted: each observation contributes its prior\n",
        "  weight as case/control mass; cases/controls are unweighted counts.)\n")
    } else {
      paste0(
        "  (AUC is design-weighted: each observation contributes its sampling\n",
        "  weight; cases/controls are unweighted counts.)\n")
    }
    cat(pal$muted(msg))
  }
  invisible(x)
}

# ---- internal helpers -------------------------------------------------------

# Resolve the family name ("binomial"/"gaussian"/"poisson"/...) of a
# maihda_model, tolerating either a stored family object/string or falling back
# to the fitted model's family. Every path is canonicalised via
# maihda_normalize_family_name() so engine-specific labels (e.g. a fixed-theta
# MASS::negative.binomial(2) family object stored as "Negative Binomial(2)")
# compare against fixed names.
maihda_model_family_name <- function(model) {
  fam <- model$family
  if (is.list(fam) && !is.null(fam$family)) {
    return(maihda_normalize_family_name(fam$family))
  }
  if (is.character(fam) && length(fam) == 1) {
    return(maihda_normalize_family_name(fam))
  }
  ff <- tryCatch(maihda_family(model$model), error = function(e) NULL)
  if (!is.null(ff) && !is.null(ff$family)) {
    return(ff$family)
  }
  NA_character_
}

# Resolve the link name ("logit"/"probit"/...) of a maihda_model, preferring the
# stored family object and falling back to the fitted model's family. Used to gate
# the Median Odds Ratio, which is defined only for the logit link.
maihda_model_link_name <- function(model) {
  fam <- model$family
  if (is.list(fam) && !is.null(fam$link)) {
    return(fam$link)
  }
  ff <- tryCatch(maihda_family(model$model), error = function(e) NULL)
  if (!is.null(ff) && !is.null(ff$link)) {
    return(ff$link)
  }
  NA_character_
}

# Observed 0/1 response aligned with predict_maihda()'s individual predictions.
# For lme4, getME(, "y") returns the numeric 0/1 response used in fitting (the
# approach used in the binary_outcomes vignette). For other engines, fall back to
# the response column of the model frame, coerced to 0/1.
maihda_da_observed_response <- function(model) {
  if (identical(model$engine, "lme4")) {
    return(as.numeric(lme4::getME(model$model, "y")))
  }

  resp <- all.vars(model$formula)[1]
  y <- model$data[[resp]]
  if (is.logical(y)) {
    return(as.integer(y))
  }
  if (is.factor(y)) {
    return(as.integer(y) - 1L)
  }
  as.numeric(y)
}

# Per-row success / trial counts for an aggregated-binomial lme4 fit, or NULL when
# the fit is not aggregated. An lme4 cbind(success, failure) response is stored in the
# fitted model frame (model$data) as a two-column [successes, failures] matrix, so it
# carries a `dim` regardless of the per-row proportions -- the STRUCTURAL signal of
# aggregation, used the same way maihda_prior_weights() recognises it. Reading the
# counts from the matrix is exact (no proportion x trials rounding). A fallback
# recovers the trial counts from the prior weights (glmer stores cbind() trials there)
# when the model frame did not retain the matrix response; it is gated on the fit
# being unweighted so a design/precision-weighted Bernoulli fit -- whose prior weights
# are the sampling weights, not trial counts -- is never misread as aggregated.
maihda_da_aggregated_counts <- function(model) {
  resp <- tryCatch(stats::model.response(model$data), error = function(e) NULL)
  if (!is.null(resp) && !is.null(dim(resp)) && ncol(resp) == 2L) {
    successes <- as.numeric(resp[, 1])
    trials <- successes + as.numeric(resp[, 2])
    if (all(is.finite(trials))) {
      return(list(successes = successes, trials = trials))
    }
  }
  # Fallback: a genuine Bernoulli fit has unit prior weights, so INTEGER trial
  # counts > 1 anywhere mark an aggregated fit; pair them with the success
  # proportions from getME(, "y"). Non-integer prior weights are NOT trial
  # counts -- they are precision weights on a Bernoulli fit, and rounding
  # y * weight into pseudo-counts against fractional trials produced negative
  # implied failures (round(1.5) = 2 successes out of 1.5 trials) and an AUC
  # above 1; such fits return NULL here and take the precision-weighted
  # concordance path in maihda_discriminatory_accuracy(). Skipped when sampling
  # weights are in play (their prior weights are not trial counts either).
  if (!is.null(model$sampling_weights)) {
    return(NULL)
  }
  y <- tryCatch(as.numeric(lme4::getME(model$model, "y")), error = function(e) NULL)
  w <- tryCatch(as.numeric(stats::weights(model$model, type = "prior")),
                error = function(e) NULL)
  if (!is.null(y) && !is.null(w) && length(y) == length(w) &&
      all(is.finite(w)) && any(w > 1) &&
      all(abs(w - round(w)) < 1e-8)) {
    return(list(successes = round(y * w), trials = round(w)))
  }
  NULL
}

# Per-row success / trial counts for a brms `y | trials(n)` aggregated-binomial fit,
# or NULL when the fit is not a brms aggregated binomial. The lme4 counterpart above
# reads trials from the matrix response / prior weights; brms exposes no
# weights.brmsfit, so the trials come from maihda_brms_trial_counts() (the same path
# the prediction weighting uses), which parses the trials() addition term off the
# stored formula. The successes are the response column -- the per-row success counts an
# aggregated binomial models. NULL for a brms Bernoulli fit (no trials() term), which
# then takes the ordinary rank-based AUC path.
maihda_da_brms_aggregated_counts <- function(model) {
  trials <- maihda_brms_trial_counts(model)
  if (is.null(trials)) {
    return(NULL)
  }
  successes <- tryCatch(maihda_da_observed_response(model), error = function(e) NULL)
  if (is.null(successes) || length(successes) != length(trials) ||
      !all(is.finite(successes)) || !all(is.finite(trials))) {
    return(NULL)
  }
  list(successes = as.numeric(successes), trials = as.numeric(trials))
}

# Split the model's random-effect groupings into the intersectional partition --
# the stratum interaction plus, for a crossed-dimensions fit ($cc_info), the
# additive dimension REs -- and everything else (a contextual (1 | school) from
# fit_maihda(context = ), or an explicit extra grouping like (1 | site)).
# Grouping names are compared in deparsed form, so a non-syntactic dimension
# name (which deparses backtick-quoted) matches either way.
maihda_da_re_scopes <- function(model) {
  bars <- reformulas::findbars(model$formula)
  groups <- unique(vapply(bars, function(b) {
    paste(deparse(b[[3]]), collapse = "")
  }, character(1)))
  intersectional <- unique(c("stratum",
                             model$cc_info$interaction_group,
                             unname(model$cc_info$dim_groups)))
  quoted <- vapply(intersectional, maihda_quote_name, character(1))
  keep <- groups %in% c(intersectional, quoted)
  list(intersectional = groups[keep], other = groups[!keep])
}

# Link-scale scores from the fixed effects plus ONLY the given random-effect
# groupings: the ranking the intersectional-scope AUC is computed over. The
# link is monotone, so ranking by this linear predictor equals ranking by the
# response-scale probability with the other random effects excluded (for an
# aggregated-binomial fit it ranks by the per-trial probability, since the
# trials multiplier is not part of the linear predictor). lme4/brms only -- the
# wemix/ordinal engines enforce the canonical single (1 | stratum) structure,
# so they never carry another grouping to exclude.
maihda_da_scope_scores <- function(model, keep_groups) {
  bars <- reformulas::findbars(model$formula)
  keep <- vapply(bars, function(b) {
    paste(deparse(b[[3]]), collapse = "") %in% keep_groups
  }, logical(1))
  re_txt <- vapply(bars[keep], function(b) {
    paste0("(", paste(deparse(b), collapse = " "), ")")
  }, character(1))
  re_form <- stats::as.formula(paste("~", paste(re_txt, collapse = " + ")))
  if (identical(model$engine, "lme4")) {
    as.numeric(stats::predict(model$model, re.form = re_form, type = "link"))
  } else if (identical(model$engine, "brms")) {
    maihda_brms_linpred_mean(model$model, re_formula = re_form)
  } else {
    stop("Internal error: intersectional-scope scores requested for engine '",
         model$engine, "'.", call. = FALSE)
  }
}

# Count-weighted AUC / C-statistic for an aggregated-binomial fit. Each row i carries
# a shared predicted probability prob_i with successes_i observed cases and
# (trials_i - successes_i) controls. This equals the Mann-Whitney AUC of the expanded
# individual-level 0/1 data -- P(case score > control score), ties counted as one
# half -- computed by grouping cases/controls at each distinct probability level
# rather than materialising the expansion.
maihda_auc_weighted <- function(prob, successes, trials) {
  failures <- trials - successes
  keep <- is.finite(prob) & is.finite(successes) & is.finite(failures) &
    (successes + failures) > 0
  prob <- prob[keep]
  successes <- successes[keep]
  failures <- failures[keep]

  # Negative mass in either class breaks the concordance bounds (the AUC can
  # exceed 1); every caller must pass successes in [0, trials].
  if (any(successes < 0) || any(failures < 0)) {
    stop("Internal error: negative case/control mass in the weighted AUC ",
         "(successes must lie in [0, trials]).", call. = FALSE)
  }

  # as.double guards against 32-bit integer overflow in n1 * n0 below when the
  # counts are large integers (see maihda_auc); a count-weighted call already
  # passes doubles, but an integer-count caller would otherwise silently get NA.
  n1 <- as.double(sum(successes))
  n0 <- as.double(sum(failures))
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  ord <- order(prob)
  prob <- prob[ord]
  successes <- successes[ord]
  failures <- failures[ord]

  # Group rows sharing a probability so ties (all individuals in a row, and rows with
  # equal fitted probabilities) are handled together.
  level <- cumsum(c(TRUE, diff(prob) > 0))
  c_k <- as.numeric(tapply(successes, level, sum))
  d_k <- as.numeric(tapply(failures, level, sum))

  # Controls strictly below each probability level (concordant with a case there),
  # plus half the same-level controls (ties).
  controls_below <- cumsum(c(0, d_k[-length(d_k)]))
  concordant <- sum(c_k * controls_below) + 0.5 * sum(c_k * d_k)
  concordant / (n1 * n0)
}
