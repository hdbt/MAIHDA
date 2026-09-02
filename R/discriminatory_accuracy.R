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
#' \strong{Scope.} \code{V_A} is the \emph{between-stratum} variance. Contextual
#' (\code{context = }) and other non-stratum random effects are never included:
#' the MOR compares two individuals from randomly chosen strata \emph{within
#' the same context}.
#'
#' \strong{Crossed-dimensions fits use a different calculation.} The closed form
#' above assumes the two strata's random effects are \emph{independent}, which is
#' what makes their difference \eqn{N(0, 2 V_A)}. In a crossed-dimensions fit
#' (from \code{maihda(decomposition = "crossed-dimensions")}) a stratum's effect
#' is the sum of its dimension effects plus the intersection effect, so two
#' strata sharing a dimension -- say two strata that are both "female" -- share
#' that dimension's random effect and are \strong{correlated}. Their difference
#' is then a \emph{mixture} of normals, one component per pattern of shared
#' dimensions, and applying the closed form to the summed variance overstates the
#' MOR (substantially so when the variance sits mainly in the additive
#' dimensions, and not at all when it sits entirely in the interaction).
#'
#' The MOR reported for such a fit is therefore computed from that mixture, under
#' an explicit sampling scheme: \strong{two distinct strata drawn uniformly at
#' random from the strata present in the analytic sample}. Writing
#' \eqn{\tau^2_d} for dimension \eqn{d}'s variance and \eqn{\tau^2_I} for the
#' intersection variance, a pair differing on the dimension set \eqn{D^*} has
#' difference variance \eqn{v = 2(\tau^2_I + \sum_{d \in D^*} \tau^2_d)}, and the
#' MOR is \code{exp(x)} for the \code{x} solving
#' \eqn{\sum_{pairs} (2\Phi(x/\sqrt{v}) - 1) / n_{pairs} = 0.5}. For a canonical
#' single-stratum fit the two calculations coincide, and that closed form is
#' used.
#'
#' @param model A \code{maihda_model} from \code{\link{fit_maihda}} fitted with a
#'   \code{binomial} (lme4), \code{bernoulli} (brms), or \code{cumulative}
#'   (ordinal) family and a \strong{logit} link.
#'
#' @return A single number (the MOR, \eqn{\ge 1}), or \code{NA_real_} if the
#'   between-stratum variance is unavailable -- which for a crossed-dimensions fit
#'   also covers the case where the stratum grid needed for the mixture cannot be
#'   resolved (fewer than two strata, absent dimension columns, or more than 12
#'   dimensions).
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

  # A crossed-dimensions fit's strata are NOT independent draws -- two strata that
  # share a dimension share that dimension's random effect -- so the closed form
  # below does not apply to them; their difference is a mixture (see Details).
  if (!is.null(model$cc_info)) {
    return(maihda_mor_crossed(model))
  }

  v_a <- tryCatch(maihda_mor_between_variance(model), error = function(e) NA_real_)
  if (!is.numeric(v_a) || length(v_a) != 1 || !is.finite(v_a) || v_a < 0) {
    return(NA_real_)
  }

  exp(sqrt(2 * v_a) * stats::qnorm(0.75))
}

# MOR of a crossed-dimensions fit: the median |difference| between two DISTINCT
# strata drawn uniformly at random from those present in the analytic sample.
#
# Stratum (a1, ..., aD) carries the effect sum_d u_d[a_d] + w[cell], every term
# independent. Two distinct strata therefore differ by
#   sum_{d : a_d != a'_d} (u_d[a_d] - u_d[a'_d])  +  (w[cell] - w[cell']),
# whose variance is 2 * (tau_I^2 + sum over the DIFFERING dimensions of tau_d^2):
# a dimension the two strata share contributes NOTHING, because its random effect
# cancels. Summing every variance for every pair -- the independent-strata closed
# form -- charges the shared dimensions to all pairs and overstates the MOR.
#
# The pattern of differing dimensions varies across pairs, so |difference| follows
# a MIXTURE of half-normals and its median has no closed form; solve the mixture
# CDF for 0.5 instead. Returns NA_real_ when the variances or the stratum grid
# cannot be resolved.
maihda_mor_crossed <- function(model) {
  parts <- tryCatch(maihda_mor_crossed_parts(model), error = function(e) NULL)
  if (is.null(parts)) {
    return(NA_real_)
  }
  maihda_mor_crossed_from_parts(parts)
}

# Solve the pair-mixture for the crossed MOR, given the parts (per-pattern contrast
# variance and pair weights) from maihda_mor_crossed_parts(). Split out as a pure
# function so the mixture-median logic -- in particular the zero-variance-mass
# boundary -- can be tested directly without a fitted model. Returns the MOR, or
# NA_real_ when the median cannot be resolved.
maihda_mor_crossed_from_parts <- function(parts) {
  v_pair <- 2 * (parts$interaction + as.numeric(parts$pattern %*% parts$dims))
  keep <- is.finite(v_pair) & v_pair > 0
  # Every pair has variance 0 (all components at the boundary): no heterogeneity,
  # so the two strata never differ and the median odds ratio is exactly 1.
  if (!any(keep)) {
    return(1)
  }
  w <- parts$weight
  w <- w / sum(w)
  # The zero-variance pairs put a point mass at 0. When that mass alone is at
  # least half the total weight, the median |difference| is exactly 0 (the CDF
  # already reaches 0.5 at x = 0), so the MOR is exp(0) = 1. Handle it here: the
  # root search below cannot bracket this root -- both endpoints of [eps, hi]
  # give cdf >= 0.5, so uniroot() errors ("values at end points not of opposite
  # sign") and the crossed MOR would be reported NA instead of 1.
  if (sum(w[!keep]) >= 0.5) {
    return(1)
  }
  # P(|difference| <= x) marginalised over the pairs. Pairs whose variance is 0
  # contribute a point mass at 0, i.e. probability 1 for every x > 0.
  cdf <- function(x) {
    sum(w[keep] * (2 * stats::pnorm(x / sqrt(v_pair[keep])) - 1)) + sum(w[!keep])
  }
  hi <- 10 * sqrt(max(v_pair[keep]))
  if (cdf(hi) < 0.5) {
    return(NA_real_)
  }
  root <- tryCatch(
    stats::uniroot(function(x) cdf(x) - 0.5, c(.Machine$double.eps, hi),
                   tol = .Machine$double.eps^0.5)$root,
    error = function(e) NA_real_)
  if (!is.finite(root)) {
    return(NA_real_)
  }
  exp(root)
}

# The pieces maihda_mor_crossed() needs: the per-dimension variances, the
# interaction variance, and -- over the DISTINCT unordered pairs of observed
# strata -- which dimensions differ.
#
# Pairs are grouped by their differing-dimension PATTERN rather than enumerated,
# so the cost is O(2^D * n_strata) instead of O(n_strata^2) -- which matters
# because a MAIHDA has few dimensions but can have many levels within them.
# `pattern` is one row per distinct pattern and `weight` how many pairs share it.
# Counting uses inclusion-exclusion over the subset lattice: the number of ordered
# pairs matching on every dimension in S is a sum of squared group sizes, and the
# exact-agreement counts follow by Mobius inversion.
maihda_mor_crossed_parts <- function(model) {
  cc <- model$cc_info
  dim_groups <- unname(cc$dim_groups)
  if (length(dim_groups) < 1L || is.null(cc$interaction_group)) {
    return(NULL)
  }
  v <- maihda_cc_variances(model)
  tau_dim <- unname(v$dims)
  tau_int <- v$interaction
  if (!all(is.finite(c(tau_dim, tau_int))) || any(c(tau_dim, tau_int) < 0)) {
    return(NULL)
  }

  dat <- model$data
  if (is.null(dat) || !all(dim_groups %in% names(dat))) {
    return(NULL)
  }
  cells <- unique(dat[stats::complete.cases(dat[, dim_groups, drop = FALSE]),
                      dim_groups, drop = FALSE])
  n <- nrow(cells)
  if (n < 2L) {
    return(NULL)
  }
  D <- length(dim_groups)
  # One grouping pass per subset, so the work is 2^D tabulations. A MAIHDA has a
  # handful of dimensions -- 12 already implies at least 4096 strata -- so cap it
  # rather than grind, and let the caller report NA.
  if (D > 12L) {
    return(NULL)
  }

  subsets <- seq_len(2^D) - 1L
  bits <- bitwShiftL(1L, seq_len(D) - 1L)
  # counts[S] = ordered pairs (i, j), i == j included, agreeing on every dimension
  # in S. Each is a sum of squared group sizes, hence an exact integer.
  counts <- vapply(subsets, function(s) {
    dims <- which(bitwAnd(s, bits) > 0L)
    if (length(dims) == 0L) {
      return(n^2)
    }
    key <- do.call(paste, c(lapply(cells[, dims, drop = FALSE], as.character),
                            sep = "\r"))
    sum(as.numeric(table(key))^2)
  }, numeric(1))
  # Superset Mobius inversion, in place: counts[S] becomes the number of pairs
  # agreeing on EXACTLY S. One vectorised sweep per dimension (D * 2^D) rather
  # than the 4^D of an explicit alternating sum over every superset.
  for (b in bits) {
    lo <- subsets[bitwAnd(subsets, b) == 0L]
    counts[lo + 1L] <- counts[lo + 1L] - counts[bitwOr(lo, b) + 1L]
  }
  agree_exactly <- counts

  full <- 2^D - 1L
  # Drop S = all dimensions: those "pairs" are a stratum with itself.
  keep <- subsets != full & agree_exactly > 0.5
  if (!any(keep)) {
    return(NULL)
  }
  # One row per surviving pattern, flagging the dimensions the pair DIFFERS on.
  pattern <- t(vapply(subsets[keep], function(s) as.numeric(bitwAnd(s, bits) == 0L),
                      numeric(D)))
  # vapply drops to a vector when D == 1; restore the matrix shape.
  if (D == 1L) {
    pattern <- matrix(pattern, ncol = 1L)
  }
  list(dims = tau_dim, interaction = tau_int, pattern = pattern,
       weight = agree_exactly[keep] / 2)   # ordered -> unordered pairs
}

# The intersectional random-effect variances of a crossed-dimensions fit: one per
# additive dimension (in cc_info$dim_groups order) plus the interaction. Shared by
# the total-variance accessor below and the MOR mixture above, so the two cannot
# drift apart on engine support or on which effects count as intersectional.
maihda_cc_variances <- function(model) {
  cc <- model$cc_info
  dim_groups <- unname(cc$dim_groups)
  inter_group <- cc$interaction_group
  var_named <- if (identical(model$engine, "lme4")) {
    maihda_random_variances_lme4(model$model)
  } else if (identical(model$engine, "brms")) {
    maihda_random_variances_brms(model$model)
  } else {
    stop("Crossed-dimensions MOR is supported for the lme4 and brms engines ",
         "only.", call. = FALSE)
  }
  missing_re <- setdiff(unique(c(inter_group, dim_groups)), names(var_named))
  if (length(missing_re) > 0) {
    stop("Crossed-dimensions MOR is missing the random effect(s): ",
         paste(missing_re, collapse = ", "), ".", call. = FALSE)
  }
  list(dims = stats::setNames(as.numeric(var_named[dim_groups]), dim_groups),
       interaction = as.numeric(var_named[inter_group]))
}

# TOTAL between-stratum (intersectional) variance. For the canonical fit this is
# the stratum intercept variance. For a crossed-dimensions fit ($cc_info) the
# between-stratum effect is the SUM of the independent crossed REs -- each
# additive dimension effect plus the interaction -- so its variance is their
# total; reading only the "stratum" (interaction) component would understate the
# between-stratum heterogeneity. Contextual and other non-intersectional REs are
# never included.
#
# NOTE this total is the variance of ONE stratum's effect. It is NOT half the
# variance of the DIFFERENCE between two strata, because crossed strata sharing a
# dimension are correlated -- which is why maihda_mor() routes a crossed fit
# through maihda_mor_crossed() instead of through the 2 * V_A closed form.
maihda_mor_between_variance <- function(model) {
  if (is.null(model$cc_info)) {
    return(extract_between_variance(model))
  }
  v <- maihda_cc_variances(model)
  sum(c(v$dims, v$interaction))
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
#' Aggregated-binomial fits are supported in all three forms R accepts -- an lme4
#' \code{cbind(success, failure)} response, an lme4 response in [0, 1] whose trial
#' counts are supplied as \code{weights =} (see \code{\link[stats]{glm}}), and a brms
#' \code{y | trials(n)} response: the AUC is the count-weighted C-statistic over the
#' implied individual-level 0/1 data, and \code{n_case} / \code{n_control} are the
#' total successes / failures. All three spellings of the same data give the same AUC.
#' This includes individual 0/1 records collapsed to frequency cells (every row
#' all-success or all-failure), whose weights are trial counts like any others; see
#' \code{binomial_weights} to override that reading.
#'
#' \strong{Scope.} When the model carries random effects \emph{beyond} the
#' intersectional partition -- a contextual \code{(1 | school)} from
#' \code{fit_maihda(context = )} or an explicit extra grouping such as
#' \code{(1 | site)} -- the headline \code{auc} is the \emph{intersectional-scope}
#' concordance: it \strong{excludes} those other random effects but keeps the fixed
#' effects plus the stratum random effect (and, for a crossed-dimensions fit, the
#' additive dimension effects). The concordance of the full model including the other
#' random effects is reported separately as \code{auc_full}
#' (\code{auc_scope = "intersectional"}). \strong{Caveat -- this is not strata-only
#' discrimination.} The intersectional-scope score retains the \emph{entire}
#' fixed-effects predictor, so when the model is adjusted for individual-level
#' covariates (e.g. \code{age}, \code{income} that vary \emph{within} strata) those
#' covariates enter this AUC too; it is the concordance of the adjusted fixed effects
#' plus the intersectional random effect(s), not of the strata alone, and it then
#' matches the between-stratum MOR's scope only when the fixed part is intercept-only.
#' For a strata-only discriminatory accuracy, score the null (strata-only) model. For
#' the canonical single-\code{(1 | stratum)} model with no other random effects, the
#' full and intersectional scopes coincide and only \code{auc} is reported
#' (\code{auc_scope = "model"}, \code{auc_full} absent).
#'
#' @param model A \code{maihda_model} from \code{\link{fit_maihda}} fitted with a
#'   \code{binomial} family -- including an aggregated response (an lme4
#'   \code{cbind(success, failure)} or a brms \code{y | trials(n)}) -- or the
#'   \code{bernoulli} family a binary 0/1 outcome is fit with under
#'   \code{engine = "brms"}.
#' @param binomial_weights How to read non-unit \code{weights =} on an \strong{lme4}
#'   binomial fit when computing the AUC. \code{"auto"} (default) reads integral
#'   weights as \strong{trial counts}, which is what R documents a binomial prior
#'   weight to be (\code{\link[stats]{glm}}: "For a binomial GLM prior weights are used
#'   to give the number of trials"), and leaves non-integral weights -- which cannot be
#'   counts -- on the observation-level path, flagged
#'   \code{precision_weights_ignored}. \code{"trials"} forces the trial-count reading
#'   even for non-integral weights (the fractional case/control mass is carried as it
#'   stands, with a warning). \code{"analytic"} forces the ordinary observation-level
#'   concordance, ignoring the weights and setting
#'   \code{precision_weights_ignored = TRUE}; it is an error for a response carrying
#'   values strictly between 0 and 1, which has no observation-level case/control
#'   reading. Has no effect on a \code{cbind(success, failure)} response (its
#'   denominator is structural), on a brms \code{y | trials(n)} fit, or on a
#'   design-weighted (\code{sampling_weights}) fit.
#'
#' @return An object of class \code{maihda_da}: a list with \code{auc},
#'   \code{auc_scope}, \code{auc_full}, \code{mor},
#'   \code{n_case}, \code{n_control}, \code{family}, \code{link}, \code{engine},
#'   \code{weighted}, \code{weight_type}, \code{precision_weights_ignored} and
#'   \code{apparent} (always \code{TRUE} -- the AUC is in-sample; see Description).
#'   \code{mor} is \code{NA} for a non-logit binomial link, where the AUC is still
#'   reported. For an aggregated-binomial fit \code{n_case} / \code{n_control} are the
#'   total successes / failures. \code{weighted} is \code{TRUE} only when the AUC is a
#'   genuinely weighted (population-mass) concordance, with \code{weight_type}
#'   \code{"sampling"} for a design-weighted fit (each observation contributes its
#'   sampling weight as case/control mass, estimating the population discriminatory
#'   accuracy); \code{n_case} / \code{n_control} stay unweighted observation counts.
#'   \code{weight_type} is \code{NULL} for an unweighted AUC. An
#'   \strong{aggregated-binomial} fit is reported \code{weighted = FALSE} with
#'   \code{weight_type = NULL}: its count-weighted AUC equals the ordinary
#'   individual-level concordance over the implied 0/1 data (the trial counts are
#'   real observations, not sampling weights), so it is not a design-weighted
#'   population quantity. Aggregation is recognised from an lme4
#'   \code{cbind(success, failure)} matrix response, from non-unit integral lme4
#'   \code{weights=} on a response in [0, 1] (the trial counts, see
#'   \code{\link[stats]{glm}}), or from a brms \code{y | trials(n)} term. The
#'   weights-based form covers a response of exactly 0/1 -- individual records
#'   collapsed to frequency cells -- because for the binomial family a prior weight
#'   \emph{is} a trial count: the weighted log-likelihood of one 0/1 row is exactly
#'   that of \code{w} trials sharing its outcome, so the fit equals the row-expanded
#'   one and its AUC should too. Pass \code{binomial_weights} to override. When a
#'   response times its trial counts is \emph{not} a whole number of successes --
#'   a malformed binomial, which \code{glmer} itself warns about as "non-integer
#'   #successes" -- the fractional case/control mass is carried into the AUC as it
#'   stands, with a warning, rather than rounded into whole observations that were
#'   never collected; \code{n_case} / \code{n_control} are then fractional.
#'   \strong{Weights not read as trial counts} -- non-integral \code{weights=}, or any
#'   weights under \code{binomial_weights = "analytic"} -- are ignored by the AUC,
#'   which is then the ordinary observation-level concordance
#'   (\code{weighted = FALSE}), with \code{precision_weights_ignored = TRUE} flagging
#'   that such weights were present.
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
maihda_discriminatory_accuracy <- function(model,
                                           binomial_weights = c("auto", "trials",
                                                                "analytic")) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  binomial_weights <- match.arg(binomial_weights)
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
  # non-0/1 response to maihda_auc() (which errors). An lme4 fit whose trial counts are
  # supplied as `weights =` has a plain vector response, so it is recognised from the
  # weights instead -- see maihda_da_aggregated_counts(), which also documents why a
  # binomial prior weight IS a trial count. agg_counts is NULL for a Bernoulli fit with
  # unit (or non-integral) weights, which takes the ordinary rank-based path, and is
  # aligned to prob's rows.
  agg_counts <- if (identical(model$engine, "lme4")) {
    maihda_da_aggregated_counts(model, binomial_weights)
  } else if (identical(model$engine, "brms")) {
    maihda_da_brms_aggregated_counts(model)
  } else {
    NULL
  }
  aggregated <- !is.null(agg_counts) && length(agg_counts$successes) == length(prob)
  # binomial_weights = "analytic" on a genuine PROPORTION response is incoherent: an
  # observation-level concordance needs a 0/1 outcome, and a proportion row only has a
  # case/control interpretation through its trial counts. Say that, rather than let
  # maihda_auc() report a generic "must be a binary 0/1 outcome".
  if (identical(binomial_weights, "analytic") && !aggregated &&
      isTRUE(any(resp > 0 & resp < 1, na.rm = TRUE))) {
    stop("binomial_weights = \"analytic\" is not available for a proportion ",
         "response: the response carries values strictly between 0 and 1, which have ",
         "a case/control interpretation only through their trial counts. Use ",
         "binomial_weights = \"trials\" (or \"auto\").", call. = FALSE)
  }
  # Design-weighted fit (sampling_weights supplied): compute the design-weighted
  # AUC, where each observation contributes its sampling weight as case (y = 1) or
  # control (y = 0) mass -- the weighted Mann-Whitney concordance, estimating the
  # POPULATION discriminatory accuracy rather than the sample's. The reported
  # case/control totals stay unweighted observation counts.
  sw <- if (!is.null(model$sampling_weights)) maihda_prior_weights(model) else NULL
  design_weighted <- !is.null(sw) && length(sw) == length(prob) &&
    any(is.finite(sw)) && !isTRUE(all(abs(sw - 1) < sqrt(.Machine$double.eps)))
  # lme4 non-unit prior weights that were NOT read as trial counts: either they are
  # non-integral (so they cannot be counts) or the caller passed
  # binomial_weights = "analytic". The AUC is then the ordinary observation-level
  # concordance, which does not fold the weights into case/control mass, and `pw`
  # flags that (precision_weights_ignored) in the result and print output so the
  # number is not mistaken for the trial-count one. A fit whose weights WERE read as
  # trial counts is `aggregated` and never reaches here.
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
  } else {
    # Ordinary observation-level AUC. This covers the unweighted case and the case
    # where non-unit weights were NOT read as trial counts (pw non-NULL): non-integral
    # weights, or binomial_weights = "analytic". Those weights are not folded into
    # case/control mass, so the fit reports the same observation-level concordance as
    # an unweighted one, flagged precision_weights_ignored.
    auc_for <- function(score) maihda_auc(score, resp)
    n_case <- sum(resp == 1, na.rm = TRUE)
    n_control <- sum(resp == 0, na.rm = TRUE)
  }

  # Scope. With non-intersectional random effects in the model (a contextual
  # (1 | school) or an explicit (1 | site)), the full-model prediction folds
  # their effects into the score, which a between-STRATUM summary should not: a
  # strong site effect can carry a high full-model AUC over a negligible stratum
  # effect. The headline AUC is therefore the intersectional-SCOPE concordance --
  # it excludes the other random effects but keeps the fixed effects plus the
  # stratum RE (and, for a crossed-dimensions fit, the additive dimension REs) --
  # with the full-model AUC reported alongside as auc_full. NOTE the scope score
  # (maihda_da_scope_scores()) retains the WHOLE fixed-effects predictor, so when
  # the model is adjusted for individual-level covariates this AUC includes them
  # too: it is an adjusted intersectional concordance, NOT strata-only, and it
  # matches the MOR's between-stratum scope only when the fixed part is
  # intercept-only. Labelled auc_scope = "intersectional" (not "strata") to avoid
  # over-claiming. For the canonical single-(1 | stratum) model the full and
  # intersectional scopes coincide and only auc is reported.
  scopes <- maihda_da_re_scopes(model)
  auc_full <- NULL
  auc_scope <- "model"
  if (length(scopes$other) > 0) {
    score <- tryCatch(maihda_da_scope_scores(model, scopes$intersectional),
                      error = function(e) NULL)
    if (!is.null(score) && length(score) == length(prob)) {
      auc_full <- auc_for(prob)
      auc <- auc_for(score)
      auc_scope <- "intersectional"
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
      # `weighted` marks a design/sampling-weighted (population-representative)
      # concordance: one that folds sampling weights into case/control mass to
      # estimate a POPULATION AUC that differs from the naive sample concordance.
      # An aggregated trial-count AUC is deliberately NOT flagged here: its trial
      # counts are real individual observations, so the count-weighted concordance
      # EQUALS the ordinary observation-level concordance over the implied 0/1 data
      # (see the aggregated branch above and its tests), reported with the true
      # success/failure totals as n_case/n_control -- exactly analogous to the
      # precision-weighted case, which is likewise the ordinary observation-level
      # concordance (precision_weights_ignored flags it). So both the aggregated and
      # the precision-weighted fit report weighted = FALSE and weight_type = NULL;
      # only a genuine sampling-weighted fit reports weighted = TRUE. (Note that
      # print.maihda_da() renders weighted = TRUE as "design-weighted ... sampling
      # weight", which would misdescribe an aggregated cbind()/trials() fit.)
      weighted = design_weighted,
      weight_type = if (design_weighted) "sampling" else NULL,
      precision_weights_ignored = !is.null(pw),
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
              if (identical(x$auc_scope, "intersectional")) " [intersectional scope]" else ""))
  if (!is.null(x$auc_full) && is.finite(x$auc_full)) {
    cat(sprintf("  AUC, full model:   %s\n",
                pal$accent(sprintf("%.3f", x$auc_full))))
    cat(pal$muted(paste0(
      "  (The model carries non-stratum random effects. The headline AUC excludes\n",
      "  those (e.g. a contextual/site effect); it is the concordance of the fixed\n",
      "  effects plus the intersectional random effect(s), so it includes any\n",
      "  individual-level covariates in the model and is strata-only only when the\n",
      "  fixed part is intercept-only. The full-model AUC adds the other random effects.)\n")))
  }
  mor_str <- if (is.finite(x$mor)) {
    pal$accent(sprintf("%.3f", x$mor))
  } else if (!is.null(x$link) && !identical(x$link, "logit")) {
    sprintf("NA (requires the logit link; model uses '%s')", x$link)
  } else {
    "NA"
  }
  cat(sprintf("  Median Odds Ratio: %s\n", mor_str))
  # Not "%d": an aggregated-binomial fit whose proportion response does not resolve
  # to whole successes reports FRACTIONAL case/control mass (see
  # maihda_da_proportion_successes()), and sprintf("%d", 1.5) is an error.
  cat(sprintf("  Cases / controls:  %s / %s\n",
              maihda_format_mass(x$n_case), maihda_format_mass(x$n_control)))
  if (!isFALSE(x$apparent)) {
    cat(pal$muted(paste0(
      "  (AUC is apparent / in-sample: scored on the same rows used to fit the\n",
      "  model, so it is optimistically biased -- more so with sparse strata. It\n",
      "  is a descriptive measure, not cross-validated out-of-sample discrimination.)\n")))
  }
  if (isTRUE(x$weighted)) {
    cat(pal$muted(paste0(
      "  (AUC is design-weighted: each observation contributes its sampling\n",
      "  weight; cases/controls are unweighted counts.)\n")))
  }
  if (isTRUE(x$precision_weights_ignored)) {
    cat(pal$muted(paste0(
      "  (The fit carries non-unit lme4 weights that were not read as trial counts\n",
      "  (they are not whole numbers, or binomial_weights = \"analytic\"); the AUC\n",
      "  ignores them, so it is the ordinary observation-level concordance over the\n",
      "  model's rows, not over the trials those weights would imply.)\n")))
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
# counts from the matrix is exact (no proportion x trials rounding).
#
# A fallback covers R's OTHER aggregated-binomial idiom -- a PROPORTION response with
# the trial counts supplied as prior weights (?glm: "For a binomial GLM prior weights
# are used to give the number of trials when the response is the proportion of
# successes"), i.e. successes/trials ~ ... , weights = trials. That response is a plain
# vector, so the structural test above cannot see it.
#
# The fallback keys on the WEIGHTS: non-unit, integral prior weights on a binomial
# fit are trial counts, which is what R itself documents them to mean. This is not a
# heuristic about the response -- for the binomial family a prior weight IS a trial
# count, because the weighted log-likelihood w * [y log p + (1 - y) log(1 - p)] is
# exactly the log-likelihood of w trials all sharing that outcome. Verified rather
# than assumed: a 600-row Bernoulli fit with weights 1:5 and the row-EXPANDED fit of
# the same data agree to 8e-14 in the fixed effects, to 13 digits in tau^2, and to
# 3e-13 in the fitted probabilities. There is no separate "precision weight" model to
# preserve -- unlike a Gaussian fit, a binomial GLM has no dispersion parameter for a
# weight to rescale, so a weight of w on a 0/1 row is w observations or it is nothing.
#
# An earlier revision keyed on the RESPONSE instead -- aggregated only if some value
# lay strictly inside (0, 1), which a Bernoulli 0/1 response never does -- to stop an
# estimand turning on the integrality of a fitting control. That test is wrong on the
# single most common shape of frequency-weighted binary data: individual 0/1 records
# collapsed to (covariate pattern x outcome) cells with a count, where EVERY row is
# all-success or all-failure by construction. Such a fit took the observation-level
# path, and because each cell contributes one case and one control at an identical
# score, every pair tied and the AUC was EXACTLY 0.5 -- no information at all, printed
# with a confident note explaining why it was right. Measured on 3000 individuals in
# 12 strata: cbind(), the coarse proportion+weights spelling, and the individual 0/1
# data all give AUC 0.786929 with 1427 / 1573 cases / controls, while the same
# individuals collapsed one level finer gave 0.500000 with 24 / 24 -- from a fit whose
# fixed effects match to 4e-10. So the response test made the AUC depend on collapse
# granularity, which is not a property of the model at all.
#
# Integrality remains required under "auto" (below): a non-integral weight cannot be a
# trial count, so those fits keep the observation-level path and are flagged
# precision_weights_ignored. `binomial_weights` lets the caller override both ways.
#
# Skipped when sampling weights are in play (their prior weights are not trial counts
# either). A non-aggregated fit returns NULL and takes the ordinary observation-level
# AUC path in maihda_discriminatory_accuracy(), where the weights are ignored and
# flagged precision_weights_ignored.
maihda_da_aggregated_counts <- function(model, binomial_weights = "auto") {
  # The cbind(success, failure) matrix response is STRUCTURAL evidence of aggregation
  # and carries its own denominator, so `binomial_weights` -- which only ever governs
  # how PRIOR WEIGHTS are read -- does not apply to it.
  resp <- tryCatch(stats::model.response(model$data), error = function(e) NULL)
  if (!is.null(resp) && !is.null(dim(resp)) && ncol(resp) == 2L) {
    successes <- as.numeric(resp[, 1])
    trials <- successes + as.numeric(resp[, 2])
    if (all(is.finite(trials))) {
      return(list(successes = successes, trials = trials))
    }
  }
  if (identical(binomial_weights, "analytic")) {
    return(NULL)
  }
  if (!is.null(model$sampling_weights)) {
    return(NULL)
  }
  y <- tryCatch(as.numeric(lme4::getME(model$model, "y")), error = function(e) NULL)
  w <- tryCatch(as.numeric(stats::weights(model$model, type = "prior")),
                error = function(e) NULL)
  if (is.null(y) || is.null(w) || length(y) != length(w) ||
      !all(is.finite(y)) || !all(is.finite(w)) || !all(w > 0)) {
    return(NULL)
  }
  # `all(y >= 0 & y <= 1)` is defensive: glmer's binomial initialize already refuses a
  # response outside [0, 1], but fractional mass built from such a row could exceed the
  # trials and trip maihda_auc_weighted()'s negative-mass guard.
  if (!all(y >= 0 & y <= 1)) {
    return(NULL)
  }
  forced <- identical(binomial_weights, "trials")
  # `any(w > 1)` keeps an ALL-UNIT-weight fit out under "auto". For a 0/1 response
  # that is an ordinary Bernoulli fit; for a proportion response it is a malformed
  # binomial (glmer warns "non-integer #successes"), and reading its proportions as
  # successes out of 1 trial would invent data -- it falls through to maihda_auc(),
  # which rejects the non-0/1 response outright.
  integral <- all(abs(w - round(w)) < 1e-8)
  if (!forced && !(any(w > 1) && integral)) {
    return(NULL)
  }
  # Under an explicit binomial_weights = "trials" the weights are taken at face value
  # even when non-integral: maihda_da_proportion_successes() snaps float noise, keeps
  # genuinely fractional mass as-is with a warning, and maihda_auc_weighted() accepts
  # real-valued mass. Rounding fractional trials would imply negative failures.
  trials <- if (integral) round(w) else w
  list(successes = maihda_da_proportion_successes(y, trials), trials = trials)
}

# Per-row success mass implied by a PROPORTION response y with trial counts w.
#
# For the documented idiom -- successes/trials ~ ... , weights = trials -- y * w IS
# the success count, up to floating-point noise from the division: (k / n) * n lands
# within a few ulp of k, so those rows are SNAPPED to the integer and the well-formed
# fit is bit-identical to the cbind(successes, failures) spelling.
#
# A row whose y * w is genuinely fractional is a malformed binomial: no integer
# success count out of w trials produces that proportion, and glmer has already
# warned "non-integer #successes in a binomial glm!" at fit time. Rounding it into a
# whole number INVENTS observations: 240 rows of 3.5 successes out of 7 reconstructed
# to 960 cases where the data carry 840. That moves the AUC as well as the reported
# totals -- on a 12-stratum fit with half-integral successes the totals went 1140/540
# to 1160/520 and the AUC 0.698 to 0.662.
#
# The fractional mass is therefore kept as-is: maihda_auc_weighted() takes real-valued
# case/control mass (the design-weighted branch already passes it some), so the
# concordance is computed under exactly the binomial weighting the model was fitted
# with, and nothing is fabricated. Warned rather than silent, because a fractional
# success count may mean the response is not successes/trials at all -- a proportion
# outcome fitted with integer PRECISION weights reaches this branch too, and its
# "trials" are not trial counts. (It may equally be a genuine successes/trials table
# whose proportions were rounded before they reached R; there the fractional mass
# stays within ~1e-4 AUC of the exact counts, so the warning is the honest signal
# either way.)
maihda_da_proportion_successes <- function(y, w) {
  successes <- y * w
  # Relative tolerance: the representation error of (k / n) * n grows with k, so a
  # fixed absolute epsilon would stop snapping large counts. It stays far below the
  # 0.5 that a genuinely fractional count sits at.
  whole <- abs(successes - round(successes)) <= 1e-8 * pmax(1, abs(successes))
  successes[whole] <- round(successes[whole])
  if (!all(whole)) {
    warning("Aggregated-binomial AUC: the proportion response times the prior ",
            "weights is not a whole number of successes in ", sum(!whole), " of ",
            length(successes), " rows (e.g. ",
            format(successes[!whole][1], digits = 6), " successes out of ",
            format(w[!whole][1], digits = 6), " trials), so these are not ",
            "successes/trials counts. The case/control mass is kept fractional ",
            "rather than rounded -- rounding would invent observations. Supply ",
            "cbind(successes, failures) if the outcome is an aggregated binomial.",
            call. = FALSE)
  }
  successes
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

# Render a case/control total. Whole numbers print as counts ("840"); the fractional
# mass an ill-formed aggregated-binomial fit produces prints with enough digits to be
# recognisable as fractional rather than being silently rounded back into a count.
maihda_format_mass <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  if (abs(x - round(x)) <= 1e-8 * max(1, abs(x))) {
    return(format(round(x), scientific = FALSE))
  }
  format(round(x, 2), nsmall = 2, scientific = FALSE)
}
