# Flag which intersectional strata carry a credibly non-zero interaction.
#
# The scientific payoff of MAIHDA is locating *genuine intersectionality*: strata
# whose outcome departs from what the additive main effects of their defining
# dimensions predict. In the adjusted model that departure IS the stratum random
# effect (BLUP) -- the same quantity plot_effect_decomposition() treats as the
# intersectional component. This file packages those per-stratum interaction BLUPs
# (with their conditional SE / posterior tail) into a flag of which strata are
# credibly != 0, with an honest multiple-comparison story.

#' Flag strata with credibly non-zero intersectional interaction
#'
#' @description
#' Reports, for each intersectional stratum, the \strong{interaction} component of
#' its outcome -- the stratum random effect (BLUP) of an \emph{adjusted} MAIHDA
#' model, i.e. how far the stratum departs from the additive main-effects
#' prediction of its defining dimensions -- and \strong{flags} the strata whose
#' interaction is credibly different from zero. This is the heart of "where is
#' there genuine intersectionality": a flagged stratum is one whose outcome
#' departs credibly from what the additive parts predict -- a descriptive
#' statement about the stratum's outcome, not a causal claim about identity.
#'
#' @details
#' \strong{It must be read off the adjusted model.} Only when the dimensions'
#' additive main effects are in the model (the \emph{adjusted} model of the
#' two-model decomposition, or the crossed-dimensions model) does the stratum
#' random effect isolate the \emph{pure interaction}. On a null model the stratum
#' random effect is the total between-stratum deviation (additive + interaction),
#' so passing one is flagged with a warning. The opposite mis-specification is
#' flagged too: a model that adds a \emph{fixed} interaction among the dimensions
#' (e.g. \code{var1 * var2}) absorbs the intersectional effect into fixed cell means,
#' so the stratum random effect is no longer the pure interaction. Passing a
#' \code{\link{maihda}} result uses the right model automatically.
#'
#' \strong{Frequentist vs. Bayesian evidence.} For the frequentist engines
#' (\code{lme4}, \code{wemix}, \code{ordinal}) the flag comes from the BLUP's
#' conditional standard error: a Wald interval at \code{conf_level} and a two-sided
#' p-value, with an optional multiplicity correction (\code{adjust}). For
#' \code{brms} the full posterior is already available, so the \emph{exact}
#' posterior tail is used -- a credible interval at \code{conf_level} and the
#' probability of direction \code{pd = max(P(BLUP > 0), P(BLUP < 0))} (in
#' \code{[0.5, 1]}; the sign is in \code{direction}) -- and \code{adjust} is not
#' applied: the posterior is already partially pooled, the hierarchical-shrinkage
#' answer to multiplicity (Gelman, Hill & Yajima 2012). That is a per-stratum
#' answer; the marginal intervals carry no formal joint error-rate guarantee for
#' the collective claim "an interaction exists \emph{somewhere}" (the disjunction
#' reading below applies to Bayesian flags too).
#'
#' \strong{Multiplicity: partial pooling and a correction are different things, and
#' the experts disagree.}
#' \itemize{
#'   \item \emph{Shrinkage (magnitude/sign).} The stratum BLUP is partially pooled,
#'     so extreme values are regularised toward the grand mean, attenuating
#'     exaggerated-magnitude and wrong-sign (Type M/S) error (Gelman & Carlin 2014).
#'     Gelman, Hill & Yajima (2012) argue this shrinkage \emph{usually substitutes}
#'     for a classical multiple-comparisons correction (the problem can "disappear
#'     entirely" in the hierarchical model); on that view the flag/no-flag step
#'     itself is what to avoid -- the null of an \emph{exactly} zero interaction is
#'     rarely the question (McShane, Gelman et al. 2019) -- so report the estimate
#'     and its interval.
#'   \item \emph{Whether to correct.} If you do want an error-rate screen, whether a
#'     correction is warranted depends on the \emph{inferential structure} of the
#'     claim. Each
#'     stratum as its own pre-specified hypothesis ("does \emph{this} stratum
#'     interact?") is \emph{individual} testing and needs none -- \strong{only} if you
#'     do not also read the flags collectively. Once the question is "is there an
#'     interaction \emph{somewhere}?" -- which an automated all-strata scan
#'     effectively is -- it is \emph{disjunction} testing and a correction applies.
#' }
#' \strong{\code{adjust = "BH"} is the default}: fitting and flagging every stratum
#' in one call is the disjunction/screening case, where controlling the expected
#' \emph{proportion} of false discoveries (FDR) is the appropriate goal. Pass
#' \code{adjust = "none"} only when each stratum is a genuine, pre-specified
#' individual hypothesis. The FDR choice (over family-wise
#' \code{"bonferroni"}/\code{"holm"}) is this package's, matching that screening
#' goal. The flag itself is a Wald test on a shrunken BLUP whose
#' conditional SE treats the variance components as known, so it (and any
#' \code{adjust} on it) is an explicit, \emph{conservative} screen -- strict, not
#' liberal. Partial pooling deflates a truly-null
#' stratum's BLUP more than its conditional SE, so the null z-statistic is
#' sub-normal (variance about the shrinkage fraction, below 1) and the screen
#' under-flags rather than over-flags, degenerating to no flags at the
#' singular/zero-variance boundary. Lead with the interval (and, for
#' \code{brms}, the probability of direction); the substantive question is often not
#' whether an interaction differs from zero but whether it exceeds a smallest
#' interaction of interest (an equivalence reading; Schuirmann 1987; Kruschke 2018),
#' read from the interval.
#'
#' The interaction is reported on the model's link (latent) scale -- a log-odds
#' deviation for a logistic model, etc. -- because the additive/interaction split
#' is only exact there. On a non-identity link that makes it a \emph{multiplicative}
#' departure, which is a different quantity from a shift in the outcome itself; see
#' the section below before reporting it.
#'
#' @section What the interaction means on a non-identity link:
#' For a Gaussian (identity-link) fit the stratum BLUP \eqn{u_j} \emph{is} the part
#' of the stratum's mean outcome attributable to interaction, in the outcome's own
#' units. For a logistic fit it is not: \eqn{u_j} is a deviation in
#' \strong{log-odds}, so what this function reports is the \emph{multiplicative},
#' odds-scale interaction. (What follows is written for the logistic case, the one
#' the literature works in. The same logic holds on any non-identity link, in its own
#' units: a Poisson fit's \eqn{u_j} is a log-rate departure and \code{scale =
#' "response"} gives it as a difference in expected counts, and a cumulative
#' (ordinal) fit's is a latent shift reported as a difference in expected category
#' score. The printed column is named for whichever applies.) Evans et al. (2024, section 2.5.1)
#' put it directly -- in a logistic MAIHDA one can "no longer directly interpret
#' \eqn{u_j} as the change in mean outcomes (i.e., shift in probabilities)
#' attributable to interaction effects".
#'
#' \strong{Three quantities, of which this function reports two.} Keeping them apart
#' is the whole of the difficulty:
#' \enumerate{
#'   \item \strong{The multiplicative interaction, \eqn{u_j}} -- the departure from
#'     additivity \emph{in log-odds}, i.e. whether a dimension multiplies the odds by
#'     the same factor at every level of the others. This is what the default
#'     \code{scale = "link"} reports and what \code{flagged} is about.
#'   \item \strong{The same departure in outcome units, \eqn{\pi^B_j}} -- what
#'     \code{scale = "response"} returns. Evans et al. (2024) define it as
#'     \deqn{\pi^B_j = \pi_j - \pi^A_j, \qquad
#'           \pi_j = \mathrm{logit}^{-1}(x_j'\beta + u_j), \qquad
#'           \pi^A_j = \mathrm{logit}^{-1}(x_j'\beta),}
#'     the gap between a stratum's total predicted probability and the probability
#'     implied by the additive main effects alone, and rank-plot \eqn{\pi^B_j} where
#'     the linear case plots \eqn{u_j}. This is quantity 1 \emph{re-expressed}, not a
#'     second finding: it is zero exactly when \eqn{u_j} is zero and it flags the same
#'     strata, so it changes the units you report, not what you may conclude.
#'   \item \strong{The additive (risk-difference) interaction} -- whether a dimension
#'     adds the same number of percentage points of risk at every level of the others.
#'     \strong{Neither of the above reports this}, and it is generally non-zero even
#'     where \eqn{u_j} is exactly zero, because the logistic curve is steeper in the
#'     middle than in the tails. With a \eqn{-2} baseline and \eqn{+0.7} for each of
#'     two dimensions and no interaction at all (\eqn{u_j = 0} throughout), the second
#'     dimension still adds 9.5 percentage points of risk at one level of the first
#'     and 14.0 at the other. That excess is a property of the link, not a finding.
#' }
#'
#' So "no strata flagged" supports "no credible \strong{multiplicative} interaction",
#' not "no interaction": additivity in log-odds does not carry over to probabilities,
#' so quantity 3 is generally non-zero regardless, and it is often the one a policy
#' audience cares about. (\emph{Generally}, not always. The risk-difference
#' interaction is positive below the curve's midpoint and negative above it, so it
#' passes through zero: with two dimensions of effect \eqn{a} and \eqn{b} it vanishes
#' exactly at an intercept of \eqn{-(a + b)/2}, where the two comparisons straddle the
#' midpoint symmetrically. That is a single configuration, not the general case, and
#' the flags say nothing about which one you are in either way.)
#'
#' \strong{Both scales flag the same strata.} Writing \eqn{g(u)} for the map from a
#' stratum's BLUP to its \eqn{\pi^B_j}, \eqn{g} is strictly increasing with
#' \eqn{g(0) = 0}, so estimate, interval endpoints and zero all carry across
#' together. \code{flagged}, \code{direction}, \code{p_value} / \code{p_adjusted} and
#' \code{pd} are therefore identical under either \code{scale}, and the response-scale
#' interval is the \emph{exact} image of the link-scale one rather than a delta-method
#' approximation -- inheriting its conditionality (fixed effects and variance
#' components held at their point estimates), which is why no simulation step is
#' needed. What does change is the size and the ranking: the same log-odds departure
#' is worth more probability near \eqn{\pi = 0.5} than in the tail, so strata are
#' ordered by the quantity actually reported.
#'
#' \strong{If quantity 3 is what you need}, no argument here will give it to you --
#' fit a linear probability model instead: \code{fit_maihda(..., family = "gaussian")}
#' on a 0/1 outcome, which \code{\link{fit_maihda}} signposts when it auto-detects a
#' binary outcome. There the BLUP \emph{is} the risk-difference interaction by
#' construction, at the usual costs (predictions outside \code{[0, 1]},
#' heteroskedastic residuals).
#'
#' @param object A \code{maihda_analysis} from \code{\link{maihda}} (preferred --
#'   its adjusted / crossed-dimensions model is used automatically) or a
#'   \code{maihda_model} from \code{\link{fit_maihda}} (which should be the
#'   \emph{adjusted} model; a null model is accepted but warned about).
#' @param conf_level Confidence / credible level for the interval and the flag.
#'   Default 0.95.
#' @param adjust Multiple-comparison adjustment for the per-stratum p-values
#'   (frequentist engines only): \code{"BH"} (default; false-discovery rate) or any
#'   method accepted by \code{\link[stats]{p.adjust}}, including \code{"none"} for
#'   the uncorrected, per-stratum individual-testing view. Ignored for \code{brms}
#'   (which uses the posterior tail directly; a message is shown only if you set it
#'   explicitly).
#' @param scale Scale the interaction is reported on. \code{"link"} (default) gives
#'   the stratum BLUP \eqn{u_j} on the model's link scale -- for a Gaussian fit the
#'   outcome's own units, for a logistic fit a log-odds departure.
#'   \code{"response"} instead gives Evans et al.'s (2024, sec. 2.5.1)
#'   \eqn{\pi^B_j = \pi_j - \pi^A_j}: the stratum's total predicted outcome minus the
#'   outcome implied by its additive main effects alone, so a logistic fit reports a
#'   difference in \emph{probability}, a count fit a difference in expected count, and
#'   a cumulative (ordinal) fit a difference in expected category score -- the printed
#'   column is named for whichever it is.
#'   The two scales flag the same strata (see the section below); an identity-link
#'   fit returns the same numbers either way. For a
#'   \code{decomposition = "crossed-dimensions"} model the additive baseline is the
#'   dimension random effects, matching that decomposition.
#' @param rope Optional equivalence region (a "smallest interaction of interest")
#'   for an "is the interaction \emph{negligible}?" reading (Schuirmann 1987;
#'   Kruschke 2018), read on the requested \code{scale} -- log-odds under the
#'   default, probability points under \code{scale = "response"} for a logistic
#'   fit. \code{NULL} (default) gives
#'   only the usual zero-centred flag. A single positive number \code{d} means the
#'   symmetric region \code{c(-d, d)}; or supply \code{c(lower, upper)}. When set, the
#'   result gains a \code{decision} column classifying each stratum from its
#'   \code{conf_level} interval relative to the region: \code{"relevant"} (interval
#'   entirely outside it), \code{"negligible"} (entirely inside it), or
#'   \code{"inconclusive"} (straddling a bound).
#' @param ... Currently unused.
#'
#' @return An object of class \code{maihda_interactions} (a data frame), one row
#'   per stratum, sorted flagged-first then by \code{abs(interaction)}. Columns
#'   common to every engine: \code{stratum}, \code{label}, \code{n} (stratum size),
#'   \code{interaction} (the BLUP, under the default \code{scale}),
#'   \code{lower}/\code{upper} (the interval),
#'   \code{flagged} (logical), and \code{direction} (\code{"above"}/\code{"below"}
#'   the additive expectation). Under \code{scale = "response"} the
#'   \code{interaction}/\code{lower}/\code{upper} columns hold \eqn{\pi^B_j} and its
#'   interval, and \code{se} is dropped: the conditional standard error is a
#'   link-scale quantity, and the response-scale interval is deliberately not
#'   symmetric about the estimate, so no single SE would reproduce it. Frequentist
#'   fits add \code{se} (link scale only) and \code{p_value}
#'   (and \code{p_adjusted} when \code{adjust != "none"}). \code{p_value} is a
#'   \emph{conditional} screening statistic -- a Wald tail on the
#'   shrunken BLUP's conditional SE, with the variance components treated as known
#'   -- \strong{not} a calibrated frequentist p-value; it is \emph{conservative}
#'   (stochastically large under a true null), so the BH-adjusted flag under-flags
#'   truly-null strata rather than delivering an exact error-rate guarantee (see
#'   Details). \code{brms} instead adds
#'   \code{pd} (probability of direction, \code{max(P(>0), P(<0))} in
#'   \code{[0.5, 1]}). When \code{rope} is set, a
#'   \code{decision} column (\code{"relevant"}/\code{"negligible"}/\code{"inconclusive"})
#'   is added. Attributes record \code{conf_level}, \code{adjust}, \code{rope},
#'   \code{engine}, \code{model_type}, \code{n_strata}, \code{n_flagged},
#'   \code{scale}, \code{link} (the model's link name, so a caller can tell a
#'   log-odds departure from one in the outcome's own units), \code{response_kind}
#'   (what the model's response scale is -- \code{"probability"}, \code{"count"},
#'   \code{"score"} or \code{"response"} -- which the link alone does not settle,
#'   since a cumulative fit is logit-linked but scores categories), \code{singular},
#'   and -- on a non-identity link under the default \code{scale} --
#'   \code{response_interaction}, the outcome-scale estimates keyed by stratum. That
#'   last attribute exists so \code{print()} can show the outcome-scale size beside
#'   the link-scale one; the columns are the same either way, and
#'   \code{scale = "response"} is how to obtain those numbers, with their interval,
#'   as data.
#'
#' @references
#' Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V. (2018). A
#' multilevel approach to modeling health inequalities at the intersection of
#' multiple social identities. \emph{Social Science & Medicine}, 203, 64-73.
#'
#' Merlo, J. (2018). Multilevel analysis of individual heterogeneity and
#' discriminatory accuracy (MAIHDA) within an intersectional framework.
#' \emph{Social Science & Medicine}, 203, 74-80.
#'
#' Evans, C. R., Leckie, G., Subramanian, S. V., Bell, A., & Merlo, J. (2024). A
#' tutorial for conducting intersectional multilevel analysis of individual
#' heterogeneity and discriminatory accuracy (MAIHDA).
#' \emph{SSM - Population Health}, 26, 101664. \doi{10.1016/j.ssmph.2024.101664}
#'
#' Gelman, A., Hill, J., & Yajima, M. (2012). Why we (usually) don't have to worry
#' about multiple comparisons. \emph{Journal of Research on Educational
#' Effectiveness}, 5(2), 189-211.
#'
#' Gelman, A., & Carlin, J. (2014). Beyond power calculations: assessing Type S
#' (sign) and Type M (magnitude) errors. \emph{Perspectives on Psychological
#' Science}, 9(6), 641-651.
#'
#' Rubin, M. (2021). When to adjust alpha during multiple testing: a consideration
#' of disjunction, conjunction, and individual testing. \emph{Synthese}, 199(3-4),
#' 10969-11000. \doi{10.1007/s11229-021-03276-4}
#'
#' McShane, B. B., Gal, D., Gelman, A., Robert, C., & Tackett, J. L. (2019). Abandon
#' statistical significance. \emph{The American Statistician}, 73(sup1), 235-245.
#'
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and
#' the power approach for assessing the equivalence of average bioavailability.
#' \emph{Journal of Pharmacokinetics and Biopharmaceutics}, 15(6), 657-680.
#'
#' Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian
#' estimation. \emph{Advances in Methods and Practices in Psychological Science},
#' 1(2), 270-280.
#'
#' @seealso \code{\link{maihda}}, \code{\link{calculate_pcv}},
#'   \code{\link{summary.maihda_model}}; and \code{plot(\dots,
#'   highlight_interactions = TRUE)} to mark the flagged strata on the
#'   effect-decomposition / predicted / shrinkage plots.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#' a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race),
#'             data = maihda_health_data)
#' maihda_interactions(a)                  # FDR-screened (default adjust = "BH")
#' maihda_interactions(a, adjust = "none") # uncorrected per-stratum individual view
#' maihda_interactions(a, rope = 0.1)      # equivalence: |interaction| within 0.1?
#'
#' # On a logistic fit the BLUP is a log-odds departure; scale = "response" reports
#' # Evans et al.'s (2024) pi_j - pi^A_j, the same interaction in probability points.
#' b <- maihda(Obese ~ Gender + Race + (1 | Gender:Race),
#'             data = maihda_health_data, family = "binomial")
#' maihda_interactions(b, scale = "response")
#' }
#'
#' @export
#' @importFrom stats qnorm pnorm quantile median p.adjust terms p.adjust.methods
#' @importFrom reformulas nobars
maihda_interactions <- function(object, conf_level = 0.95, adjust = "BH",
                                rope = NULL, scale = c("link", "response"), ...) {
  adjust_was_set <- !missing(adjust)
  resolved <- maihda_resolve_interaction_model(object)
  model <- resolved$model
  summary_obj <- resolved$summary
  model_type <- resolved$model_type

  conf_level <- maihda_validate_conf_level(conf_level)
  adjust <- match.arg(adjust, c("none", stats::p.adjust.methods))
  scale <- match.arg(scale)
  rope <- maihda_normalize_rope(rope)

  se_tab <- summary_obj$stratum_estimates
  if (is.null(se_tab) || nrow(se_tab) == 0) {
    stop("No stratum estimates are available to assess interaction.", call. = FALSE)
  }

  strata <- as.character(se_tab$stratum)
  label <- if ("label" %in% names(se_tab)) as.character(se_tab$label) else strata
  est <- as.numeric(se_tab$random_effect)
  n <- maihda_interaction_strata_n(model, strata)
  engine <- model$engine
  alpha <- 1 - conf_level

  if (identical(engine, "brms")) {
    if (adjust_was_set && !identical(adjust, "none")) {
      message("maihda_interactions(): 'adjust' is ignored for brms models -- ",
              "the posterior is already partially pooled, and no additional ",
              "p-value correction is applied to its tail probabilities.")
    }
    group <- maihda_interaction_group(model)
    tail <- maihda_interaction_brms_tail(model$model, group, conf_level)
    idx <- match(strata, tail$stratum)
    interaction <- tail$interaction[idx]
    lower <- tail$lower[idx]
    upper <- tail$upper[idx]
    pd <- tail$pd[idx]
    flagged <- !is.na(lower) & !is.na(upper) & ((lower > 0) | (upper < 0))
    out <- data.frame(
      stratum = strata, label = label, n = n,
      interaction = interaction, lower = lower, upper = upper, pd = pd,
      flagged = flagged,
      direction = ifelse(interaction >= 0, "above", "below"),
      stringsAsFactors = FALSE
    )
  } else {
    se <- as.numeric(se_tab$se)
    z <- stats::qnorm((1 + conf_level) / 2)
    lower <- est - z * se
    upper <- est + z * se
    # Wald two-sided p; undefined where the BLUP SE is zero/NA (singular/boundary
    # fit), where "no flag" is not evidence of no interaction.
    p_value <- 2 * stats::pnorm(-abs(est / se))
    p_value[is.na(se) | se <= 0] <- NA_real_
    p_adjusted <- if (identical(adjust, "none")) {
      p_value
    } else {
      stats::p.adjust(p_value, method = adjust)
    }
    flagged <- !is.na(p_adjusted) & p_adjusted < alpha
    out <- data.frame(
      stratum = strata, label = label, n = n,
      interaction = est, se = se, lower = lower, upper = upper,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
    if (!identical(adjust, "none")) {
      out$p_adjusted <- p_adjusted
    }
    out$flagged <- flagged
    out$direction <- ifelse(est >= 0, "above", "below")
  }

  # Carry the interaction onto the response scale (Evans et al. 2024, sec. 2.5.1)
  # BEFORE the ROPE block and the sort, so the equivalence region is read on the
  # scale the caller asked for and the ranking follows the reported quantity. The
  # map is monotone through zero, so flagged / direction / p_value / pd are
  # untouched; `se` is dropped because the conditional SE is a link-scale quantity
  # and the response-scale interval is deliberately not symmetric about the estimate,
  # so no single SE reproduces it.
  link <- maihda_model_link_name(model)
  # What the response scale IS for this family -- not inferable from the link alone
  # (a cumulative fit is logit-linked but scores categories, not probabilities).
  response_kind <- maihda_interaction_response_kind(model)
  response_map <- NULL
  if (identical(scale, "response")) {
    resp <- maihda_interaction_response_values(
      model, summary_obj, out$stratum, out$interaction, out$lower, out$upper)
    out$interaction <- resp$interaction
    out$lower <- resp$lower
    out$upper <- resp$upper
    out$se <- NULL
  } else if (!is.na(link) && !identical(link, "identity")) {
    # Link scale on a non-identity link: also carry the outcome-scale estimate, which
    # print() shows beside the log-odds departure so the reader sees the size in
    # probability points without re-running. Kept as an attribute rather than a
    # column so the returned shape is unchanged; scale = "response" is the way to
    # get it (with its interval) as data. Never fatal -- an engine that cannot
    # produce per-stratum predictions just prints as before.
    response_map <- tryCatch(
      stats::setNames(
        maihda_interaction_response_values(model, summary_obj, out$stratum,
                                           out$interaction, out$lower,
                                           out$upper)$interaction,
        as.character(out$stratum)),
      error = function(e) NULL)
  }

  # Equivalence / ROPE reading (Schuirmann 1987; Kruschke 2018): an
  # "is the interaction negligible?" classification from each stratum's interval
  # relative to the smallest-interaction-of-interest region, separate from the
  # zero-centred flag. Uses the conf_level interval (lower/upper), so it is valid
  # across engines without re-deriving anything.
  if (!is.null(rope)) {
    lo <- rope[1]; hi <- rope[2]
    have <- !is.na(out$lower) & !is.na(out$upper)
    inside  <- have & out$lower >= lo & out$upper <= hi
    outside <- have & (out$lower > hi | out$upper < lo)
    out$decision <- ifelse(!have, NA_character_,
                    ifelse(inside, "negligible",
                    ifelse(outside, "relevant", "inconclusive")))
  }

  # Flagged strata first, then by interaction magnitude (most extreme first).
  ord <- order(!out$flagged, -abs(out$interaction))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  singular <- isTRUE(summary_obj$diagnostics$singular) ||
    isTRUE(model$diagnostics$singular)

  attr(out, "conf_level") <- conf_level
  attr(out, "adjust") <- adjust
  attr(out, "rope") <- rope
  attr(out, "engine") <- engine
  attr(out, "model_type") <- model_type
  attr(out, "n_strata") <- nrow(out)
  attr(out, "n_flagged") <- sum(out$flagged)
  attr(out, "scale") <- scale
  # The link name travels with the result so print() -- and any caller -- can tell a
  # log-odds departure from one in the outcome's own units. NA for an engine that
  # exposes no family (nothing is claimed then).
  attr(out, "link") <- link
  attr(out, "response_kind") <- response_kind
  attr(out, "response_interaction") <- response_map
  attr(out, "singular") <- singular
  class(out) <- c("maihda_interactions", "data.frame")
  out
}

#' Print a MAIHDA interaction diagnostic
#'
#' @param x A \code{maihda_interactions} object from \code{\link{maihda_interactions}}.
#' @param ... Additional arguments (not used).
#' @return No return value, called for side effects.
#' @export
#' @importFrom utils head
print.maihda_interactions <- function(x, ...) {
  # Subsetting a maihda_interactions object -- head(), `[`, or a dplyr verb -- keeps
  # the class but drops these attributes, so read each defensively and fall back to a
  # value recomputed from the data frame. Without this, a plain head(x) errors (e.g.
  # `if (n_flagged > 0)` on a NULL attribute is "argument of length zero"). The full
  # object returned by maihda_interactions() carries every attribute and prints
  # exactly as before.
  get_attr <- function(name, default) {
    v <- attr(x, name)
    if (length(v) == 0) default else v
  }
  has_flagged <- "flagged" %in% names(x)
  conf       <- get_attr("conf_level", 0.95)
  adjust     <- get_attr("adjust", NA_character_)
  engine     <- get_attr("engine", NA_character_)
  model_type <- get_attr("model_type", NA_character_)
  n_strata   <- get_attr("n_strata", nrow(x))
  n_flagged  <- attr(x, "n_flagged")
  if (length(n_flagged) != 1L || is.na(n_flagged)) {
    n_flagged <- if (has_flagged) sum(x$flagged, na.rm = TRUE) else NA_integer_
  }

  # Semantic colouring via the shared cli palette (auto-degrades to plain text when
  # the destination has no ANSI support: knitr/vignettes, R CMD check, testthat,
  # NO_COLOR). Colour encodes emphasis, not valence.
  pal <- maihda_palette()
  count_style <- if (isTRUE(n_flagged > 0)) function(s) pal$bold(pal$accent(s)) else pal$muted

  cat(cli::rule(left = pal$bold("Intersectional interactions")), "\n", sep = "")

  conf_pct <- conf * 100
  evidence <- if (identical(engine, "brms")) {
    sprintf("%.0f%% credible interval; probability of direction", conf_pct)
  } else if (is.na(adjust)) {
    sprintf("%.0f%% interval", conf_pct)
  } else if (identical(adjust, "none")) {
    sprintf("%.0f%% interval; no multiplicity correction", conf_pct)
  } else {
    # "conservative": the null p-value is a Wald tail on a shrunken BLUP's
    # conditional SE (variance components treated as known), stochastically large
    # under a true null -- so the BH flag under-flags, it does not over-flag; not a
    # calibrated frequentist p-value that BH could repair. See ?maihda_interactions.
    sprintf("%.0f%% interval; %s-adjusted conservative p-values", conf_pct, adjust)
  }
  count_label <- if (is.na(n_flagged)) "?" else as.character(n_flagged)
  cat(sprintf("%s of %d strata flagged (%s).\n",
              count_style(count_label), n_strata, pal$muted(evidence)))
  scale_attr <- get_attr("scale", "link")
  scale_word <- if (identical(scale_attr, "response")) {
    "response (outcome)"
  } else {
    "link (latent)"
  }
  if (!is.na(model_type)) {
    cat(pal$muted(sprintf("Model: %s; interaction on the %s scale.\n",
                          model_type, scale_word)))
  } else {
    cat(pal$muted(sprintf("Interaction reported on the %s scale.\n", scale_word)))
  }
  # Show the outcome-scale interaction next to the link-scale one, so the reader at
  # the console sees the departure in probability points (or expected counts) without
  # having to re-run. Display only: the object's columns are unchanged, and the
  # footer says how to get these numbers as data.
  resp_map <- get_attr("response_interaction", NULL)
  kind_attr <- get_attr("response_kind", NULL)
  show_resp <- function(df) {
    maihda_interactions_add_response_col(df, resp_map, kind_attr)
  }

  # On a non-identity link the BLUP is a departure in LINK units, so a flag (or its
  # absence) is a statement about MULTIPLICATIVE interaction only -- the one thing a
  # reader cannot infer from the numbers on screen. ?maihda_interactions carries the
  # rest, including the probability-scale quantity this is not.
  link_note <- maihda_interaction_link_note(get_attr("link", NA_character_), scale_attr,
                                            has_response = !is.null(resp_map),
                                            kind = kind_attr)
  if (!is.null(link_note)) {
    cat(pal$muted(link_note))
  }

  rope <- attr(x, "rope")
  if (!is.null(rope) && "decision" %in% names(x)) {
    dec <- factor(x$decision, levels = c("relevant", "negligible", "inconclusive"))
    tab <- table(dec)
    # relevant -> accent ("look here"); negligible -> the neutral data colour (a
    # firm "practically zero", NOT green/good); inconclusive -> muted (undecided).
    cat(sprintf("Equivalence vs ROPE [%g, %g]: %s | %s | %s.\n", rope[1], rope[2],
                pal$accent(sprintf("%d relevant", tab[["relevant"]])),
                pal$secondary(sprintf("%d negligible", tab[["negligible"]])),
                pal$muted(sprintf("%d inconclusive", tab[["inconclusive"]]))))
  }

  if (isTRUE(attr(x, "singular"))) {
    cat(pal$warn(paste0(
      "\n", cli::symbol$warning,
      " singular/boundary fit (between-stratum variance ~ 0); the BLUP SEs\n",
      "  collapse toward zero, so 'no flag' is not evidence of no interaction.\n")))
  }

  cat("\n")
  if (!has_flagged) {
    # The 'flagged' column was selected away in a subset; show the rows as-is.
    print(show_resp(utils::head(as.data.frame(x), 10)), row.names = FALSE, digits = 4)
  } else {
    flagged_rows <- x[x$flagged %in% TRUE, , drop = FALSE]
    if (nrow(flagged_rows) == 0) {
      cat(pal$muted(
        "No strata show interaction credibly different from zero at this level.\n"))
    } else {
      print(show_resp(utils::head(as.data.frame(flagged_rows), 10)),
            row.names = FALSE, digits = 4)
      if (nrow(flagged_rows) > 10) {
        cat(pal$muted(sprintf("  ... and %d more flagged strata\n",
                              nrow(flagged_rows) - 10)))
      }
    }
  }

  footer <- if (!identical(engine, "brms") && identical(adjust, "none") && n_strata > 1) {
    paste0("\nFlagging many strata inflates false positives; for a screening error-rate\n",
           "  story use adjust = \"BH\" (FDR). Interaction BLUPs are shrunken estimates,\n",
           "  so correction is optional -- see ?maihda_interactions.\n")
  } else {
    paste0("\nInteraction BLUPs are shrunken (partially pooled) estimates; treat flags as\n",
           "  exploratory. See ?maihda_interactions.\n")
  }
  cat(pal$muted(footer))

  invisible(x)
}

# ---- internal helpers -------------------------------------------------------

# The response-scale interaction of Evans et al. (2024, sec. 2.5.1):
# pi^B_j = pi_j - pi^A_j, a stratum's total predicted outcome minus the outcome
# implied by the additive main effects alone. On the link scale the interaction IS
# the BLUP; on the response scale it is this difference -- for a logistic fit a
# difference in PROBABILITY rather than in log-odds.
#
# Computed by reusing the engine's own per-stratum prediction helper (the one behind
# plot(type = "predicted")), so offsets, prior/sampling/trial weights,
# aggregated-binomial responses and the crossed-dimensions additive baseline are all
# handled in one place -- and the pi_j here is the same pi_j that plot draws. Those
# helpers read their three link-scale inputs from random_effect / lower_95 /
# upper_95, so the conf_level-specific interval this function was given is
# substituted in before the call. Where covariates vary within a stratum the two
# predictions are each averaged over that stratum's own rows, which reduces exactly
# to Evans et al.'s definition when the fixed part is constant within a stratum (the
# canonical adjusted model, whose fixed effects are the dimension main effects).
#
# g(u) = weighted mean over a stratum's rows of [linkinv(eta_i + u) - linkinv(eta_i)]
# is strictly increasing with g(0) = 0, so estimate, interval endpoints and zero all
# map across together: the flag, the direction and the p-value / pd are IDENTICAL on
# both scales, and the returned interval is the exact image of the link-scale one,
# not a delta-method approximation. It inherits that interval's conditionality (fixed
# effects and variance components held at their point estimates).
maihda_interaction_response_values <- function(model, summary_obj, strata,
                                               est, lower, upper) {
  se_tab <- summary_obj$stratum_estimates
  idx <- match(as.character(se_tab$stratum), as.character(strata))
  se_tab$random_effect <- est[idx]
  se_tab$lower_95 <- lower[idx]
  se_tab$upper_95 <- upper[idx]
  summary_obj$stratum_estimates <- se_tab

  pred <- switch(
    model$engine,
    lme4 = maihda_stratum_predictions_lme4(model, summary_obj, scale = "response"),
    brms = maihda_stratum_predictions_brms(model, summary_obj, scale = "response"),
    wemix = maihda_stratum_predictions_wemix(model, summary_obj, scale = "response"),
    ordinal = maihda_stratum_predictions_ordinal(model, summary_obj, scale = "response"),
    stop("scale = \"response\" is not available for engine \"", model$engine, "\".",
         call. = FALSE))

  k <- match(as.character(strata), as.character(pred$stratum))
  # A stratum whose link-scale input was NA (a singular fit's undefined BLUP SE)
  # aggregates to NaN; report it as the NA it came from rather than as NaN.
  clean <- function(x) {
    x[!is.finite(x)] <- NA_real_
    x
  }
  list(interaction = clean(pred$predicted_row[k] - pred$fixed_row[k]),
       lower = clean(pred$lower_row[k] - pred$fixed_row[k]),
       upper = clean(pred$upper_row[k] - pred$fixed_row[k]))
}

# Column header for the printed outcome-scale interaction, named for what it is on
# this link so the header itself says which units the number is in.
# What the RESPONSE scale of a model actually is. The link alone does not settle it:
# a cumulative (ordinal) fit has a LOGIT link, but its response scale is the expected
# category score, not a probability -- so switching on the link labels an ordinal
# fit's interaction "prob_diff ... probability points", which is simply wrong. Drives
# the printed header and its gloss; the DEPARTURE units still come from the link.
maihda_response_kind_from_link <- function(link) {
  switch(as.character(link),
         logit = "probability",
         probit = "probability",
         cloglog = "probability",
         log = "count",
         "response")
}

maihda_interaction_response_kind <- function(model) {
  is_ord <- isTRUE(tryCatch(maihda_family_is_ordinal(model$family),
                            error = function(e) FALSE))
  if (is_ord) {
    return("score")
  }
  maihda_response_kind_from_link(maihda_model_link_name(model))
}

maihda_interaction_response_header <- function(kind) {
  switch(as.character(kind),
         probability = "prob_diff",
         count = "count_diff",
         score = "score_diff",
         "resp_diff")
}

# Insert the outcome-scale interaction into a frame about to be PRINTED, directly
# after the link-scale one. Display only -- maihda_interactions() returns the same
# columns it always did, and scale = "response" is how these numbers are obtained as
# data. `resp_map` is keyed by stratum, so a reordered or subset frame still lines
# up, and anything the map does not cover is left out rather than silently NA-filled.
maihda_interactions_add_response_col <- function(df, resp_map, kind) {
  if (is.null(resp_map) || !is.data.frame(df) || nrow(df) == 0 ||
      !"stratum" %in% names(df)) {
    return(df)
  }
  vals <- unname(resp_map[as.character(df$stratum)])
  if (all(is.na(vals))) {
    return(df)
  }
  header <- maihda_interaction_response_header(kind)
  if (header %in% names(df)) {
    return(df)
  }
  out <- df
  out[[header]] <- vals
  pos <- match("interaction", names(df))
  if (!is.na(pos)) {
    # Sit next to the quantity it restates, not at the far right of a wide table.
    out <- out[, append(names(df), header, after = pos), drop = FALSE]
  }
  out
}

# One printed line naming what a non-identity-link interaction actually is. The
# additive/interaction split is exact only on the link scale, so for a logit fit the
# BLUP is a LOG-ODDS departure. NULL for an identity link (where the BLUP is already
# in the outcome's units) or an unknown link (claim nothing). `kind` names the
# response scale (see maihda_interaction_response_kind); it defaults to the link's
# own reading, which is right for every family except the cumulative one.
maihda_interaction_link_note <- function(link, scale = "link", has_response = FALSE,
                                         kind = NULL) {
  if (length(link) != 1L || is.na(link) || identical(link, "identity")) {
    return(NULL)
  }
  if (is.null(kind) || is.na(kind)) {
    kind <- maihda_response_kind_from_link(link)
  }
  units <- switch(link,
                  logit = "log-odds",
                  probit = "probit (latent)",
                  cloglog = "complementary log-log",
                  log = "log-rate",
                  sprintf("%s-scale", link))
  if (identical(scale, "response")) {
    # On this scale the number already IS the outcome-scale interaction, so the
    # caveat below does not apply; what a reader needs instead is which quantity it
    # is and that the flags did not change with the scale. The pi notation is
    # Evans et al.'s and is probability notation, so it appears only where the
    # response scale really is a probability.
    lead <- switch(as.character(kind),
                   probability = "A probability difference (pi_j - pi^A_j, Evans et al. 2024)",
                   count = "An expected-count difference",
                   score = "An expected-score difference",
                   "An outcome-scale difference")
    return(sprintf(paste0("%s: the
",
                          "  interaction carried onto the outcome scale; flags match the
",
                          "  %s scale.
"),
                   lead, units))
  }
  # Name the two scales and which one the flags are about. Two wordings to avoid:
  # "odds combine additively" (when log-odds add, odds MULTIPLY), and "the
  # probability scale" for what is ruled out -- that phrase would collide with the
  # probability-POINTS restatement below and reintroduce the conflation between the
  # risk-scale departure and pi^B_j that this note exists to prevent. What the flags
  # do and do not license is in ?maihda_interactions; the line stays to one sentence.
  ruled_out <- switch(as.character(kind),
                      probability = "risks",
                      count = "counts",
                      score = "expected scores",
                      "outcomes")
  note <- sprintf(paste0("A %s departure (%s link): the flags test additivity on the link
",
                         "  scale, not in %s.
"),
                  units, link, ruled_out)
  if (isTRUE(has_response)) {
    # A printed column nobody asked for needs one line saying which question it
    # answers and where to get it as data; the rest is in ?maihda_interactions.
    resp_units <- switch(as.character(kind),
                         probability = "probability points",
                         count = "expected counts",
                         score = "expected category score",
                         "outcome units")
    formula_txt <- if (identical(as.character(kind), "probability")) {
      " (pi_j - pi^A_j)"
    } else {
      ""
    }
    note <- paste0(note,
                   sprintf(paste0("  %s is that same departure in %s%s;
",
                                  "  scale = \"response\" returns it with its interval.
"),
                           maihda_interaction_response_header(kind), resp_units,
                           formula_txt))
  }
  note
}

# Resolve the model + summary the interaction diagnostic should read, and label
# the model type. A maihda_analysis routes to its adjusted (two-model) or single
# crossed-dimensions model -- both carry the additive part, so the stratum random
# effect is the pure interaction. A bare maihda_model is used directly, with a
# guardrail warning when it looks like a null model.
maihda_resolve_interaction_model <- function(object) {
  # A longitudinal MAIHDA's per-stratum interaction is a TRAJECTORY (random
  # intercept + slope(s)), so the scalar BLUP diagnostic is undefined: collapsing it
  # to one number drops the slope and silently returns a cross-sectional-looking
  # value. The automatic attachment path skips longitudinal objects; this guards the
  # direct maihda_interactions() call that bypasses it (analysis mode
  # "longitudinal", or a bare longitudinal maihda_model carrying longitudinal_info).
  is_longitudinal <-
    (inherits(object, "maihda_analysis") && identical(object$mode, "longitudinal")) ||
    (inherits(object, "maihda_model") && !is.null(object$longitudinal_info))
  if (is_longitudinal) {
    maihda_stop_longitudinal_scalar(
      "A scalar per-stratum interaction diagnostic",
      stratum_slope = maihda_object_stratum_slope(object))
  }

  if (inherits(object, "maihda_analysis")) {
    if (identical(object$mode, "two-model")) {
      if (is.null(object$model_adjusted)) {
        stop("This analysis has no adjusted model, so the pure-interaction BLUPs ",
             "are unavailable.", call. = FALSE)
      }
      summ <- object$summary_adjusted
      if (is.null(summ)) summ <- summary(object$model_adjusted)
      return(list(model = object$model_adjusted, summary = summ,
                  model_type = "adjusted (two-model)"))
    }
    # crossed-dimensions: the single model's interaction RE is the interaction.
    summ <- object$summary
    if (is.null(summ)) summ <- summary(object$model)
    return(list(model = object$model, summary = summ,
                model_type = "crossed-dimensions"))
  }

  if (inherits(object, "maihda_model")) {
    model_type <- if (!is.null(object$cc_info)) "crossed-dimensions" else "adjusted"
    if (is.null(object$cc_info)) {
      maihda_warn_if_not_adjusted(object)
    }
    return(list(model = object, summary = summary(object), model_type = model_type))
  }

  stop("'object' must be a maihda_analysis (from maihda()) or a maihda_model ",
       "(from fit_maihda()).", call. = FALSE)
}

# Warn when a bare maihda_model looks like a NULL model -- the dimensions' additive
# main effects are not all in the fixed part -- so its stratum random effects are
# the total between-stratum deviation (additive + interaction), not the pure
# interaction the diagnostic claims. Mirrors maihda()'s dimension-present logic.
maihda_warn_if_not_adjusted <- function(model) {
  sv <- model$strata_vars
  if (is.null(sv) || length(sv) < 2) {
    return(invisible(NULL))
  }
  expected <- tryCatch(
    maihda_adjusted_terms(sv, model$strata_autobin_info, model$original_data)$terms,
    error = function(e) sv)
  fixed_terms <- tryCatch(
    attr(stats::terms(maihda_nobars(model$formula)), "term.labels"),
    error = function(e) character(0))
  # terms() backtick-quotes non-syntactic labels (e.g. `gender var`), so the
  # expected dimension main effects must be compared in their quoted form -- a
  # raw-name intersect would miss them and falsely warn that a fully-specified
  # adjusted model is a null model. Mirrors maihda_workflow.R's dim_terms_quoted.
  expected_quoted <- vapply(expected, maihda_quote_name, character(1))
  # A fixed interaction among the dimensions (e.g. gender * race, which the
  # main-effects check above passes because both main effects ARE present) is a
  # different corruption: the fixed interaction absorbs the intersectional effect, so
  # the stratum random effect is no longer the pure interaction this diagnostic
  # reports. Flag it first -- it is the more specific, more serious problem.
  flagged_int <- tryCatch(
    maihda_dimension_interaction_terms(model$formula, sv, expected),
    error = function(e) character(0))
  if (length(flagged_int) > 0) {
    warning("maihda_interactions(): the model's fixed part contains the interaction ",
            "term(s) ", paste(flagged_int, collapse = ", "), " among the stratum ",
            "dimensions (", paste(sv, collapse = ", "), "). That fixed interaction ",
            "absorbs the intersectional effect, so the stratum random effect is no ",
            "longer the PURE interaction this diagnostic reports (its BLUPs are ",
            "driven toward a singular boundary). Fit the adjusted model with only the ",
            "dimensions' ADDITIVE main effects (e.g. ", paste(sv, collapse = " + "),
            ", not ", paste(sv, collapse = " * "), "), or pass a maihda() result.",
            call. = FALSE)
  } else if (length(flagged_tf <- tryCatch(
      maihda_transformed_dimension_terms(model$formula, sv, expected),
      error = function(e) character(0))) > 0) {
    # A dimension written in transformed form (e.g. factor(gender)) fails the bare-name
    # presence test above, so without this branch the model would be reported as a NULL
    # model -- true in effect but misleading about the cause and about the remedy.
    tf_dim <- maihda_transformed_dimension_vars(flagged_tf, sv, expected)[1]
    warning("maihda_interactions(): the model's fixed part contains ",
            paste(flagged_tf, collapse = ", "), " -- a transformed appearance of the ",
            "stratum dimension(s) ",
            paste(maihda_transformed_dimension_vars(flagged_tf, sv, expected),
                  collapse = ", "), ". Only a bare ",
            "column name counts as a dimension's additive main effect, so this model ",
            "is treated as a NULL model: its stratum random effects are the TOTAL ",
            "between-stratum deviation, not the pure interaction. Transform the column ",
            "in the data (e.g. data$", tf_dim, " <- factor(data$", tf_dim, ")) and write ",
            "the bare name, or pass a maihda() result.", call. = FALSE)
  } else if (!all(expected_quoted %in% fixed_terms)) {
    warning("maihda_interactions(): this looks like a null model -- the stratum ",
            "dimensions' additive main effects (", paste(sv, collapse = ", "),
            ") are not all in the fixed part, so the stratum random effects capture ",
            "the TOTAL between-stratum deviation (additive + interaction), not the ",
            "pure interaction. Pass the adjusted model, or a maihda() result, for ",
            "the interaction diagnostic.", call. = FALSE)
  }
  invisible(NULL)
}

# The grouping factor whose random effect is the interaction: the named
# interaction group of a crossed-dimensions fit, or "stratum" otherwise.
maihda_interaction_group <- function(model) {
  if (!is.null(model$cc_info) && !is.null(model$cc_info$interaction_group)) {
    model$cc_info$interaction_group
  } else {
    "stratum"
  }
}

# Per-stratum sample sizes aligned to `strata`, from the model's refreshed
# strata_info$n; NA when unavailable.
maihda_interaction_strata_n <- function(model, strata) {
  info <- model$strata_info
  if (is.null(info) || !"stratum" %in% names(info) || !"n" %in% names(info)) {
    return(rep(NA_integer_, length(strata)))
  }
  as.integer(info$n[match(strata, as.character(info$stratum))])
}

# Exact posterior tail of the stratum interaction for a brms fit: per stratum the
# posterior median, a central credible interval at conf_level, and the probability
# of direction pd = max(P(BLUP > 0), P(BLUP < 0)) -- the share of the posterior on
# the side of zero the estimate actually favours, so pd is in [0.5, 1] and the SIGN
# is reported separately (the `direction` column). Uses the full random-effect draws
# (summary = FALSE) rather than a Gaussian approximation to the posterior SD.
maihda_interaction_brms_tail <- function(brmsfit, group, conf_level) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required to work with brms models. Please install it ",
         "with: install.packages('brms')", call. = FALSE)
  }
  re <- brms::ranef(brmsfit, summary = FALSE, groups = group)
  if (!group %in% names(re)) {
    stop("No '", group, "' random effects found in the brms model.", call. = FALSE)
  }
  arr <- re[[group]]  # [draws, levels, effects]
  effect_names <- dimnames(arr)[[3]]
  idx <- match(TRUE, effect_names %in% c("(Intercept)", "Intercept"))
  if (is.na(idx)) {
    idx <- if (length(effect_names) == 1) 1L else
      stop("The '", group, "' random effect must include an intercept.", call. = FALSE)
  }
  levels <- dimnames(arr)[[2]]
  draws_mat <- arr[, , idx]
  if (is.null(dim(draws_mat))) {
    draws_mat <- matrix(draws_mat, ncol = length(levels))
  }
  a <- 1 - conf_level
  data.frame(
    stratum = levels,
    interaction = apply(draws_mat, 2, stats::median),
    lower = apply(draws_mat, 2, stats::quantile, probs = a / 2, names = FALSE),
    upper = apply(draws_mat, 2, stats::quantile, probs = 1 - a / 2, names = FALSE),
    pd = apply(draws_mat, 2, maihda_pd),
    stringsAsFactors = FALSE
  )
}

# Conventional probability of direction for a posterior draw vector (cf.
# bayestestR::p_direction): the larger of the two one-sided tail masses, in
# [0.5, 1]. NOT P(d > 0), which would read ~0 for a strong NEGATIVE interaction
# even though its direction is near-certain -- the sign is carried separately by
# the `direction` column. NA for an all-non-finite vector.
maihda_pd <- function(d) {
  d <- d[is.finite(d)]
  if (length(d) == 0) return(NA_real_)
  max(mean(d > 0), mean(d < 0))
}

# Validate and normalise the equivalence/ROPE argument of maihda_interactions():
# NULL stays NULL; a single positive d -> c(-d, d); a length-2 c(lower, upper) with
# lower < upper passes through. Anything else is an error.
maihda_normalize_rope <- function(rope) {
  if (is.null(rope)) return(NULL)
  if (!is.numeric(rope) || any(!is.finite(rope))) {
    stop("'rope' must be NULL, a single positive number, or a finite c(lower, upper).",
         call. = FALSE)
  }
  if (length(rope) == 1L) {
    if (rope <= 0) stop("A single-number 'rope' must be a positive half-width.", call. = FALSE)
    return(c(-abs(rope), abs(rope)))
  }
  if (length(rope) == 2L) {
    if (rope[1] >= rope[2]) stop("'rope = c(lower, upper)' must have lower < upper.", call. = FALSE)
    return(rope)
  }
  stop("'rope' must have length 1 (a half-width) or 2 (lower, upper).", call. = FALSE)
}

# Compute and attach the per-stratum interaction diagnostic to a fitted object
# (a maihda_analysis or a maihda_model), honouring the `interactions` request:
# FALSE/NULL skips; TRUE uses maihda_interactions()'s own default correction (BH);
# a character p.adjust method name uses that correction. The longitudinal interaction
# is a trajectory (intercept + slope), for which the scalar per-stratum diagnostic is
# undefined, so it is skipped. A genuine error degrades to NULL (rather than
# breaking the fit) but is re-emitted as a warning carrying the original message,
# so a requested-but-failed diagnostic is never dropped silently. (The "looks like
# a null model" warning on the opt-in fit_maihda path is informative and left to
# surface; maihda() never triggers it.)
maihda_attach_interactions <- function(object, interactions, conf_level = 0.95) {
  is_longitudinal <- identical(object$mode, "longitudinal") ||
    !is.null(object$longitudinal_info)
  if (is.null(interactions) || isFALSE(interactions) || is_longitudinal) {
    object$interactions <- NULL
    return(object)
  }
  if (!isTRUE(interactions)) {
    if (!is.character(interactions) || length(interactions) != 1L) {
      stop("'interactions' must be TRUE, FALSE, or a multiple-comparison method ",
           "name accepted by p.adjust (e.g. \"BH\", \"none\").", call. = FALSE)
    }
    interactions <- match.arg(interactions, c("none", stats::p.adjust.methods))
  }
  object$interactions <- maihda_try_optional(
    if (isTRUE(interactions)) {
      maihda_interactions(object, conf_level = conf_level)
    } else {
      maihda_interactions(object, conf_level = conf_level, adjust = interactions)
    },
    "Interaction diagnostics")
  object
}

# One-line interaction summary for a print method, naming the multiplicity stance
# actually used so an uncorrected scan is never silently read as corrected.
maihda_print_interactions_line <- function(ints, indent = "") {
  if (is.null(ints) || !inherits(ints, "maihda_interactions")) return(invisible(NULL))
  n_flag <- attr(ints, "n_flagged"); n_str <- attr(ints, "n_strata")
  adjust <- attr(ints, "adjust"); conf <- attr(ints, "conf_level")
  engine <- attr(ints, "engine")
  # Subsetting drops these attributes (see print.maihda_interactions); fall back to
  # values recomputed from the data frame so the one-line summary never errors.
  if (length(n_flag) != 1L || is.na(n_flag)) {
    n_flag <- if ("flagged" %in% names(ints)) sum(ints$flagged, na.rm = TRUE) else NA_integer_
  }
  if (length(n_str) != 1L) n_str <- nrow(ints)
  conf_pct <- if (is.null(conf)) 95 else conf * 100
  basis <- if (identical(engine, "brms")) {
    sprintf("%.0f%% credible interval", conf_pct)
  } else if (is.null(adjust) || identical(adjust, "none")) {
    sprintf("%.0f%% interval, no multiplicity correction", conf_pct)
  } else {
    sprintf("%.0f%% interval, %s-adjusted", conf_pct, adjust)
  }
  pal <- maihda_palette()
  count <- if (isTRUE(n_flag > 0)) {
    pal$bold(pal$accent(as.character(n_flag)))
  } else {
    pal$muted(as.character(n_flag))
  }
  cat(sprintf("%s%s %s of %d strata flagged (%s)\n", indent,
              pal$bold("Intersectional interactions:"), count, n_str,
              pal$muted(basis)))
  if (isTRUE(n_flag > 0)) {
    top <- ints[ints$flagged %in% TRUE, , drop = FALSE]
    top <- top[order(-abs(top$interaction)), , drop = FALSE][1, ]
    cat(sprintf("%s  strongest: %s (%s, %s)\n", indent, top$label,
                pal$accent(sprintf("%+.3f", top$interaction)), top$direction))
  }
  if ((is.null(adjust) || identical(adjust, "none")) &&
      !identical(engine, "brms") && isTRUE(n_str > 1)) {
    cat(pal$muted(sprintf(
      "%s  uncorrected across %d strata; maihda_interactions(x, adjust = \"BH\") for an FDR screen\n",
      indent, n_str)))
  }
  invisible(NULL)
}
