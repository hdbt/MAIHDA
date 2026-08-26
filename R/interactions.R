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
#' is only exact there.
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
#' @param rope Optional equivalence region (a "smallest interaction of interest")
#'   for an "is the interaction \emph{negligible}?" reading (Schuirmann 1987;
#'   Kruschke 2018), on the link (latent) scale. \code{NULL} (default) gives
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
#'   \code{interaction} (the BLUP), \code{lower}/\code{upper} (the interval),
#'   \code{flagged} (logical), and \code{direction} (\code{"above"}/\code{"below"}
#'   the additive expectation). Frequentist fits add \code{se} and \code{p_value}
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
#'   \code{scale} and \code{singular}.
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
#' }
#'
#' @export
#' @importFrom stats qnorm pnorm quantile median p.adjust terms p.adjust.methods
#' @importFrom reformulas nobars
maihda_interactions <- function(object, conf_level = 0.95, adjust = "BH",
                                rope = NULL, ...) {
  adjust_was_set <- !missing(adjust)
  resolved <- maihda_resolve_interaction_model(object)
  model <- resolved$model
  summary_obj <- resolved$summary
  model_type <- resolved$model_type

  conf_level <- maihda_validate_conf_level(conf_level)
  adjust <- match.arg(adjust, c("none", stats::p.adjust.methods))
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
  attr(out, "scale") <- "link"
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
  if (!is.na(model_type)) {
    cat(pal$muted(sprintf("Model: %s; interaction on the link (latent) scale.\n",
                          model_type)))
  } else {
    cat(pal$muted("Interaction reported on the link (latent) scale.\n"))
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
    print(utils::head(as.data.frame(x), 10), row.names = FALSE, digits = 4)
  } else {
    flagged_rows <- x[x$flagged %in% TRUE, , drop = FALSE]
    if (nrow(flagged_rows) == 0) {
      cat(pal$muted(
        "No strata show interaction credibly different from zero at this level.\n"))
    } else {
      print(utils::head(as.data.frame(flagged_rows), 10), row.names = FALSE, digits = 4)
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
