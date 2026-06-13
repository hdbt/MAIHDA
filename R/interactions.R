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
#' there genuine intersectionality": a flagged stratum is one whose joint identity
#' produces an outcome the additive parts do not.
#'
#' @details
#' \strong{It must be read off the adjusted model.} Only when the dimensions'
#' additive main effects are in the model (the \emph{adjusted} model of the
#' two-model decomposition, or the crossed-dimensions model) does the stratum
#' random effect isolate the \emph{pure interaction}. On a null model the stratum
#' random effect is the total between-stratum deviation (additive + interaction),
#' so passing one is flagged with a warning. Passing a \code{\link{maihda}} result
#' uses the right model automatically.
#'
#' \strong{Frequentist vs. Bayesian evidence.} For the frequentist engines
#' (\code{lme4}, \code{wemix}, \code{ordinal}) the flag comes from the BLUP's
#' conditional standard error: a Wald interval at \code{conf_level} and a two-sided
#' p-value, with an optional multiplicity correction (\code{adjust}). For
#' \code{brms} the full posterior is already available, so the \emph{exact}
#' posterior tail is used -- a credible interval at \code{conf_level} and the
#' probability of direction \code{pd = P(BLUP > 0)} -- and \code{adjust} is not
#' applied (the Bayesian answer is multiplicity-free).
#'
#' \strong{On multiplicity (why \code{adjust = "none"} is the default).} The
#' stratum BLUP is already a shrunken (partially pooled) estimate, which provides
#' built-in protection against spurious extremes (Gelman, Hill & Yajima 2012); the
#' relevant risk is sign/magnitude (Type S/M) error rather than family-wise Type I
#' error (Gelman & Carlin 2014). Multiplicity corrections are therefore optional
#' and, applied to shrunken estimates, conservative -- \code{"none"} is the
#' default. When an explicit error-rate story is required for screening many
#' strata, \strong{FDR (\code{adjust = "BH"})} matches the goal (discovery), not
#' the family-wise methods (\code{"bonferroni"}/\code{"holm"}), which remain
#' available but answer a question this exploratory scan does not ask. Whatever the
#' choice, treat the flags as exploratory rather than confirmatory.
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
#'   (frequentist engines only): \code{"none"} (default) or any method accepted by
#'   \code{\link[stats]{p.adjust}} (e.g. \code{"BH"}, \code{"holm"},
#'   \code{"bonferroni"}). Ignored for \code{brms} (with a message), which uses the
#'   posterior tail directly.
#' @param ... Currently unused.
#'
#' @return An object of class \code{maihda_interactions} (a data frame), one row
#'   per stratum, sorted flagged-first then by \code{abs(interaction)}. Columns
#'   common to every engine: \code{stratum}, \code{label}, \code{n} (stratum size),
#'   \code{interaction} (the BLUP), \code{lower}/\code{upper} (the interval),
#'   \code{flagged} (logical), and \code{direction} (\code{"above"}/\code{"below"}
#'   the additive expectation). Frequentist fits add \code{se} and \code{p_value}
#'   (and \code{p_adjusted} when \code{adjust != "none"}); \code{brms} adds
#'   \code{pd} (probability of direction). Attributes record \code{conf_level},
#'   \code{adjust}, \code{engine}, \code{model_type}, \code{n_strata},
#'   \code{n_flagged}, \code{scale} and \code{singular}.
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
#' @seealso \code{\link{maihda}}, \code{\link{calculate_pvc}},
#'   \code{\link{summary.maihda_model}}; and \code{plot(\dots,
#'   highlight_interactions = TRUE)} to mark the flagged strata on the
#'   effect-decomposition / predicted / shrinkage plots.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#' a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race),
#'             data = maihda_health_data)
#' maihda_interactions(a)                 # which strata interact (95%, no correction)
#' maihda_interactions(a, adjust = "BH")  # FDR-controlled screen
#' }
#'
#' @export
#' @importFrom stats qnorm pnorm quantile median p.adjust terms p.adjust.methods
#' @importFrom reformulas nobars
maihda_interactions <- function(object, conf_level = 0.95, adjust = "none", ...) {
  resolved <- maihda_resolve_interaction_model(object)
  model <- resolved$model
  summary_obj <- resolved$summary
  model_type <- resolved$model_type

  conf_level <- maihda_validate_conf_level(conf_level)
  adjust <- match.arg(adjust, c("none", stats::p.adjust.methods))

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
    if (!identical(adjust, "none")) {
      message("maihda_interactions(): the Bayesian posterior tail is ",
              "multiplicity-free; 'adjust' is ignored for brms models.")
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

  # Flagged strata first, then by interaction magnitude (most extreme first).
  ord <- order(!out$flagged, -abs(out$interaction))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  singular <- isTRUE(summary_obj$diagnostics$singular) ||
    isTRUE(model$diagnostics$singular)

  attr(out, "conf_level") <- conf_level
  attr(out, "adjust") <- adjust
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
  conf <- attr(x, "conf_level")
  adjust <- attr(x, "adjust")
  engine <- attr(x, "engine")
  model_type <- attr(x, "model_type")
  n_strata <- attr(x, "n_strata")
  n_flagged <- attr(x, "n_flagged")

  cat("Strata with credibly non-zero intersectional interaction\n")
  cat("========================================================\n\n")

  evidence <- if (identical(engine, "brms")) {
    sprintf("%.0f%% credible interval; probability of direction", conf * 100)
  } else if (identical(adjust, "none")) {
    sprintf("%.0f%% interval; no multiplicity correction", conf * 100)
  } else {
    sprintf("%.0f%% interval; %s-adjusted p-values", conf * 100, adjust)
  }
  cat(sprintf("%d of %d strata flagged (%s).\n", n_flagged, n_strata, evidence))
  cat(sprintf("Model: %s; interaction on the link (latent) scale.\n", model_type))

  if (isTRUE(attr(x, "singular"))) {
    cat("\nNote: singular/boundary fit (between-stratum variance ~ 0); the BLUP SEs\n",
        "  collapse toward zero, so 'no flag' is not evidence of no interaction.\n",
        sep = "")
  }

  cat("\n")
  flagged_rows <- x[x$flagged %in% TRUE, , drop = FALSE]
  if (nrow(flagged_rows) == 0) {
    cat("No strata show interaction credibly different from zero at this level.\n")
  } else {
    print(utils::head(as.data.frame(flagged_rows), 10), row.names = FALSE, digits = 4)
    if (nrow(flagged_rows) > 10) {
      cat(sprintf("  ... and %d more flagged strata\n", nrow(flagged_rows) - 10))
    }
  }

  if (!identical(engine, "brms") && identical(adjust, "none") && n_strata > 1) {
    cat("\nFlagging many strata inflates false positives; for a screening error-rate\n",
        "  story use adjust = \"BH\" (FDR). Interaction BLUPs are shrunken estimates,\n",
        "  so correction is optional -- see ?maihda_interactions.\n", sep = "")
  } else {
    cat("\nInteraction BLUPs are shrunken (partially pooled) estimates; treat flags as\n",
        "  exploratory. See ?maihda_interactions.\n", sep = "")
  }

  invisible(x)
}

# ---- internal helpers -------------------------------------------------------

# Resolve the model + summary the interaction diagnostic should read, and label
# the model type. A maihda_analysis routes to its adjusted (two-model) or single
# crossed-dimensions model -- both carry the additive part, so the stratum random
# effect is the pure interaction. A bare maihda_model is used directly, with a
# guardrail warning when it looks like a null model.
maihda_resolve_interaction_model <- function(object) {
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
    attr(stats::terms(reformulas::nobars(model$formula)), "term.labels"),
    error = function(e) character(0))
  if (!all(expected %in% fixed_terms)) {
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
# of direction pd = P(BLUP > 0). Uses the full random-effect draws (summary =
# FALSE) rather than a Gaussian approximation to the posterior SD.
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
    pd = apply(draws_mat, 2, function(d) mean(d > 0)),
    stringsAsFactors = FALSE
  )
}
