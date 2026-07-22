# Information criteria for MAIHDA models, for comparing model structures.
#
# The VPC/PCV answer "how much intersectional inequality is there", but choosing
# between alternative model STRUCTURES (which covariates, which strata dimensions,
# Gaussian vs ordinal, ...) is a relative-fit question that the variance summaries
# do not address. This file surfaces the standard fit criteria for each engine:
# AIC / BIC for the likelihood engines (lme4, ordinal::clmm) and the Bayesian
# WAIC / LOOIC for brms. They are wired into compare_maihda() and available
# directly via maihda_ic().

#' Information criteria for MAIHDA models
#'
#' @description
#' Reports the relative-fit information criteria for one or more MAIHDA models, to
#' help choose between model \emph{structures} (different covariate sets, strata
#' definitions, or families) -- a question the VPC/ICC and PCV do not address. The
#' criteria reported depend on the engine: \strong{AIC} and \strong{BIC} for the
#' likelihood engines (\code{lme4}, and \code{ordinal::clmm}), and the Bayesian
#' \strong{WAIC} and \strong{LOOIC} (leave-one-out information criterion) for
#' \code{brms}. Lower is better for all four.
#'
#' @details
#' \strong{REML vs ML.} \code{lmer} fits Gaussian models by REML by default, and a
#' REML log-likelihood (hence its AIC/BIC) is \emph{not} comparable across models
#' with different fixed effects -- exactly the canonical MAIHDA null-vs-adjusted
#' comparison. When more than one model is supplied, \code{maihda_ic()} therefore
#' refits any REML \code{lmer} model with maximum likelihood
#' (\code{\link[lme4]{refitML}}) before computing AIC/BIC, matching the behaviour of
#' \code{anova()} on \code{lme4} models; the \code{estimator} column records when
#' this happened. For a single model the criterion is reported as fitted (the
#' \code{estimator} column then reads \code{"REML"}).
#'
#' \strong{Comparability.} Like the VPC, information criteria are only comparable
#' across models fitted to the \emph{same} analytic sample (same rows and outcome)
#' with the \emph{same} weights -- prior (precision) weights and sampling
#' (design) weights each change which likelihood, or pseudo-likelihood, is being
#' maximised, so the criteria of a weighted and an unweighted fit of the identical
#' model are not on a common scale. AIC/BIC additionally require the same response
#' distribution -- they are not comparable across families (e.g. a Gaussian vs a
#' Poisson fit), nor between the likelihood engines and \code{brms} (AIC/BIC vs
#' WAIC/LOOIC are different scales). When the supplied models differ in any of
#' these respects \code{maihda_ic()} warns and omits the \code{delta} column,
#' still reporting each model's own criteria; \code{\link{compare_maihda}} warns
#' on the same grounds.
#'
#' \strong{Predictive target of the Bayesian criteria.} \code{brms::waic()} and
#' \code{brms::loo()} are computed from pointwise log-likelihoods
#' \emph{conditional on the fitted random effects}, so WAIC/LOOIC assess
#' prediction of new observations \emph{within the strata (and persons or
#' contexts) already represented in the data} -- not performance for a new,
#' unseen intersectional stratum. Choosing between strata definitions on LOOIC
#' therefore compares conditional predictive fit; generalisation to new strata
#' is a leave-one-group-out cross-validation question (e.g.
#' \code{brms::kfold()} with \code{group = "stratum"}), which this package does
#' not wrap. AIC/BIC for the likelihood engines are instead computed from the
#' \emph{marginal} likelihood (random effects integrated out) -- a further
#' reason the likelihood and Bayesian criteria are never comparable with each
#' other.
#'
#' \strong{Design-weighted fits.} For the \code{wemix} (design-weighted) engine the
#' criteria are reported as \code{NA}: a pseudo-likelihood with sampling weights does
#' not define a standard AIC/BIC. A \code{brms} fit with \code{sampling_weights} is
#' treated the same way: its sampling weights enter as likelihood weights, giving a
#' pseudo-posterior whose weighted pointwise log-likelihoods are not log predictive
#' densities, so WAIC/LOOIC are likewise reported as \code{NA} (the
#' \code{estimator} column reads \code{"Bayesian (weighted pseudo-posterior)"}).
#'
#' @param ... One or more \code{maihda_model} objects (from \code{\link{fit_maihda}})
#'   or \code{maihda_analysis} objects (from \code{\link{maihda}}). A
#'   \code{maihda_analysis} contributes its null model and, when present, its
#'   adjusted model as separate rows.
#' @param model_names Optional character vector of names, one per \code{...}
#'   argument. A \code{maihda_analysis} argument's null/adjusted rows are suffixed
#'   from its name.
#'
#' @return A \code{data.frame} of class \code{maihda_ic} with one row per model and
#'   the columns that apply: \code{model}, \code{n} (analytic sample size),
#'   \code{estimator}, \code{df} (number of parameters; likelihood engines),
#'   \code{logLik}, \code{AIC}, \code{BIC} (likelihood engines), \code{WAIC},
#'   \code{LOOIC} (brms), and -- when more than one model is supplied -- \code{delta}
#'   (the difference from the best model on the primary criterion: AIC for the
#'   likelihood engines, LOOIC for brms). Columns that are entirely \code{NA} across
#'   the supplied models are dropped.
#'
#' @seealso \code{\link{compare_maihda}}, which reports these criteria alongside the
#'   VPC/ICC, and \code{\link{calculate_pcv}} for the variance decomposition.
#'
#' @examples
#' \donttest{
#' strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#' null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
#' adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)
#'
#' # AIC/BIC for two nested structures (REML lmer fits are ML-refitted first)
#' maihda_ic(null_model, adj_model, model_names = c("Null", "Adjusted"))
#'
#' # Or straight from a one-call maihda() analysis (null + adjusted rows)
#' a <- maihda(health_outcome ~ age + gender + race + (1 | gender:race),
#'             data = maihda_sim_data)
#' maihda_ic(a)
#' }
#'
#' @export
#' @importFrom stats AIC BIC logLik
maihda_ic <- function(..., model_names = NULL) {
  args <- list(...)
  if (length(args) == 0) {
    stop("maihda_ic() needs at least one maihda_model or maihda_analysis object.",
         call. = FALSE)
  }

  # Default per-argument names; a maihda_analysis expands to Null/Adjusted rows
  # below, suffixed from the argument's name.
  if (is.null(model_names)) {
    model_names <- paste0("Model", seq_along(args))
  } else if (length(model_names) != length(args)) {
    stop("Length of 'model_names' (", length(model_names), ") must match the number ",
         "of model arguments (", length(args), ").", call. = FALSE)
  }

  # Flatten the arguments into a list of (name, maihda_model) pairs.
  named_models <- list()
  for (i in seq_along(args)) {
    a <- args[[i]]
    base <- model_names[i]
    if (inherits(a, "maihda_analysis")) {
      named_models[[length(named_models) + 1L]] <-
        list(name = paste0(base, " (Null)"), model = a$model)
      if (!is.null(a$model_adjusted)) {
        named_models[[length(named_models) + 1L]] <-
          list(name = paste0(base, " (Adjusted)"), model = a$model_adjusted)
      }
    } else if (inherits(a, "maihda_model")) {
      named_models[[length(named_models) + 1L]] <- list(name = base, model = a)
    } else {
      stop("maihda_ic() arguments must be maihda_model or maihda_analysis objects; ",
           "argument ", i, " is of class '", paste(class(a), collapse = "/"), "'.",
           call. = FALSE)
    }
  }

  # With more than one model the table IS a comparison, so refit REML lmer fits
  # with ML (see Details) before reading AIC/BIC.
  ml <- length(named_models) > 1L

  rows <- lapply(named_models, function(nm) {
    one <- maihda_ic_one(nm$model, ml = ml)
    data.frame(model = nm$name, one, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  # Primary criterion for the delta column: AIC where the likelihood engines are
  # used, otherwise the Bayesian LOOIC (then WAIC). Only meaningful with >1 row.
  primary <- maihda_ic_primary(out)
  if (nrow(out) > 1L && !is.na(primary)) {
    # A delta is only meaningful across models fitted to the SAME analytic sample
    # and, for AIC/BIC, the SAME response distribution (see Details) -- and it is
    # never meaningful across the likelihood/Bayesian divide (AIC/BIC and WAIC/LOOIC
    # are different scales). maihda_ic() does not otherwise enforce this, so a direct
    # call could rank a Gaussian against a Poisson fit, fits on different rows, or an
    # lme4 fit against a brms fit -- each with a seemingly meaningful delta (the
    # cross-scale case picking AIC as the primary criterion and leaving the Bayesian
    # row a bare NA). Withhold the delta (and say why) when the supplied models are
    # not mutually comparable; the per-model criteria are still reported. The
    # canonical null-vs-adjusted comparison (same outcome/sample/family, differing
    # only in fixed effects) is comparable, so its delta is unaffected.
    delta_issues <- maihda_ic_delta_issues(lapply(named_models, function(x) x$model))
    # Whether the POPULATED criterion columns span both scales -- the same guard
    # compare_maihda() applies, which the outcome/family/sample check above does NOT
    # catch (a same-family lme4-vs-brms comparison agrees on all three). Test only the
    # columns that carry a finite value: `out` still has all four IC columns here,
    # most of them all-NA, so an all-lme4 table is not misread as mixed.
    populated_ic <- Filter(function(col) any(is.finite(out[[col]])),
                           intersect(c("AIC", "BIC", "WAIC", "LOOIC"), names(out)))
    if (maihda_ic_spans_scales(populated_ic)) {
      delta_issues <- c(delta_issues,
                        "scale (likelihood AIC/BIC vs Bayesian WAIC/LOOIC)")
    }
    # An intended ML refit (ml = TRUE whenever there is >1 model) that FELL BACK to
    # REML leaves that row's estimator at "REML" -- a successful refit reads "ML
    # (refit from REML)", a glmer fit "ML", so a "REML" row here means refitML() failed.
    # AIC/BIC then mix REML and ML bases, and a delta across the differing fixed effects
    # is not comparable; withhold it rather than rank on a mixed basis.
    if (isTRUE(ml) && any(out$estimator == "REML", na.rm = TRUE)) {
      delta_issues <- c(delta_issues,
                        "estimation basis (an ML refit failed, leaving a REML fit)")
    }
    if (length(delta_issues) > 0) {
      warning("maihda_ic(): the models differ in ",
              paste(delta_issues, collapse = " and "),
              ", so a delta is not meaningful and is omitted -- information ",
              "criteria are only comparable across models fitted to the same ",
              "analytic sample with the same weights and, for AIC/BIC, the same ",
              "response distribution, and are never comparable across the ",
              "likelihood/Bayesian divide. The per-model criteria are still ",
              "reported.", call. = FALSE)
    } else {
      vals <- out[[primary]]
      if (any(is.finite(vals))) {
        out$delta <- vals - min(vals, na.rm = TRUE)
      }
    }
  }

  # Drop value columns that are entirely NA (e.g. no WAIC column for an all-lme4
  # comparison), keeping the bookkeeping columns.
  keep_always <- c("model", "n", "estimator", "delta")
  value_cols <- setdiff(names(out), keep_always)
  for (col in value_cols) {
    if (all(is.na(out[[col]]))) out[[col]] <- NULL
  }
  rownames(out) <- NULL

  attr(out, "ic_primary") <- primary
  class(out) <- c("maihda_ic", "data.frame")
  out
}

# Evaluate a brms IC call (waic/loo) keeping brms's progress MESSAGES quiet, but
# CAPTURING its warnings and re-emitting them as a single branded summary. brms's
# own reliability warnings -- high Pareto-k for PSIS-LOO, low effective sample
# size, large p_waic -- are exactly the signal that a criterion (and therefore any
# model ranking built on it) may be untrustworthy, so they must be surfaced, not
# swallowed. `what` names the criterion for the message.
maihda_ic_quiet_but_warn <- function(expr, what) {
  warns <- character(0)
  val <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warns) > 0) {
    warning(sprintf("%s reliability diagnostics flagged this fit: %s", what,
                    paste(unique(trimws(warns)), collapse = " | ")),
            call. = FALSE)
  }
  val
}

#' Information criteria for a single MAIHDA model
#'
#' Internal worker for \code{\link{maihda_ic}}: returns a one-row data frame of the
#' fit criteria for one \code{maihda_model}, dispatched on the fitted object's class
#' (mirroring \code{maihda_fit_diagnostics}).
#'
#' @param model A \code{maihda_model}.
#' @param ml Logical; for a REML \code{lmer} fit, refit with ML via
#'   \code{\link[lme4]{refitML}} before reading AIC/BIC (used when comparing
#'   models that may differ in fixed effects).
#' @return A one-row data frame with \code{n}, \code{estimator}, \code{df},
#'   \code{logLik}, \code{AIC}, \code{BIC}, \code{WAIC}, \code{LOOIC} (NA where not
#'   applicable to the engine).
#' @keywords internal
maihda_ic_one <- function(model, ml = FALSE) {
  if (!inherits(model, "maihda_model")) {
    stop("'model' must be a maihda_model object from fit_maihda().", call. = FALSE)
  }
  fm <- model$model
  na_real <- NA_real_
  row <- list(
    n = NA_integer_, estimator = NA_character_, df = na_real,
    logLik = na_real, AIC = na_real, BIC = na_real,
    WAIC = na_real, LOOIC = na_real
  )
  # maihda_wrapper_nobs() falls back to nrow(model$data) for engines whose fitted
  # object has no nobs() method (WeMixResults), so the IC table reports the
  # analytic n rather than NA (glance() already reports it via the same frame).
  n <- maihda_wrapper_nobs(model)
  row$n <- if (is.finite(n)) as.integer(n) else NA_integer_

  if (inherits(fm, "merMod")) {
    ml_used <- FALSE
    reml <- tryCatch(isTRUE(lme4::isREML(fm)), error = function(e) FALSE)
    fit_for_ic <- fm
    if (isTRUE(ml) && reml) {
      fit_for_ic <- tryCatch(lme4::refitML(fm), error = function(e) fm)
      ml_used <- !identical(fit_for_ic, fm)
    }
    ll <- tryCatch(stats::logLik(fit_for_ic), error = function(e) NULL)
    if (!is.null(ll)) {
      row$logLik <- as.numeric(ll)
      row$df <- attr(ll, "df")
    }
    row$AIC <- tryCatch(as.numeric(stats::AIC(fit_for_ic)), error = function(e) na_real)
    row$BIC <- tryCatch(as.numeric(stats::BIC(fit_for_ic)), error = function(e) na_real)
    row$estimator <- if (ml_used) "ML (refit from REML)" else if (reml) "REML" else "ML"

  } else if (inherits(fm, "clmm")) {
    # ordinal::clmm is maximum-likelihood; AIC/BIC dispatch through the stats
    # generics on its logLik method.
    ll <- tryCatch(stats::logLik(fm), error = function(e) NULL)
    if (!is.null(ll)) {
      row$logLik <- as.numeric(ll)
      row$df <- attr(ll, "df")
    }
    row$AIC <- tryCatch(as.numeric(stats::AIC(fm)), error = function(e) na_real)
    row$BIC <- tryCatch(as.numeric(stats::BIC(fm)), error = function(e) na_real)
    row$estimator <- "ML"

  } else if (inherits(fm, "brmsfit")) {
    if (!is.null(model$sampling_weights)) {
      # Sampling weights enter brms as likelihood weights, giving a
      # PSEUDO-posterior (see fit_maihda): the weighted pointwise
      # log-likelihood terms are not log predictive densities for any
      # observation process, so they define no standard WAIC/LOOIC -- exactly
      # as the wemix pseudo-likelihood defines no standard AIC/BIC. Report NA
      # (silently, matching the wemix branch below; fit_maihda already
      # messaged the pseudo-posterior caveat at fit time).
      row$estimator <- "Bayesian (weighted pseudo-posterior)"
    } else {
      # Bayesian analogues: WAIC and the leave-one-out IC (LOOIC). brms emits
      # progress messages (kept quiet) alongside genuine reliability warnings (high
      # Pareto-k, low ESS); the latter are surfaced via maihda_ic_quiet_but_warn
      # rather than suppressed, so a criterion is never reported as if reliable when
      # its own diagnostics say otherwise.
      row$estimator <- "Bayesian"
      if (requireNamespace("brms", quietly = TRUE)) {
        row$WAIC <- tryCatch({
          w <- maihda_ic_quiet_but_warn(brms::waic(fm), "WAIC")
          as.numeric(w$estimates["waic", "Estimate"])
        }, error = function(e) na_real)
        row$LOOIC <- tryCatch({
          l <- maihda_ic_quiet_but_warn(brms::loo(fm), "PSIS-LOO")
          as.numeric(l$estimates["looic", "Estimate"])
        }, error = function(e) na_real)
      }
    }

  } else if (inherits(fm, "WeMixResults")) {
    # Design-weighted pseudo-maximum-likelihood does not define a standard AIC/BIC;
    # report NA (silently -- a wemix comparison must not add a warning).
    row$estimator <- "pseudo-ML (weighted)"
  }

  as.data.frame(row, stringsAsFactors = FALSE)
}

# Choose the criterion the delta column is computed on: AIC for the likelihood
# engines, else the Bayesian LOOIC, else WAIC, else BIC. Returns NA when no
# criterion column is populated.
maihda_ic_primary <- function(df) {
  has <- function(col) col %in% names(df) && any(is.finite(df[[col]]))
  if (has("AIC")) return("AIC")
  if (has("LOOIC")) return("LOOIC")
  if (has("WAIC")) return("WAIC")
  if (has("BIC")) return("BIC")
  NA_character_
}

# The ways a set of models differ that make a delta between their information
# criteria meaningless: a differing outcome, family/link, analytic sample, or set
# of weights. A delta is a difference of criteria, so it inherits exactly the
# comparability requirements the criteria themselves have (same sample; for
# AIC/BIC same response distribution) -- the caveat spelled out in the maihda_ic
# Details and enforced by compare_maihda()'s warning. Returns the human-readable
# list of differences (empty when the models are mutually comparable). Uses the
# same response/family/sample/weight fingerprints as compare_maihda() so the two
# agree, and deliberately does NOT compare fixed effects: the canonical
# null-vs-adjusted comparison differs only there and must keep its delta.
maihda_ic_delta_issues <- function(models) {
  if (length(models) < 2L) {
    return(character(0))
  }
  responses <- vapply(models, function(m) {
    paste(deparse(m$formula[[2]]), collapse = "")
  }, character(1))
  fam_keys <- vapply(models, maihda_model_family_key, character(1))
  nobs_vec <- vapply(models, function(m) {
    n <- maihda_wrapper_nobs(m)
    if (is.finite(n)) as.integer(n) else NA_integer_
  }, integer(1))
  # Sort the id set so a reordered-but-identical sample is not read as a different
  # one; the response fingerprint is likewise row-order independent.
  row_keys <- vapply(models, function(m) {
    rid <- maihda_wrapper_row_ids(m)
    if (is.null(rid)) NA_character_ else paste(sort(rid), collapse = "\r")
  }, character(1))
  response_keys <- vapply(models, maihda_wrapper_response_fingerprint, character(1))
  # Weights define WHICH likelihood is being maximised, so criteria built on
  # different weights are not on a common scale at all -- an unweighted fit and a
  # weighted fit of the identical model report different AICs, and their
  # difference measures the weighting, not the models. The same two fingerprints
  # calculate_pcv() hard-errors on and compare_maihda() warns about; PRIOR
  # (precision) weights live on the fitted object, SAMPLING (design) weights on
  # the maihda_model wrapper.
  weight_keys <- vapply(models, function(m) maihda_weight_fingerprint(m$model),
                        character(1))
  sampling_keys <- vapply(models, maihda_sampling_weight_fingerprint, character(1))

  issues <- character(0)
  if (length(unique(responses)) > 1L) {
    issues <- c(issues, paste0("outcomes (", paste(unique(responses), collapse = ", "), ")"))
  }
  if (length(unique(fam_keys)) > 1L) {
    issues <- c(issues, paste0("families/links (", paste(unique(fam_keys), collapse = ", "), ")"))
  }
  if (length(unique(stats::na.omit(weight_keys))) > 1L) {
    issues <- c(issues, "prior weights")
  }
  if (length(unique(stats::na.omit(sampling_keys))) > 1L) {
    issues <- c(issues, "sampling weights")
  }
  sample_differs <- length(unique(stats::na.omit(nobs_vec))) > 1L ||
    length(unique(stats::na.omit(row_keys))) > 1L ||
    length(unique(stats::na.omit(response_keys))) > 1L
  if (sample_differs) {
    issues <- c(issues, paste0("analytic sample (n = ", paste(nobs_vec, collapse = ", "), ")"))
  }
  issues
}

# TRUE when a set of information-criterion column names spans BOTH the likelihood
# scale (AIC/BIC, from lme4/ordinal) and the Bayesian scale (WAIC/LOOIC, from
# brms). Such a table puts criteria from different scales side by side, which are
# not comparable to each other (see the maihda_ic Details). This happens for a
# same-family lme4-vs-brms comparison -- a case the family/link check in
# compare_maihda() does not flag -- so compare_maihda() uses this to warn.
maihda_ic_spans_scales <- function(cols) {
  any(c("AIC", "BIC") %in% cols) && any(c("WAIC", "LOOIC") %in% cols)
}

#' Print MAIHDA information criteria
#'
#' @param x A \code{maihda_ic} object from \code{\link{maihda_ic}}.
#' @param ... Additional arguments (not used).
#' @return No return value, called for side effects.
#' @export
print.maihda_ic <- function(x, ...) {
  pal <- maihda_palette()
  cat(pal$bold("MAIHDA Information Criteria"), "\n", sep = "")
  cat("===========================\n\n")
  print(as.data.frame(x), row.names = FALSE, digits = 4)

  # `primary` is NULL when a subset dropped the ic_primary attribute (older objects,
  # or a `[` that did not preserve it); guard the length so the `if` never tests a
  # zero-length `!is.na()` (which errors after the table has already printed).
  primary <- attr(x, "ic_primary")
  if ("delta" %in% names(x) && length(primary) == 1L && !is.na(primary)) {
    cat(pal$muted(sprintf("\ndelta = difference from the best model on %s (lower is better).\n",
                primary)))
  }
  if (!is.null(x$estimator) && any(grepl("refit from REML", x$estimator))) {
    cat(pal$muted(paste0("REML lmer fit(s) were refitted with ML so AIC/BIC are comparable across ",
        "different fixed effects.\n")))
  }
  if (any(c("AIC", "BIC", "WAIC", "LOOIC") %in% names(x))) {
    cat(pal$muted(paste0("Information criteria are only comparable across models fitted to the same ",
        "analytic sample with the same weights (and, for AIC/BIC, the same family).\n")))
  }
  invisible(x)
}

#' Subset MAIHDA information criteria
#'
#' Indexing method for \code{\link{maihda_ic}} results that preserves the
#' \code{ic_primary} metadata attribute (which names the criterion the
#' \code{delta} column is computed on) after a row/column subset. Plain
#' \code{[.data.frame} keeps the \code{maihda_ic} class but drops that attribute,
#' which would leave \code{print.maihda_ic()} testing a dropped (\code{NULL})
#' attribute as a scalar condition and erroring after the table prints.
#'
#' @param x A \code{maihda_ic} object.
#' @param ... Row/column indices forwarded to \code{[.data.frame}.
#' @return A \code{maihda_ic} with \code{ic_primary} carried over, or a bare vector
#'   when the selection drops to a single column.
#' @keywords internal
#' @export
`[.maihda_ic` <- function(x, ...) {
  out <- NextMethod()
  # A single-column selection with drop = TRUE returns a bare vector; leave it
  # unclassed. Otherwise restore the class and the ic_primary attribute the print
  # method relies on.
  if (is.data.frame(out)) {
    attr(out, "ic_primary") <- attr(x, "ic_primary")
    class(out) <- class(x)
  }
  out
}
