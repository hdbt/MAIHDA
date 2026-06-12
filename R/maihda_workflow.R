#' Run a Complete MAIHDA Analysis
#'
#' A single high-level entry point that runs the standard two-model MAIHDA workflow
#' and returns one bundled object. It fits the \strong{null} model (covariates plus
#' the intersectional random intercept, \emph{excluding} the stratum dimensions' main
#' effects) and the \strong{adjusted} model (the null plus the additive main effects of
#' the stratum-defining dimensions), summarises the variance partition (VPC/ICC) of the
#' null model, and reports the \strong{PCV} -- the proportional change in
#' between-stratum variance from the null to the adjusted model, i.e. the additive
#' share of the intersectional inequality. When a higher-level grouping variable is
#' supplied it also compares this decomposition across that variable's levels.
#'
#' \strong{Binomial companions.} For a binary outcome the model summaries also carry
#' the discriminatory accuracy (AUC / C-statistic and Median Odds Ratio) -- the "DA"
#' in MAIHDA -- automatically, so the null model's strata-only AUC sits alongside its
#' VPC; set \code{response_vpc = TRUE} to add the (simulation-based) response-scale
#' VPC as well. These are read from \code{summary(x)} and the attached
#' \code{summary_adjusted}, and the headline \code{print()} shows the null-model AUC.
#'
#' This is a convenience wrapper around \code{\link{fit_maihda}},
#' \code{\link{calculate_pvc}}, \code{\link{summary.maihda_model}} and
#' \code{\link{compare_maihda_groups}}. It is \emph{intrinsically} a two-model
#' decomposition and has no single-model mode -- for a single fit (e.g. just the
#' null-model VPC / discriminatory accuracy), call \code{\link{fit_maihda}} directly.
#'
#' \strong{The dimensions' additive main effects.} You may write them in the formula --
#' the fully-specified, lme4-native adjusted model
#' \code{outcome ~ covars + var1 + var2 + (1 | var1:var2)} -- or omit them. Either way
#' the null \emph{excludes} the dimension main effects and the adjusted \emph{includes}
#' them: when the formula already lists them it is taken as the adjusted model and the
#' null is derived by dropping them; when they are missing \code{maihda()} adds them to
#' the adjusted model and emits a \code{message()} so the decomposition stays explicit.
#' The dimensions themselves are read from the random term: the shorthand
#' \code{(1 | var1:var2)} and \code{\link{make_strata}} both record them, and a numeric
#' dimension that \code{make_strata()} auto-binned enters the adjusted model as its
#' reconstructed tertile factor (matching the strata), not as a linear term. Because
#' \code{maihda()} is intrinsically a decomposition, it \strong{errors} (rather than
#' returning a null-only result) when it cannot build the adjusted model -- when the
#' dimensions cannot be recovered (a hand-built \code{stratum} column records none) or
#' there is only one dimension (no intersection to decompose). Use \code{\link{fit_maihda}}
#' for those single-model fits.
#'
#' @param formula A model formula, using either the intersectional shorthand
#'   \code{outcome ~ covars + (1 | var1:var2)} or \code{... + (1 | stratum)} when
#'   \code{data} already has a \code{stratum} column from \code{\link{make_strata}}.
#'   The dimensions' additive main effects may be listed (the fully-specified adjusted
#'   model) or omitted (added automatically, with a message); see Details.
#' @param data A data frame with the model variables (and the \code{group}
#'   variable, if used).
#' @param group Optional character string naming a higher-level grouping variable
#'   (e.g. \code{"country"}). When supplied, \code{\link{compare_maihda_groups}}
#'   is run and attached to the result. \code{group} runs a \emph{stratified}
#'   analysis (one independent model per level); to instead model a higher level
#'   \emph{jointly}, crossed with the strata, use \code{context}. The two are
#'   different designs, so supplying both errors.
#' @param context Optional character vector naming higher-level \emph{context}
#'   column(s) in \code{data} (e.g. \code{"school"}, \code{"hospital"},
#'   \code{"region"}), forwarded to \code{\link{fit_maihda}}. Each context enters
#'   every fitted model as a crossed random intercept --
#'   \code{outcome ~ covars + (1 | stratum) + (1 | context)} -- the
#'   \emph{contextual cross-classified MAIHDA} of the literature. The summaries
#'   then partition the unexplained variance into between-stratum vs.
#'   between-context vs. residual, the headline VPC/ICC becomes the between-stratum
#'   share \emph{net of} the context, and the PCV decomposition is computed with the
#'   context partialled out (the context random intercept is carried by both the
#'   null and the adjusted model). Cannot be combined with \code{group}; a context
#'   with few levels weakly identifies its variance (consider \code{engine = "brms"}).
#' @param engine Modeling engine, "lme4" (default) or "brms".
#' @param family Model family. Default "gaussian". As in \code{\link{fit_maihda}},
#'   a binary outcome is auto-detected when \code{family} is left at the default,
#'   and the same resolved family is then used for the group comparison so all
#'   models agree.
#' @param decomposition How to decompose the intersectional inequality into additive
#'   and interaction parts. \code{"two-model"} (default) is the classic MAIHDA
#'   approach: a null model and an adjusted model (the dimensions' additive main
#'   effects as \emph{fixed} effects), with the additive share read from the PCV
#'   (proportional change in between-stratum variance). \code{"crossed-dimensions"}
#'   fits a \strong{single} model that enters each dimension's additive main effect as
#'   a \emph{random} intercept -- \code{outcome ~ covars + (1 | dim1) + ... +
#'   (1 | stratum)} -- so each dimension's RE variance is its additive contribution and
#'   the intersection (\code{stratum}) RE variance is the interaction beyond additive;
#'   the additive and interaction \emph{shares} of the total between-strata variance are
#'   read directly from that one fit. The two modes target the same scientific question
#'   with different estimators, so their additive shares are conceptually parallel but
#'   not numerically identical. The crossed-dimensions additive share is a partial-pooling
#'   estimate: dimensions with few levels (e.g. a binary sex variable, whose variance is
#'   estimated from two groups) are poorly identified and often give a singular lme4 fit
#'   -- the \code{brms} engine handles this better. See Details.
#'   (\code{"cross-classified"} is accepted as a deprecated alias for
#'   \code{"crossed-dimensions"}, with a warning: in the MAIHDA literature
#'   \dQuote{cross-classified} refers to the contextual stratum-by-place model,
#'   which this package fits via \code{context}.)
#' @param autobin Logical passed to \code{\link{make_strata}}; tertile-bins numeric
#'   grouping variables. Default TRUE.
#' @param shared_strata Logical, forwarded to \code{\link{compare_maihda_groups}}
#'   when \code{group} is supplied: build strata once on the full data so VPCs are
#'   comparable across groups (TRUE, default) or rebuild them within each group.
#' @param min_group_n Minimum group size for the per-group comparison, forwarded
#'   to \code{\link{compare_maihda_groups}}. Default 30.
#' @param bootstrap Logical; compute parametric-bootstrap VPC confidence intervals
#'   (lme4 only) for both the overall summary and the per-group comparison.
#'   Default FALSE.
#' @param n_boot Number of bootstrap samples when \code{bootstrap = TRUE}.
#' @param conf_level Confidence level for bootstrap intervals. Default 0.95.
#' @param response_vpc Logical; for a binomial (lme4) outcome, also attach the
#'   response-scale VPC (\code{\link{maihda_vpc_response}}) to the model summaries.
#'   It is estimated by simulation, so it is opt-in (default \code{FALSE}) and uses
#'   \code{seed}. The discriminatory accuracy (AUC + MOR) is attached automatically
#'   for a binomial/Bernoulli outcome regardless of this flag (see Details).
#' @param seed Optional integer seed for the response-scale VPC simulation.
#' @param ... Additional arguments passed to \code{\link{fit_maihda}} (and on to
#'   \code{lmer}/\code{glmer}).
#'
#' @return An object of class \code{maihda_analysis}: a list with
#'   \item{model}{the fitted \code{maihda_model} (see \code{\link{fit_maihda}}); the
#'     \strong{null} model in \code{"two-model"} mode, or the single
#'     \strong{crossed-dimensions} model in \code{"crossed-dimensions"} mode}
#'   \item{summary}{the model's \code{maihda_summary} (VPC/ICC, variance components,
#'     stratum estimates; plus the additive/interaction \code{decomposition} in
#'     crossed-dimensions mode, and the stratum-vs-context \code{context} partition
#'     when \code{context} is supplied)}
#'   \item{model_adjusted}{the fitted \strong{adjusted} \code{maihda_model}
#'     (\code{"two-model"} mode only; \code{NULL} otherwise)}
#'   \item{summary_adjusted}{the adjusted model's \code{maihda_summary}, or \code{NULL}}
#'   \item{pcv}{the \code{pvc_result} from \code{\link{calculate_pvc}}
#'     (\code{"two-model"} mode only; \code{NULL} otherwise)}
#'   \item{decomposition}{the additive/interaction partition (additive and interaction
#'     variances and shares, per-dimension variances; \code{"crossed-dimensions"} mode
#'     only, \code{NULL} otherwise)}
#'   \item{groups}{a \code{maihda_group_comparison} when \code{group} is supplied,
#'     otherwise \code{NULL}}
#'   \item{mode}{\code{"two-model"} or \code{"crossed-dimensions"}}
#'   \item{context_vars}{the context variable name(s) when \code{context} was
#'     supplied, otherwise \code{NULL}}
#'   \item{formula, adjusted_formula, group_var, call}{bookkeeping for printing}
#'
#' @seealso \code{\link{fit_maihda}} for the single-model fitter,
#'   \code{\link{compare_maihda_groups}} for the group comparison, and
#'   \code{\link{summary.maihda_model}} for the variance summary.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#'
#' # One call: null + adjusted fit, VPC summary, and PCV decomposition. Writing the
#' # dimensions' additive main effects (Gender + Race) gives the fully-specified
#' # adjusted model; maihda() derives the null by dropping them.
#' a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race), data = maihda_health_data)
#' a                              # VPC (null) and PCV (null -> adjusted)
#' a$pcv                          # proportional change in between-stratum variance
#' a$formula                      # null:     BMI ~ Age + (1 | stratum)
#' a$adjusted_formula             # adjusted: null + Gender + Race main effects
#'
#' # Omitting them is equivalent -- maihda() adds them to the adjusted model and
#' # emits a message; the null and PCV are identical to the explicit form above.
#' a0 <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data)
#'
#' plot(a, type = "vpc")          # null model
#' plot(a, type = "effect_decomp")# adjusted model (additive vs intersectional)
#'
#' # Crossed-dimensions decomposition: one model, the dimensions' main effects entered
#' # as RANDOM intercepts. The additive and interaction shares of the between-strata
#' # variance are read directly from the single fit (no null/adjusted pair).
#' cc <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data,
#'              decomposition = "crossed-dimensions")
#' cc                                    # VPC and additive/interaction shares
#' cc$decomposition$additive_share       # crossed-dimensions analogue of the PCV
#' cc$formula                            # BMI ~ Age + (1|Gender) + (1|Race) + (1|stratum)
#'
#' # Add a higher-level grouping variable to also compare across its levels.
#' # maihda_country_data has a real country grouping (PISA achievement data):
#' data(maihda_country_data)
#' a2 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
#'              group = "country")
#' a2
#' plot(a2, type = "group_vpc")
#' plot(a2, type = "group_pcv")
#'
#' # Contextual cross-classified MAIHDA: instead of one model per country (group=),
#' # model the strata CROSSED with country in a single fit. The summary partitions
#' # the unexplained variance into between-stratum vs. between-country vs. residual,
#' # and the PCV is computed with country partialled out.
#' a3 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
#'              context = "country")
#' a3
#' a3$summary$context$vpc_context_total  # the country (general contextual) share
#' plot(a3, type = "context_vpc")
#' }
#'
#' @export
maihda <- function(formula, data, group = NULL, context = NULL, engine = "lme4",
                   family = "gaussian",
                   decomposition = c("two-model", "crossed-dimensions"),
                   autobin = TRUE, shared_strata = TRUE,
                   min_group_n = 30, bootstrap = FALSE,
                   n_boot = 1000, conf_level = 0.95,
                   response_vpc = FALSE, seed = NULL, ...) {
  call <- match.call()
  decomposition <- maihda_resolve_decomposition(decomposition)

  # group= (stratified: one independent model per level) and context= (joint:
  # strata crossed with the higher level in ONE model) are different designs for
  # bringing in a higher level; combining them is ambiguous, so error early.
  if (!is.null(group) && !is.null(context)) {
    stop("Supply either 'group' (stratified comparison: one independent model per ",
         "level) or 'context' (contextual cross-classified model: strata crossed ",
         "with the higher level in one fit), not both.", call. = FALSE)
  }

  # Fit the supplied formula first -- this resolves the stratum dimensions (and the
  # stratum column) and the family. Depending on whether the formula already lists the
  # dimensions' additive main effects, it serves as either the null or the adjusted
  # model below. When the user leaves 'family' at the default we omit it so fit_maihda()
  # can auto-detect a binary outcome, then reuse the resolved family for every model so
  # they agree.
  if (missing(family)) {
    model <- fit_maihda(formula, data, engine = engine, autobin = autobin,
                        context = context, ...)
  } else {
    model <- fit_maihda(formula, data, engine = engine, family = family,
                        autobin = autobin, context = context, ...)
  }
  family_used <- model$family

  # maihda() IS the two-model decomposition -- it has no single-model mode. Building
  # the adjusted model requires at least two identifiable stratum dimensions, so error
  # (rather than return a null-only result) when that is impossible, pointing the user
  # at the shorthand / make_strata() or at fit_maihda() for a single fit. The usual
  # shorthand and make_strata() paths both record the dimensions; only a hand-built
  # 'stratum' column from a custom grouping records none.
  strata_vars <- model$strata_vars
  if (is.null(strata_vars) || length(strata_vars) == 0) {
    stop("maihda() builds an adjusted model from the stratum-defining variables to ",
         "compute the PCV, but they could not be identified (a pre-built 'stratum' ",
         "column from a custom grouping records none). Use the shorthand ",
         "(1 | var1:var2) or run make_strata() so the dimensions are recorded -- or ",
         "call fit_maihda() directly for a single-model fit.", call. = FALSE)
  }
  if (length(strata_vars) < 2) {
    stop("maihda() decomposes the intersectional inequality into additive and ",
         "interaction parts, which needs at least two stratum dimensions, but the ",
         "strata are defined by a single dimension (", paste(strata_vars, collapse = ", "),
         "). With one dimension there is no intersection to decompose; use ",
         "fit_maihda() for a single-dimension random-intercept fit.", call. = FALSE)
  }

  # --- Null and adjusted models ------------------------------------------------
  # The stratum dimensions' additive main effects belong in the ADJUSTED model. They
  # may be written explicitly in the formula (the fully-specified, lme4-native form)
  # or omitted; either way the null model excludes them and the adjusted includes
  # them. maihda_adjusted_terms() gives the correct term per dimension -- an
  # auto-binned numeric dimension enters as its reconstructed tertile factor
  # (.maihda_dim_*), matching the strata, not the raw column.
  dim_terms <- maihda_adjusted_terms(strata_vars, model$strata_autobin_info,
                                     model$original_data)$terms
  supplied_fixed <- attr(stats::terms(reformulas::nobars(model$formula)), "term.labels")
  present_terms <- intersect(dim_terms, supplied_fixed)
  missing_vars <- strata_vars[!(dim_terms %in% supplied_fixed)]

  remove_terms <- function(f, terms) {
    stats::update(f, stats::as.formula(
      paste(". ~ . -", paste(sprintf("`%s`", terms), collapse = " - "))))
  }

  # --- Crossed-dimensions decomposition -----------------------------------------
  # A single model that estimates each dimension's additive main effect as a RANDOM
  # intercept (its variance is that dimension's additive contribution) and the
  # intersection as the interaction RE. This is the crossed-dimensions alternative
  # to the two-model fixed-effects PCV below. Any dimension main effects the user
  # wrote as fixed terms are stripped first -- here they enter as random effects.
  if (decomposition == "crossed-dimensions") {
    base_formula <- if (length(present_terms) > 0) {
      remove_terms(model$formula, present_terms)
    } else {
      model$formula
    }
    cc <- maihda_cross_classified_formula(base_formula, strata_vars,
                                          model$strata_autobin_info,
                                          model$original_data)
    # The formula builder strips ALL random effects before adding the dimension +
    # intersection intercepts; passing `context` again re-appends (and re-tags) any
    # contextual random intercept, so the two structures compose in one fit.
    cc_model <- fit_maihda(cc$formula, cc$data, engine = engine,
                           family = family_used, context = context, ...)
    # Tag the fit so summary()/plot() read the additive-vs-interaction partition.
    cc_model$cc_info <- list(dim_groups = cc$dim_groups,
                             interaction_group = cc$interaction_group,
                             dim_labels = strata_vars)
    summary_obj <- summary(cc_model, bootstrap = bootstrap, n_boot = n_boot,
                           conf_level = conf_level, response_vpc = response_vpc,
                           seed = seed)

    groups <- NULL
    if (!is.null(group)) {
      group_formula <- if (length(present_terms) > 0) {
        remove_terms(formula, present_terms)
      } else {
        formula
      }
      groups <- compare_maihda_groups(
        group_formula, data, group = group, engine = engine, family = family_used,
        shared_strata = shared_strata, min_group_n = min_group_n,
        autobin = autobin, bootstrap = bootstrap, n_boot = n_boot,
        conf_level = conf_level, decomposition = "crossed-dimensions", ...
      )
    }

    return(structure(
      list(
        model = cc_model,
        summary = summary_obj,
        model_adjusted = NULL,
        summary_adjusted = NULL,
        pcv = NULL,
        decomposition = summary_obj$decomposition,
        groups = groups,
        formula = cc_model$formula,
        adjusted_formula = NULL,
        group_var = group,
        mode = "crossed-dimensions",
        context_vars = context,
        call = call
      ),
      class = "maihda_analysis"
    ))
  }

  # In the branches below every derived formula already carries any contextual
  # random intercept (it is part of model$formula and survives remove_terms /
  # maihda_adjusted_formula); passing `context` again is idempotent in
  # fit_maihda() -- it only re-tags the fit so the summaries report the
  # stratum-vs-context partition. The null and adjusted models therefore share the
  # same context structure and the PCV is computed with the context partialled out.
  if (length(present_terms) == 0) {
    # No dimension main effects in the formula: the supplied fit IS the null; build
    # the adjusted by adding the dimensions' additive main effects.
    null_model <- model
    af <- maihda_adjusted_formula(null_model$formula, strata_vars,
                                  null_model$strata_autobin_info,
                                  null_model$original_data)
    adjusted_model <- fit_maihda(af$formula, af$data, engine = engine,
                                 family = family_used, context = context, ...)
    adjusted_formula <- af$formula
  } else if (length(missing_vars) == 0) {
    # Every dimension main effect is already present (the fully-specified adjusted
    # model): the supplied fit IS the adjusted; derive the null by removing them.
    adjusted_model <- model
    adjusted_formula <- model$formula
    null_model <- fit_maihda(remove_terms(model$formula, present_terms),
                             model$original_data, engine = engine,
                             family = family_used, context = context, ...)
  } else {
    # Some dimension main effects present, some missing: fit a clean null (covariates
    # only) and a clean adjusted (all dimension main effects).
    null_model <- fit_maihda(remove_terms(model$formula, present_terms),
                             model$original_data, engine = engine,
                             family = family_used, context = context, ...)
    af <- maihda_adjusted_formula(null_model$formula, strata_vars,
                                  null_model$strata_autobin_info,
                                  null_model$original_data)
    adjusted_model <- fit_maihda(af$formula, af$data, engine = engine,
                                 family = family_used, context = context, ...)
    adjusted_formula <- af$formula
  }

  # When maihda() supplied any dimension main effects itself, say so -- the
  # decomposition stays explicit rather than silent. Listing them in the formula (the
  # lme4-native form) makes the adjusted model self-documenting and suppresses this.
  if (length(missing_vars) > 0) {
    message("maihda(): added the additive main effect(s) of the stratum dimension(s) ",
            paste(missing_vars, collapse = ", "),
            " to the adjusted model; the null model excludes them. List them in the ",
            "formula to specify the adjusted model explicitly.")
  }

  summary_adj <- summary(adjusted_model, bootstrap = bootstrap, n_boot = n_boot,
                         conf_level = conf_level, response_vpc = response_vpc,
                         seed = seed)
  # A successfully fitted pair can still leave the PCV undefined when the null model
  # has zero between-stratum variance (a boundary/singular fit); keep both models and
  # warn rather than aborting in that numerical edge case.
  pcv <- tryCatch(
    calculate_pvc(null_model, adjusted_model, bootstrap = bootstrap, n_boot = n_boot,
                  conf_level = conf_level),
    error = function(e) {
      warning("maihda(): the PCV could not be computed (", conditionMessage(e),
              "). Returning the fitted null and adjusted models without a PCV.",
              call. = FALSE)
      NULL
    })

  summary_obj <- summary(null_model, bootstrap = bootstrap, n_boot = n_boot,
                         conf_level = conf_level, response_vpc = response_vpc,
                         seed = seed)

  groups <- NULL
  if (!is.null(group)) {
    # Strip any dimension main effects the user wrote so the per-group decomposition
    # matches the overall one: each group's null excludes them and its adjusted adds
    # them. (compare_maihda_groups() builds the per-group strata and adjusted itself.)
    group_formula <- if (length(present_terms) > 0) remove_terms(formula, present_terms) else formula
    groups <- compare_maihda_groups(
      group_formula, data, group = group, engine = engine, family = family_used,
      shared_strata = shared_strata, min_group_n = min_group_n,
      autobin = autobin, bootstrap = bootstrap, n_boot = n_boot,
      conf_level = conf_level, ...
    )
  }

  structure(
    list(
      model = null_model,
      summary = summary_obj,
      model_adjusted = adjusted_model,
      summary_adjusted = summary_adj,
      pcv = pcv,
      decomposition = NULL,
      groups = groups,
      formula = null_model$formula,
      adjusted_formula = adjusted_formula,
      group_var = group,
      mode = "two-model",
      context_vars = context,
      call = call
    ),
    class = "maihda_analysis"
  )
}

# Compact discriminatory-accuracy line for the analysis headline: the null /
# strata-only model's AUC + MOR (Merlo's central DA quantity), the adjusted model's
# AUC for comparison, and the response-scale VPC when present. A no-op for a
# non-binomial outcome, where the summaries carry no `discriminatory_accuracy`.
maihda_print_analysis_da <- function(summary_obj, summary_adjusted = NULL) {
  da <- summary_obj$discriminatory_accuracy
  if (is.null(da)) return(invisible(NULL))
  fmt <- function(v, d = 3) {
    if (isTRUE(is.finite(v))) formatC(v, format = "f", digits = d) else "NA"
  }
  cat("\nDiscriminatory accuracy (null model):\n")
  cat(sprintf("  AUC: %s | MOR: %s | cases/controls: %d/%d\n",
              fmt(da$auc), fmt(da$mor), da$n_case, da$n_control))
  da_adj <- if (!is.null(summary_adjusted)) summary_adjusted$discriminatory_accuracy else NULL
  if (!is.null(da_adj) && isTRUE(is.finite(da_adj$auc))) {
    cat(sprintf("  Adjusted-model AUC: %s\n", fmt(da_adj$auc)))
  }
  vr <- summary_obj$vpc_response
  if (!is.null(vr) && isTRUE(is.finite(vr$estimate))) {
    cat(sprintf("  Response-scale VPC (null): %s\n", fmt(vr$estimate, 4)))
  }
  invisible(NULL)
}

#' Print a MAIHDA Analysis
#'
#' @param x A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param ... Additional arguments (not used).
#' @return No return value, called for side effects.
#' @export
print.maihda_analysis <- function(x, ...) {
  cat("MAIHDA Analysis\n")
  cat("===============\n\n")

  # Crossed-dimensions mode: one model, additive (dimension REs) vs interaction
  # (intersection RE) read directly. Distinct layout from the two-model PCV below.
  if (identical(x$mode, "crossed-dimensions")) {
    cat("Decomposition:   crossed-dimensions (single model)\n")
    cat("Formula:         ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
    cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n", sep = "")
    if (!is.null(x$context_vars)) {
      cat("Context:         ", paste(x$context_vars, collapse = ", "),
          " (crossed contextual random intercept)\n", sep = "")
    }
    maihda_print_fit_diagnostics(x$model$diagnostics)

    vpc <- x$summary$vpc
    if (maihda_vpc_has_interval(vpc)) {
      cat(sprintf("VPC/ICC: %.4f [%.4f, %.4f]\n", vpc$estimate, vpc$ci_lower, vpc$ci_upper))
    } else {
      cat(sprintf("VPC/ICC: %.4f\n", vpc$estimate))
    }
    if (!is.null(x$decomposition)) {
      cat("\n")
      maihda_print_cc_decomposition(x$decomposition)
    }
    if (!is.null(x$summary$context)) {
      maihda_print_context_partition(x$summary$context)
    }
    maihda_print_analysis_da(x$summary, x$summary_adjusted)
    if (!is.null(x$summary$stratum_estimates)) {
      cat("Strata: ", nrow(x$summary$stratum_estimates), "\n", sep = "")
    }
    if (!is.null(x$groups)) {
      cat(sprintf("\nGroup comparison by '%s':\n", x$group_var))
      print(x$groups)
    }
    cat("\nUse summary() for variance components and plot(type = ...) for figures.\n")
    return(invisible(x))
  }

  cat("Null formula:    ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  if (!is.null(x$adjusted_formula)) {
    cat("Adjusted formula:", paste(deparse(x$adjusted_formula), collapse = " "), "\n", sep = "")
  }
  cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n", sep = "")
  if (!is.null(x$context_vars)) {
    cat("Context: ", paste(x$context_vars, collapse = ", "),
        " (crossed contextual random intercept in the null and adjusted models)\n",
        sep = "")
  }

  vpc <- x$summary$vpc
  if (maihda_vpc_has_interval(vpc)) {
    cat(sprintf("VPC/ICC (null): %.4f [%.4f, %.4f]\n", vpc$estimate, vpc$ci_lower, vpc$ci_upper))
  } else {
    cat(sprintf("VPC/ICC (null): %.4f\n", vpc$estimate))
  }
  if (!is.null(x$summary$context)) {
    ctx <- x$summary$context
    cat(sprintf("Context share (null): %.4f (between-%s share of unexplained variance)\n",
                ctx$vpc_context_total, paste(ctx$context_vars, collapse = "/")))
  }

  if (!is.null(x$pcv)) {
    pcv <- x$pcv
    if (isTRUE(pcv$bootstrap) && !is.null(pcv$ci_lower)) {
      cat(sprintf("PCV (null -> adjusted): %.4f [%.4f, %.4f]\n",
                  pcv$pvc, pcv$ci_lower, pcv$ci_upper))
    } else {
      cat(sprintf("PCV (null -> adjusted): %.4f\n", pcv$pvc))
    }
    cat(sprintf("Between-stratum variance: %.4f (null) -> %.4f (adjusted)\n",
                pcv$var_model1, pcv$var_model2))
    if (pcv$pvc >= 0) {
      cat(sprintf(paste0("  ~%.1f%% of the between-stratum variance is additive (the ",
                         "dimensions' main\n  effects); the remainder is the between-stratum ",
                         "variance remaining after the\n  additive main effects -- a ",
                         "model-dependent quantity, often interpreted as\n  intersectional ",
                         "interaction, but interpret it cautiously.\n"),
                  pcv$pvc * 100))
    } else {
      cat("  PCV < 0: the additive main effects do not account for the between-stratum\n",
          "  variance (possible suppression/rescaling) -- interpret cautiously.\n", sep = "")
    }
  }

  maihda_print_analysis_da(x$summary, x$summary_adjusted)

  if (!is.null(x$summary$stratum_estimates)) {
    cat("Strata: ", nrow(x$summary$stratum_estimates), "\n", sep = "")
  }

  if (!is.null(x$groups)) {
    cat(sprintf("\nGroup comparison by '%s':\n", x$group_var))
    print(x$groups)
  }

  cat("\nUse summary() for variance components and plot(type = ...) for figures.\n")
  invisible(x)
}

#' Summarize a MAIHDA Analysis
#'
#' Returns the variance summary (VPC/ICC, variance components, stratum estimates)
#' of the fitted model. The per-group comparison, when present, is attached as the
#' \code{"groups"} attribute.
#'
#' @param object A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param ... Additional arguments (not used).
#' @return The \code{maihda_summary} for the fitted model.
#' @export
summary.maihda_analysis <- function(object, ...) {
  out <- object$summary
  attr(out, "groups") <- object$groups
  attr(out, "pcv") <- object$pcv
  attr(out, "adjusted") <- object$summary_adjusted
  out
}

#' Plot a MAIHDA Analysis
#'
#' Dispatches each \code{type} to the model it is valid on. The VPC and shrinkage
#' views (\code{"vpc"}, \code{"obs_vs_shrunken"}, \code{"predicted"}) use the
#' \strong{null} model. The additive-vs-intersectional views (\code{"risk_vs_effect"},
#' \code{"effect_decomp"}, \code{"ternary"}, \code{"prediction_deviation"}) use the
#' \strong{adjusted} model, whose fixed effects carry the dimensions' additive part so
#' the stratum random effect is the pure interaction; with fewer than two dimensions
#' (no adjusted model) they fall back to the null model. Group types
#' (\code{"group_vpc"}, \code{"group_components"}, \code{"group_between_variance"},
#' \code{"group_pcv"}) use the group comparison when \code{\link{maihda}} was called
#' with a \code{group}.
#'
#' @param x A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param type One of the model types ("all", "vpc", "obs_vs_shrunken", "predicted",
#'   "risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation"), the
#'   contextual type ("context_vpc", a stratum-vs-context variance bar; requires
#'   \code{maihda(context = )}), or a group type ("group_vpc", "group_components",
#'   "group_between_variance", "group_pcv"). Default "all".
#' @param ... Additional arguments passed to the underlying plot method.
#' @return A ggplot2 object, or (for \code{type = "all"}) an invisible list of them.
#' @export
plot.maihda_analysis <- function(x, type = "all", ...) {
  type <- match.arg(type, c(
    "all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
    "effect_decomp", "ternary", "prediction_deviation", "context_vpc",
    "group_vpc", "group_components", "group_between_variance", "group_pcv",
    "group_additive_share"
  ))

  if (type == "context_vpc" && is.null(x$context_vars)) {
    stop("No contextual partition is available. Call maihda() with a 'context' ",
         "argument (the contextual cross-classified model).", call. = FALSE)
  }

  group_types <- c("group_vpc", "group_components", "group_between_variance",
                   "group_pcv", "group_additive_share")
  if (type %in% group_types) {
    if (is.null(x$groups)) {
      stop("No group comparison is available. Call maihda() with a 'group' argument.",
           call. = FALSE)
    }
    gtype <- sub("^group_", "", type)
    return(plot(x$groups, type = gtype))
  }

  # The additive-vs-intersectional views are only interpretable on the adjusted
  # model -- its fixed effects (two-model) or dimension random effects
  # (crossed-dimensions) carry the dimensions' additive part, so the stratum random
  # effect is the pure interaction. The VPC, shrinkage, and context views belong to
  # the null / crossed-dimensions model. In crossed-dimensions mode there is no
  # separate adjusted model, so all views use x$model (which already carries the
  # additive structure).
  adjusted_types <- c("risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation")
  adj_model <- if (!is.null(x$model_adjusted)) x$model_adjusted else x$model
  adj_summary <- if (!is.null(x$model_adjusted)) x$summary_adjusted else x$summary

  if (type == "all") {
    null_plots <- list(vpc = plot(x$model, type = "vpc", summary_obj = x$summary, ...))
    null_plots$obs_vs_shrunken <- tryCatch(
      plot(x$model, type = "obs_vs_shrunken", summary_obj = x$summary, ...),
      error = function(e) NULL)
    null_plots$predicted <- tryCatch(
      plot(x$model, type = "predicted", summary_obj = x$summary, ...),
      error = function(e) NULL)
    if (!is.null(x$context_vars)) {
      null_plots$context_vpc <- tryCatch(
        plot(x$model, type = "context_vpc", summary_obj = x$summary, ...),
        error = function(e) NULL)
    }

    adj_plots <- list()
    for (t in adjusted_types) {
      adj_plots[[t]] <- tryCatch(
        plot(adj_model, type = t, summary_obj = adj_summary, ...),
        error = function(e) NULL)
    }

    model_plots <- c(null_plots, adj_plots)
    if (!is.null(x$groups)) {
      group_plots <- list(
        group_vpc = plot(x$groups, type = "vpc"),
        group_components = plot(x$groups, type = "components"),
        group_between_variance = plot(x$groups, type = "between_variance")
      )
      group_plots$group_pcv <- tryCatch(plot(x$groups, type = "pcv"),
                                         error = function(e) NULL)
      model_plots <- c(model_plots, group_plots)
    }
    for (p in model_plots[!vapply(model_plots, is.null, logical(1))]) print(p)
    return(invisible(model_plots))
  }

  if (type %in% adjusted_types) {
    return(plot(adj_model, type = type, summary_obj = adj_summary, ...))
  }

  # vpc, obs_vs_shrunken, predicted -> null model
  plot(x$model, type = type, summary_obj = x$summary, ...)
}
