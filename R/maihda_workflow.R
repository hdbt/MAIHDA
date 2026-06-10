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
#'   is run and attached to the result.
#' @param engine Modeling engine, "lme4" (default) or "brms".
#' @param family Model family. Default "gaussian". As in \code{\link{fit_maihda}},
#'   a binary outcome is auto-detected when \code{family} is left at the default,
#'   and the same resolved family is then used for the group comparison so all
#'   models agree.
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
#' @param ... Additional arguments passed to \code{\link{fit_maihda}} (and on to
#'   \code{lmer}/\code{glmer}).
#'
#' @return An object of class \code{maihda_analysis}: a list with
#'   \item{model}{the fitted \strong{null} \code{maihda_model} (see
#'     \code{\link{fit_maihda}}); the VPC/ICC and discriminatory-accuracy source}
#'   \item{summary}{the null model's \code{maihda_summary} (VPC/ICC, variance
#'     components, stratum estimates)}
#'   \item{model_adjusted}{the fitted \strong{adjusted} \code{maihda_model} (null plus
#'     the dimensions' additive main effects), or \code{NULL} with fewer than two
#'     dimensions}
#'   \item{summary_adjusted}{the adjusted model's \code{maihda_summary} (its residual
#'     VPC), or \code{NULL}}
#'   \item{pcv}{the \code{pvc_result} from \code{\link{calculate_pvc}} (null vs
#'     adjusted between-stratum variance), or \code{NULL}}
#'   \item{groups}{a \code{maihda_group_comparison} when \code{group} is supplied,
#'     otherwise \code{NULL}}
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
#' # Add a higher-level grouping variable to also compare across its levels.
#' # maihda_country_data has a real country grouping (PISA achievement data):
#' data(maihda_country_data)
#' a2 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
#'              group = "country")
#' a2
#' plot(a2, type = "group_vpc")
#' plot(a2, type = "group_pcv")
#' }
#'
#' @export
maihda <- function(formula, data, group = NULL, engine = "lme4",
                   family = "gaussian", autobin = TRUE, shared_strata = TRUE,
                   min_group_n = 30, bootstrap = FALSE,
                   n_boot = 1000, conf_level = 0.95, ...) {
  call <- match.call()

  # Fit the supplied formula first -- this resolves the stratum dimensions (and the
  # stratum column) and the family. Depending on whether the formula already lists the
  # dimensions' additive main effects, it serves as either the null or the adjusted
  # model below. When the user leaves 'family' at the default we omit it so fit_maihda()
  # can auto-detect a binary outcome, then reuse the resolved family for every model so
  # they agree.
  if (missing(family)) {
    model <- fit_maihda(formula, data, engine = engine, autobin = autobin, ...)
  } else {
    model <- fit_maihda(formula, data, engine = engine, family = family,
                        autobin = autobin, ...)
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

  if (length(present_terms) == 0) {
    # No dimension main effects in the formula: the supplied fit IS the null; build
    # the adjusted by adding the dimensions' additive main effects.
    null_model <- model
    af <- maihda_adjusted_formula(null_model$formula, strata_vars,
                                  null_model$strata_autobin_info,
                                  null_model$original_data)
    adjusted_model <- fit_maihda(af$formula, af$data, engine = engine,
                                 family = family_used, ...)
    adjusted_formula <- af$formula
  } else if (length(missing_vars) == 0) {
    # Every dimension main effect is already present (the fully-specified adjusted
    # model): the supplied fit IS the adjusted; derive the null by removing them.
    adjusted_model <- model
    adjusted_formula <- model$formula
    null_model <- fit_maihda(remove_terms(model$formula, present_terms),
                             model$original_data, engine = engine,
                             family = family_used, ...)
  } else {
    # Some dimension main effects present, some missing: fit a clean null (covariates
    # only) and a clean adjusted (all dimension main effects).
    null_model <- fit_maihda(remove_terms(model$formula, present_terms),
                             model$original_data, engine = engine,
                             family = family_used, ...)
    af <- maihda_adjusted_formula(null_model$formula, strata_vars,
                                  null_model$strata_autobin_info,
                                  null_model$original_data)
    adjusted_model <- fit_maihda(af$formula, af$data, engine = engine,
                                 family = family_used, ...)
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
                         conf_level = conf_level)
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
                         conf_level = conf_level)

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
      groups = groups,
      formula = null_model$formula,
      adjusted_formula = adjusted_formula,
      group_var = group,
      call = call
    ),
    class = "maihda_analysis"
  )
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
  cat("Null formula:    ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  if (!is.null(x$adjusted_formula)) {
    cat("Adjusted formula:", paste(deparse(x$adjusted_formula), collapse = " "), "\n", sep = "")
  }
  cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n", sep = "")

  vpc <- x$summary$vpc
  if (maihda_vpc_has_interval(vpc)) {
    cat(sprintf("VPC/ICC (null): %.4f [%.4f, %.4f]\n", vpc$estimate, vpc$ci_lower, vpc$ci_upper))
  } else {
    cat(sprintf("VPC/ICC (null): %.4f\n", vpc$estimate))
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
#'   "risk_vs_effect", "effect_decomp", "ternary", "prediction_deviation") or a group
#'   type ("group_vpc", "group_components", "group_between_variance", "group_pcv").
#'   Default "all".
#' @param ... Additional arguments passed to the underlying plot method.
#' @return A ggplot2 object, or (for \code{type = "all"}) an invisible list of them.
#' @export
plot.maihda_analysis <- function(x, type = "all", ...) {
  type <- match.arg(type, c(
    "all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
    "effect_decomp", "ternary", "prediction_deviation",
    "group_vpc", "group_components", "group_between_variance", "group_pcv"
  ))

  group_types <- c("group_vpc", "group_components", "group_between_variance", "group_pcv")
  if (type %in% group_types) {
    if (is.null(x$groups)) {
      stop("No group comparison is available. Call maihda() with a 'group' argument.",
           call. = FALSE)
    }
    gtype <- sub("^group_", "", type)
    return(plot(x$groups, type = gtype))
  }

  # The additive-vs-intersectional views are only interpretable on the adjusted
  # model -- its fixed effects carry the dimensions' additive part, so the stratum
  # random effect is the pure interaction. The VPC and shrinkage views belong to the
  # null model. Fall back to the null model (with its built-in caveat captions) when
  # no adjusted model is available (< 2 dimensions).
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
