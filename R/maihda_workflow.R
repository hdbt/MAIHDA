#' Run a Complete MAIHDA Analysis
#'
#' A single high-level entry point that runs the standard MAIHDA workflow and
#' returns one bundled object: it fits the multilevel model, summarises the
#' variance partition (VPC/ICC) and components, and -- when a higher-level
#' grouping variable is supplied -- also compares intersectional inequality across
#' that variable's levels.
#'
#' This is a convenience wrapper around \code{\link{fit_maihda}},
#' \code{\link{summary.maihda_model}} and \code{\link{compare_maihda_groups}}; it
#' always returns the same \code{maihda_analysis} structure (the \code{groups}
#' slot is simply \code{NULL} when \code{group} is not given), so downstream code
#' never has to branch on the return type.
#'
#' @param formula A model formula, using either the intersectional shorthand
#'   \code{outcome ~ covars + (1 | var1:var2)} or \code{... + (1 | stratum)} when
#'   \code{data} already has a \code{stratum} column from \code{\link{make_strata}}.
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
#'   \item{model}{the fitted \code{maihda_model} (see \code{\link{fit_maihda}})}
#'   \item{summary}{the \code{maihda_summary} (VPC/ICC, variance components,
#'     stratum estimates)}
#'   \item{groups}{a \code{maihda_group_comparison} when \code{group} is supplied,
#'     otherwise \code{NULL}}
#'   \item{formula, group_var, call}{bookkeeping for printing}
#'
#' @seealso \code{\link{fit_maihda}} for the single-model fitter,
#'   \code{\link{compare_maihda_groups}} for the group comparison, and
#'   \code{\link{summary.maihda_model}} for the variance summary.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#'
#' # One call: fit + VPC summary
#' a <- maihda(BMI ~ Age + (1 | Gender:Race), data = maihda_health_data)
#' a
#' plot(a, type = "vpc")
#'
#' # Add a higher-level grouping variable to also compare across its levels.
#' # maihda_country_data has a real country grouping (PISA achievement data):
#' data(maihda_country_data)
#' a2 <- maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
#'              group = "country")
#' a2
#' plot(a2, type = "group_vpc")
#' }
#'
#' @export
maihda <- function(formula, data, group = NULL, engine = "lme4",
                   family = "gaussian", autobin = TRUE, shared_strata = TRUE,
                   min_group_n = 30, bootstrap = FALSE,
                   n_boot = 1000, conf_level = 0.95, ...) {
  call <- match.call()

  # Fit the overall model. When the user leaves 'family' at the default we omit
  # it so fit_maihda() can auto-detect a binary outcome; we then reuse whatever
  # family it resolved for the group comparison so every model agrees.
  if (missing(family)) {
    model <- fit_maihda(formula, data, engine = engine, autobin = autobin, ...)
  } else {
    model <- fit_maihda(formula, data, engine = engine, family = family,
                        autobin = autobin, ...)
  }
  family_used <- model$family

  summary_obj <- summary(model, bootstrap = bootstrap, n_boot = n_boot,
                         conf_level = conf_level)

  groups <- NULL
  if (!is.null(group)) {
    groups <- compare_maihda_groups(
      formula, data, group = group, engine = engine, family = family_used,
      shared_strata = shared_strata, min_group_n = min_group_n,
      autobin = autobin, bootstrap = bootstrap, n_boot = n_boot,
      conf_level = conf_level, ...
    )
  }

  structure(
    list(
      model = model,
      summary = summary_obj,
      groups = groups,
      formula = model$formula,
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
  cat("Formula:", paste(deparse(x$formula), collapse = " "), "\n")
  cat("Engine: ", x$model$engine, " | Family: ", x$model$family$family, "\n", sep = "")

  vpc <- x$summary$vpc
  if (maihda_vpc_has_interval(vpc)) {
    cat(sprintf("VPC/ICC: %.4f [%.4f, %.4f]\n", vpc$estimate, vpc$ci_lower, vpc$ci_upper))
  } else {
    cat(sprintf("VPC/ICC: %.4f\n", vpc$estimate))
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
  out
}

#' Plot a MAIHDA Analysis
#'
#' Dispatches to the model's plots (see \code{\link{plot.maihda_model}}) for the
#' model-level \code{type}s, and to the group comparison for \code{"group_vpc"},
#' \code{"group_components"}, and \code{"group_between_variance"} when
#' \code{\link{maihda}} was called with a \code{group}.
#'
#' @param x A \code{maihda_analysis} object from \code{\link{maihda}}.
#' @param type One of the \code{\link{plot.maihda_model}} types ("all", "vpc",
#'   "obs_vs_shrunken", "predicted", "risk_vs_effect", "effect_decomp", "ternary",
#'   "prediction_deviation") or a group type ("group_vpc", "group_components",
#'   "group_between_variance"). Default "all".
#' @param ... Additional arguments passed to the underlying plot method.
#' @return A ggplot2 object, or (for \code{type = "all"}) an invisible list of them.
#' @export
plot.maihda_analysis <- function(x, type = "all", ...) {
  type <- match.arg(type, c(
    "all", "vpc", "obs_vs_shrunken", "predicted", "risk_vs_effect",
    "effect_decomp", "ternary", "prediction_deviation",
    "group_vpc", "group_components", "group_between_variance"
  ))

  if (type %in% c("group_vpc", "group_components", "group_between_variance")) {
    if (is.null(x$groups)) {
      stop("No group comparison is available. Call maihda() with a 'group' argument.",
           call. = FALSE)
    }
    gtype <- sub("^group_", "", type)
    return(plot(x$groups, type = gtype))
  }

  if (type == "all") {
    model_plots <- plot(x$model, type = "all", summary_obj = x$summary, ...)
    if (!is.null(x$groups)) {
      group_plots <- list(
        group_vpc = plot(x$groups, type = "vpc"),
        group_components = plot(x$groups, type = "components"),
        group_between_variance = plot(x$groups, type = "between_variance")
      )
      for (p in group_plots) print(p)
      model_plots <- c(model_plots, group_plots)
    }
    return(invisible(model_plots))
  }

  plot(x$model, type = type, summary_obj = x$summary, ...)
}
