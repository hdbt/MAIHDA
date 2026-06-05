#' Compare MAIHDA Models
#'
#' Compares variance partition coefficients (VPC/ICC) across multiple MAIHDA models,
#' with optional bootstrap confidence intervals.
#'
#' @param ... Multiple maihda_model objects to compare.
#' @param model_names Optional character vector of names for the models.
#' @param bootstrap Logical; for \strong{lme4} models, compute parametric-bootstrap
#'   VPC confidence intervals. Default FALSE. It does not apply to \strong{brms}
#'   models, which always return a posterior credible interval (so passing
#'   \code{bootstrap = TRUE} with brms models errors) -- their interval is included
#'   regardless.
#' @param n_boot Number of bootstrap samples if bootstrap = TRUE. Default is 1000.
#' @param conf_level Confidence level for the VPC interval (lme4 bootstrap CI or
#'   brms credible interval). Default is 0.95.
#'
#' @return A \code{maihda_comparison} data frame of VPC/ICC by model. Interval
#'   columns (\code{ci_lower}/\code{ci_upper}) are included when any model supplies
#'   an interval -- an lme4 bootstrap CI or a brms posterior credible interval.
#'
#' @details
#' VPCs are only directly comparable when the models share an outcome,
#' family/link, analytic sample, and strata -- the canonical use is nested models
#' (e.g. null vs covariate-adjusted) on the \emph{same} data and strata, to show
#' how the VPC attenuates. If the supplied models differ in any of these,
#' \code{compare_maihda()} still returns the table but issues a single warning,
#' because the VPCs are then not directly comparable.
#'
#' @examples
#' \donttest{
#' # Canonical use: nested models on the SAME data and strata (null vs adjusted)
#' strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#'
#' null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
#' adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)
#'
#' # Compare without bootstrap
#' comparison <- compare_maihda(null_model, adj_model,
#'                             model_names = c("Null", "Adjusted"))
#'
#' # Compare with bootstrap CI
#' comparison_boot <- compare_maihda(null_model, adj_model,
#'                                  model_names = c("Null", "Adjusted"),
#'                                  bootstrap = TRUE, n_boot = 500)
#' }
#'
#' @export
compare_maihda <- function(..., model_names = NULL, bootstrap = FALSE,
                          n_boot = 1000, conf_level = 0.95) {
  models <- list(...)

  # Validate inputs
  if (length(models) == 0) {
    stop("At least one model must be provided")
  }

  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "maihda_model")) {
      stop("All arguments must be maihda_model objects")
    }
  }

  # VPCs are only comparable on a shared outcome and family/link. Warn (rather
  # than error, to keep the function flexible) when they differ.
  if (length(models) > 1) {
    responses <- vapply(models, function(m) {
      paste(deparse(m$formula[[2]]), collapse = "")
    }, character(1))
    fam_keys <- vapply(models, function(m) {
      fam <- maihda_family(m$model)
      paste0(if (!is.null(fam$family)) fam$family else NA_character_, "(",
             if (!is.null(fam$link)) fam$link else NA_character_, ")")
    }, character(1))
    nobs_vec <- vapply(models, function(m) {
      n <- maihda_nobs(m$model)
      if (is.finite(n)) as.integer(n) else NA_integer_
    }, integer(1))
    # Row identity of the analytic sample and the row-wise stratum assignment, so
    # disjoint samples of the SAME size and strata are still flagged (the same
    # checks calculate_pvc() makes).
    row_keys <- vapply(models, function(m) {
      rid <- maihda_row_ids(m$model)
      if (is.null(rid)) NA_character_ else paste(rid, collapse = "\r")
    }, character(1))
    row_stratum_keys <- vapply(models, function(m) {
      if (!is.null(m$data) && "stratum" %in% names(m$data)) {
        paste(as.character(m$data$stratum), collapse = "\r")
      } else {
        NA_character_
      }
    }, character(1))
    # Content fingerprint of the analytic response, so unrelated datasets that
    # happen to share n, default 1:n row names and stratum ids are still caught.
    response_keys <- vapply(models, function(m) {
      maihda_response_fingerprint(m$model)
    }, character(1))
    # Key the strata by their DEFINITIONS (the grouping variables and the stratum
    # labels), not just the integer IDs: two models can both number their strata
    # 1..k while defining them from different variables (e.g. a:b vs a:c).
    strata_keys <- vapply(models, function(m) {
      vars_key <- if (!is.null(m$strata_vars)) {
        paste(m$strata_vars, collapse = ",")
      } else {
        ""
      }
      info <- m$strata_info
      if (!is.null(info) && "label" %in% names(info)) {
        paste0(vars_key, "::",
               paste(sort(unique(as.character(info$label))), collapse = "|"))
      } else if (!is.null(m$data) && "stratum" %in% names(m$data)) {
        paste0(vars_key, "::",
               paste(sort(unique(as.character(stats::na.omit(m$data$stratum)))), collapse = "|"))
      } else {
        NA_character_
      }
    }, character(1))

    issues <- character(0)
    if (length(unique(responses)) > 1) {
      issues <- c(issues, paste0("outcomes (", paste(unique(responses), collapse = ", "), ")"))
    }
    if (length(unique(fam_keys)) > 1) {
      issues <- c(issues, paste0("families/links (", paste(unique(fam_keys), collapse = ", "), ")"))
    }
    sample_differs <- length(unique(stats::na.omit(nobs_vec))) > 1 ||
      length(unique(stats::na.omit(row_keys))) > 1 ||
      length(unique(stats::na.omit(row_stratum_keys))) > 1 ||
      length(unique(stats::na.omit(response_keys))) > 1
    if (sample_differs) {
      issues <- c(issues, paste0("analytic sample (n = ", paste(nobs_vec, collapse = ", "), ")"))
    }
    if (length(unique(stats::na.omit(strata_keys))) > 1) {
      issues <- c(issues, "stratum definitions")
    }
    if (length(issues) > 0) {
      # Single aggregated warning even when several aspects differ.
      warning("compare_maihda(): models differ in ", paste(issues, collapse = " and "),
              ". VPCs are only directly comparable across models that share an ",
              "outcome, family/link, analytic sample, and strata.", call. = FALSE)
    }
  }

  # Create model names if not provided
  if (is.null(model_names)) {
    model_names <- paste0("Model", seq_along(models))
  } else {
    if (length(model_names) != length(models)) {
      stop("Length of model_names must match number of models")
    }
  }

  # Calculate VPC for each model. Include an interval whenever the summary
  # provides one -- an lme4 bootstrap CI (bootstrap = TRUE) or a brms posterior
  # credible interval (always available) -- rather than keying off the bootstrap
  # flag, which dropped brms intervals.
  comparison_list <- lapply(seq_along(models), function(i) {
    summary_obj <- summary(models[[i]], bootstrap = bootstrap,
                                 n_boot = n_boot, conf_level = conf_level)
    vpc <- summary_obj$vpc
    has_ci <- maihda_vpc_has_interval(vpc)
    data.frame(
      model = model_names[i],
      vpc = vpc$estimate,
      ci_lower = if (has_ci) vpc$ci_lower else NA_real_,
      ci_upper = if (has_ci) vpc$ci_upper else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  # Combine results
  comparison_df <- do.call(rbind, comparison_list)
  # Drop the interval columns only if no model supplied one, preserving the
  # plain two-column output for unbootstrapped lme4 comparisons.
  if (all(is.na(comparison_df$ci_lower)) && all(is.na(comparison_df$ci_upper))) {
    comparison_df$ci_lower <- NULL
    comparison_df$ci_upper <- NULL
  }
  # Class the result so plot() dispatches to plot.maihda_comparison(). It remains
  # a data.frame, so existing column access and printing are unaffected.
  class(comparison_df) <- c("maihda_comparison", "data.frame")
  return(comparison_df)
}

#' Plot a MAIHDA Model Comparison
#'
#' Plots VPC/ICC across the models compared by \code{\link{compare_maihda}}.
#' Dispatched via \code{plot()} on the classed result.
#'
#' @param x A \code{maihda_comparison} object from \code{\link{compare_maihda}}.
#' @param ... Additional arguments (not used).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \donttest{
#' strata <- make_strata(maihda_sim_data, vars = c("gender", "race"))
#'
#' null_model <- fit_maihda(health_outcome ~ 1 + (1 | stratum), data = strata$data)
#' adj_model  <- fit_maihda(health_outcome ~ age + (1 | stratum), data = strata$data)
#'
#' comparison <- compare_maihda(null_model, adj_model, bootstrap = TRUE)
#' plot(comparison)
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom rlang .data
plot.maihda_comparison <- function(x, ...) {
  required_cols <- c("model", "vpc")
  if (!is.data.frame(x) || !all(required_cols %in% names(x))) {
    stop("A maihda_comparison must be a data frame with 'model' and 'vpc' columns.",
         call. = FALSE)
  }

  has_ci <- all(c("ci_lower", "ci_upper") %in% names(x))

  p <- ggplot(x, aes(x = .data$model, y = .data$vpc)) +
    geom_point(size = 4, color = "#0072B2") +
    labs(
      title = "Comparison of Variance Partition Coefficients",
      x = "Model",
      y = "VPC (ICC)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_cartesian(ylim = c(0, 1))

  if (has_ci) {
    p <- p + geom_errorbar(aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                          width = 0.2, color = "#0072B2")
  }

  return(p)
}

#' Plot Model Comparison (deprecated)
#'
#' Deprecated. Use \code{plot()} on the \code{\link{compare_maihda}} result
#' instead, e.g. \code{plot(compare_maihda(...))}.
#'
#' @param comparison_df A data frame from \code{compare_maihda()}.
#' @return A ggplot2 object.
#' @keywords internal
#' @export
plot_comparison <- function(comparison_df) {
  .Deprecated("plot", msg = paste(
    "'plot_comparison()' is deprecated.",
    "Use plot() on the compare_maihda() result, e.g. plot(compare_maihda(...))."
  ))
  if (is.data.frame(comparison_df) && !inherits(comparison_df, "maihda_comparison")) {
    class(comparison_df) <- c("maihda_comparison", class(comparison_df))
  }
  plot(comparison_df)
}

#' Compare MAIHDA Metrics Across Levels of a Grouping Variable
#'
#' Fits a separate random-intercept MAIHDA model (intercept-only \emph{random}
#' effects; any fixed-effect covariates in \code{formula} are still used) within
#' each level of a higher-level grouping variable (for example country, region, or
#' survey wave) and reports how the variance partition coefficient (VPC/ICC) and
#' the between-/within-stratum variance components differ across those groups.
#'
#' It estimates one VPC per group as a stratified analysis: each group is modelled
#' independently. It is \emph{not} a cross-classified model and does not adjust the
#' strata for the grouping variable.
#'
#' The VPC is the \emph{share} of the unexplained variance that lies between strata,
#' not the absolute magnitude of intersectional inequality. Because it is a ratio,
#' a group's VPC can differ from another's because the between-stratum variance
#' differs, because the within-stratum (residual) variance differs, or both -- two
#' groups with the same between-stratum variance can have very different VPCs. To
#' compare the absolute amount of between-stratum (intersectional) variation across
#' groups, read the returned \code{var_between} column alongside the VPC rather than
#' treating a higher VPC as "more inequality".
#'
#' It is \strong{descriptive}: it reports each group's VPC (with an interval when
#' available -- an lme4 bootstrap CI or a brms credible interval) for side-by-side
#' comparison, but does not test whether the VPCs differ between groups. The
#' per-group intervals describe each group's own uncertainty; whether two intervals
#' overlap is \emph{not} a valid test of the difference between their VPCs, which
#' would require modelling that difference directly.
#'
#' @param formula A model formula. Either the shorthand intersectional form
#'   \code{outcome ~ covars + (1 | var1:var2)} (strata are built automatically)
#'   or \code{outcome ~ covars + (1 | stratum)} when \code{data} already contains
#'   a \code{stratum} column from \code{\link{make_strata}}.
#' @param data A data frame containing the variables in \code{formula} and the
#'   grouping variable.
#' @param group Character string naming the grouping variable in \code{data}
#'   (e.g. \code{"country"}). A separate model is fitted for each non-missing
#'   level.
#' @param engine Modeling engine, "lme4" (default) or "brms".
#' @param family Model family. Default "gaussian". As in \code{\link{fit_maihda}},
#'   a binary outcome is auto-detected once on the full data and switched to
#'   "binomial" (with a warning) so every group uses the same family.
#' @param shared_strata Logical. When TRUE (default) intersectional strata are
#'   defined once on the full data so that a stratum denotes the same combination
#'   in every group; this makes the stratum \emph{definitions} comparable across
#'   groups. Note that a group may still not contain every stratum, so two groups'
#'   VPCs can be estimated over different sets of populated strata -- they are then
#'   not strictly directly comparable, and the function warns when this happens.
#'   When FALSE, strata are rebuilt independently within each group (stratum
#'   identities are then not comparable across groups at all).
#' @param min_group_n Minimum size of the \emph{analytic} sample a group must have
#'   -- the rows that survive the model frame (covariate transformations applied,
#'   rows with a missing outcome/covariate dropped) -- to be modelled. Groups with
#'   a smaller usable sample are skipped with a warning, even if they have more raw
#'   rows. Default 30.
#' @param bootstrap Logical; compute per-group parametric-bootstrap VPC
#'   confidence intervals. lme4 engine only. Default FALSE.
#' @param n_boot Number of bootstrap samples when \code{bootstrap = TRUE}.
#'   Default 1000.
#' @param conf_level Confidence level for bootstrap intervals. Default 0.95.
#' @param autobin Logical passed to \code{\link{make_strata}} controlling tertile
#'   binning of numeric grouping variables. Default TRUE.
#' @param ... Additional arguments passed to \code{\link{fit_maihda}} (and on to
#'   \code{lmer}/\code{glmer}).
#'
#' @return A \code{data.frame} of class \code{maihda_group_comparison} with one
#'   row per group and columns \code{group}, \code{n}, \code{n_strata},
#'   \code{vpc}, \code{var_between}, \code{var_other}, \code{var_residual},
#'   \code{status} (and \code{ci_lower}/\code{ci_upper} when
#'   \code{bootstrap = TRUE}). \code{n} is the analytic sample size used by the
#'   model (after dropping rows with a missing outcome/covariate) for both fitted
#'   and skipped groups, falling back to the raw row count only when the model
#'   frame cannot be built. \code{var_other} is the variance of any additional
#'   random effects and is 0 for the canonical single-stratum model. Groups that
#'   were skipped or failed have \code{NA} metrics and an explanatory
#'   \code{status}.
#'
#' @details
#' Robustness: a group whose \emph{analytic} sample (rows surviving the model
#' frame) has fewer than \code{min_group_n} observations is always skipped with a
#' warning. A group with fewer than two populated strata is also skipped
#' (VPC is undefined with a single stratum) when the stratum membership is known
#' before fitting -- that is, when \code{shared_strata = TRUE} or \code{data}
#' already carries a \code{stratum} column. Under \code{shared_strata = FALSE}
#' strata are rebuilt inside each group, so a degenerate single-stratum group is
#' instead reported with a "fit failed" status rather than a pre-fit skip. A
#' singular fit yields a VPC of 0 rather than an error (unlike
#' \code{\link{calculate_pvc}}). A hard fit failure in one group records \code{NA}
#' and a status note without aborting the whole comparison.
#'
#' Fit-quality diagnostics: for the \code{lme4} engine, groups whose model is
#' singular or fails to converge keep a \code{status} of \code{"ok"} (the fit did
#' complete) but are named in a single aggregated warning, because their VPC/ICC
#' may be unreliable -- a singular fit usually pins the between-stratum variance at
#' the boundary, giving a VPC of 0.
#'
#' @examples
#' \donttest{
#' data(maihda_country_data)
#' # How does gender x SES inequality in PISA math scores differ across countries?
#' cmp <- compare_maihda_groups(
#'   math ~ 1 + (1 | gender:ses),
#'   data = maihda_country_data,
#'   group = "country"
#' )
#' print(cmp)
#' plot(cmp, type = "vpc")
#' }
#'
#' @seealso \code{\link{compare_maihda}} for comparing different models on the
#'   same data; \code{\link{plot.maihda_group_comparison}} for visualising the result.
#' @export
#' @importFrom reformulas findbars nobars
#' @importFrom stats update na.omit
compare_maihda_groups <- function(formula, data, group, engine = "lme4",
                                  family = "gaussian", shared_strata = TRUE,
                                  min_group_n = 30, bootstrap = FALSE,
                                  n_boot = 1000, conf_level = 0.95,
                                  autobin = TRUE, ...) {
  # ---- input validation ----
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1 || is.na(group)) {
    stop("'group' must be a single column name.", call. = FALSE)
  }
  if (!group %in% names(data)) {
    stop("Group variable not found in data: ", group, call. = FALSE)
  }
  if (!is.character(engine) || length(engine) != 1 || !engine %in% c("lme4", "brms")) {
    stop("'engine' should be one of: lme4, brms", call. = FALSE)
  }
  if (!is.logical(shared_strata) || length(shared_strata) != 1 || is.na(shared_strata)) {
    stop("'shared_strata' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(min_group_n) || length(min_group_n) != 1 ||
      is.na(min_group_n) || !is.finite(min_group_n) || min_group_n < 1) {
    stop("'min_group_n' must be a single positive number.", call. = FALSE)
  }
  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be TRUE or FALSE.", call. = FALSE)
  }
  if (bootstrap) {
    if (engine != "lme4") {
      stop("Bootstrap VPC confidence intervals are only supported for the lme4 engine.",
           call. = FALSE)
    }
    bootstrap_args <- maihda_validate_bootstrap_args(n_boot, conf_level)
    n_boot <- bootstrap_args$n_boot
    conf_level <- bootstrap_args$conf_level
  }

  # Evaluate the forwarded engine arguments (e.g. weights=/subset=) once, with the
  # data mask, as fit_maihda() does. Capturing them as quosures keeps them working
  # when compare_maihda_groups() is itself called through maihda(). The resulting
  # full-length values are then SLICED per group before both the min_group_n guard
  # and the per-group fit, so an external `weights = w` / `subset = keep` vector
  # (one not stored as a data column) lines up with each group's rows instead of
  # being passed at full length (which fails with a length mismatch or, for a
  # subset, silently recycles onto the wrong rows).
  dot_vals <- lapply(rlang::enquos(...), function(q) rlang::eval_tidy(q, data = data))
  subset_value <- dot_vals[["subset"]]
  weights_value <- dot_vals[["weights"]]
  n_full <- nrow(data)
  # Per-group slice of a single full-length value (used for the row-count guard).
  slice_full <- function(val, idx) {
    if (is.null(val) || length(val) != n_full) return(NULL)
    val[idx]
  }
  # Per-group view of all forwarded dots: full-length atomic vectors (weights,
  # subset, offset, ...) are sliced to the group's rows; scalars and objects
  # (control, REML, nAGQ, ...) are passed through unchanged.
  slice_dots_for_group <- function(idx) {
    lapply(dot_vals, function(v) {
      if (is.atomic(v) && is.null(dim(v)) && length(v) == n_full) v[idx] else v
    })
  }

  # ---- resolve family once on the full data (mirrors fit_maihda) ----
  # Detect on the analytic sample (transformations, subset and weight-NA applied),
  # not the raw outcome column.
  if (missing(family)) {
    is_binary <- tryCatch(
      maihda_response_is_binary(formula, data, subset = subset_value,
                                weights = weights_value),
      error = function(e) FALSE)
    if (isTRUE(is_binary)) {
      warning("The outcome variable appears to be binary. Using family = 'binomial' ",
              "for every group. To fit linear probability models, specify ",
              "family = 'gaussian' explicitly.", call. = FALSE)
      family <- "binomial"
    }
  }

  # ---- build shared strata (or defer to per-group) and the fitting formula ----
  prepared <- maihda_prepare_group_strata(formula, data, shared_strata, autobin)
  data <- prepared$data
  fit_formula <- prepared$formula

  strata_attr_names <- c("strata_info", "strata_vars", "strata_sep", "strata_autobin_info")
  carried_attrs <- stats::setNames(
    lapply(strata_attr_names, function(a) attr(data, a)),
    strata_attr_names
  )

  group_values <- as.character(data[[group]])
  group_levels <- sort(unique(group_values[!is.na(group_values)]))
  if (length(group_levels) == 0) {
    stop("Group variable '", group, "' has no non-missing levels.", call. = FALSE)
  }

  rows <- vector("list", length(group_levels))
  # Collect groups whose lme4 fit was singular or failed to converge, so a single
  # aggregated warning can be raised at the end (these make a group's VPC
  # unreliable). brms fits report NA here and are diagnosed via Rhat elsewhere.
  singular_groups <- character(0)
  nonconverged_groups <- character(0)
  # Set of strata actually populated in each successfully fitted group, used to
  # warn when shared strata still leave groups with different stratum support
  # (their VPCs are then estimated over different level sets).
  populated_strata <- list()

  for (gi in seq_along(group_levels)) {
    g <- group_levels[gi]
    idx <- which(as.character(data[[group]]) == g)
    sub <- data[idx, , drop = FALSE]
    # base `[` drops custom attributes; re-attach so fit_maihda sees the strata
    for (a in strata_attr_names) {
      attr(sub, a) <- carried_attrs[[a]]
    }

    # Size of the analytic sample the model will actually fit: the rows surviving
    # the model frame (covariate transformations applied; rows missing, excluded by
    # `subset`, or carrying a missing weight dropped), not the raw group row count.
    # min_group_n guards this so a group with enough raw rows but a tiny usable
    # sample is still skipped.
    analytic_fr <- maihda_analytic_model_frame(
      fit_formula, sub,
      subset = slice_full(subset_value, idx),
      weights = slice_full(weights_value, idx)
    )
    n_g <- if (is.null(analytic_fr)) nrow(sub) else nrow(analytic_fr)
    row <- data.frame(
      group = g, n = n_g, n_strata = NA_integer_,
      vpc = NA_real_, var_between = NA_real_, var_other = NA_real_,
      var_residual = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      status = NA_character_, stringsAsFactors = FALSE
    )

    # pre-fit guards (only when a stratum column is already present, i.e. shared)
    if ("stratum" %in% names(sub)) {
      row$n_strata <- length(unique(stats::na.omit(sub$stratum)))
    }
    if (n_g < min_group_n) {
      row$status <- sprintf("skipped: analytic n = %d < min_group_n = %g", n_g, min_group_n)
      warning(sprintf("Group '%s' has %d analytic rows (< min_group_n = %g); skipped.",
                      g, n_g, min_group_n), call. = FALSE)
      rows[[gi]] <- row
      next
    }
    if (!is.na(row$n_strata) && row$n_strata < 2) {
      row$status <- sprintf("skipped: %d populated stratum (VPC undefined)", row$n_strata)
      warning(sprintf("Group '%s' has %d populated stratum; VPC is undefined, skipped.",
                      g, row$n_strata), call. = FALSE)
      rows[[gi]] <- row
      next
    }

    fit_obj <- tryCatch(
      {
        # Pass the per-group SLICED dot values (not the raw full-length `...`) so
        # weights/subset/offset align with this group's rows.
        model <- do.call(
          fit_maihda,
          c(list(fit_formula, sub, engine = engine, family = family),
            slice_dots_for_group(idx))
        )
        summ <- summary(model, bootstrap = bootstrap, n_boot = n_boot,
                        conf_level = conf_level)
        list(model = model, summ = summ, error = NULL)
      },
      error = function(e) list(model = NULL, summ = NULL, error = conditionMessage(e))
    )

    if (!is.null(fit_obj$error)) {
      row$status <- paste0("fit failed: ", fit_obj$error)
      warning(sprintf("Group '%s' model fit failed: %s", g, fit_obj$error), call. = FALSE)
      rows[[gi]] <- row
      next
    }

    vc <- fit_obj$summ$variance_components
    row$vpc <- fit_obj$summ$vpc$estimate
    row$var_between <- vc$variance[vc$component == "Between-stratum (random)"][1]
    row$var_residual <- vc$variance[vc$component == "Within-stratum (residual)"][1]
    # Variance of any additional random effects (0 for the canonical single
    # stratum model). Captured so the VPC, the variance columns, and the
    # components plot stay mutually consistent when extra random effects exist.
    other_var <- vc$variance[vc$component == "Other random effects"]
    row$var_other <- if (length(other_var) > 0) sum(other_var) else 0
    # Report the analytic sample size actually used by the model (lme4 drops rows
    # with missing outcome/covariates), not the raw group row count.
    analytic_n <- maihda_nobs(fit_obj$model$model)
    if (is.finite(analytic_n)) row$n <- as.integer(analytic_n)
    group_strata <- sort(unique(stats::na.omit(as.character(fit_obj$model$data$stratum))))
    row$n_strata <- length(group_strata)
    populated_strata[[g]] <- group_strata
    row$status <- "ok"

    # Flag fit-quality problems (singular / non-converged) for an aggregated
    # warning. The fit still "succeeded" (status stays "ok"), but the VPC may be
    # unreliable -- a singular fit typically pins the stratum variance at 0.
    diag <- fit_obj$model$diagnostics
    if (isTRUE(diag$singular)) singular_groups <- c(singular_groups, g)
    if (isFALSE(diag$converged)) nonconverged_groups <- c(nonconverged_groups, g)

    # Keep whatever interval the summary provides -- an lme4 bootstrap CI or a
    # brms posterior credible interval -- rather than only the bootstrap flag,
    # which dropped brms intervals.
    if (maihda_vpc_has_interval(fit_obj$summ$vpc)) {
      row$ci_lower <- fit_obj$summ$vpc$ci_lower
      row$ci_upper <- fit_obj$summ$vpc$ci_upper
    }

    rows[[gi]] <- row
  }

  # Even with shared/global strata, groups can end up with different *populated*
  # strata (some combinations are simply absent in a given group). Each group's VPC
  # is then estimated over a different set of strata, so the VPCs are not strictly
  # directly comparable; warn (this does not apply to shared_strata = FALSE, where
  # the strata are explicitly per-group and already documented as non-comparable).
  if (isTRUE(shared_strata) && length(populated_strata) >= 2) {
    distinct_supports <- unique(populated_strata)
    if (length(distinct_supports) > 1) {
      support_sizes <- vapply(populated_strata, length, integer(1))
      warning("compare_maihda_groups(): groups have different populated strata (",
              paste(sprintf("%s: %d", names(populated_strata), support_sizes),
                    collapse = ", "),
              "). Even with shared_strata = TRUE each group's VPC is estimated over ",
              "the strata present in that group, so VPCs based on different stratum ",
              "support are not strictly directly comparable.", call. = FALSE)
    }
  }

  # Single aggregated warning when any group's lme4 fit was singular or did not
  # converge, naming the affected groups so their VPCs can be read with caution.
  diag_notes <- character(0)
  if (length(singular_groups) > 0) {
    diag_notes <- c(diag_notes,
                    paste0("singular fit: ", paste(singular_groups, collapse = ", ")))
  }
  if (length(nonconverged_groups) > 0) {
    diag_notes <- c(diag_notes,
                    paste0("did not converge: ", paste(nonconverged_groups, collapse = ", ")))
  }
  if (length(diag_notes) > 0) {
    warning("compare_maihda_groups(): lme4 reported fit problems (",
            paste(diag_notes, collapse = "; "),
            "). A singular fit means the between-stratum variance is at the boundary ",
            "(often VPC = 0); a non-converged fit may be unreliable. Interpret these ",
            "groups' VPC/ICC values with caution.", call. = FALSE)
  }

  out <- do.call(rbind, rows)
  # Drop the interval columns only when no group supplied one (e.g. unbootstrapped
  # lme4); brms groups carry a posterior credible interval without bootstrap.
  if (all(is.na(out$ci_lower)) && all(is.na(out$ci_upper))) {
    out$ci_lower <- NULL
    out$ci_upper <- NULL
  }
  rownames(out) <- NULL

  attr(out, "group_var") <- group
  attr(out, "engine") <- engine
  attr(out, "family") <- if (is.character(family)) family else family$family
  attr(out, "shared_strata") <- shared_strata
  class(out) <- c("maihda_group_comparison", "data.frame")
  out
}

#' Resolve shared strata and the fitting formula for group comparison
#'
#' @param formula User formula.
#' @param data Full data frame.
#' @param shared_strata Logical; build strata once on the full data.
#' @param autobin Logical passed to make_strata.
#' @return list(data, formula, strata_vars).
#' @keywords internal
#' @importFrom reformulas findbars nobars
maihda_prepare_group_strata <- function(formula, data, shared_strata, autobin = TRUE) {
  re_terms <- reformulas::findbars(formula)
  if (length(re_terms) == 0) {
    stop("'formula' must include a random effect such as (1 | gender:race) or (1 | stratum).",
         call. = FALSE)
  }

  grouping_vars_by_term <- lapply(re_terms, function(x) all.vars(x[[3]]))
  has_stratum_group <- any(vapply(grouping_vars_by_term, function(v) {
    identical(v, "stratum")
  }, logical(1)))

  # Case 1: the formula explicitly references (1 | stratum) -> reuse the existing
  # column. The decision is driven by the formula, not by the incidental presence
  # of a 'stratum' column, so a shorthand formula is always handled in Case 2.
  if (has_stratum_group) {
    if (!"stratum" %in% names(data)) {
      stop("Formula references (1 | stratum) but 'data' has no 'stratum' column. ",
           "Run make_strata() first, or use the shorthand (1 | var1:var2).",
           call. = FALSE)
    }
    strata_vars <- attr(data, "strata_vars")
    if (is.null(strata_vars)) {
      strata_vars <- maihda_infer_strata_vars(attr(data, "strata_info"))
    }
    return(list(data = data, formula = formula, strata_vars = strata_vars))
  }

  # Case 2: shorthand intersectional term -> need a single intercept-only RE.
  if (length(re_terms) != 1) {
    stop("Group comparison with automatic strata supports a single intercept-only ",
         "random effect such as (1 | gender:race). Call make_strata() first for ",
         "more complex structures.", call. = FALSE)
  }
  random_lhs <- paste(deparse(re_terms[[1]][[2]]), collapse = " ")
  if (random_lhs != "1") {
    stop("Group comparison supports intercept-only random effects such as ",
         "(1 | gender:race).", call. = FALSE)
  }
  if (!maihda_is_colon_interaction(re_terms[[1]][[3]])) {
    stop("Automatic strata creation supports a single variable or a colon ",
         "interaction such as (1 | gender:race). For other grouping expressions ",
         "(e.g. interaction(), paste(), cut()), call make_strata() first and use ",
         "(1 | stratum).", call. = FALSE)
  }

  strata_vars <- grouping_vars_by_term[[1]]
  missing_vars <- setdiff(strata_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Grouping variables not found in data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  if ("stratum" %in% names(data)) {
    stop("'data' already has a 'stratum' column but the formula uses the shorthand ",
         "(1 | ", paste(strata_vars, collapse = ":"), "). Use (1 | stratum) to reuse ",
         "the existing column, or remove it to rebuild strata from these variables.",
         call. = FALSE)
  }

  fixed_formula <- reformulas::nobars(formula)
  fit_formula <- stats::update(fixed_formula, . ~ . + (1 | stratum))

  if (shared_strata) {
    strata_result <- make_strata(data, vars = strata_vars, autobin = autobin)
    return(list(data = strata_result$data, formula = fit_formula,
                strata_vars = strata_vars))
  }

  # Per-group strata: keep the shorthand so fit_maihda rebuilds strata in each
  # subset. Return data unchanged (no stratum column yet).
  list(data = data, formula = formula, strata_vars = strata_vars)
}

#' Print method for MAIHDA group comparisons
#'
#' @param x A maihda_group_comparison object.
#' @param ... Additional arguments (not used).
#' @return No return value, called for side effects.
#' @export
print.maihda_group_comparison <- function(x, ...) {
  cat("MAIHDA Group Comparison\n")
  cat("=======================\n\n")
  cat("Group variable:", attr(x, "group_var"), "\n")
  cat("Engine:", attr(x, "engine"),
      " | Family:", attr(x, "family"),
      " | Strata:", if (isTRUE(attr(x, "shared_strata"))) "shared/global" else "per-group",
      "\n\n")
  print(as.data.frame(x), row.names = FALSE, digits = 4)
  invisible(x)
}

#' Plot a MAIHDA Group Comparison
#'
#' Visualises the output of \code{\link{compare_maihda_groups}} either as a
#' point/forest plot of the VPC/ICC by group, or as stacked variance-composition
#' bars (between- vs within-stratum share) by group. Dispatched via \code{plot()}
#' on the classed result.
#'
#' @param x A \code{maihda_group_comparison} object from
#'   \code{\link{compare_maihda_groups}}.
#' @param type Either "vpc" (default) for VPC by group with optional bootstrap
#'   confidence intervals, or "components" for stacked variance proportions.
#' @param ... Additional arguments (not used).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#' cmp <- compare_maihda_groups(BMI ~ Age + (1 | Gender:Race),
#'                              data = maihda_health_data, group = "Education")
#' plot(cmp, type = "vpc")
#' plot(cmp, type = "components")
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom rlang .data
plot.maihda_group_comparison <- function(x, type = c("vpc", "components"), ...) {
  if (!inherits(x, "maihda_group_comparison")) {
    stop("'x' must be a maihda_group_comparison object from compare_maihda_groups().",
         call. = FALSE)
  }
  type <- match.arg(type)

  df <- as.data.frame(x)
  df <- df[!is.na(df$vpc), , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No groups with an estimable VPC to plot.", call. = FALSE)
  }
  group_var <- attr(x, "group_var")

  if (type == "vpc") {
    df <- df[order(df$vpc), , drop = FALSE]
    df$group <- factor(df$group, levels = df$group)
    has_ci <- all(c("ci_lower", "ci_upper") %in% names(df)) &&
      any(is.finite(df$ci_lower) & is.finite(df$ci_upper))

    p <- ggplot(df, aes(x = .data$group, y = .data$vpc)) +
      geom_point(size = 4, color = "#0072B2") +
      labs(
        title = "Intersectional VPC/ICC by Group",
        x = group_var,
        y = "VPC (ICC)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      coord_cartesian(ylim = c(0, 1))

    if (has_ci) {
      p <- p + geom_errorbar(aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                             width = 0.2, color = "#0072B2")
    }
    return(p)
  }

  # type == "components": stacked variance proportions per group. Include any
  # "Other random effects" variance so the slices match the VPC denominator
  # (between + other + residual). var_other is 0 for the canonical single
  # stratum model and the slice is omitted when no group has additional REs.
  var_other <- if ("var_other" %in% names(df)) df$var_other else rep(0, nrow(df))
  var_other[is.na(var_other)] <- 0
  totals <- df$var_between + var_other + df$var_residual

  comp_blocks <- list(
    data.frame(group = df$group, component = "Between-stratum (random)",
               variance = df$var_between, stringsAsFactors = FALSE)
  )
  if (any(var_other > sqrt(.Machine$double.eps))) {
    comp_blocks <- c(comp_blocks, list(
      data.frame(group = df$group, component = "Other random effects",
                 variance = var_other, stringsAsFactors = FALSE)
    ))
  }
  comp_blocks <- c(comp_blocks, list(
    data.frame(group = df$group, component = "Within-stratum (residual)",
               variance = df$var_residual, stringsAsFactors = FALSE)
  ))
  comp <- do.call(rbind, comp_blocks)

  total_map <- stats::setNames(totals, as.character(df$group))
  comp$proportion <- comp$variance / total_map[as.character(comp$group)]
  comp$group <- factor(comp$group, levels = df$group[order(df$vpc)])
  comp$component <- factor(
    comp$component,
    levels = c("Between-stratum (random)", "Other random effects",
               "Within-stratum (residual)")
  )
  component_colors <- c(
    "Between-stratum (random)" = "#E69F00",
    "Other random effects" = "#009E73",
    "Within-stratum (residual)" = "#56B4E9"
  )

  ggplot(comp, aes(x = .data$group, y = .data$proportion, fill = .data$component)) +
    geom_bar(stat = "identity", color = "white") +
    scale_fill_manual(values = component_colors) +
    labs(
      title = "Variance Composition by Group",
      x = group_var,
      y = "Proportion of Variance",
      fill = "Component"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Plot a MAIHDA Group Comparison (deprecated)
#'
#' Deprecated. Use \code{plot()} on the \code{\link{compare_maihda_groups}}
#' result instead, e.g. \code{plot(cmp, type = "vpc")}.
#'
#' @param x A \code{maihda_group_comparison} object.
#' @param type Either "vpc" (default) or "components".
#' @return A ggplot2 object.
#' @keywords internal
#' @export
plot_group_comparison <- function(x, type = c("vpc", "components")) {
  .Deprecated("plot", msg = paste(
    "'plot_group_comparison()' is deprecated.",
    "Use plot() on the compare_maihda_groups() result, e.g. plot(cmp, type = 'vpc')."
  ))
  plot(x, type = match.arg(type))
}
