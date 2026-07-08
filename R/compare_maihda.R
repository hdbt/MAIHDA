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
#' @param ic Logical; append relative-fit information criteria to the table for
#'   comparing model \emph{structures}: \code{AIC}/\code{BIC} for the likelihood
#'   engines (lme4, ordinal) and \code{WAIC}/\code{LOOIC} for brms (see
#'   \code{\link{maihda_ic}}). Default TRUE. REML \code{lmer} fits are refitted with
#'   ML so AIC/BIC are comparable across different fixed effects. Set FALSE for the
#'   lean VPC-only table.
#'
#' @return A \code{maihda_comparison} data frame of VPC/ICC by model. Interval
#'   columns (\code{ci_lower}/\code{ci_upper}) are included when any model supplies
#'   an interval -- an lme4 bootstrap CI or a brms posterior credible interval. When
#'   \code{ic = TRUE}, information-criteria columns (\code{AIC}/\code{BIC} or
#'   \code{WAIC}/\code{LOOIC}, whichever apply) are appended.
#'
#' @details
#' VPCs are only directly comparable when the models share an outcome,
#' family/link, analytic sample, and strata -- the canonical use is nested models
#' (e.g. null vs covariate-adjusted) on the \emph{same} data and strata, to show
#' how the VPC attenuates. If the supplied models differ in any of these,
#' \code{compare_maihda()} still returns the table but issues a single warning,
#' because the VPCs are then not directly comparable. The same comparability caveat
#' applies to the appended information criteria (see \code{\link{maihda_ic}}). In
#' addition, when the appended criteria mix scales -- likelihood \code{AIC}/\code{BIC}
#' (lme4/ordinal) shown alongside Bayesian \code{WAIC}/\code{LOOIC} (brms), which can
#' happen for a same-family lme4-vs-brms comparison that the family/link check does
#' not flag -- \code{compare_maihda()} warns, because those criteria are on different
#' scales and are not comparable to each other.
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
                          n_boot = 1000, conf_level = 0.95, ic = TRUE) {
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

  # A longitudinal fit has a time-varying VPC; compare_maihda() reports the single
  # reference-time (baseline) value, which is not the whole picture. Warn so the
  # comparison is not read as if the VPC were a scalar.
  if (any(vapply(models, function(m) !is.null(m$longitudinal_info), logical(1)))) {
    warning("One or more models are longitudinal (growth-curve) fits whose VPC is ",
            "time-varying; compare_maihda() compares only the reference-time ",
            "(baseline) VPC. Use summary()/plot(type = \"vpc_trajectory\") for the ",
            "full trajectory.", call. = FALSE)
  }

  # VPCs are only comparable on a shared outcome and family/link. Warn (rather
  # than error, to keep the function flexible) when they differ.
  if (length(models) > 1) {
    responses <- vapply(models, function(m) {
      paste(deparse(m$formula[[2]]), collapse = "")
    }, character(1))
    # Canonical "family(link)" keys; the helper falls back to the wrapper-recorded
    # family for engines where stats::family() is undefined (previously a wemix
    # model compared as "NA(NA)") and normalises engine-specific labels such as
    # lme4's theta-embedding "Negative Binomial(<theta>)".
    fam_keys <- vapply(models, maihda_model_family_key, character(1))
    # Wrapper helpers fall back to the stored analytic $data for engines whose
    # fitted object exposes no model frame (WeMixResults), so the n, row-identity
    # and response checks below apply to those engines too instead of degrading to
    # a silent NA pass (the same checks calculate_pcv() makes).
    nobs_vec <- vapply(models, function(m) {
      n <- maihda_wrapper_nobs(m)
      if (is.finite(n)) as.integer(n) else NA_integer_
    }, integer(1))
    # Row identity of the analytic sample and the row-wise stratum assignment, so
    # disjoint samples of the SAME size and strata are still flagged.
    row_keys <- vapply(models, function(m) {
      rid <- maihda_wrapper_row_ids(m)
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
      maihda_wrapper_response_fingerprint(m)
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
    # Prior weights change the variance estimates, so VPCs from differently
    # weighted fits are not comparable even on the same rows/strata. (Unweighted
    # and explicit unit-weight fits map to the same "unit" key.)
    weight_keys <- vapply(models, function(m) {
      maihda_weight_fingerprint(m$model)
    }, character(1))
    # Likewise for SAMPLING weights (design-weighted fits): the prior-weight
    # fingerprint above cannot see them (it degrades to "unit" for wemix/brms),
    # so differing design weights -- or a weighted vs. unweighted mix -- get
    # their own key, mirroring the calculate_pcv() guard.
    sampling_keys <- vapply(models, function(m) {
      maihda_sampling_weight_fingerprint(m)
    }, character(1))

    issues <- character(0)
    if (length(unique(responses)) > 1) {
      issues <- c(issues, paste0("outcomes (", paste(unique(responses), collapse = ", "), ")"))
    }
    if (length(unique(stats::na.omit(weight_keys))) > 1) {
      issues <- c(issues, "prior weights")
    }
    if (length(unique(stats::na.omit(sampling_keys))) > 1) {
      issues <- c(issues, "sampling weights")
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

  # Append relative-fit information criteria (AIC/BIC for the likelihood engines,
  # WAIC/LOOIC for brms) for comparing model structures. ml = TRUE refits any REML
  # lmer fit with ML so AIC/BIC are comparable across different fixed effects. The
  # whole block is guarded so an IC failure on an exotic fit never breaks the VPC
  # comparison, and only the populated criterion columns are appended.
  if (!is.logical(ic) || length(ic) != 1 || is.na(ic)) {
    stop("'ic' must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(ic)) {
    ic_cols <- tryCatch({
      ic_rows <- lapply(models, function(m) maihda_ic_one(m, ml = TRUE))
      ic_df <- do.call(rbind, ic_rows)
      ic_df[, c("AIC", "BIC", "WAIC", "LOOIC"), drop = FALSE]
    }, error = function(e) NULL)
    if (!is.null(ic_cols)) {
      ic_cols <- ic_cols[, vapply(ic_cols, function(col) !all(is.na(col)),
                                  logical(1)), drop = FALSE]
      if (ncol(ic_cols) > 0) {
        # The likelihood engines (lme4/ordinal) report AIC/BIC and brms reports
        # WAIC/LOOIC; when both kinds of column survive (a same-family lme4 vs
        # brms comparison, which the family/link check above does NOT flag) the
        # appended columns place criteria from different scales side by side.
        # AIC/BIC are not comparable to WAIC/LOOIC (see maihda_ic Details), so
        # warn -- the differences block does not catch this, since the models can
        # agree on outcome, family, sample and strata.
        if (maihda_ic_spans_scales(names(ic_cols))) {
          warning("compare_maihda(): the appended information criteria mix scales ",
                  "-- likelihood AIC/BIC (lme4/ordinal) and Bayesian WAIC/LOOIC ",
                  "(brms) are on different scales and are not comparable to each ",
                  "other. Compare a model to another on the same criterion (AIC with ",
                  "AIC, WAIC/LOOIC within brms), not across the likelihood/Bayesian ",
                  "divide.", call. = FALSE)
        }
        comparison_df <- cbind(comparison_df, ic_cols, stringsAsFactors = FALSE)
      }
    }
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
    theme_maihda() +
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
#' the between-/within-stratum variance components differ across those groups. When
#' the strata are defined by at least two dimensions it also fits the adjusted model
#' (the dimensions' additive main effects) within each group and reports the per-group
#' \code{pcv} -- the proportional change in between-stratum variance, i.e. the additive
#' share of that group's intersectional inequality.
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
#' @param engine Modeling engine, "lme4" (default), "brms", "wemix" (the
#'   design-weighted fit; requires \code{sampling_weights} and is selected
#'   automatically when they are supplied with the default engine), or "ordinal"
#'   (cumulative link mixed model via \code{ordinal::clmm()}; selected
#'   automatically for an ordinal family or an ordered-factor outcome).
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
#' @param decomposition Per-group additive-vs-interaction decomposition: the two-model
#'   null -> adjusted PCV (\code{"two-model"}, default) or the single crossed-dimensions
#'   model (\code{"crossed-dimensions"}; \code{"cross-classified"} is a deprecated
#'   alias that warns). The crossed-dimensions form requires
#'   \code{shared_strata = TRUE} and at least two stratum dimensions, and adds the
#'   \code{var_additive}, \code{var_interaction}, \code{additive_share} and
#'   \code{interaction_share} columns (in place of \code{pcv} /
#'   \code{var_between_adjusted}); \code{var_between} is then the total between-strata
#'   variance (additive + interaction). See \code{\link{maihda}} for the underlying
#'   model and its caveats.
#' @param context Optional character vector naming higher-level \emph{context}
#'   column(s) in \code{data} (e.g. \code{"school"}, \code{"region"}). Each
#'   per-group fit then becomes a contextual cross-classified model --
#'   \code{outcome ~ covars + (1 | stratum) + (1 | context)} -- so within every
#'   group the stratum random intercept is crossed with the context random
#'   intercept(s). \code{vpc}/\code{var_between} are then the between-stratum
#'   quantities \emph{net of} the context, and two columns report the per-group
#'   context partition: \code{var_context} (the between-context variance, summed
#'   over contexts) and \code{vpc_context} (the contexts' share of the group's
#'   unexplained variance). The per-context split is kept on the
#'   \code{"context_per"} attribute. A per-group subset shrinks each context's
#'   level count, so groups with too few context levels (< 10) to identify the
#'   context variance are named in a single warning. Forwarded to
#'   \code{\link{fit_maihda}}; not available for the \code{wemix}/\code{ordinal}
#'   engines (they fit no crossed random effect), and the \code{context} may not
#'   name the \code{group} variable itself.
#' @param sampling_weights Optional name of a sampling-weight column in
#'   \code{data} for design-weighted per-group fits; see \code{\link{fit_maihda}}.
#'   The column is sliced with each group's rows, so every group is fitted with
#'   its own members' weights. Not compatible with \code{engine = "lme4"},
#'   \code{bootstrap = TRUE}, or (under the wemix engine)
#'   \code{decomposition = "crossed-dimensions"}.
#' @param ... Additional arguments passed to \code{\link{fit_maihda}} (and on to
#'   \code{lmer}/\code{glmer}).
#'
#' @return A \code{data.frame} of class \code{maihda_group_comparison} with one
#'   row per group and columns \code{group}, \code{n}, \code{n_strata},
#'   \code{vpc}, \code{var_between}, \code{var_other}, \code{var_residual},
#'   \code{status} (and \code{ci_lower}/\code{ci_upper} when
#'   \code{bootstrap = TRUE}). When the strata are defined by at least two
#'   dimensions, two further columns report the per-group null -> adjusted
#'   decomposition: \code{pcv} (proportional change in between-stratum variance when
#'   the dimensions' additive main effects are added; computed on the
#'   maximum-likelihood scale -- see \code{\link{calculate_pcv}} -- because REML
#'   variances are not comparable across the null vs. adjusted fixed effects),
#'   \code{var_between_adjusted} (a \emph{derived} coherence quantity, reported as
#'   \code{var_between * (1 - pcv)} so it shares the scale of the REML
#'   \code{var_between}/\code{vpc} and the table satisfies
#'   \code{pcv = (var_between - var_between_adjusted) / var_between} exactly -- it is
#'   \strong{not} the adjusted fit's own variance), and
#'   \code{var_between_adjusted_ml} (the adjusted model's \emph{actual}
#'   between-stratum variance, read straight off the adjusted fit on the same
#'   maximum-likelihood scale as the PCV; it differs from \code{var_between_adjusted}
#'   only by the small REML-vs-ML gap in the null variance). All three are
#'   \code{NA} for a group whose adjusted fit failed, and the columns are
#'   omitted entirely when the strata have a single dimension. When \code{context}
#'   is supplied, two further columns report each group's contextual partition:
#'   \code{var_context} (the between-context variance, summed over contexts) and
#'   \code{vpc_context} (the contexts' share of the group's unexplained variance);
#'   the per-context split is on the \code{"context_per"} attribute and the
#'   context name(s) on \code{"context_var"}. These are dropped when no context is
#'   supplied. \code{n} is the analytic sample size used by the
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
#' \code{\link{calculate_pcv}}). A hard fit failure in one group records \code{NA}
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
                                  autobin = TRUE,
                                  decomposition = c("two-model", "crossed-dimensions"),
                                  context = NULL,
                                  sampling_weights = NULL,
                                  warn_linear = TRUE,
                                  ...) {
  decomposition <- maihda_resolve_decomposition(decomposition)
  if (identical(decomposition, "longitudinal")) {
    stop("compare_maihda_groups() does not support decomposition = ",
         "\"longitudinal\" (a per-group longitudinal comparison is out of scope). ",
         "Use maihda(decomposition = \"longitudinal\") on the pooled data.",
         call. = FALSE)
  }
  # ---- input validation ----
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object", call. = FALSE)
  }
  if (decomposition == "crossed-dimensions" && !isTRUE(shared_strata)) {
    stop("decomposition = \"crossed-dimensions\" requires shared_strata = TRUE so the ",
         "stratum dimensions are recorded once on the full data and the crossed-dimensions ",
         "model can be built per group.", call. = FALSE)
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

  # Contextual cross-classified comparison: each per-group fit carries the
  # higher-level context random intercept(s). Validate the names/collisions once
  # on the full data (mirrors fit_maihda()); the columns then ride along in each
  # per-group slice and only the name(s) are forwarded to each fit. The
  # wemix/ordinal engines fit no crossed context random effect, so they reject it
  # up front (below) instead of failing every per-group fit.
  context <- maihda_validate_context(context, data)
  if (!is.null(context) && any(context %in% group)) {
    stop("'context' cannot include the group variable '", group, "': it defines ",
         "the stratified comparison (one independent fit per level), not a ",
         "within-group crossed context. Drop it from 'context'.", call. = FALSE)
  }

  # Sampling weights select the design-weighted engine, mirroring fit_maihda().
  # The weights are a data COLUMN, so the per-group slicing below carries them
  # automatically; only the name is forwarded to each group's fit.
  if (!is.null(sampling_weights)) {
    sampling_weights <- maihda_validate_sampling_weights(sampling_weights, data)
    if (missing(engine)) {
      engine <- "wemix"
      message("compare_maihda_groups(): 'sampling_weights' supplied; using ",
              "engine = \"wemix\" (design-weighted pseudo-maximum-likelihood via ",
              "WeMix). Set 'engine' explicitly to silence this message or to ",
              "choose engine = \"brms\".")
    } else if (identical(engine, "lme4")) {
      stop("Sampling weights are not supported by engine = \"lme4\" (lme4's ",
           "weights are precision weights, not sampling weights). Use ",
           "engine = \"wemix\" or \"brms\".", call. = FALSE)
    }
  }

  # Evaluate the forwarded engine arguments (e.g. weights=/subset=) once, with the
  # data mask, as fit_maihda() does. Capturing them as quosures keeps them working
  # when compare_maihda_groups() is itself called through maihda(). The resulting
  # full-length values are then SLICED per group before both the min_group_n guard
  # and the per-group fit, so an external `weights = w` / `subset = keep` vector
  # (one not stored as a data column) lines up with each group's rows instead of
  # being passed at full length (which fails with a length mismatch or, for a
  # subset, silently recycles onto the wrong rows). Resolved here, BEFORE the
  # engine handshake, so the ordinal pre-check detects the cumulative outcome on
  # the same analytic sample the per-group family resolution uses.
  dot_vals <- lapply(rlang::enquos(...), function(q) rlang::eval_tidy(q, data = data))
  n_full <- nrow(data)
  # A numeric `subset` (e.g. subset = c(1:10, 31:40)) holds GLOBAL row indices into
  # `data`, and a recycled logical mask is likewise positional over the full data.
  # Expand both to a full-length logical mask BEFORE the per-group slicing below;
  # otherwise the same vector is reinterpreted relative to each subgroup (numeric
  # index k picks the k-th row of the group, not global row k), silently fitting
  # the wrong rows. Store the normalized value back into dot_vals so the per-group
  # fit (slice_dots_for_group) slices it too, not just the row-count guard.
  subset_value <- maihda_normalize_subset(dot_vals[["subset"]], n_full)
  dot_vals[["subset"]] <- subset_value
  weights_value <- dot_vals[["weights"]]
  # Family detection and the analytic row-count must see the same sample the engine
  # fits. Sampling weights drop non-finite/<= 0 rows (mapped to NA here so the shared
  # row mask drops them); precision weights only drop NA. Mirrors fit_maihda().
  detect_weights <- if (!is.null(sampling_weights)) {
    maihda_sampling_weight_mask(data[[sampling_weights]])
  } else {
    weights_value
  }

  # Ordinal (cumulative) family <-> engine handshake, mirroring fit_maihda():
  # the per-group fits receive 'engine' explicitly, so fit_maihda()'s own
  # missing(engine) auto-switch could never fire through them. An ordered-factor
  # outcome under all-default family/engine likewise selects the ordinal engine
  # here (the per-group family resolution below then picks up the family). Detect
  # on the analytic sample (subset/weights applied), not the raw outcome column,
  # so the engine choice cannot contradict the analytic family (e.g. an ordered
  # factor subset to two observed levels resolves binomial, not ordinal, and
  # would otherwise fail every group with an engine/family contradiction).
  if (missing(family) && missing(engine) && is.null(sampling_weights) &&
      isTRUE(tryCatch(maihda_response_is_ordinal(formula, data,
                                                 subset = subset_value,
                                                 weights = detect_weights),
                      error = function(e) FALSE))) {
    engine <- "ordinal"
    message("compare_maihda_groups(): the outcome is an ordered factor; using ",
            "the cumulative (ordinal) model with engine = \"ordinal\" ",
            "(ordinal::clmm). Specify 'family'/'engine' explicitly to override.")
  }
  if (maihda_family_is_ordinal(
    if (is.function(family)) tryCatch(family(), error = function(e) NULL) else family
  )) {
    if (missing(engine) && is.null(sampling_weights)) {
      engine <- "ordinal"
      message("compare_maihda_groups(): ordinal (cumulative) family; using ",
              "engine = \"ordinal\" (ordinal::clmm). Set 'engine' explicitly to ",
              "silence this message or to choose engine = \"brms\".")
    } else if (identical(engine, "lme4")) {
      stop("lme4 cannot fit a cumulative (ordinal) model. Use engine = ",
           "\"ordinal\" (ordinal::clmm, the default for this family) or ",
           "engine = \"brms\" (brms::cumulative).", call. = FALSE)
    }
  }

  if (!is.character(engine) || length(engine) != 1 ||
      !engine %in% c("lme4", "brms", "wemix", "ordinal")) {
    stop("'engine' should be one of: lme4, brms, wemix, ordinal", call. = FALSE)
  }
  if (identical(engine, "wemix")) {
    if (decomposition == "crossed-dimensions") {
      stop("decomposition = \"crossed-dimensions\" needs crossed random effects, ",
           "which WeMix does not fit. Use the default two-model decomposition with ",
           "engine = \"wemix\", or engine = \"brms\" for the crossed-dimensions form.",
           call. = FALSE)
    }
    if (!is.null(context)) {
      stop("engine = \"wemix\" does not support 'context' (WeMix fits no crossed ",
           "random effects). Use engine = \"lme4\" or \"brms\" for a contextual ",
           "cross-classified model.", call. = FALSE)
    }
  }
  if (identical(engine, "ordinal")) {
    if (decomposition == "crossed-dimensions") {
      stop("decomposition = \"crossed-dimensions\" needs crossed random effects, ",
           "which the ordinal (clmm) engine does not fit. Use the default ",
           "two-model decomposition, or engine = \"brms\" for the ",
           "crossed-dimensions form.", call. = FALSE)
    }
    if (isTRUE(bootstrap)) {
      stop("Bootstrap intervals are not available for engine = \"ordinal\" ",
           "(ordinal::clmm has no simulate()/refit() machinery). Set ",
           "bootstrap = FALSE, or use engine = \"brms\" for posterior credible ",
           "intervals.", call. = FALSE)
    }
    if (!is.null(context)) {
      stop("engine = \"ordinal\" does not support 'context' (the clmm path fits ",
           "the canonical single (1 | stratum) structure only). Use engine = ",
           "\"brms\" for a contextual cross-classified cumulative model.",
           call. = FALSE)
    }
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

  # `dot_vals`, `n_full`, `subset_value`, `weights_value` and `detect_weights`
  # were resolved above (before the engine handshake) so the ordinal pre-check and
  # the family resolution below both see the same analytic sample. Define the
  # per-group slicers that consume them.
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
                                weights = detect_weights),
      error = function(e) FALSE)
    if (isTRUE(is_binary)) {
      warning("The outcome variable appears to be binary. Using family = 'binomial' ",
              "for every group. To fit linear probability models, specify ",
              "family = 'gaussian' explicitly.", call. = FALSE)
      family <- "binomial"
    } else if (isTRUE(tryCatch(
      maihda_response_is_ordinal(formula, data, subset = subset_value,
                                 weights = detect_weights),
      error = function(e) FALSE))) {
      warning("The outcome variable is an ordered factor. Using the cumulative ",
              "(ordinal) model, family = 'ordinal', for every group. Specify a ",
              "family explicitly to override.", call. = FALSE)
      family <- "ordinal"
    }
  }

  # `family` may be supplied as a (possibly uncalled) family function such as
  # stats::gaussian -- a documented form that fit_maihda accepts. Evaluate it to a
  # family object now so the per-group fits and the family metadata recorded on the
  # result (attr "family", read below as family$family) treat it identically to a
  # family object, rather than erroring on the closure.
  if (is.function(family)) {
    family <- family()
  }

  # ---- build shared strata (or defer to per-group) and the fitting formula ----
  prepared <- maihda_prepare_group_strata(formula, data, shared_strata, autobin)
  data <- prepared$data
  fit_formula <- prepared$formula

  # Whether to run a per-group additive-vs-interaction decomposition. Both forms need
  # at least two stratum dimensions (with one dimension there is no intersection to
  # decompose). do_decomp = the two-model null -> adjusted PCV; do_cc = the single
  # crossed-dimensions model. When neither runs the decomposition columns are dropped.
  decomp_vars <- prepared$strata_vars
  has_two_dims <- !is.null(decomp_vars) && length(decomp_vars) >= 2
  if (decomposition == "crossed-dimensions" && !has_two_dims) {
    stop("decomposition = \"crossed-dimensions\" needs at least two stratum dimensions, ",
         "but the strata are defined by a single dimension (",
         paste(decomp_vars, collapse = ", "), ").", call. = FALSE)
  }
  do_cc <- decomposition == "crossed-dimensions" && has_two_dims
  do_decomp <- decomposition == "two-model" && has_two_dims

  strata_attr_names <- c("strata_info", "strata_vars", "strata_sep", "strata_autobin_info")
  carried_attrs <- stats::setNames(
    lapply(strata_attr_names, function(a) attr(data, a)),
    strata_attr_names
  )

  # When a per-group decomposition will run, the BASE model fitted below must be the
  # NULL (covariates only): the two-model PCV compares it against an adjusted model
  # that ADDS the dimensions' additive main effects, and the crossed-dimensions form
  # enters those dimensions as random effects. If the supplied formula already lists
  # the dimension main effects (e.g. math ~ gender + ses + (1 | gender:ses)), strip
  # them here -- exactly as maihda() does before it delegates -- so the adjusted model
  # does not re-add the same terms (collapsing PCV to 0/NA and pinning the stratum
  # variance to a singular boundary) and the crossed-dimensions model does not enter a
  # dimension as both a fixed and a random effect. Idempotent when maihda() already
  # stripped them, and a no-op for the bare shorthand (1 | a:b) that carries no main
  # effects. An auto-binned numeric dimension is reconstructed as its tertile factor
  # (.maihda_dim_*) for the comparison, which a user formula never names, so it is
  # left in place -- matching maihda()'s quoted-term presence test.
  if (do_decomp || do_cc) {
    dim_terms <- maihda_adjusted_terms(decomp_vars,
                                       carried_attrs[["strata_autobin_info"]],
                                       data)$terms
    # Reject a fixed interaction among the stratum dimensions (e.g. `gender * ses`).
    # Stripping only the main effects below would leave the fixed gender:ses term in
    # place, duplicating the (1 | stratum) random intercept and driving every group's
    # PCV to NA. (maihda() catches this earlier; this guards a direct call.)
    maihda_check_no_dimension_interaction(fit_formula, decomp_vars, dim_terms,
                                          fn = "compare_maihda_groups")
    # Warn once about any numeric stratum dimension entering the per-group adjusted
    # models as a raw linear term. Suppressed (warn_linear = FALSE) when maihda()
    # delegates here, because it already warned on the same overall data.
    if (isTRUE(warn_linear)) {
      maihda_warn_linear_strata_dims(decomp_vars,
                                     carried_attrs[["strata_autobin_info"]], data,
                                     fn = "compare_maihda_groups")
    }
    supplied_fixed <- attr(stats::terms(reformulas::nobars(fit_formula)),
                           "term.labels")
    present_terms <- dim_terms[
      vapply(dim_terms, maihda_quote_name, character(1)) %in% supplied_fixed]
    if (length(present_terms) > 0) {
      # Quote through maihda_quote_name() (the same path used to DETECT the terms
      # against supplied_fixed) rather than a manual sprintf("`%s`", .): a bare
      # backtick wrap mis-parses a name that itself contains a backtick, so the two
      # paths would disagree for that rare-but-legal case.
      quoted_terms <- vapply(present_terms, maihda_quote_name, character(1))
      fit_formula <- stats::update(fit_formula, stats::as.formula(
        paste(". ~ . -", paste(quoted_terms, collapse = " - "))))
    }
  }

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
  # For a contextual comparison: the per-context between-context variances of each
  # fitted group (a named numeric vector), stashed on the "context_per" attribute
  # so the summed table columns keep the per-context split available. And the
  # group x context cells whose context level count is too small to identify the
  # context variance, aggregated into a single weak-identification warning below.
  context_per <- list()
  weak_context_cells <- character(0)

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
      # Rows with a non-finite or <= 0 sampling weight are dropped by the weighted
      # engines, so count the analytic sample the same way: mask those weights to NA
      # (the weight column is part of `sub`, already sliced to this group's rows).
      weights = if (!is.null(sampling_weights)) {
        maihda_sampling_weight_mask(sub[[sampling_weights]])
      } else {
        slice_full(weights_value, idx)
      }
    )
    n_g <- if (is.null(analytic_fr)) nrow(sub) else nrow(analytic_fr)
    row <- data.frame(
      group = g, n = n_g, n_strata = NA_integer_,
      vpc = NA_real_, var_between = NA_real_, var_other = NA_real_,
      var_residual = NA_real_,
      pcv = NA_real_, var_between_adjusted = NA_real_,
      var_between_adjusted_ml = NA_real_,
      var_additive = NA_real_, var_interaction = NA_real_,
      additive_share = NA_real_, interaction_share = NA_real_,
      var_context = NA_real_, vpc_context = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      status = NA_character_, stringsAsFactors = FALSE
    )

    # pre-fit guards (only when a stratum column is already present, i.e. shared).
    # Count populated strata on the ANALYTIC frame (the rows the model actually
    # fits) when available, not the raw subgroup: a group with >= 2 raw strata but
    # only one stratum left after missing-row removal must be skipped as
    # VPC-undefined here, otherwise lme4 fails later with "grouping factors must
    # have > 1 sampled level".
    if ("stratum" %in% names(sub)) {
      stratum_src <- if (!is.null(analytic_fr) && "stratum" %in% names(analytic_fr)) {
        analytic_fr$stratum
      } else {
        sub$stratum
      }
      row$n_strata <- length(unique(stats::na.omit(stratum_src)))
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

    # Fit the per-group model. In crossed-dimensions mode this is the single
    # crossed-dimensions model (dimension REs + intersection RE); otherwise the
    # canonical single-stratum model. Both share the fit/summary/error plumbing.
    fit_obj <- tryCatch(
      {
        # Pass the per-group SLICED dot values (not the raw full-length `...`) so
        # weights/subset/offset align with this group's rows.
        if (do_cc) {
          cc <- maihda_cross_classified_formula(
            fit_formula, decomp_vars, carried_attrs[["strata_autobin_info"]], sub)
          model <- do.call(
            fit_maihda,
            c(list(cc$formula, cc$data, engine = engine, family = family,
                   context = context, sampling_weights = sampling_weights),
              slice_dots_for_group(idx))
          )
          model$cc_info <- list(dim_groups = cc$dim_groups,
                                interaction_group = cc$interaction_group,
                                dim_labels = decomp_vars)
        } else {
          model <- do.call(
            fit_maihda,
            c(list(fit_formula, sub, engine = engine, family = family,
                   context = context, sampling_weights = sampling_weights),
              slice_dots_for_group(idx))
          )
        }
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

    row$vpc <- fit_obj$summ$vpc$estimate
    # Report the analytic sample size actually used by the model (lme4 drops rows
    # with missing outcome/covariates), not the raw group row count.
    analytic_n <- maihda_nobs(fit_obj$model$model)
    if (is.finite(analytic_n)) row$n <- as.integer(analytic_n)
    group_strata <- sort(unique(stats::na.omit(as.character(fit_obj$model$data$stratum))))
    row$n_strata <- length(group_strata)
    populated_strata[[g]] <- group_strata
    row$status <- "ok"

    # Contextual partition (when context = ): summary()$context carries the summed
    # between-context variance and the contexts' total share of the unexplained
    # variance; row$vpc above is already the between-stratum share NET of the
    # context. Read the same way in both the two-model and crossed-dimensions
    # branches (all four summary builders expose these fields). The per-context
    # split is stashed for the "context_per" attribute so the summed table columns
    # do not lose it. A per-group subset shrinks each context's level count, so
    # flag the group x context cells whose analytic sample holds too few context
    # levels (< 10, the documented weak-identification threshold) for one
    # aggregated warning after the loop.
    if (!is.null(context) && !is.null(fit_obj$summ$context)) {
      cs <- fit_obj$summ$context
      row$var_context <- cs$context_var_total
      row$vpc_context <- cs$vpc_context_total
      if (!is.null(cs$per_context)) context_per[[g]] <- cs$per_context
      for (cv in context) {
        n_lvl <- length(unique(stats::na.omit(fit_obj$model$data[[cv]])))
        if (is.finite(n_lvl) && n_lvl < 10) {
          weak_context_cells <- c(weak_context_cells,
                                  sprintf("%s x %s (%d level%s)", g, cv, n_lvl,
                                          if (n_lvl == 1) "" else "s"))
        }
      }
    }

    if (do_cc) {
      # Crossed-dimensions partition read from the single fit: var_between is the total
      # between-strata variance (additive + interaction); the additive share is the
      # crossed-dimensions analogue of the PCV.
      dcmp <- fit_obj$summ$decomposition
      row$var_between <- dcmp$between_var
      row$var_residual <- dcmp$within_var
      row$var_other <- 0
      row$var_additive <- dcmp$additive_var
      row$var_interaction <- dcmp$interaction_var
      row$additive_share <- dcmp$additive_share
      row$interaction_share <- dcmp$interaction_share
    } else {
      vc <- fit_obj$summ$variance_components
      row$var_between <- vc$variance[vc$component == "Between-stratum (random)"][1]
      row$var_residual <- vc$variance[vc$component == "Within-stratum (residual)"][1]
      # Variance of any additional random effects (0 for the canonical single
      # stratum model). Captured so the VPC, the variance columns, and the
      # components plot stay mutually consistent when extra random effects exist.
      other_var <- vc$variance[vc$component == "Other random effects"]
      row$var_other <- if (length(other_var) > 0) sum(other_var) else 0

      # Per-group PCV decomposition: fit the adjusted model (null + the dimensions'
      # additive main effects) on this group's same sample/strata and record the
      # proportional change in between-stratum variance. The adjusted formula uses the
      # global (shared) auto-bin recipe carried on the fitted model, so binned factors
      # match across groups. A failure here (a singular adjusted fit, or zero null
      # between-variance) leaves pcv as NA without aborting the comparison.
      if (do_decomp) {
        af <- maihda_adjusted_formula(fit_obj$model$formula, fit_obj$model$strata_vars,
                                      fit_obj$model$strata_autobin_info,
                                      fit_obj$model$original_data)
        if (!is.null(af)) {
          pcv_obj <- tryCatch({
            adj_model <- do.call(
              fit_maihda,
              c(list(af$formula, af$data, engine = engine, family = family,
                     context = context, sampling_weights = sampling_weights),
                slice_dots_for_group(idx))
            )
            calculate_pcv(fit_obj$model, adj_model)
          }, error = function(e) NULL)
          if (!is.null(pcv_obj)) {
            row$pcv <- pcv_obj$pcv
            # Report the adjusted between-stratum variance on the SAME scale as
            # var_between / vpc (the REML single-model VPC variance), so the table stays
            # internally coherent: PCV = (var_between - var_between_adjusted) /
            # var_between. The PCV itself is calculate_pcv()'s ML-refit value, since REML
            # variances are not comparable across the null vs. adjusted fixed effects;
            # the small REML-vs-ML gap in the null variance is absorbed here rather than
            # left as an apparent inconsistency between the columns. NOTE: this is a
            # derived bookkeeping value, NOT the adjusted fit's own variance.
            row$var_between_adjusted <- row$var_between * (1 - pcv_obj$pcv)
            # The adjusted MODEL's actual between-stratum variance, on the ML scale the
            # PCV is computed on (var_model2 from calculate_pcv()). Unlike
            # var_between_adjusted above, this is read straight off the adjusted fit, so
            # it is the literal "adjusted between-stratum variance" -- it differs from
            # var_between_adjusted only by the REML-vs-ML gap in the null variance.
            row$var_between_adjusted_ml <- pcv_obj$var_model2
          }
        }
      }
    }

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

  # Contextual weak identification: stratifying by group shrinks each group x
  # context cell, so a context that identifies its variance on the pooled data can
  # have too few levels within a group. Name the offending cells in one warning so
  # the per-group context variances there are read with caution (the brms engine's
  # weakly-informative priors regularise this better than lme4).
  if (length(weak_context_cells) > 0) {
    warning("compare_maihda_groups(): few context levels within group(s) (",
            paste(weak_context_cells, collapse = "; "),
            "). A context with < 10 levels weakly identifies its variance (often a ",
            "singular lme4 fit, VPC_context near 0); interpret those groups' context ",
            "partition with caution or use engine = \"brms\".", call. = FALSE)
  }

  out <- do.call(rbind, rows)
  # Drop the interval columns only when no group supplied one (e.g. unbootstrapped
  # lme4); brms groups carry a posterior credible interval without bootstrap.
  if (all(is.na(out$ci_lower)) && all(is.na(out$ci_upper))) {
    out$ci_lower <- NULL
    out$ci_upper <- NULL
  }
  # Keep only the decomposition columns for the mode that ran: the two-model PCV
  # columns (pcv, var_between_adjusted, var_between_adjusted_ml) in "two-model" mode,
  # the crossed-dimensions columns (var_additive, var_interaction, additive_share,
  # interaction_share) in "crossed-dimensions" mode; drop both when no decomposition
  # ran (single dimension).
  if (!do_decomp) {
    out$pcv <- NULL
    out$var_between_adjusted <- NULL
    out$var_between_adjusted_ml <- NULL
  }
  if (!do_cc) {
    out$var_additive <- NULL
    out$var_interaction <- NULL
    out$additive_share <- NULL
    out$interaction_share <- NULL
  }
  # The contextual columns only apply to a contextual comparison (context = );
  # drop them entirely otherwise, mirroring the decomposition columns above.
  if (is.null(context)) {
    out$var_context <- NULL
    out$vpc_context <- NULL
  }
  rownames(out) <- NULL

  attr(out, "group_var") <- group
  attr(out, "engine") <- engine
  attr(out, "family") <- if (is.character(family)) family else family$family
  attr(out, "shared_strata") <- shared_strata
  attr(out, "decomposition") <- decomposition
  # The context variable name(s) and the per-group, per-context between-context
  # variances (a named list keyed by group). NULL when no context was supplied.
  attr(out, "context_var") <- context
  attr(out, "context_per") <- if (is.null(context)) NULL else context_per
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
  cat(maihda_palette()$bold("MAIHDA Group Comparison"), "\n", sep = "")
  cat("=======================\n\n")
  # `[.maihda_group_comparison` re-attaches the metadata attributes after a row/column
  # subset, but degrade gracefully if they are missing anyway (e.g. an object built by
  # some other path): print the bare data frame rather than blank header fields.
  group_var <- attr(x, "group_var")
  if (!is.null(group_var)) {
    cat("Group variable:", group_var, "\n")
    cat("Engine:", attr(x, "engine"),
        " | Family:", attr(x, "family"),
        " | Strata:", if (isTRUE(attr(x, "shared_strata"))) "shared/global" else "per-group",
        "\n")
    # Contextual comparison: each per-group fit is a contextual cross-classified
    # model, so vpc/var_between are the between-stratum share NET of the context and
    # var_context/vpc_context are the (summed) between-context partition.
    context_var <- attr(x, "context_var")
    if (!is.null(context_var)) {
      cat("Context: ", paste(context_var, collapse = ", "),
          " (per-group crossed contextual random intercept; vpc/var_between are ",
          "net of context)\n", sep = "")
    }
    cat("\n")
  }
  print(as.data.frame(x), row.names = FALSE, digits = 4)
  invisible(x)
}

#' Subset a MAIHDA group comparison
#'
#' Indexing method for \code{\link{compare_maihda_groups}} results that preserves the
#' comparison's metadata attributes (group variable, engine, family, ...) so the
#' \code{print()} header stays populated after a row/column subset. Plain
#' \code{[.data.frame} keeps the S3 class but drops those attributes, which would
#' otherwise blank the printed header.
#'
#' @param x A \code{maihda_group_comparison} object.
#' @param ... Row/column indices forwarded to \code{[.data.frame}.
#' @return A \code{maihda_group_comparison} with the metadata attributes carried over,
#'   or a bare vector when the selection drops to a single column.
#' @keywords internal
#' @export
`[.maihda_group_comparison` <- function(x, ...) {
  out <- NextMethod()
  # A single-column selection with drop = TRUE returns a bare vector; leave it
  # unclassed. Otherwise restore the class and the metadata attributes the print
  # method relies on.
  if (is.data.frame(out)) {
    for (a in c("group_var", "engine", "family", "shared_strata", "decomposition",
                "context_var", "context_per")) {
      attr(out, a) <- attr(x, a)
    }
    class(out) <- class(x)
  }
  out
}

#' Join non-empty caption lines for a plot
#'
#' Pastes its arguments into a single newline-separated caption, dropping any that
#' are NULL, NA, or empty. Returns NULL when nothing remains, so a plot with no
#' caption content keeps a clean (caption-free) look.
#'
#' @param ... Character scalars (or NULL).
#' @return A single newline-joined string, or NULL.
#' @keywords internal
maihda_compose_caption <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0) {
    return(NULL)
  }
  paste(parts, collapse = "\n")
}

#' Plot a MAIHDA Group Comparison
#'
#' Visualises the output of \code{\link{compare_maihda_groups}} as a point/forest
#' plot of the VPC/ICC by group, as stacked variance-composition bars (between- vs
#' within-stratum share) by group, as bars of the absolute between-stratum
#' (intersectional) variance by group, or as bars of the additive share (PCV) by
#' group. Dispatched via \code{plot()} on the classed result.
#'
#' @param x A \code{maihda_group_comparison} object from
#'   \code{\link{compare_maihda_groups}}.
#' @param type One of "vpc" (default) for VPC by group with optional bootstrap
#'   confidence intervals, "components" for stacked variance proportions (additive /
#'   interaction / residual for a crossed-dimensions comparison, between / other /
#'   residual otherwise, with a separate context slice for a contextual comparison),
#'   "between_variance" for the absolute between-stratum variance
#'   by group, "pcv" for the two-model additive share (null -> adjusted proportional
#'   change in between-stratum variance) by group, or "additive_share" for the
#'   crossed-dimensions additive share by group. The VPC is a \emph{share} of the
#'   unexplained variance;
#'   "between_variance" shows the \emph{magnitude} the ratio cannot convey (two groups
#'   with very different VPCs can share a between-stratum variance, and vice versa);
#'   "pcv" requires strata defined by at least two dimensions.
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
#' plot(cmp, type = "between_variance")
#' }
#'
#' @export
#' @import ggplot2
#' @importFrom rlang .data
plot.maihda_group_comparison <- function(x, type = c("vpc", "components", "between_variance", "pcv", "additive_share"), ...) {
  if (!inherits(x, "maihda_group_comparison")) {
    stop("'x' must be a maihda_group_comparison object from compare_maihda_groups().",
         call. = FALSE)
  }
  type <- match.arg(type)

  df_all <- as.data.frame(x)
  df <- df_all[!is.na(df_all$vpc), , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No groups with an estimable VPC to plot.", call. = FALSE)
  }
  group_var <- attr(x, "group_var")

  # Groups skipped (small analytic n, a single populated stratum) or whose fit failed
  # carry an NA vpc and are dropped from every view. Note them on the plot rather than
  # letting them disappear silently (long lists are truncated past five names).
  omitted <- as.character(df_all$group[is.na(df_all$vpc)])
  omit_note <- if (length(omitted) > 0) {
    shown <- omitted
    if (length(shown) > 5) {
      shown <- c(shown[seq_len(5)], sprintf("(+%d more)", length(omitted) - 5))
    }
    sprintf("%d group(s) omitted (no estimable VPC): %s.",
            length(omitted), paste(shown, collapse = ", "))
  } else {
    NULL
  }

  if (type == "vpc") {
    df <- df[order(df$vpc), , drop = FALSE]
    df$group <- factor(df$group, levels = df$group)
    has_ci <- all(c("ci_lower", "ci_upper") %in% names(df)) &&
      any(is.finite(df$ci_lower) & is.finite(df$ci_upper))

    vpc_caption <- maihda_compose_caption(omit_note)

    p <- ggplot(df, aes(x = .data$group, y = .data$vpc)) +
      geom_point(size = 4, color = "#0072B2") +
      labs(
        title = "Intersectional VPC/ICC by Group",
        x = group_var,
        y = "VPC (ICC)",
        caption = vpc_caption
      ) +
      theme_maihda() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      coord_cartesian(ylim = c(0, 1))

    if (has_ci) {
      p <- p + geom_errorbar(aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                             width = 0.2, color = "#0072B2")
    }
    return(p)
  }

  if (type == "between_variance") {
    # Absolute between-stratum variance -- the magnitude the VPC ratio cannot convey:
    # two groups with very different VPCs can share a var_between, and vice versa.
    # Singular groups (var_between pinned at the boundary) stay as an informative zero
    # bar; only NA (skipped/failed) groups are dropped, matching the "vpc" view.
    df <- df[order(df$var_between), , drop = FALSE]
    df$group <- factor(df$group, levels = df$group)

    bv_caption <- maihda_compose_caption(omit_note)

    return(
      ggplot(df, aes(x = .data$group, y = .data$var_between)) +
        geom_col(fill = "#E69F00") +
        labs(
          title = "Between-stratum (intersectional) variance by group",
          x = group_var,
          y = "Between-stratum variance",
          caption = bv_caption
        ) +
        theme_maihda() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
  }

  if (type == "pcv") {
    # Additive share by group: the proportional change in between-stratum variance
    # from the null to the adjusted (dimension main effects) model.
    if (!"pcv" %in% names(df)) {
      stop("No PCV available to plot: this comparison has fewer than two stratum ",
           "dimensions, so there is no additive-vs-intersectional decomposition.",
           call. = FALSE)
    }
    df_pcv <- df[is.finite(df$pcv), , drop = FALSE]
    if (nrow(df_pcv) == 0) {
      stop("No groups with an estimable PCV to plot.", call. = FALSE)
    }
    df_pcv <- df_pcv[order(df_pcv$pcv), , drop = FALSE]
    df_pcv$group <- factor(df_pcv$group, levels = df_pcv$group)

    pcv_caption <- maihda_compose_caption(omit_note)

    return(
      ggplot(df_pcv, aes(x = .data$group, y = .data$pcv)) +
        geom_col(fill = "#CC79A7") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        labs(
          title = "Additive share (PCV) by group",
          x = group_var,
          y = "PCV (null -> adjusted)",
          caption = pcv_caption
        ) +
        theme_maihda() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
  }

  if (type == "additive_share") {
    # Crossed-dimensions additive share by group (parallel to the two-model "pcv").
    if (!"additive_share" %in% names(df)) {
      stop("No additive share to plot: this comparison was not run with ",
           "decomposition = \"crossed-dimensions\".", call. = FALSE)
    }
    df_as <- df[is.finite(df$additive_share), , drop = FALSE]
    if (nrow(df_as) == 0) {
      stop("No groups with an estimable additive share to plot.", call. = FALSE)
    }
    df_as <- df_as[order(df_as$additive_share), , drop = FALSE]
    df_as$group <- factor(df_as$group, levels = df_as$group)

    as_caption <- maihda_compose_caption(omit_note)

    return(
      ggplot(df_as, aes(x = .data$group, y = .data$additive_share)) +
        geom_col(fill = "#CC79A7") +
        labs(
          title = "Additive share by group (crossed-dimensions)",
          x = group_var,
          y = "Additive share of between-strata variance",
          caption = as_caption
        ) +
        theme_maihda() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        coord_cartesian(ylim = c(0, 1))
    )
  }

  # type == "components": stacked variance proportions per group. In crossed-dimensions
  # mode the between-strata variance is itself split into the additive (dimension main
  # effects) and intersectional interaction parts; otherwise it is the single
  # between-stratum component plus any "Other random effects" (0 for the canonical
  # single-stratum model). Either way the slices sum to the VPC denominator.
  if ("var_additive" %in% names(df) && "var_interaction" %in% names(df)) {
    var_add <- df$var_additive; var_add[is.na(var_add)] <- 0
    var_int <- df$var_interaction; var_int[is.na(var_int)] <- 0
    # A contextual crossed-dimensions comparison also carries a between-context
    # slice, which must enter the denominator alongside additive + interaction +
    # residual so the shares still sum to 1. 0/absent otherwise.
    var_context <- if ("var_context" %in% names(df)) df$var_context else rep(0, nrow(df))
    var_context[is.na(var_context)] <- 0
    has_context <- any(var_context > sqrt(.Machine$double.eps))
    totals <- var_add + var_int + var_context + df$var_residual
    comp_blocks <- list(
      data.frame(group = df$group, component = "Additive (dimension main effects)",
                 variance = var_add, stringsAsFactors = FALSE),
      data.frame(group = df$group, component = "Intersectional interaction",
                 variance = var_int, stringsAsFactors = FALSE)
    )
    if (has_context) {
      comp_blocks <- c(comp_blocks, list(
        data.frame(group = df$group, component = "Context (higher level)",
                   variance = var_context, stringsAsFactors = FALSE)))
    }
    comp_blocks <- c(comp_blocks, list(
      data.frame(group = df$group, component = "Within-stratum (residual)",
                 variance = df$var_residual, stringsAsFactors = FALSE)))
    comp <- do.call(rbind, comp_blocks)
    total_map <- stats::setNames(totals, as.character(df$group))
    comp$proportion <- comp$variance / total_map[as.character(comp$group)]
    comp$group <- factor(comp$group, levels = df$group[order(df$vpc)])
    comp$component <- factor(
      comp$component,
      levels = c("Additive (dimension main effects)", "Intersectional interaction",
                 "Context (higher level)", "Within-stratum (residual)")
    )
    cc_colors <- c(
      "Additive (dimension main effects)" = "#CC79A7",
      "Intersectional interaction" = "#E69F00",
      "Context (higher level)" = "#009E73",
      "Within-stratum (residual)" = "#56B4E9"
    )
    return(
      ggplot(comp, aes(x = .data$group, y = .data$proportion, fill = .data$component)) +
        geom_bar(stat = "identity", color = "white") +
        scale_fill_manual(values = cc_colors) +
        labs(
          title = "Variance Composition by Group (crossed-dimensions)",
          x = group_var,
          y = "Proportion of Variance",
          fill = "Component",
          caption = maihda_compose_caption(omit_note)
        ) +
        theme_maihda() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
  }

  var_other <- if ("var_other" %in% names(df)) df$var_other else rep(0, nrow(df))
  var_other[is.na(var_other)] <- 0
  # Contextual comparison: the (summed) between-context variance is its own slice
  # (the "general contextual effect"), so it must enter the VPC denominator too --
  # otherwise the between-stratum share would be overstated. 0/absent otherwise.
  var_context <- if ("var_context" %in% names(df)) df$var_context else rep(0, nrow(df))
  var_context[is.na(var_context)] <- 0
  has_context <- any(var_context > sqrt(.Machine$double.eps))
  totals <- df$var_between + var_context + var_other + df$var_residual

  comp_blocks <- list(
    data.frame(group = df$group, component = "Between-stratum (random)",
               variance = df$var_between, stringsAsFactors = FALSE)
  )
  if (has_context) {
    comp_blocks <- c(comp_blocks, list(
      data.frame(group = df$group, component = "Context (higher level)",
                 variance = var_context, stringsAsFactors = FALSE)
    ))
  }
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
    levels = c("Between-stratum (random)", "Context (higher level)",
               "Other random effects", "Within-stratum (residual)")
  )
  component_colors <- c(
    "Between-stratum (random)" = "#E69F00",
    "Context (higher level)" = "#CC79A7",
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
      fill = "Component",
      caption = maihda_compose_caption(omit_note)
    ) +
    theme_maihda() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.caption = element_text(hjust = 0.5, face = "italic", size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Plot a MAIHDA Group Comparison (deprecated)
#'
#' Deprecated. Use \code{plot()} on the \code{\link{compare_maihda_groups}}
#' result instead, e.g. \code{plot(cmp, type = "vpc")}.
#'
#' @param x A \code{maihda_group_comparison} object.
#' @param type One of "vpc" (default), "components", "between_variance", "pcv", or
#'   "additive_share".
#' @return A ggplot2 object.
#' @keywords internal
#' @export
plot_group_comparison <- function(x, type = c("vpc", "components", "between_variance", "pcv", "additive_share")) {
  .Deprecated("plot", msg = paste(
    "'plot_group_comparison()' is deprecated.",
    "Use plot() on the compare_maihda_groups() result, e.g. plot(cmp, type = 'vpc')."
  ))
  plot(x, type = match.arg(type))
}
