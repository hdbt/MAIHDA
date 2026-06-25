# Publication-ready MAIHDA results tables.
#
# The canonical MAIHDA write-up (e.g. Evans et al. 2024's tutorial) reports two
# deliverables: (a) a model-comparison table -- the null (Model 1) and adjusted
# (Model 2) fits side by side, with the intercept, between-stratum SD, VPC/ICC,
# PCV, and (for a binary outcome) the AUC; and (b) a ranked table of the strata
# with the highest and lowest predicted outcomes (their Table 4). The package
# already computes every piece -- this file assembles them into one object so a
# single maihda_table() call produces both, ready to export or print.

#' Canonical MAIHDA results table and ranked-strata table
#'
#' @description
#' Assembles the two standard MAIHDA write-up deliverables from a fitted analysis
#' in one call: (a) a \strong{model-results table} contrasting the null and adjusted
#' models (intercept, between-stratum variance and SD, VPC/ICC, the PCV, and -- for a
#' binary outcome -- the AUC and Median Odds Ratio), and (b) a \strong{ranked-strata
#' table} ordering the intersectional strata by their predicted outcome, so the
#' best- and worst-off strata can be read directly. It introduces no new estimator:
#' the model-results table reuses the quantities from \code{summary()} (calling
#' \code{summary()} itself for a bare \code{\link{fit_maihda}} model), and the
#' ranked-strata table reuses the same stratum predictions as
#' \code{plot(type = "predicted")}, so the table agrees exactly with \code{summary()}
#' and \code{plot()}.
#'
#' @details
#' The model-results table is mostly numeric and export-ready (e.g.
#' \code{write.csv(maihda_table(a)$models, ...)} or pass it to \code{knitr::kable()}):
#' statistics are rows, models are columns, and each estimate has accompanying
#' \code{*_lower}/\code{*_upper} columns that hold the confidence/credible interval
#' when one is available (the VPC bootstrap or posterior interval, and the bootstrap
#' PCV interval) and \code{NA} otherwise. The intercept and the variance/SD rows are
#' point estimates. The \code{print()} method renders the same table in the
#' familiar \dQuote{estimate [low, high]} layout.
#'
#' For a \code{"crossed-dimensions"} analysis (one model, no null/adjusted pair) the
#' results table has a single estimate column and gains \dQuote{Additive share} /
#' \dQuote{Interaction share} rows instead of the PCV. For a contextual
#' cross-classified analysis (\code{maihda(context = )}) it gains a
#' \dQuote{Context share (VPC)} row. A bare \code{\link{fit_maihda}} model is also
#' accepted and yields a single-model table (no PCV).
#'
#' The ranked-strata table ranks every stratum by its model-predicted outcome (on the
#' \code{scale} requested), using the same stratum predictions as
#' \code{plot(type = "predicted")}: the predicted value carries the conditional
#' (random-effect) interval, and the stratum random effect (BLUP) is reported
#' alongside it. By default the ranking uses the \strong{null} model -- the headline
#' intersectional inequality (which strata fare best/worst overall); set
#' \code{which = "adjusted"} to rank by the adjusted model instead. The full ranked
#' table is returned in \code{$strata}; \code{print()} shows the top and bottom
#' \code{n_strata}.
#'
#' @param x A \code{maihda_analysis} from \code{\link{maihda}} (the usual input), or a
#'   single \code{maihda_model} from \code{\link{fit_maihda}}.
#' @param n_strata Number of strata to show at each end (top and bottom) in the
#'   printed ranked-strata table. The returned \code{$strata} holds all strata
#'   (\code{NULL} for a longitudinal fit, whose strata are trajectories rather than
#'   single ranked values -- see \code{$strata} in Value). Default 10.
#' @param scale Scale for the predicted stratum values: \code{"response"} (default)
#'   or \code{"link"}. For a cumulative (ordinal) model the response scale is the
#'   expected category score.
#' @param which For a two-model analysis, which model's predictions to rank the
#'   strata by: \code{"null"} (default) or \code{"adjusted"}. Ignored for a
#'   crossed-dimensions analysis or a single model.
#' @param digits Number of decimal places for the \code{print()} method. Default 3.
#' @param ... For a \code{maihda_model} input, additional arguments passed to
#'   \code{\link{summary.maihda_model}} (e.g. \code{bootstrap = TRUE}); ignored for a
#'   \code{maihda_analysis} input, whose summaries are already computed.
#'
#' @return An object of class \code{maihda_table}: a list with
#'   \item{models}{a data frame of the model-results table (statistics in rows; one
#'     estimate column per model, each with \code{*_lower}/\code{*_upper} interval
#'     columns)}
#'   \item{strata}{a data frame of all strata ranked by predicted outcome, with
#'     \code{rank}, \code{stratum}, \code{label}, \code{n}, the predicted value and
#'     its conditional interval, and the stratum random effect and its interval.
#'     \code{NULL} when a single ranked value per stratum is not defined -- in
#'     particular for a longitudinal (growth-curve) fit, whose strata are
#'     trajectories (see \code{strata_note}) -- or if the stratum predictions could
#'     not otherwise be computed}
#'   \item{strata_note}{a character note explaining why \code{strata} is \code{NULL}
#'     (e.g. a longitudinal fit), or \code{NULL} when a ranked-strata table was
#'     produced}
#'   \item{model_keys, model_labels}{the estimate-column keys and their display labels}
#'   \item{family, engine, mode, scale, ranked_by, n_obs, n_strata_total, context_vars}{
#'     metadata used by \code{print()}}
#'
#' @seealso \code{\link{maihda}}, \code{\link{summary.maihda_model}},
#'   \code{\link{calculate_pvc}}, \code{\link{maihda_discriminatory_accuracy}}.
#'
#' @examples
#' \donttest{
#' data(maihda_health_data)
#' a <- maihda(BMI ~ Age + Gender + Race + (1 | Gender:Race), data = maihda_health_data)
#'
#' tab <- maihda_table(a)
#' tab                 # printed: model-results table + top/bottom strata
#' tab$models          # the numeric, export-ready results table
#' tab$strata          # all strata ranked by predicted BMI
#'
#' # write.csv(tab$models, "results.csv", row.names = FALSE)
#' }
#'
#' @export
#' @importFrom utils head tail
maihda_table <- function(x, n_strata = 10L, scale = c("response", "link"),
                         which = c("null", "adjusted"), digits = 3, ...) {
  scale <- match.arg(scale)
  which <- match.arg(which)

  # Resolve the inputs to: the per-model summaries to tabulate, the PCV (if any),
  # the model/summary the strata are ranked from, and a little metadata. Two shapes
  # are accepted -- the usual maihda_analysis bundle, and a bare fitted model.
  if (inherits(x, "maihda_analysis")) {
    mode <- x$mode
    null_summary <- x$summary
    adj_summary <- x$summary_adjusted
    pcv <- x$pcv
    context_vars <- x$context_vars
    if (identical(mode, "two-model") && !is.null(adj_summary)) {
      model_keys <- c("null", "adjusted")
      model_labels <- c(null = "Null (Model 1)", adjusted = "Adjusted (Model 2)")
      model_stats <- list(null = maihda_collect_model_stats(null_summary),
                          adjusted = maihda_collect_model_stats(adj_summary))
    } else {
      # crossed-dimensions (single fit) or any two-model fit lacking its adjusted
      # summary: one estimate column.
      model_keys <- "estimate"
      model_labels <- c(estimate = if (identical(mode, "crossed-dimensions"))
        "Crossed-dimensions" else "Estimate")
      model_stats <- list(estimate = maihda_collect_model_stats(null_summary))
      pcv <- NULL
    }
    primary_model <- x$model
    primary_summary <- x$summary
    if (identical(which, "adjusted") && !is.null(x$model_adjusted)) {
      rank_model <- x$model_adjusted
      rank_summary <- x$summary_adjusted
      ranked_by <- "adjusted"
    } else {
      rank_model <- x$model
      rank_summary <- x$summary
      ranked_by <- if (identical(mode, "two-model")) "null" else "model"
    }
  } else if (inherits(x, "maihda_model")) {
    mode <- "single"
    primary_summary <- summary(x, ...)
    null_summary <- primary_summary
    adj_summary <- NULL
    pcv <- NULL
    context_vars <- x$context_vars
    model_keys <- "estimate"
    model_labels <- c(estimate = "Estimate")
    model_stats <- list(estimate = maihda_collect_model_stats(primary_summary))
    primary_model <- x
    rank_model <- x
    rank_summary <- primary_summary
    ranked_by <- "model"
  } else {
    stop("'x' must be a maihda_analysis (from maihda()) or a maihda_model ",
         "(from fit_maihda()).", call. = FALSE)
  }

  models_df <- maihda_build_results_table(model_stats, pcv)

  # Ranked-strata table -- a bonus deliverable; never let it break the results
  # table (e.g. a brms predict failure on an exotic family). A longitudinal fit
  # has no single ranked value per stratum (each stratum is a trajectory), so the
  # ranking is intentionally omitted with a note pointing to the trajectory tools.
  strata_note <- if (!is.null(rank_model$longitudinal_info)) {
    paste0("Strata are trajectories (growth-curve model): a single ranked ",
           "value per stratum is not defined. Use predict(type = \"strata\") ",
           "for per-stratum intercept/slope, or plot(type = \"trajectories\").")
  } else {
    NULL
  }
  strata_df <- tryCatch(
    maihda_strata_ranking(rank_model, rank_summary, scale = scale),
    error = function(e) NULL)

  fam_name <- tryCatch(maihda_model_family_name(primary_model),
                       error = function(e) NA_character_)
  n_obs <- tryCatch(nrow(primary_model$data), error = function(e) NA_integer_)
  n_strata_total <- if (!is.null(strata_df)) nrow(strata_df) else {
    se <- primary_summary$stratum_estimates
    if (!is.null(se)) nrow(se) else NA_integer_
  }

  structure(
    list(
      models = models_df,
      strata = strata_df,
      model_keys = model_keys,
      model_labels = model_labels,
      family = fam_name,
      engine = primary_model$engine,
      mode = mode,
      scale = scale,
      ranked_by = ranked_by,
      strata_note = strata_note,
      n_obs = n_obs,
      n_strata_total = n_strata_total,
      context_vars = context_vars,
      n_strata = as.integer(n_strata),
      digits = digits,
      call = match.call()
    ),
    class = "maihda_table"
  )
}

# Canonical order of the statistics in the model-results table. A statistic only
# becomes a row if at least one model supplies it (e.g. AUC only for a binary
# outcome, the share rows only for a crossed-dimensions fit).
maihda_table_stat_order <- c(
  "Intercept",
  "Between-stratum variance",
  "Between-stratum SD",
  "VPC/ICC",
  "Context share (VPC)",
  "Additive share",
  "Interaction share",
  "AUC",
  "MOR"
)

# Pull the statistics for one model out of its maihda_summary, as a named list of
# length-3 numeric vectors c(estimate, lower, upper) (lower/upper NA when the
# package does not supply an interval for that quantity).
maihda_collect_model_stats <- function(s) {
  triple <- function(est, lo = NA_real_, hi = NA_real_) {
    c(as.numeric(est)[1], as.numeric(lo)[1], as.numeric(hi)[1])
  }
  st <- list()
  st[["Intercept"]] <- maihda_extract_intercept(s)

  bv <- maihda_between_stratum_variance(s)
  st[["Between-stratum variance"]] <- triple(bv)
  st[["Between-stratum SD"]] <- triple(if (is.finite(bv) && bv >= 0) sqrt(bv) else NA_real_)

  v <- s$vpc
  if (!is.null(v)) {
    st[["VPC/ICC"]] <- triple(v$estimate,
                              if (!is.null(v$ci_lower)) v$ci_lower else NA_real_,
                              if (!is.null(v$ci_upper)) v$ci_upper else NA_real_)
  }

  if (!is.null(s$context)) {
    cx <- s$context
    ci <- cx$vpc_context_total_ci
    st[["Context share (VPC)"]] <- triple(
      cx$vpc_context_total,
      if (!is.null(ci) && length(ci) == 2) ci[1] else NA_real_,
      if (!is.null(ci) && length(ci) == 2) ci[2] else NA_real_)
  }

  if (!is.null(s$decomposition)) {
    d <- s$decomposition
    aci <- d$additive_share_ci
    ici <- d$interaction_share_ci
    st[["Additive share"]] <- triple(
      d$additive_share,
      if (!is.null(aci) && length(aci) == 2) aci[1] else NA_real_,
      if (!is.null(aci) && length(aci) == 2) aci[2] else NA_real_)
    st[["Interaction share"]] <- triple(
      d$interaction_share,
      if (!is.null(ici) && length(ici) == 2) ici[1] else NA_real_,
      if (!is.null(ici) && length(ici) == 2) ici[2] else NA_real_)
  }

  da <- s$discriminatory_accuracy
  if (!is.null(da)) {
    st[["AUC"]] <- triple(da$auc)
    st[["MOR"]] <- triple(da$mor)
  }
  st
}

# The intercept (grand mean, beta0) from a summary's fixed_effects, as a point
# estimate c(est, NA, NA). Handles the lme4/wemix/ordinal data-frame form
# (term/estimate) and the brms matrix form (rownames + an "Estimate" column). A
# cumulative (ordinal) model has thresholds rather than a single intercept, so it
# returns NA -- the thresholds are reported by summary()$thresholds.
maihda_extract_intercept <- function(s) {
  na3 <- c(NA_real_, NA_real_, NA_real_)
  fe <- s$fixed_effects
  if (is.null(fe)) return(na3)

  if (is.matrix(fe)) {
    rn <- rownames(fe)
    i <- which(rn %in% c("Intercept", "(Intercept)"))
    if (length(i) == 0) return(na3)
    est_col <- if ("Estimate" %in% colnames(fe)) "Estimate" else 1L
    return(c(as.numeric(fe[i[1], est_col]), NA_real_, NA_real_))
  }

  if (!all(c("term", "estimate") %in% names(fe))) return(na3)
  i <- which(fe$term %in% c("(Intercept)", "Intercept"))
  if (length(i) == 0) return(na3)
  c(as.numeric(fe$estimate[i[1]]), NA_real_, NA_real_)
}

# The between-stratum (intersectional) variance from a summary, regardless of fit
# type: the total between-strata variance for a crossed-dimensions fit
# (additive + interaction), the stratum component for a contextual fit, the
# baseline (ref_time) between-stratum variance for a longitudinal fit (it is
# time-varying, so there is no single scalar -- report it at the reference time
# the VPC is anchored to), otherwise the "Between-stratum (random)" row of the
# variance-components table.
maihda_between_stratum_variance <- function(s) {
  if (!is.null(s$decomposition) && !is.null(s$decomposition$between_var)) {
    return(as.numeric(s$decomposition$between_var))
  }
  if (!is.null(s$context) && !is.null(s$context$var_stratum)) {
    return(as.numeric(s$context$var_stratum))
  }
  if (!is.null(s$longitudinal)) {
    lng <- s$longitudinal
    if (!is.null(lng$Sigma_stratum) && !is.null(lng$ref_time)) {
      return(as.numeric(maihda_var_at_time(lng$Sigma_stratum, lng$ref_time)))
    }
    return(NA_real_)
  }
  vc <- s$variance_components
  if (!is.null(vc) && "component" %in% names(vc)) {
    row <- which(vc$component == "Between-stratum (random)")
    if (length(row) >= 1) return(as.numeric(vc$variance[row[1]]))
  }
  NA_real_
}

# Assemble the wide model-results data frame from per-model statistics. Columns:
# statistic, then for each model key k the triple k / k_lower / k_upper. The PCV
# (a null -> adjusted quantity) is added as its own row under the adjusted column
# when a two-model PCV is supplied.
maihda_build_results_table <- function(model_stats, pcv = NULL) {
  keys <- names(model_stats)
  present <- character(0)
  for (k in keys) present <- union(present, names(model_stats[[k]]))
  stats_present <- maihda_table_stat_order[maihda_table_stat_order %in% present]

  empty <- c(NA_real_, NA_real_, NA_real_)
  rows <- lapply(stats_present, function(stat) {
    row <- list(statistic = stat)
    for (k in keys) {
      tri <- model_stats[[k]][[stat]]
      if (is.null(tri)) tri <- empty
      row[[k]] <- tri[1]
      row[[paste0(k, "_lower")]] <- tri[2]
      row[[paste0(k, "_upper")]] <- tri[3]
    }
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)

  # PCV row: only meaningful for a two-model fit (null + adjusted). It sits under
  # the adjusted column (it is the null -> adjusted change); the null column is NA.
  if (!is.null(pcv) && all(c("null", "adjusted") %in% keys)) {
    has_ci <- isTRUE(pcv$bootstrap) && !is.null(pcv$ci_lower) && !is.null(pcv$ci_upper)
    pcv_row <- data.frame(statistic = "PCV (null -> adjusted)",
                          null = NA_real_, null_lower = NA_real_, null_upper = NA_real_,
                          adjusted = as.numeric(pcv$pvc),
                          adjusted_lower = if (has_ci) as.numeric(pcv$ci_lower) else NA_real_,
                          adjusted_upper = if (has_ci) as.numeric(pcv$ci_upper) else NA_real_,
                          stringsAsFactors = FALSE)
    pos <- which(df$statistic == "VPC/ICC")
    if (length(pos) == 1 && pos < nrow(df)) {
      df <- rbind(df[seq_len(pos), , drop = FALSE], pcv_row,
                  df[(pos + 1):nrow(df), , drop = FALSE])
    } else {
      df <- rbind(df, pcv_row)
    }
  }

  rownames(df) <- NULL
  df
}

# Rank every stratum by its model-predicted outcome, reusing the same per-stratum
# predictions as plot(type = "predicted") so the table matches the figure. Returns
# a data frame ordered from highest to lowest predicted value, with the predicted
# value (and its conditional interval), the stratum size, and the stratum random
# effect (BLUP, and its interval).
maihda_strata_ranking <- function(object, summary_obj, scale = c("response", "link")) {
  scale <- match.arg(scale)
  if (!is.null(object$longitudinal_info)) {
    maihda_stop_longitudinal_scalar("A ranked-strata table")
  }
  pred <- switch(object$engine,
    lme4 = maihda_stratum_predictions_lme4(object, summary_obj, scale = scale),
    brms = maihda_stratum_predictions_brms(object, summary_obj, scale = scale),
    wemix = maihda_stratum_predictions_wemix(object, summary_obj, scale = scale),
    ordinal = maihda_stratum_predictions_ordinal(object, summary_obj, scale = scale),
    stop("Unsupported engine: ", object$engine, call. = FALSE))

  se <- summary_obj$stratum_estimates
  if (is.null(se) || nrow(se) == 0) {
    stop("No stratum estimates available to rank.", call. = FALSE)
  }

  idx <- match(as.character(se$stratum), as.character(pred$stratum))
  out <- data.frame(
    stratum = as.character(se$stratum),
    label = if ("label" %in% names(se)) as.character(se$label) else NA_character_,
    n = pred$n[idx],
    predicted = pred$predicted_row[idx],
    predicted_lower = pred$lower_row[idx],
    predicted_upper = pred$upper_row[idx],
    random_effect = se$random_effect,
    re_lower = if ("lower_95" %in% names(se)) se$lower_95 else NA_real_,
    re_upper = if ("upper_95" %in% names(se)) se$upper_95 else NA_real_,
    stringsAsFactors = FALSE
  )

  out <- out[order(-out$predicted), , drop = FALSE]
  out <- data.frame(rank = seq_len(nrow(out)), out, row.names = NULL,
                    stringsAsFactors = FALSE)
  out
}

# Format one estimate (+ optional interval) for printing.
maihda_table_fmt <- function(est, lo = NA_real_, hi = NA_real_, digits = 3) {
  if (length(est) == 0 || is.na(est)) return("")
  s <- formatC(est, format = "f", digits = digits)
  if (!is.na(lo) && !is.na(hi)) {
    s <- sprintf("%s [%s, %s]", s,
                 formatC(lo, format = "f", digits = digits),
                 formatC(hi, format = "f", digits = digits))
  }
  s
}

#' Print a MAIHDA results table
#'
#' @param x A \code{maihda_table} object from \code{\link{maihda_table}}.
#' @param digits Decimal places (defaults to the value stored on \code{x}).
#' @param ... Additional arguments (not used).
#' @return Invisibly, \code{x}.
#' @export
print.maihda_table <- function(x, digits = x$digits, ...) {
  pal <- maihda_palette()
  cat(pal$bold("MAIHDA Results Table"), "\n", sep = "")
  cat("====================\n\n")

  cat(sprintf("Engine: %s | Family: %s | Mode: %s\n",
              x$engine, x$family, x$mode))
  cat(sprintf("Observations: %s | Strata: %s\n",
              format(x$n_obs), format(x$n_strata_total)))
  if (!is.null(x$context_vars)) {
    cat(sprintf("Context: %s (crossed contextual random intercept)\n",
                paste(x$context_vars, collapse = ", ")))
  }

  # --- Model-results table ---------------------------------------------------
  disp <- data.frame(Statistic = x$models$statistic,
                     check.names = FALSE, stringsAsFactors = FALSE)
  for (k in x$model_keys) {
    vals <- vapply(seq_len(nrow(x$models)), function(i) {
      maihda_table_fmt(x$models[[k]][i], x$models[[paste0(k, "_lower")]][i],
                       x$models[[paste0(k, "_upper")]][i], digits = digits)
    }, character(1))
    disp[[x$model_labels[[k]]]] <- vals
  }
  cat("\n", pal$bold("Model results:"), "\n", sep = "")
  print(disp, row.names = FALSE)

  # --- Ranked-strata table ---------------------------------------------------
  if (!is.null(x$strata) && nrow(x$strata) > 0) {
    st <- x$strata
    scale_lab <- if (identical(x$scale, "response")) "predicted value" else "predicted (link scale)"
    by_lab <- switch(x$ranked_by, null = "null model", adjusted = "adjusted model", "model")
    cat("\n", pal$bold(sprintf("Strata ranked by %s (%s):", scale_lab, by_lab)), "\n", sep = "")

    fmt_row <- function(rows) {
      data.frame(
        Rank = as.character(rows$rank),
        Stratum = ifelse(!is.na(rows$label), rows$label, rows$stratum),
        N = as.character(rows$n),
        Predicted = vapply(seq_len(nrow(rows)), function(i)
          maihda_table_fmt(rows$predicted[i], rows$predicted_lower[i],
                           rows$predicted_upper[i], digits = digits), character(1)),
        `Stratum RE` = vapply(seq_len(nrow(rows)), function(i)
          maihda_table_fmt(rows$random_effect[i], rows$re_lower[i],
                           rows$re_upper[i], digits = digits), character(1)),
        check.names = FALSE, stringsAsFactors = FALSE)
    }

    n <- x$n_strata
    if (nrow(st) > 2 * n) {
      top <- fmt_row(utils::head(st, n))
      bottom <- fmt_row(utils::tail(st, n))
      sep <- top[1, ]
      sep[] <- "..."
      combined <- rbind(top, sep, bottom)
      print(combined, row.names = FALSE)
      cat(sprintf("  (%d strata between ranks %d and %d not shown; see $strata)\n",
                  nrow(st) - 2 * n, n + 1, nrow(st) - n))
    } else {
      print(fmt_row(st), row.names = FALSE)
    }
    cat("  Predicted intervals are conditional (random-effect) only; ",
        "Stratum RE is the stratum BLUP.\n", sep = "")
  } else if (!is.null(x$strata_note)) {
    cat("\nRanked strata: ", x$strata_note, "\n", sep = "")
  }

  cat("\nEstimates are point values unless a [low, high] interval is shown ",
      "(VPC/PCV).\n", sep = "")
  invisible(x)
}
