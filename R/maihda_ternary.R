#' Compute Ternary Data for MAIHDA Models
#'
#' @param model A fitted MAIHDA model object from `fit_maihda()`.
#' @param scale Character, either "link" or "response".
#' @param reference_values List or data.frame of reference values for covariates.
#' @param uncertainty_method Character indicating how to extract uncertainty.
#'   "auto" uses conditional standard errors for lme4 models and posterior
#'   standard deviations for brms models. "ci_width" uses the 95\% interval width.
#' @param include_na_strata Logical, whether to include strata with missing data.
#' @param verbose Logical, whether to print messages.
#'
#' @return A tidy tibble with ternary coordinates.
#' @export
#'
#' @importFrom dplyr select mutate left_join bind_rows distinct all_of across where as_tibble
#' @importFrom stats predict
compute_maihda_ternary_data <- function(
    model,
    scale = c("link", "response"),
    reference_values = NULL,
    uncertainty_method = c("auto", "se", "ci_width", "posterior_sd"),
    include_na_strata = FALSE,
    verbose = TRUE
) {
  scale <- match.arg(scale)
  uncertainty_method <- match.arg(uncertainty_method)

  if (!inherits(model, "maihda_model")) {
    stop("model must be a maihda_model object.")
  }

  engine <- model$engine
  if (is.null(engine)) engine <- "unknown"

  if (identical(engine, "ordinal")) {
    stop("The ternary diagnostic is not yet supported for the ordinal (clmm) ",
         "engine. Use plot(type = \"predicted\"), \"risk_vs_effect\", ",
         "\"effect_decomp\", or plot_prediction_deviation_panels() for a ",
         "cumulative MAIHDA model.", call. = FALSE)
  }

  if (verbose && scale == "response") {
    warning("Ternary decomposition is most coherent on the link scale.")
  }

  fitted_mod <- model$model
  if (is.null(fitted_mod)) stop("Could not find fitted model within maihda_model object.")

  # Retrieve strata info: prefer the stored strata_info table, otherwise rebuild
  # minimal stratum counts from the model data below.
  if (!is.null(model$strata_info)) {
      strata_info_df <- model$strata_info
  } else if (!is.null(model$data) && "stratum" %in% names(model$data)) {
      strata_counts <- table(model$data$stratum)
      strata_info_df <- data.frame(
          stratum = names(strata_counts),
          n = as.numeric(strata_counts),
          stringsAsFactors = FALSE
      )
  } else {
      stop("Cannot derive strata info. Missing base data with 'stratum' column.")
  }

  u_j_raw <- NULL
  u_j_se <- NULL

  if (engine == "lme4" || inherits(fitted_mod, "merMod")) {
    re_stratum <- maihda_stratum_ranef_lme4(fitted_mod)
    re_df <- data.frame(
      stratum = re_stratum$stratum,
      u_j = re_stratum$random_effect,
      lower_95 = re_stratum$lower_95,
      upper_95 = re_stratum$upper_95,
      uncertainty = maihda_ternary_uncertainty(re_stratum, uncertainty_method, "lme4"),
      stringsAsFactors = FALSE
    )

  } else if (engine == "brms" || inherits(fitted_mod, "brmsfit")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms package is required for brms engine models.")
    }
    re_stratum <- maihda_stratum_ranef_brms(fitted_mod)
    re_df <- data.frame(
      stratum = re_stratum$stratum,
      u_j = re_stratum$random_effect,
      lower_95 = re_stratum$lower_95,
      upper_95 = re_stratum$upper_95,
      uncertainty = maihda_ternary_uncertainty(re_stratum, uncertainty_method, "brms"),
      stringsAsFactors = FALSE
    )
  } else {
    stop(sprintf("Engine '%s' is not fully supported for ternary plots yet.", engine))
  }

  pred_data <- model$data
  stratum_keep <- !is.na(pred_data$stratum)
  pred_data <- pred_data[stratum_keep, , drop = FALSE]
  # Prior/precision weights aligned to pred_data, so the per-stratum aggregation
  # and the reference centre below are weighted for a weighted fit (consistent with
  # the weighted VPC and the other stratum-level plots). Unit weights reduce both
  # to the previous plain/size-weighted means. These are lme4 prior/precision
  # weights, not a complex survey design (no design-based variance is computed).
  prior_w_full <- maihda_prior_weights(model)
  pred_w <- prior_w_full[stratum_keep]

  if (!is.null(reference_values)) {
    strata_levels <- re_df$stratum
    first_idx <- match(strata_levels, as.character(pred_data$stratum))
    pred_data <- pred_data[first_idx, , drop = FALSE]
    pred_w <- pred_w[first_idx]

    if (is.data.frame(reference_values)) {
      ref <- reference_values
      if (nrow(ref) == 1) {
        for (nm in names(ref)) {
          pred_data[[nm]] <- ref[[nm]][1]
        }
      } else if ("stratum" %in% names(ref)) {
        ref_idx <- match(strata_levels, as.character(ref$stratum))
        for (nm in setdiff(names(ref), "stratum")) {
          pred_data[[nm]] <- ref[[nm]][ref_idx]
        }
      } else if (nrow(ref) == nrow(pred_data)) {
        for (nm in names(ref)) {
          pred_data[[nm]] <- ref[[nm]]
        }
      } else {
        stop("reference_values must have one row, one row per stratum, or a 'stratum' column.")
      }
    } else if (is.list(reference_values)) {
      for (nm in names(reference_values)) {
        pred_data[[nm]] <- reference_values[[nm]][1]
      }
    } else {
      stop("reference_values must be a list or data frame.")
    }
  }

  fam <- maihda_family(fitted_mod)
  linkinv <- maihda_linkinv(fam)

  if (engine == "lme4" || inherits(fitted_mod, "merMod")) {
      fe_link <- stats::predict(fitted_mod, newdata = pred_data, re.form = NA, type = "link")
  } else if (engine == "brms" || inherits(fitted_mod, "brmsfit")) {
      fe_link <- brms::posterior_linpred(fitted_mod, newdata = pred_data, re_formula = NA, summary = TRUE)[, "Estimate"]
  }

  re_idx <- match(as.character(pred_data$stratum), as.character(re_df$stratum))
  u_by_row <- re_df$u_j[re_idx]

  # Cross-classified model: fold the dimension random effects into the additive
  # baseline so the "additive signal" reflects the dimensions' (random) main effects,
  # not just the fixed covariates. fe_link becomes fixed + dimension REs =
  # eta(all REs) - u_stratum; the interaction signal (u_j, the stratum RE) is unchanged.
  if (!is.null(model$cc_info)) {
    total_link <- if (engine == "brms" || inherits(fitted_mod, "brmsfit")) {
      brms::posterior_linpred(fitted_mod, newdata = pred_data, summary = TRUE)[, "Estimate"]
    } else {
      stats::predict(fitted_mod, newdata = pred_data, type = "link")
    }
    fe_link <- as.numeric(total_link) - u_by_row
  }

  if (scale == "response") {
    additive_values <- linkinv(fe_link)
    full_values <- linkinv(fe_link + u_by_row)

    if (uncertainty_method == "ci_width") {
      uncertainty_values <- abs(
        linkinv(fe_link + re_df$upper_95[re_idx]) -
          linkinv(fe_link + re_df$lower_95[re_idx])
      )
    } else {
      uncertainty_values <- abs(
        linkinv(fe_link + u_by_row + re_df$uncertainty[re_idx]) -
          linkinv(fe_link + u_by_row)
      )
    }
  } else {
    additive_values <- as.numeric(fe_link)
    full_values <- as.numeric(fe_link + u_by_row)
    uncertainty_values <- re_df$uncertainty[re_idx]
  }

  # Prior-weight-weighted per-stratum means, so a weighted fit's stratum
  # additive/full predictions match the weighted VPC; reduces to plain means when
  # the fit is unweighted.
  additive_by_stratum <- maihda_weighted_stratum_aggregate(
    data.frame(
      stratum = as.character(pred_data$stratum),
      additive_only = as.numeric(additive_values),
      full_prediction = as.numeric(full_values),
      uncertainty = as.numeric(uncertainty_values),
      weight = pred_w,
      stringsAsFactors = FALSE
    ),
    c("additive_only", "full_prediction", "uncertainty")
  )
  re_df$additive_only <- additive_by_stratum$additive_only[
    match(as.character(re_df$stratum), additive_by_stratum$stratum)
  ]
  re_df$full_prediction <- additive_by_stratum$full_prediction[
    match(as.character(re_df$stratum), additive_by_stratum$stratum)
  ]
  if (scale == "response") {
    re_df$uncertainty <- additive_by_stratum$uncertainty[
      match(as.character(re_df$stratum), additive_by_stratum$stratum)
    ]
  }

  # Reference is the mean additive prediction across strata, weighted so rare and
  # common strata are not given equal weight. For a weighted fit use each stratum's
  # summed prior weights (so the centre is on the same weighted footing as the
  # per-stratum aggregation above); otherwise use the stratum size. Both reduce to
  # the plain mean when neither is available.
  is_weighted <- !isTRUE(all.equal(prior_w_full, rep(1, length(prior_w_full))))
  strata_w <- if (is_weighted) {
    pop_w <- tapply(prior_w_full[stratum_keep],
                    as.character(model$data$stratum[stratum_keep]),
                    sum, na.rm = TRUE)
    as.numeric(pop_w[as.character(re_df$stratum)])
  } else if (!is.null(strata_info_df) &&
             all(c("stratum", "n") %in% names(strata_info_df))) {
    strata_info_df$n[match(as.character(re_df$stratum),
                           as.character(strata_info_df$stratum))]
  } else {
    NULL
  }
  grand_mean_additive <- if (!is.null(strata_w) && any(is.finite(strata_w))) {
    stats::weighted.mean(re_df$additive_only, strata_w, na.rm = TRUE)
  } else {
    mean(re_df$additive_only, na.rm = TRUE)
  }

  res <- re_df
  res$grand_mean_additive <- grand_mean_additive
  res$additive_signal <- abs(res$additive_only - grand_mean_additive)
  res$interaction_signal <- if (scale == "response") {
    abs(res$full_prediction - res$additive_only)
  } else {
    abs(res$u_j)
  }
  res$u_sign <- ifelse(res$u_j >= 0, "Positive", "Negative")

  if (!include_na_strata) {
    res <- res[!is.na(res$u_j) & !is.na(res$uncertainty), ]
  }

  row_sums <- res$additive_signal + res$interaction_signal + res$uncertainty
  res$additive_prop <- ifelse(row_sums > 0, res$additive_signal / row_sums, NA_real_)
  res$interaction_prop <- ifelse(row_sums > 0, res$interaction_signal / row_sums, NA_real_)
  res$uncertainty_prop <- ifelse(row_sums > 0, res$uncertainty / row_sums, NA_real_)

  if (!is.null(strata_info_df) && "n" %in% names(strata_info_df)) {
      res <- maihda_add_strata_columns(res, strata_info_df)
      if (!"label" %in% names(res) || all(is.na(res$label))) {
        strat_vars <- setdiff(names(strata_info_df), c("stratum", "n", "label"))
        if (length(strat_vars) > 0) {
            res$label <- apply(res[, strat_vars, drop = FALSE], 1, paste, collapse = "\n")
        } else {
            res$label <- as.character(res$stratum)
        }
      }
  } else {
      res$n <- 1
      res$label <- as.character(res$stratum)
  }

  res <- res[order(res$interaction_signal, decreasing = TRUE), ]
  rownames(res) <- NULL

  out <- tibble::as_tibble(res)
  # Class the result so plot() dispatches to plot.maihda_ternary(). It remains a
  # tibble/data.frame, so existing column access is unaffected.
  class(out) <- c("maihda_ternary", class(out))
  out
}

maihda_ternary_uncertainty <- function(re_stratum, uncertainty_method, engine) {
  method <- uncertainty_method
  if (method == "auto") {
    method <- if (engine == "brms") "posterior_sd" else "se"
  }

  if (method == "se") {
    return(re_stratum$se)
  }
  if (method == "ci_width") {
    return(abs(re_stratum$upper_95 - re_stratum$lower_95))
  }
  if (method == "posterior_sd") {
    if (engine != "brms") {
      stop("uncertainty_method = 'posterior_sd' is only available for brms models.",
           call. = FALSE)
    }
    return(re_stratum$se)
  }

  stop("Unsupported uncertainty_method: ", uncertainty_method, call. = FALSE)
}

#' Plot MAIHDA Ternary Diagram
#'
#' Renders the ternary decomposition produced by
#' \code{\link{compute_maihda_ternary_data}}. Dispatched via \code{plot()} on the
#' classed result.
#'
#' @note This method \strong{attaches the \pkg{ggtern} package to the search
#'   path} (as if by \code{library(ggtern)}) if it is not already attached. This
#'   is a deliberate, unavoidable side effect: \pkg{ggtern} replaces several
#'   \pkg{ggplot2} build/print internals at attach time, and without it the
#'   ternary coordinate system and themes do not render (you get a blank or
#'   distorted plot). The attachment persists after the call so the returned
#'   object can still be printed later in the session; it is not detached on
#'   exit. If you need a pristine search path, attach \pkg{ggtern} yourself
#'   before plotting and manage its lifecycle, or run plotting in a separate
#'   session.
#'
#' @param x A \code{maihda_ternary} object from \code{compute_maihda_ternary_data}.
#' @param size_var Column name for point sizing.
#' @param color_var Column name for point colors.
#' @param label_top_n Number of top strata to label.
#' @param label_by Variable used to determine top strata.
#' @param alpha Point transparency.
#' @param ... Additional arguments (not used).
#'
#' @return A plot object.
#' @export
plot.maihda_ternary <- function(
    x,
    size_var = "n",
    color_var = "label",
    label_top_n = 5,
    label_by = c("interaction_signal", "uncertainty", "n"),
    alpha = 0.7,
    ...
) {
  ternary_data <- x
  label_by <- match.arg(label_by)

  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("ggtern package is not installed. Please install ggtern.")
  }

  # Crucial to attach ggtern dynamically so ternary coordinates and themes plot correctly instead of resulting in a blank plot
  suppressPackageStartupMessages(requireNamespace("ggtern", quietly = TRUE))
  if(!"ggtern" %in% .packages()) {
    suppressPackageStartupMessages(attachNamespace("ggtern"))
  }

  p <- ggtern::ggtern(data = ternary_data, ggtern::aes(x = .data$additive_prop, y = .data$interaction_prop, z = .data$uncertainty_prop)) +
    ggplot2::geom_point(ggplot2::aes(size = .data[[size_var]], color = .data[[color_var]]), alpha = alpha) +
    ggtern::theme_bw() +
    ggtern::theme_showgrid() +
    ggtern::theme_showarrows() +
    ggtern::theme_clockwise() +
    ggplot2::labs(
      x = "Additive signal",
      y = "Intersection-specific signal",
      z = "Uncertainty",
      color = "Stratum",
      size = "Sample Size (N)",
      title = "MAIHDA Strata Effects Decomposition",
      caption = ""
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", size = 11, color = "#333333"),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 40)
    )

  if (label_top_n > 0) {
    top_idx <- seq_len(min(label_top_n, nrow(ternary_data)))
    top_data <- ternary_data[order(ternary_data[[label_by]], decreasing = TRUE), ][top_idx, ]
    p <- p + ggplot2::geom_text(data = top_data, ggplot2::aes(label = .data[["label"]]), size = 2.5, vjust = -1, color = "#222222")
  }

  p <- p + ggplot2::guides(color = "none") # Hide the large color legend, but keep the size legend
  p <- p + ggtern::theme_legend_position("middleright")
  return(p)
}

#' Plot MAIHDA Ternary Diagram (deprecated)
#'
#' Deprecated. Use \code{plot()} on the \code{\link{compute_maihda_ternary_data}}
#' result instead, e.g. \code{plot(compute_maihda_ternary_data(model))}.
#'
#' @param ternary_data Data output from \code{compute_maihda_ternary_data}.
#' @param ... Further arguments passed to \code{plot()} (e.g. \code{size_var}).
#' @return A plot object.
#' @keywords internal
#' @export
plot_maihda_ternary <- function(ternary_data, ...) {
  .Deprecated("plot", msg = paste(
    "'plot_maihda_ternary()' is deprecated.",
    "Use plot() on the compute_maihda_ternary_data() result, e.g.",
    "plot(compute_maihda_ternary_data(model))."
  ))
  if (!inherits(ternary_data, "maihda_ternary")) {
    class(ternary_data) <- c("maihda_ternary", class(ternary_data))
  }
  plot(ternary_data, ...)
}

#' Generate Ternary Plot from MAIHDA Model
#'
#' @param model A fitted MAIHDA model.
#' @param ... Additional arguments passed to \code{compute_maihda_ternary_data} and \code{\link{plot.maihda_ternary}}.
#'
#' @return A list containing \code{data} and \code{plot}.
#' @export
maihda_ternary_plot <- function(model, ...) {
  args <- list(...)

  compute_args <- args[names(args) %in% names(formals(compute_maihda_ternary_data))]
  plot_args <- args[names(args) %in% names(formals(plot.maihda_ternary))]

  compute_args$model <- model

  ternary_data <- do.call(compute_maihda_ternary_data, compute_args)

  plot_args$x <- ternary_data
  p <- do.call(plot.maihda_ternary, plot_args)

  list(
    data = ternary_data,
    plot = p
  )
}
