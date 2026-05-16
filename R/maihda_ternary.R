#' Compute Ternary Data for MAIHDA Models
#'
#' @param model A fitted MAIHDA model object from `fit_maihda()`.
#' @param summary_obj Optional output from `summary()`.
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
    summary_obj = NULL,
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

  if (verbose && scale == "response") {
    warning("Ternary decomposition is most coherent on the link scale.")
  }

  fitted_mod <- model$model
  if (is.null(fitted_mod)) stop("Could not find fitted model within maihda_model object.")

  # Try to retrieve strata info
  has_strata_info <- exists("strata_info", where = model) || !is.null(model$strata_info) || exists("data", where = model)

  if (has_strata_info && !is.null(model$strata_info)) {
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
  pred_data <- pred_data[!is.na(pred_data$stratum), , drop = FALSE]

  if (!is.null(reference_values)) {
    strata_levels <- re_df$stratum
    first_idx <- match(strata_levels, as.character(pred_data$stratum))
    pred_data <- pred_data[first_idx, , drop = FALSE]

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

  additive_by_stratum <- stats::aggregate(
    x = list(
      additive_only = as.numeric(additive_values),
      full_prediction = as.numeric(full_values),
      uncertainty = as.numeric(uncertainty_values)
    ),
    by = list(stratum = as.character(pred_data$stratum)),
    FUN = mean,
    na.rm = TRUE
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

  grand_mean_additive <- mean(re_df$additive_only, na.rm = TRUE)

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

  tibble::as_tibble(res)
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
#' @param ternary_data Data output from \code{compute_maihda_ternary_data}.
#' @param size_var Column name for point sizing.
#' @param color_var Column name for point colors.
#' @param label_top_n Number of top strata to label.
#' @param label_by Variable used to determine top strata.
#' @param alpha Point transparency.
#'
#' @return A plot object.
#' @export
plot_maihda_ternary <- function(
    ternary_data,
    size_var = "n",
    color_var = "label",
    label_top_n = 5,
    label_by = c("interaction_signal", "uncertainty", "n"),
    alpha = 0.7
) {
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
      caption = "Each point is a stratum. Proximity to a corner indicates higher proportion of that component.\nAdditive: Fixed effects only. Intersection: Random effect magnitude. Uncertainty: selected uncertainty metric."
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

#' Generate Ternary Plot from MAIHDA Model
#'
#' @param model A fitted MAIHDA model.
#' @param summary_obj Optional output from \code{summary_maihda}.
#' @param ... Additional arguments passed to \code{compute_maihda_ternary_data} and \code{plot_maihda_ternary}.
#'
#' @return A list containing \code{data} and \code{plot}.
#' @export
maihda_ternary_plot <- function(model, summary_obj = NULL, ...) {
  args <- list(...)

  compute_args <- args[names(args) %in% names(formals(compute_maihda_ternary_data))]
  plot_args <- args[names(args) %in% names(formals(plot_maihda_ternary))]

  compute_args$model <- model
  compute_args$summary_obj <- summary_obj

  ternary_data <- do.call(compute_maihda_ternary_data, compute_args)

  plot_args$ternary_data <- ternary_data
  p <- do.call(plot_maihda_ternary, plot_args)

  list(
    data = ternary_data,
    plot = p
  )
}
