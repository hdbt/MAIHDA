#' Compute Ternary Data for MAIHDA Models
#'
#' @param model A fitted MAIHDA model object from `fit_maihda()`.
#' @param summary_obj Optional output from `summary_maihda()`.
#' @param scale Character, either "link" or "response".
#' @param reference_values List or data.frame of reference values for covariates.
#' @param uncertainty_method Character indicating how to extract uncertainty.
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
    re_terms <- lme4::ranef(fitted_mod, condVar = TRUE)
    if (!"stratum" %in% names(re_terms)) {
      stop("Could not find 'stratum' random effects in the lme4 model.")
    }
    re_stratum <- re_terms[["stratum"]]

    cond_var <- attr(re_stratum, "postVar")
    strata_names <- rownames(re_stratum)
    u_j_raw <- re_stratum[, 1]

    if (!is.null(cond_var)) {
      u_j_se <- sqrt(cond_var[1, 1, ])
    } else {
      u_j_se <- rep(NA_real_, length(u_j_raw))
    }

    re_df <- data.frame(
      stratum = strata_names,
      u_j = u_j_raw,
      uncertainty = u_j_se,
      stringsAsFactors = FALSE
    )

  } else if (engine == "brms" || inherits(fitted_mod, "brmsfit")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("brms package is required for brms engine models.")
    }
    re_summary <- brms::ranef(fitted_mod, summary = TRUE)
    if (!"stratum" %in% names(re_summary)) {
      stop("Could not find 'stratum' random effects in the brms model.")
    }
    brms_re <- re_summary[["stratum"]]
    strata_names <- dimnames(brms_re)[[1]]

    u_j_raw <- brms_re[, "Estimate", 1]
    u_j_se <- brms_re[, "Est.Error", 1]

    re_df <- data.frame(
      stratum = strata_names,
      u_j = u_j_raw,
      uncertainty = u_j_se,
      stringsAsFactors = FALSE
    )
  } else {
    stop(sprintf("Engine '%s' is not fully supported for ternary plots yet.", engine))
  }

  if (engine == "lme4") {
      pred_data <- model$data
      pred_data <- pred_data[!duplicated(pred_data$stratum), ]

      fe_preds <- stats::predict(fitted_mod, newdata = pred_data, re.form = NA, type = scale)
      re_df$additive_only <- fe_preds[match(re_df$stratum, pred_data$stratum)]
  } else if (engine == "brms") {
      pred_data <- model$data
      pred_data <- pred_data[!duplicated(pred_data$stratum), ]
      fe_preds <- stats::predict(fitted_mod, newdata = pred_data, re_formula = NA)
      re_df$additive_only <- fe_preds[match(re_df$stratum, pred_data$stratum), "Estimate"]
  }

  grand_mean_additive <- mean(re_df$additive_only, na.rm = TRUE)

  res <- re_df
  res$grand_mean_additive <- grand_mean_additive
  res$additive_signal <- abs(res$additive_only - grand_mean_additive)
  res$interaction_signal <- abs(res$u_j)
  res$u_sign <- ifelse(res$u_j >= 0, "Positive", "Negative")

  if (!include_na_strata) {
    res <- res[!is.na(res$u_j) & !is.na(res$uncertainty), ]
  }

  row_sums <- res$additive_signal + res$interaction_signal + res$uncertainty
  res$additive_prop <- res$additive_signal / row_sums
  res$interaction_prop <- res$interaction_signal / row_sums
  res$uncertainty_prop <- res$uncertainty / row_sums

  if (!is.null(strata_info_df) && "n" %in% names(strata_info_df)) {
      res <- merge(res, strata_info_df, by = "stratum", all.x = TRUE)

      strat_vars <- setdiff(names(strata_info_df), c("stratum", "n"))
      if (length(strat_vars) > 0) {
          res$label <- apply(res[, strat_vars, drop = FALSE], 1, paste, collapse = "\n")
      } else {
          res$label <- as.character(res$stratum)
      }
  } else {
      res$n <- 1
      res$label <- as.character(res$stratum)
  }

  res <- res[order(res$interaction_signal, decreasing = TRUE), ]
  rownames(res) <- NULL

  tibble::as_tibble(res)
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

  p <- ggtern::ggtern(data = ternary_data, ggtern::aes(x = additive_prop, y = interaction_prop, z = uncertainty_prop)) +
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
      caption = "Each point is a stratum. Proximity to a corner indicates higher proportion of that component.\nAdditive: Fixed effects only. Intersection: Random effect magnitude. Uncertainty: Standard error."
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", size = 11, color = "#333333"),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 40)
    )

  if (label_top_n > 0) {
    top_data <- ternary_data[order(ternary_data[[label_by]], decreasing = TRUE), ][1:label_top_n, ]
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
