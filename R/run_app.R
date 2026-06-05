maihda_app_required_packages <- function() {
  c("shiny", "bslib", "DT", "future", "promises", "shinyjs", "plotly", "ggtern")
}

maihda_app_pvc_display <- function(pvc_percent) {
  pvc_percent <- suppressWarnings(as.numeric(pvc_percent)[1])
  fmt_percent <- function(x) paste0(round(x, 2), "%")

  if (!is.finite(pvc_percent)) {
    return(list(
      label = "Residual Strata Variance",
      value = "N/A",
      description = "Between-strata variance remaining after adjusting for main effects",
      remaining_value = "N/A",
      status = "unknown"
    ))
  }

  remaining_percent <- 100 - pvc_percent
  if (pvc_percent < 0) {
    return(list(
      label = "Unmasked Variance",
      value = paste0("+", round(abs(pvc_percent), 2), "%"),
      description = "Increase in between-strata variance after adjustment",
      remaining_value = fmt_percent(remaining_percent),
      status = "negative"
    ))
  }

  list(
    label = "Residual Strata Variance",
    value = fmt_percent(remaining_percent),
    description = "Between-strata variance not explained by main effects",
    remaining_value = fmt_percent(remaining_percent),
    status = "nonnegative"
  )
}

maihda_app_ternary_plotly <- function(td) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required to render the interactive ternary plot.",
         call. = FALSE)
  }

  required_cols <- c("additive_prop", "interaction_prop", "uncertainty_prop", "label", "n")
  missing_cols <- setdiff(required_cols, names(td))
  if (length(missing_cols) > 0) {
    stop("Ternary plot data is missing required columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  marker_sizes <- pmax(sqrt(td$n) * 2, 4)
  marker_colors <- as.numeric(as.factor(td$label))
  hover_text <- paste0(
    "<b>Stratum:</b> ", td$label, "<br>",
    "<b>Size (N):</b> ", td$n, "<br>",
    "<b>Additive:</b> ", round(td$additive_prop * 100, 1), "%<br>",
    "<b>Intersection-specific:</b> ", round(td$interaction_prop * 100, 1), "%<br>",
    "<b>Uncertainty:</b> ", round(td$uncertainty_prop * 100, 1), "%"
  )

  plotly::plot_ly(
    data = td,
    type = "scatterternary",
    mode = "markers",
    a = td$additive_prop,
    b = td$interaction_prop,
    c = td$uncertainty_prop,
    text = hover_text,
    hoverinfo = "text",
    marker = list(
      size = marker_sizes,
      color = marker_colors,
      colorscale = "Viridis",
      opacity = 0.8,
      line = list(color = "rgba(0,0,0,0.5)", width = 1)
    )
  ) |>
    plotly::layout(
      title = "Interactive MAIHDA Strata Effects Decomposition",
      ternary = list(
        sum = 1,
        aaxis = list(title = "Additive", min = 0, linewidth = 2, ticks = "outside", tickvals = seq(0, 1, by = 0.2)),
        baxis = list(title = "Intersection", min = 0, linewidth = 2, ticks = "outside", tickvals = seq(0, 1, by = 0.2)),
        caxis = list(title = "Uncertainty", min = 0, linewidth = 2, ticks = "outside", tickvals = seq(0, 1, by = 0.2))
      ),
      margin = list(t = 50, b = 50, l = 50, r = 50)
    )
}

maihda_app_fit_models <- function(dat, outcome_var, grouping_vars,
                                  additional_covars = character(),
                                  family = "gaussian", use_boot = FALSE,
                                  n_boot = 100, autobin = TRUE,
                                  engine = "lme4") {
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data frame.", call. = FALSE)
  }
  if (!is.character(outcome_var) || length(outcome_var) != 1 || !outcome_var %in% names(dat)) {
    stop("'outcome_var' must name one column in 'dat'.", call. = FALSE)
  }
  if (!is.character(grouping_vars) || length(grouping_vars) == 0) {
    stop("'grouping_vars' must contain at least one column name.", call. = FALSE)
  }

  additional_covars <- if (is.null(additional_covars)) character() else additional_covars
  all_required_cols <- unique(c(outcome_var, grouping_vars, additional_covars))
  missing_cols <- setdiff(all_required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Variables not found in data: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  complete_dat <- dat[stats::complete.cases(dat[, all_required_cols, drop = FALSE]), , drop = FALSE]
  if (nrow(complete_dat) == 0) {
    stop("No complete cases remaining after omitting missing values (NAs). Please select different variables.",
         call. = FALSE)
  }

  strata_dat <- make_strata(complete_dat, vars = grouping_vars, autobin = autobin)
  model_dat <- complete_dat
  model_dat$stratum <- strata_dat$data$stratum
  attr(model_dat, "strata_info") <- strata_dat$strata_info
  attr(model_dat, "strata_vars") <- strata_dat$vars
  attr(model_dat, "strata_sep") <- strata_dat$sep
  attr(model_dat, "strata_autobin_info") <- strata_dat$autobin_info

  # The app always passes a family, which suppresses fit_maihda()'s automatic
  # binary detection. Mirror that detection here for the default gaussian family so
  # a binary outcome -- whether a two-level factor/character or a numeric 0/1 -- is
  # fit as binomial rather than silently as a linear probability model. This keeps
  # the no-code app consistent with the core API and avoids surprising LPM fits. To
  # fit an LPM intentionally, call fit_maihda(..., family = "gaussian") from R.
  if (identical(family, "gaussian") &&
      maihda_is_binary_vector(complete_dat[[outcome_var]])) {
    message("maihda_app: outcome '", outcome_var,
            "' is binary; using family = 'binomial'. ",
            "For a linear probability model, fit from R with ",
            "fit_maihda(..., family = 'gaussian').")
    family <- "binomial"
  }

  adjusted_fmla <- maihda_formula_with_stratum(outcome_var, c(grouping_vars, additional_covars))
  null_fmla <- maihda_formula_with_stratum(outcome_var)

  null_model <- fit_maihda(formula = null_fmla, data = model_dat, engine = engine, family = family)
  adjusted_model <- fit_maihda(formula = adjusted_fmla, data = model_dat, engine = engine, family = family)
  pvc <- calculate_pvc(null_model, adjusted_model, bootstrap = use_boot, n_boot = n_boot)
  stepwise <- stepwise_pcv(
    model_dat,
    outcome = outcome_var,
    vars = c(grouping_vars, additional_covars),
    engine = engine,
    family = family
  )

  list(null_model = null_model, model = adjusted_model, pvc = pvc, stepwise = stepwise)
}

#' Run MAIHDA Shiny Application
#'
#' @description
#' Launches a Shiny graphical user interface that exposes core functions of the
#' MAIHDA package, allowing for visual data exploration, model fitting, and
#' performance visualization.
#'
#' @return No return value, called to launch the shiny app.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' run_maihda_app()
#' }
run_maihda_app <- function() {
  required_pkgs <- maihda_app_required_packages()
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required to run the MAIHDA Shiny app:\n",
      paste("  -", missing_pkgs, collapse = "\n"),
      "\n\nPlease install them by running:\n",
      "install.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))",
      call. = FALSE
    )
  }

  app_dir <- system.file("shiny", package = "MAIHDA")
  if (app_dir == "") {
    stop("Could not find shiny app directory. Try re-installing `MAIHDA`.", call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
