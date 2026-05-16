maihda_app_required_packages <- function() {
  c("shiny", "bslib", "DT", "future", "promises", "shinyjs", "plotly", "ggtern")
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
