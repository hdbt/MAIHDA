maihda_app_required_packages <- function() {
  c("shiny", "bslib", "DT", "future", "promises", "shinyjs", "plotly", "ggtern",
    "shinycssloaders")
}

# Sensible default outcome + grouping variables for the dashboard's variable
# pickers. Built-in datasets get curated defaults (validated against the data's
# actual columns); anything else -- uploads, or a built-in whose columns changed --
# falls back to a heuristic: the first column as the outcome and the first two
# "categorical-ish" columns (factor/character/logical, or numeric with <= 10
# distinct values) as the strata. This replaces the previous logic, which hardcoded
# the simulated dataset's lowercase names ("health_outcome", "gender", "race") and
# therefore left the NHANES dataset (Obese, Gender, Race, ...) with no usable
# defaults, disabling the Fit button.
maihda_app_default_vars <- function(dataset, data) {
  if (!is.data.frame(data) || ncol(data) == 0) {
    return(list(outcome = NULL, groups = character(0)))
  }
  cols <- names(data)

  known <- switch(
    as.character(dataset),
    pisa = list(outcome = "math", groups = c("gender", "ses")),
    # Gender x Race x Age (Age autobinned to tertiles) -> 30 well-sized strata
    # with a non-singular adjusted fit and a visible intersectional residual
    # (PCV ~72%), a far better headline MAIHDA showcase than a 2-dimension default.
    health = list(outcome = "Obese", groups = c("Gender", "Race", "Age")),
    list(outcome = NULL, groups = character(0))
  )

  outcome <- if (!is.null(known$outcome) && known$outcome %in% cols) known$outcome else NULL
  groups <- intersect(known$groups, cols)

  if (is.null(outcome)) {
    outcome <- cols[1]
  }

  # Top up the strata defaults to two from categorical-ish columns when the curated
  # names did not supply them (covers uploads and unknown datasets).
  if (length(groups) < 2) {
    candidates <- setdiff(cols, c(outcome, groups))
    is_catish <- vapply(candidates, function(cn) {
      x <- data[[cn]]
      is.factor(x) || is.character(x) || is.logical(x) ||
        (is.numeric(x) && length(unique(x[!is.na(x)])) <= 10)
    }, logical(1))
    groups <- utils::head(c(groups, candidates[is_catish]), 2)
  }

  list(outcome = outcome, groups = groups)
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

# PCV is genuinely undefined for some otherwise valid fits -- most commonly when
# the baseline (null) model has zero or negative between-stratum variance (a
# singular fit / no between-stratum variation). calculate_pvc() errors in that
# case by design, but the dashboard should still show the fitted model, VPC,
# summaries and plots rather than aborting the whole analysis. This wrapper
# returns calculate_pvc()'s result when it succeeds and, when it does not, a
# sentinel pvc_result (pvc = NA, available = FALSE) the UI can recognise. A
# bootstrap-only failure (the point PCV is fine but its CI could not be formed)
# degrades to the point estimate with a note rather than discarding the PCV.
maihda_app_calculate_pvc_safe <- function(null_model, adjusted_model,
                                          use_boot = FALSE, n_boot = 100) {
  attempt <- function(bootstrap) {
    tryCatch(
      calculate_pvc(null_model, adjusted_model, bootstrap = bootstrap, n_boot = n_boot),
      error = function(e) conditionMessage(e)
    )
  }

  boot_message <- NULL
  if (isTRUE(use_boot)) {
    res <- attempt(TRUE)
    if (inherits(res, "pvc_result")) {
      return(res)
    }
    # The bootstrap leg failed; remember why, then retry for just the point PCV so
    # a CI-only failure does not discard an otherwise valid PCV.
    boot_message <- res
  }

  res <- attempt(FALSE)
  if (inherits(res, "pvc_result")) {
    if (!is.null(boot_message)) {
      res$boot_message <- boot_message
    }
    return(res)
  }

  # The point PCV itself is undefined. Surface the variances we can still report
  # alongside an availability flag and the underlying message for the UI.
  var1 <- tryCatch(extract_between_variance(null_model), error = function(e) NA_real_)
  var2 <- tryCatch(extract_between_variance(adjusted_model), error = function(e) NA_real_)
  structure(
    list(
      pvc = NA_real_,
      var_model1 = var1,
      var_model2 = var2,
      bootstrap = FALSE,
      available = FALSE,
      message = res
    ),
    class = "pvc_result"
  )
}

# The "Compute Bootstrap CIs" control advertises uncertainty for the VPC/ICC, but
# the bootstrap is expensive. Compute the VPC/ICC intervals here, inside the
# background worker that already fits the models, so the UI session stays
# responsive; summary.maihda_model() itself is dispatched later in the main app
# session (where S3 dispatch is reliable) and these intervals are merged in via
# maihda_app_attach_vpc_ci(). Returns a list of two CI vectors (or NULL each) for
# the null and adjusted models. lme4 only -- the dashboard always fits with lme4,
# and a failed/insufficient set of refits yields NULL (no interval) rather than
# aborting the fit.
maihda_app_bootstrap_vpc_cis <- function(null_model, adjusted_model,
                                         use_boot = FALSE, n_boot = 100,
                                         conf_level = 0.95) {
  empty <- list(null = NULL, adjusted = NULL)
  if (!isTRUE(use_boot)) {
    return(empty)
  }

  args <- tryCatch(maihda_validate_bootstrap_args(n_boot, conf_level),
                   error = function(e) NULL)
  if (is.null(args)) {
    return(empty)
  }

  boot_ci <- function(model) {
    if (!identical(model$engine, "lme4")) {
      return(NULL)
    }
    tryCatch(
      bootstrap_vpc(model$model, model$data, model$formula, args$n_boot, args$conf_level),
      error = function(e) NULL
    )
  }

  list(null = boot_ci(null_model), adjusted = boot_ci(adjusted_model))
}

# Merge a worker-computed VPC/ICC bootstrap interval into a maihda_summary's vpc
# component so it carries the same fields summary.maihda_model(bootstrap = TRUE)
# would set (and maihda_vpc_has_interval()/maihda_vpc_interval_label() recognise).
# Returns the summary unchanged when no usable interval is available.
maihda_app_attach_vpc_ci <- function(summary_obj, vpc_ci, conf_level = 0.95) {
  if (is.null(summary_obj) || is.null(vpc_ci) || length(vpc_ci) < 2 ||
      !all(is.finite(vpc_ci[1:2]))) {
    return(summary_obj)
  }

  summary_obj$vpc$ci_lower <- vpc_ci[[1]]
  summary_obj$vpc$ci_upper <- vpc_ci[[2]]
  summary_obj$vpc$conf_level <- conf_level
  summary_obj$vpc$bootstrap <- TRUE
  summary_obj$vpc$method <- "bootstrap"
  summary_obj$vpc$n_boot_ok <- attr(vpc_ci, "n_ok")
  summary_obj$vpc$mc_se <- attr(vpc_ci, "mc_se")
  summary_obj
}

maihda_app_fit_models <- function(dat, outcome_var, grouping_vars,
                                  additional_covars = character(),
                                  family = "gaussian", use_boot = FALSE,
                                  n_boot = 100, autobin = TRUE,
                                  engine = "lme4", seed = NULL,
                                  decomposition = c("two-model", "cross-classified")) {
  decomposition <- match.arg(decomposition)
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

  # Reconstruct the adjusted-model main-effect terms for the stratum dimensions. A
  # numeric dimension that make_strata() auto-binned (e.g. Age tertiles) defines the
  # strata via a binned factor while the original column stays numeric, so its
  # additive main effect must be that SAME binned factor (.maihda_dim_*), not a raw
  # linear term -- mirroring core maihda()/maihda_adjusted_terms(). Covariates are not
  # stratum dimensions and enter as their raw columns (appended below).
  adj_terms <- maihda_adjusted_terms(grouping_vars, strata_dat$autobin_info, model_dat)
  model_dat <- adj_terms$data

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
  family_requested <- family
  family_autoswitched <- FALSE
  if (identical(family, "gaussian") &&
      maihda_is_binary_vector(complete_dat[[outcome_var]])) {
    message("maihda_app: outcome '", outcome_var,
            "' is binary; using family = 'binomial'. ",
            "For a linear probability model, fit from R with ",
            "fit_maihda(..., family = 'gaussian').")
    family <- "binomial"
    family_autoswitched <- TRUE
  }

  # Cross-classified decomposition: a SINGLE model that enters each stratum dimension's
  # additive main effect as a random intercept (its variance is that dimension's
  # additive contribution) plus the intersection random intercept (the interaction).
  # The additive/interaction shares come straight off this one fit. Bootstrap CIs (when
  # requested) are computed here, in the worker, via the cross-classified summary path.
  if (decomposition == "cross-classified") {
    if (length(grouping_vars) < 2) {
      stop("Cross-classified decomposition needs at least two grouping variables.",
           call. = FALSE)
    }
    cc_null_fmla <- maihda_formula_with_stratum(outcome_var, additional_covars)
    cc <- maihda_cross_classified_formula(cc_null_fmla, grouping_vars,
                                          strata_dat$autobin_info, model_dat)
    cc_model <- fit_maihda(formula = cc$formula, data = cc$data,
                           engine = engine, family = family)
    cc_model$cc_info <- list(dim_groups = cc$dim_groups,
                             interaction_group = cc$interaction_group,
                             dim_labels = grouping_vars)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    cc_summary <- summary(cc_model, bootstrap = use_boot, n_boot = n_boot)
    return(list(
      null_model = cc_model,
      model = cc_model,
      summary_obj = cc_summary,
      decomposition = cc_summary$decomposition,
      pvc = NULL,
      stepwise = NULL,
      vpc_ci_null = NULL,
      vpc_ci_adjusted = NULL,
      family_used = family,
      family_requested = family_requested,
      family_autoswitched = family_autoswitched,
      decomposition_mode = "cross-classified"
    ))
  }

  # PCV decomposition (mirrors core maihda()): the null model carries any selected
  # covariates alongside the intersectional random intercept, and the adjusted model
  # ADDS only the stratum dimensions' additive main effects. Holding the covariates in
  # BOTH models means the PCV isolates the additive share of the stratum dimensions,
  # rather than conflating covariate adjustment with the additive-vs-intersectional
  # split. With no covariates selected the null is the usual strata-only model.
  null_fmla <- maihda_formula_with_stratum(outcome_var, additional_covars)
  adjusted_fmla <- maihda_formula_with_stratum(outcome_var, c(adj_terms$terms, additional_covars))

  null_model <- fit_maihda(formula = null_fmla, data = model_dat, engine = engine, family = family)
  adjusted_model <- fit_maihda(formula = adjusted_fmla, data = model_dat, engine = engine, family = family)

  # A user-supplied seed makes the (otherwise random) bootstrap CIs reproducible
  # across runs, and lets the exported "Reproduce in R" script reproduce the same
  # PCV interval. Set it once before the bootstrap-consuming steps below; the
  # deterministic steps (strata, lme4 fits, stepwise refits) are unaffected.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # PCV degrades gracefully (a sentinel, not an error) when it is undefined; the
  # VPC/ICC bootstrap intervals are computed here in the worker and merged into
  # the main-session summaries by the app.
  pvc <- maihda_app_calculate_pvc_safe(null_model, adjusted_model, use_boot, n_boot)
  vpc_ci <- maihda_app_bootstrap_vpc_cis(null_model, adjusted_model, use_boot, n_boot)

  # Stepwise PCV: pass the raw grouping/covariate names. model_dat carries the
  # strata auto-bin recipe (strata_autobin_info), so stepwise_pcv() reconstructs an
  # auto-binned numeric dimension as the SAME tertile factor used by the adjusted
  # model (.maihda_dim_*), keeping the sequential decomposition consistent with it.
  stepwise <- stepwise_pcv(
    model_dat,
    outcome = outcome_var,
    vars = c(grouping_vars, additional_covars),
    engine = engine,
    family = family
  )

  list(
    null_model = null_model,
    model = adjusted_model,
    pvc = pvc,
    stepwise = stepwise,
    vpc_ci_null = vpc_ci$null,
    vpc_ci_adjusted = vpc_ci$adjusted,
    family_used = family,
    family_requested = family_requested,
    family_autoswitched = family_autoswitched,
    decomposition_mode = "two-model"
  )
}

# Build a runnable R script that reproduces, from the console, the analysis the
# dashboard just performed. It is a pure string builder (no Shiny, no fitting) so it
# can be unit-tested directly, and it reuses maihda_quote_name() -- the same helper
# maihda_formula_with_stratum() uses -- so the emitted formula matches what
# maihda_app_fit_models() actually fits. `family` should be the *resolved* family
# (e.g. the binomial an auto-switched binary outcome was fit with), so the script
# reproduces the real fit rather than the originally selected family.
#
# For the canonical intersectional case (>= 2 grouping variables) it delegates the
# two-model decomposition to maihda(): that puts any covariates in BOTH the null and
# adjusted models (so the PCV isolates the stratum dimensions' additive share) and
# enters an auto-binned numeric dimension as its reconstructed tertile factor --
# exactly what the dashboard fits -- without the script needing the (data-dependent)
# cut-points. With a single grouping variable there is no intersection to decompose,
# so it falls back to a single strata-only (plus covariate) fit_maihda() fit.
maihda_app_generate_code <- function(outcome_var, grouping_vars,
                                     additional_covars = character(),
                                     family = "gaussian", autobin = TRUE,
                                     use_boot = FALSE, n_boot = 100, seed = NULL,
                                     dataset = c("pisa", "health", "upload"),
                                     upload_name = NULL,
                                     decomposition = c("two-model", "cross-classified")) {
  dataset <- match.arg(dataset)
  decomposition <- match.arg(decomposition)
  additional_covars <- if (is.null(additional_covars)) character() else additional_covars

  quote_names <- function(x) vapply(x, maihda_quote_name, character(1))
  # String literals emitted into the script (column names, the upload filename) are
  # escaped with encodeString(quote = '"') so a name containing a double quote or
  # backslash produces a valid, self-contained R literal rather than breaking the
  # downloaded script or injecting code when it is sourced.
  quote_string <- function(x) encodeString(x, quote = '"')
  as_char_vec <- function(x) paste0("c(", paste0(quote_string(x), collapse = ", "), ")")
  fixed_terms <- c(grouping_vars, additional_covars)

  # The model formula's fixed part is any selected covariates (the null model); the
  # stratum dimensions' additive main effects are added by maihda() for the adjusted
  # model, so they are NOT listed here.
  covar_rhs <- if (length(additional_covars) > 0) {
    paste(quote_names(additional_covars), collapse = " + ")
  } else {
    "1"
  }
  model_fmla <- paste0(maihda_quote_name(outcome_var), " ~ ", covar_rhs, " + (1 | stratum)")

  data_line <- switch(dataset,
    pisa   = "data <- MAIHDA::maihda_country_data",
    health = "data <- MAIHDA::maihda_health_data",
    upload = sprintf('data <- read.csv(%s)  # adjust path/reader for your file',
                     quote_string(if (is.null(upload_name)) "your_data.csv" else upload_name))
  )

  model_vars <- unique(c(outcome_var, fixed_terms))

  lines <- c(
    "# Reproducible MAIHDA analysis script",
    "# Generated by MAIHDA::run_maihda_app() -- mirrors the dashboard's model fit.",
    "library(MAIHDA)",
    "",
    "# 1. Load the data",
    data_line,
    "",
    "# 2. Keep complete cases on the model variables before building strata, so the",
    "#    auto-bin cut-points match the dashboard (which fits the complete cases).",
    sprintf("model_vars <- %s", as_char_vec(model_vars)),
    "data <- data[complete.cases(data[, model_vars, drop = FALSE]), , drop = FALSE]",
    "",
    "# 3. Build intersectional strata. Keeping make_strata()'s returned data carries",
    "#    both the stratum column and the auto-bin recipe, so an auto-binned numeric",
    "#    dimension (e.g. Age tertiles) enters the adjusted model as that SAME factor.",
    sprintf("strata <- make_strata(data, vars = %s, autobin = %s)",
            as_char_vec(grouping_vars), if (isTRUE(autobin)) "TRUE" else "FALSE"),
    "data <- strata$data",
    ""
  )

  emit_stepwise <- TRUE
  if (length(grouping_vars) >= 2 && decomposition == "cross-classified") {
    # Cross-classified: a single model with each dimension's main effect as a random
    # intercept plus the intersection random intercept; the additive/interaction
    # shares are read off this one fit (no separate null/adjusted, no stepwise PCV).
    lines <- c(lines,
      "# 4. Cross-classified MAIHDA via maihda(decomposition = 'cross-classified'): a",
      "#    single model entering each stratum dimension's additive main effect as a",
      "#    random intercept plus the intersection random intercept. The additive and",
      "#    interaction shares of the between-strata variance come straight off this fit."
    )
    if (!is.null(seed)) {
      lines <- c(lines, sprintf("set.seed(%s)", seed))
    }
    maihda_call <- if (isTRUE(use_boot)) {
      sprintf('analysis <- maihda(%s, data = data, family = "%s", decomposition = "cross-classified", bootstrap = TRUE, n_boot = %s)',
              model_fmla, family, n_boot)
    } else {
      sprintf('analysis <- maihda(%s, data = data, family = "%s", decomposition = "cross-classified")',
              model_fmla, family)
    }
    lines <- c(lines,
      maihda_call,
      "analysis                              # VPC and additive/interaction shares",
      "summary(analysis$model)               # cross-classified variance components",
      "analysis$decomposition$additive_share # additive share of between-strata variance"
    )
    emit_stepwise <- FALSE
  } else if (length(grouping_vars) >= 2) {
    lines <- c(lines,
      "# 4. Two-model MAIHDA decomposition via maihda(): it fits the null model (any",
      "#    covariates + the intersectional random intercept) and the adjusted model",
      "#    (null + the stratum dimensions' additive main effects), holding the",
      "#    covariates in BOTH so the PCV isolates the dimensions' additive share."
    )
    if (!is.null(seed)) {
      lines <- c(lines, sprintf("set.seed(%s)", seed))
    }
    maihda_call <- if (isTRUE(use_boot)) {
      sprintf('analysis <- maihda(%s, data = data, family = "%s", bootstrap = TRUE, n_boot = %s)',
              model_fmla, family, n_boot)
    } else {
      sprintf('analysis <- maihda(%s, data = data, family = "%s")', model_fmla, family)
    }
    lines <- c(lines,
      maihda_call,
      "analysis                          # VPC (null) and PCV (null -> adjusted)",
      "summary(analysis$model_adjusted)  # adjusted-model variance components",
      "analysis$pcv                      # proportional change in between-stratum variance"
    )
  } else {
    lines <- c(lines,
      "# 4. A single grouping variable is an ordinary multilevel model, not an",
      "#    intersectional MAIHDA, so there is no additive-vs-intersectional",
      "#    decomposition. Fit the strata-only (plus any covariate) model and read its VPC.",
      sprintf('model <- fit_maihda(%s, data = data, family = "%s")', model_fmla, family),
      "summary(model)"
    )
  }

  if (emit_stepwise) {
    lines <- c(lines,
      "",
      "# 5. Stepwise PCV decomposition: adds each grouping dimension, then any covariates,",
      "#    one at a time (an auto-binned dimension enters as its tertile factor).",
      sprintf('stepwise_pcv(data, outcome = %s, vars = %s, family = "%s")',
              quote_string(outcome_var), as_char_vec(fixed_terms), family)
    )
  }

  paste(lines, collapse = "\n")
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
