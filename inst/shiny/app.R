library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(plotly)
library(MAIHDA)
library(future)
library(promises)

# Set up multisession for async processing
maihda_app_previous_future_plan <- future::plan()
future::plan(multisession)
shiny::onStop(function() {
  future::plan(maihda_app_previous_future_plan)
})

# Allow larger uploads than Shiny's 5 MB default so real-world CSV/DTA/SAV files
# do not silently fail before reaching reactive_data()'s reader.
options(shiny.maxRequestSize = 50 * 1024^2)

# A system-ui font stack (deliberately not font_google()): it needs no network
# call at startup, so the app never fails to launch offline, and still renders in
# a modern UI font ("Inter" when the user has it installed, otherwise the OS UI
# font). The palette keeps the existing brand navy as primary and reuses the
# plotly bar colour (#4D9DE0) as info for cross-component consistency.
maihda_fonts <- c("Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto",
                  "Helvetica Neue", "Arial", "sans-serif")
maihda_theme <- bs_theme(
  version = 5,
  primary = "#2C3E50", secondary = "#5D6D7E",
  success = "#2E7D5B", info = "#4D9DE0", warning = "#E0A458",
  base_font = maihda_fonts, heading_font = maihda_fonts,
  "border-radius" = "0.6rem"
)

ui <- page_navbar(
  id = "main_tabs",
  window_title = "MAIHDA Analysis Dashboard",
  title = tags$span(class = "maihda-brand", icon("layer-group"), "MAIHDA Dashboard"),
  theme = maihda_theme,
  fillable = FALSE,
  header = tagList(
    shinyjs::useShinyjs(),
    tags$head(tags$link(rel = "stylesheet", href = "custom.css"))
  ),

  sidebar = sidebar(
    title = "Controls",
    width = 330,
    accordion(
      id = "sidebar_steps",
      open = c("step-data", "step-model"),
      multiple = TRUE,
      accordion_panel(
        "1 · Data", value = "step-data", icon = icon("table"),
        selectInput("dataset", "Select Dataset:",
                    # NHANES first so it is the default: its Gender x Race x Age
                    # MAIHDA is a clean, non-singular showcase. PISA shines in the
                    # cross-country group-comparison tab rather than this headline.
                    choices = c("Built-in: NHANES Health Data" = "health",
                                "Built-in: PISA Country Data" = "pisa",
                                "Upload Custom Data" = "upload")),
        conditionalPanel(
          condition = "input.dataset == 'upload'",
          fileInput("upload", "Upload Data (CSV/DTA/SAV)", accept = c(".csv", ".dta", ".sav"))
        )
      ),
      accordion_panel(
        "2 · Model", value = "step-model", icon = icon("sitemap"),
        selectizeInput("outcome",
                       tagList("Outcome Variable ",
                               tooltip(icon("info-circle"),
                                       "The variable whose inequality you're analysing. A two-level outcome is fitted as logistic automatically.")),
                       choices = NULL),
        selectizeInput("group_vars",
                       tagList("Strata Grouping Variables ",
                               tooltip(icon("info-circle"),
                                       "Pick 2 or more variables (e.g. gender and race). Their combinations define the intersectional strata.")),
                       choices = NULL, multiple = TRUE),
        uiOutput("group_var_hint"),
        checkboxInput("autobin", "Auto-bin continuous strata vars (>10 unique values) into 3 groups", value = TRUE),
        selectizeInput("covariates", "Additional Covariates (Fixed Effects)", choices = NULL, multiple = TRUE)
      ),
      accordion_panel(
        "3 · Options", value = "step-options", icon = icon("sliders"),
        selectInput("family", "Family", choices = c("gaussian", "binomial", "poisson"), selected = "gaussian"),
        radioButtons("decomposition",
                     tagList("Additive vs. interaction decomposition ",
                             tooltip(icon("info-circle"),
                                     "Two-model: fit a null and an adjusted model and read the additive share from the PCV. Cross-classified: a single model entering each dimension's additive main effect as a random intercept and the intersection as the interaction, with the additive/interaction shares read off that one fit.")),
                     choices = c("Two-model (null vs adjusted, PCV)" = "two-model",
                                 "Cross-classified (single model)" = "cross-classified"),
                     selected = "two-model"),
        checkboxInput("use_boot", "Compute Bootstrap CIs (Slower)", value = FALSE),
        conditionalPanel(
          condition = "input.use_boot == true",
          numericInput("n_boot", "Bootstrap Samples", value = 100, min = 10, step = 10),
          numericInput("seed", "Random Seed (reproducible bootstrap)", value = 123, min = 1, step = 1)
        )
      )
    ),
    # Primary action, pinned to the bottom of the sidebar so it stays visible
    # regardless of how far the controls scroll. input_task_button shows a busy
    # state and disables itself while a fit is in flight (driven from the server).
    div(
      class = "maihda-fit-wrap",
      input_task_button("fit_btn", "Fit MAIHDA Model", icon = icon("play"),
                        label_busy = "Fitting…", auto_reset = FALSE)
    )
  ),

  # ---- Tabs (results-first), then a spacer and right-aligned utilities -------
  nav_panel("Overview", value = "overview", icon = icon("house"),
            uiOutput("overview_ui")),
  nav_panel("Model Summary", value = "model_summary", icon = icon("clipboard-list"),
            uiOutput("model_summary_ui")),
  nav_panel("PCV Results", value = "pcv", icon = icon("arrow-down-wide-short"),
            navset_pill(
              nav_panel("PCV summary", uiOutput("pvc_summary_ui")),
              nav_panel("Stepwise decomposition", uiOutput("stepwise_pcv_ui"))
            )),
  nav_panel("Model Comparison", value = "compare", icon = icon("code-compare"),
            MAIHDA:::mod_compare_ui("compare")),
  nav_panel("Visualizations", value = "viz", icon = icon("chart-column"),
            MAIHDA:::mod_visualizations_ui("viz")),
  nav_panel("Interactive Explorer", value = "explorer", icon = icon("compass"),
            MAIHDA:::mod_explorer_ui("explorer")),
  nav_panel("Data View", value = "data", icon = icon("table"),
            shinycssloaders::withSpinner(DTOutput("data_table"))),

  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")),
  nav_item(actionLink("help", tagList(icon("circle-question"), "Help"),
                      `aria-label` = "Open help and glossary")),
  nav_item(actionLink("show_code", tagList(icon("code"), "Reproduce in R"),
                      `aria-label` = "Show the reproducible R script for this analysis")),
  nav_item(bookmarkButton(label = "Bookmark",
                          title = "Save this analysis configuration to a shareable URL (built-in data + selections)."))
)

server <- function(input, output, session) {

  # URL bookmarking captures the analysis *configuration* (built-in dataset +
  # variable/family/seed selections), not the fitted results. Exclude the file
  # upload (an uploaded file cannot be restored from a URL) and the action controls
  # (so a restored bookmark never auto-triggers a fit). Module run-buttons are
  # excluded inside their own modules.
  setBookmarkExclude(c("upload", "fit_btn", "help", "show_code",
                       "overview_start", "overview_help"))

  # Load data: built-in PISA / NHANES datasets, or an uploaded file
  reactive_data <- reactive({
    if (input$dataset == "pisa") {
      # PISA: math/reading scores by gender x SES across countries
      return(MAIHDA::maihda_country_data)
    } else if (input$dataset == "health") {
      # Real-world NHANES health dataset
      return(MAIHDA::maihda_health_data)
    } else if (input$dataset == "upload" && !is.null(input$upload)) {
      ext <- tolower(tools::file_ext(input$upload$name))
      dat <- tryCatch({
        raw <- if (ext == "csv") {
          read.csv(input$upload$datapath)
        } else if (ext == "dta") {
          if (!requireNamespace('haven', quietly = TRUE)) stop("haven package required for DTA files")
          haven::as_factor(haven::read_dta(input$upload$datapath))
        } else if (ext == "sav") {
          if (!requireNamespace('haven', quietly = TRUE)) stop("haven package required for SAV files")
          haven::as_factor(haven::read_sav(input$upload$datapath))
        } else {
          stop("Unsupported format")
        }
        raw <- as.data.frame(raw)
        if (nrow(raw) == 0 || ncol(raw) == 0) {
          stop("The uploaded file has no usable rows or columns.")
        }
        raw
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
        NULL
      })
      return(dat)
    } else {
      # Fallback while waiting for upload
      return(NULL)
    }
  })

  observe({
    req(reactive_data())
    cols <- names(reactive_data())

    # Preserve selections if still valid
    curr_outcome <- isolate(input$outcome)
    curr_group <- isolate(input$group_vars)
    curr_covars <- isolate(input$covariates)

    # Dataset-aware defaults (fixes the NHANES dataset, whose columns are
    # capitalised); preserve the user's current selections when still valid.
    spec <- MAIHDA:::maihda_app_default_vars(input$dataset, reactive_data())
    new_outcome <- if (!is.null(curr_outcome) && curr_outcome %in% cols) curr_outcome else spec$outcome
    new_group <- if (!is.null(curr_group) && all(curr_group %in% cols)) curr_group else spec$groups

    # Calculate available covariates by excluding outcome and strata variables
    used_vars <- c(new_outcome, new_group)
    avail_covars <- setdiff(cols, used_vars)

    # Filter current covariates that might have been pushed out
    new_covars <- intersect(curr_covars, avail_covars)

    updateSelectizeInput(session, "outcome", choices = cols, selected = new_outcome, server = TRUE)
    updateSelectizeInput(session, "group_vars", choices = cols, selected = new_group, server = TRUE)
    updateSelectizeInput(session, "covariates", choices = avail_covars, selected = new_covars, server = TRUE)
  })

  observe({
    # Grey out fit button if no grouping vars are selected
    shinyjs::toggleState("fit_btn", condition = length(input$group_vars) > 0)
  })

  # Non-blocking guard: a single grouping variable is an ordinary multilevel model,
  # not an intersectional MAIHDA. Warn, but allow it (a user may want it).
  output$group_var_hint <- renderUI({
    if (length(input$group_vars) == 1) {
      div(class = "text-warning small mt-1",
          icon("triangle-exclamation"),
          " Only one grouping variable: this is an ordinary multilevel model, not an intersectional MAIHDA. Select 2+ variables to form intersectional strata.")
    } else {
      NULL
    }
  })

  # The Help / Glossary modal is reachable from both the navbar link and the
  # Overview tab's call-to-action, so build it once and show it from either.
  maihda_help_modal <- function() {
    modalDialog(
      title = "MAIHDA Dashboard -- Help & Glossary",
      easyClose = TRUE,
      size = "l",
      markdown(
        "### How to use this dashboard
1. **Select a dataset** -- a built-in example or your own CSV/DTA/SAV upload.
2. **Choose an outcome** and **two or more grouping variables**; their
   combinations form the intersectional *strata*.
3. *(Optional)* add covariates, pick a family, enable bootstrap CIs.
4. Click **Fit MAIHDA Model** -- results populate the tabs across the top.

### Glossary
- **Strata** -- the intersectional groups formed by combining the grouping
  variables (e.g. *Female x Black x Low-education*).
- **VPC / ICC** -- Variance Partition Coefficient: the share of outcome variation
  that lies *between* strata. For binary/count outcomes it is on the model's
  latent scale.
- **PCV** -- Proportional Change in Variance: how much the between-stratum variance
  shrinks once the additive main effects are added (a model-dependent comparison,
  not proof that inequality was 'explained away').
- **MAIHDA** -- Multilevel Analysis of Individual Heterogeneity and Discriminatory
  Accuracy.
- **Discriminatory accuracy** (binary outcomes):
    - **AUC / C-statistic** -- how well stratum membership predicts the individual
      outcome (0.5 = chance). A high VPC can still go with only modest AUC.
    - **MOR** (Median Odds Ratio) -- the typical change in odds between a higher-
      and a lower-risk stratum; 1 means no between-stratum heterogeneity.

*See the package vignettes for worked examples and the statistical details.*"
      ),
      footer = modalButton("Close")
    )
  }
  observeEvent(input$help, showModal(maihda_help_modal()))
  observeEvent(input$overview_help, showModal(maihda_help_modal()))

  # Welcome / Overview landing (the default tab on load).
  output$overview_ui <- renderUI({
    tagList(
      div(
        class = "maihda-hero mb-4",
        h1("MAIHDA Analysis Dashboard"),
        p(class = "maihda-hero-lead",
          "A no-code interface for Multilevel Analysis of Individual Heterogeneity
           and Discriminatory Accuracy -- build intersectional strata, fit the null
           and adjusted multilevel models, and read off the variance partition
           coefficient and its proportional change."),
        layout_columns(
          col_widths = c(4, 4, 4),
          class = "maihda-steps mt-3",
          div(class = "maihda-step",
              icon("table", class = "maihda-step-icon"),
              h6("1 · Choose data"),
              p("Pick a built-in PISA or NHANES example, or upload your own CSV/DTA/SAV.")),
          div(class = "maihda-step",
              icon("sitemap", class = "maihda-step-icon"),
              h6("2 · Define strata"),
              p("Select an outcome and 2+ grouping variables; their combinations form the intersectional strata.")),
          div(class = "maihda-step",
              icon("play", class = "maihda-step-icon"),
              h6("3 · Fit & explore"),
              p("Click Fit MAIHDA Model -- results populate the tabs above."))
        ),
        div(class = "mt-4 d-flex gap-2 flex-wrap",
            actionButton("overview_start", tagList(icon("arrow-right"), "Get started"),
                         class = "btn-light"),
            actionButton("overview_help", tagList(icon("circle-question"), "What do VPC & PCV mean?"),
                         class = "btn-outline-light"))
      ),
      card(
        class = "maihda-overview-note",
        card_body(
          markdown(
            "**MAIHDA** partitions outcome variation into a *between-strata* share
            (the **VPC / ICC**) and shows how that share changes once additive main
            effects are added (the **PCV**). A high VPC need not imply high
            individual-level *discriminatory accuracy* -- the cautionary core of the
            method. Open **Help** in the top bar for the full glossary.")
        )
      )
    )
  })

  # "Get started" jumps to the data view, where the sidebar workflow is in reach.
  observeEvent(input$overview_start, {
    nav_select("main_tabs", "data")
  })

  output$data_table <- renderDT({
    datatable(reactive_data(), options = list(pageLength = 10, scrollX = TRUE))
  })

  # Reactive values for model and results
  model_results <- reactiveVal(NULL)
  null_summary_results <- reactiveVal(NULL)
  summary_results <- reactiveVal(NULL)
  pvc_results <- reactiveVal(NULL)
  stepwise_results <- reactiveVal(NULL)
  fitted_family <- reactiveVal(NULL)   # resolved family of the last fit (for display)
  fit_params <- reactiveVal(NULL)      # inputs of the last fit (for the "Reproduce in R" dialog)
  da_results <- reactiveVal(NULL)      # discriminatory accuracy (binomial only): null + adjusted
  comparison_results <- reactiveVal(NULL)  # nested-model VPC comparison (null vs adjusted)
  decomposition_results <- reactiveVal(NULL)  # cross-classified additive/interaction partition
  decomposition_mode <- reactiveVal("two-model")  # "two-model" or "cross-classified"

  # Monotonic request token: each fit increments it, so a slower superseded future
  # can recognise it is stale and discard its result rather than overwriting a
  # newer fit. Pairs with disabling the button while a fit is in flight.
  fit_id <- reactiveVal(0)

  observeEvent(input$fit_btn, {
    dat <- reactive_data()
    req(dat)

    grouping_vars <- input$group_vars
    req(length(grouping_vars) > 0)

    additional_covars <- input$covariates
    outcome_var <- input$outcome
    eng <- "lme4"
    fam <- input$family

    use_boot <- input$use_boot
    n_boot <- input$n_boot
    autobin_opt <- input$autobin
    decomp_opt <- input$decomposition
    seed_opt <- if (isTRUE(use_boot)) input$seed else NULL

    # Claim a request token and lock the button so a second click cannot launch a
    # competing fit whose (possibly out-of-order) result would clobber this one.
    this_fit <- fit_id() + 1
    fit_id(this_fit)
    bslib::update_task_button("fit_btn", state = "busy")

    # Remember exactly what this fit used so the "Reproduce in R" dialog mirrors it,
    # independent of any later sidebar edits.
    fit_params(list(
      dataset = input$dataset,
      upload_name = if (!is.null(input$upload)) input$upload$name else NULL,
      outcome = outcome_var,
      grouping_vars = grouping_vars,
      covariates = additional_covars,
      autobin = autobin_opt,
      use_boot = use_boot,
      n_boot = n_boot,
      seed = seed_opt,
      decomposition = decomp_opt
    ))

    # Reset old results
    model_results(NULL)
    null_summary_results(NULL)
    summary_results(NULL)
    pvc_results(NULL)
    stepwise_results(NULL)
    fitted_family(NULL)
    da_results(NULL)
    comparison_results(NULL)
    decomposition_results(NULL)

    id <- showNotification("Creating strata & Fitting Models (May take a moment)...", duration = NULL, type = "message")

    future_promise({
        MAIHDA:::maihda_app_fit_models(
          dat = dat,
          outcome_var = outcome_var,
          grouping_vars = grouping_vars,
          additional_covars = additional_covars,
          family = fam,
          use_boot = use_boot,
          n_boot = n_boot,
          autobin = autobin_opt,
          engine = eng,
          seed = seed_opt,
          decomposition = decomp_opt
        )
      }, seed = TRUE) %...>% (function(res) {
        removeNotification(id)
        # A newer fit has superseded this one: drop the stale result and leave the
        # button for the newer fit to re-enable when it finishes.
        if (!identical(this_fit, fit_id())) return(invisible(NULL))
        bslib::update_task_button("fit_btn", state = "ready")

        model_results(res$model)
        fitted_family(res$family_used)
        decomposition_mode(res$decomposition_mode)

        if (identical(res$decomposition_mode, "cross-classified")) {
          # Single cross-classified model: the additive/interaction partition (with
          # bootstrap share CIs when requested) was computed in the worker. There is no
          # separate null/adjusted pair, PCV, stepwise or nested comparison.
          summary_results(res$summary_obj)
          null_summary_results(res$summary_obj)
          pvc_results(NULL)
          stepwise_results(NULL)
          decomposition_results(res$decomposition)
          comparison_results(NULL)
          if (identical(res$family_used, "binomial")) {
            da_results(list(
              null = tryCatch(MAIHDA::maihda_discriminatory_accuracy(res$model),
                              error = function(e) NULL),
              adjusted = NULL
            ))
          }
        } else {
          # Dispatch summary() in the main thread, where the MAIHDA S3 method is
          # reliably found, then merge in the VPC/ICC bootstrap intervals computed
          # in the background worker (see maihda_app_bootstrap_vpc_cis()).
          null_summary_results(MAIHDA:::maihda_app_attach_vpc_ci(summary(res$null_model), res$vpc_ci_null))
          summary_results(MAIHDA:::maihda_app_attach_vpc_ci(summary(res$model), res$vpc_ci_adjusted))
          pvc_results(res$pvc)
          stepwise_results(res$stepwise)
          decomposition_results(NULL)
          # Discriminatory accuracy is only defined for binomial fits. Compute AUC/MOR
          # for the strata-only (null) and adjusted models in the main thread (fast --
          # predict + rank, no refit).
          if (identical(res$family_used, "binomial")) {
            da_results(list(
              null = tryCatch(MAIHDA::maihda_discriminatory_accuracy(res$null_model),
                              error = function(e) NULL),
              adjusted = tryCatch(MAIHDA::maihda_discriminatory_accuracy(res$model),
                                  error = function(e) NULL)
            ))
          }
          # Nested-model VPC comparison (null vs adjusted) for the Model Comparison
          # tab -- pure (reads VPCs from the already-fitted models), so main-thread.
          comparison_results(tryCatch(
            MAIHDA::compare_maihda(res$null_model, res$model,
                                   model_names = c("Model 1: Null", "Model 2: Adjusted")),
            error = function(e) NULL
          ))
        }
        if (isTRUE(res$family_autoswitched)) {
          showNotification(
            sprintf("Outcome '%s' is binary -- fitted as 'binomial' (not the selected 'gaussian').",
                    outcome_var),
            type = "warning", duration = 12
          )
        }
        nav_select("main_tabs", "model_summary")
      }) %...!% (function(err) {
        removeNotification(id)
        if (!identical(this_fit, fit_id())) return(invisible(NULL))
        bslib::update_task_button("fit_btn", state = "ready")
        showNotification(paste("Error fitting model:", err$message), type = "error", duration = 15)
      })
  })

  output$model_summary_ui <- renderUI({
    if (is.null(summary_results())) {
      return(MAIHDA:::maihda_app_empty_state(
        "No model fitted yet",
        "**1.** Pick a dataset, **2.** choose an outcome and **2+ grouping
        variables**, then **3.** click *Fit MAIHDA Model* in the sidebar.
        Need definitions? Open **Help** in the top bar."))
    }
    res <- summary_results()

    vpc <- res$vpc
    vpc_interval <- if (MAIHDA:::maihda_vpc_has_interval(vpc)) {
      div(class = "small",
          sprintf("[%.2f%%, %.2f%%] %s",
                  vpc$ci_lower * 100, vpc$ci_upper * 100,
                  MAIHDA:::maihda_vpc_interval_label(vpc)))
    } else {
      NULL
    }

    family_line <- if (!is.null(fitted_family())) {
      div(class = "small", sprintf("Fitted with family = '%s'", fitted_family()))
    } else {
      NULL
    }

    # For binomial fits the headline VPC is latent-scale; also show the response
    # (probability) scale VPC as an interpretable complement. Seeded so the
    # simulation-based value is stable across re-renders. Skipped in cross-classified
    # mode, where the simulation helper reads only the interaction variance (not the
    # full additive+interaction between-strata variance) and would understate it.
    response_vpc_line <- if (identical(fitted_family(), "binomial") &&
                             !identical(decomposition_mode(), "cross-classified")) {
      seed_val <- if (!is.null(fit_params()) && !is.null(fit_params()$seed)) fit_params()$seed else 1L
      rv <- tryCatch(MAIHDA::maihda_vpc_response(model_results(), seed = seed_val),
                     error = function(e) NULL)
      if (!is.null(rv) && is.finite(rv$estimate)) {
        div(class = "small",
            sprintf("Response-scale VPC: %.2f%% (simulation method; latent-scale shown above)",
                    rv$estimate * 100))
      } else {
        NULL
      }
    } else {
      NULL
    }

    # Surface the fit-quality diagnostics fit_maihda() already computes (singular
    # fit / non-convergence): these silently invalidate the VPC/PCV if ignored.
    diag_lines <- MAIHDA:::maihda_format_fit_diagnostics(model_results()$diagnostics)
    diag_ui <- if (length(diag_lines) > 0) {
      div(class = "alert alert-warning",
          tags$strong("Fit diagnostics"),
          tags$ul(lapply(diag_lines, function(l) tags$li(l))))
    } else {
      div(class = "text-success small",
          icon("check-circle"), " Model converged with no singularity warnings.")
    }

    # Strata overview + small-cell warning, from the per-stratum sample sizes in
    # strata_info. Small cells make the random-effect estimates unstable.
    si <- model_results()$strata_info
    strata_ui <- if (!is.null(si) && "n" %in% names(si)) {
      small_thresh <- 10
      n_strata <- nrow(si)
      n_small <- sum(si$n < small_thresh, na.rm = TRUE)
      tagList(
        div(class = "small text-muted mt-1",
            sprintf("%d strata; sizes range %d-%d (median %d).",
                    n_strata, min(si$n, na.rm = TRUE), max(si$n, na.rm = TRUE),
                    round(stats::median(si$n, na.rm = TRUE)))),
        if (n_small > 0) {
          div(class = "alert alert-warning mt-1",
              sprintf("%d of %d strata have fewer than %d individuals. Random-effect estimates for small strata are unstable -- interpret their deviations cautiously.",
                      n_small, n_strata, small_thresh))
        } else NULL
      )
    } else {
      NULL
    }

    # Discriminatory Accuracy card (binomial only): AUC of the strata-only vs the
    # adjusted model, plus the Median Odds Ratio.
    fmt_metric <- function(x, digits = 3) {
      if (!is.null(x) && is.finite(x)) formatC(x, format = "f", digits = digits) else "NA"
    }
    da <- da_results()
    da_card <- if (!is.null(da) && (!is.null(da$null) || !is.null(da$adjusted))) {
      card(
        card_header("Discriminatory Accuracy (binary outcome)"),
        layout_columns(
          col_widths = c(4, 4, 4),
          class = "maihda-metric-row",
          value_box(
            title = "AUC — strata only",
            value = fmt_metric(if (!is.null(da$null)) da$null$auc else NA),
            showcase = icon("bullseye"), theme = "info",
            p(class = "mb-0", "C-statistic of the intersectional strata alone")
          ),
          value_box(
            title = "AUC — adjusted",
            value = fmt_metric(if (!is.null(da$adjusted)) da$adjusted$auc else NA),
            showcase = icon("bullseye"), theme = "primary",
            p(class = "mb-0", "With individual covariates added")
          ),
          value_box(
            title = "Median Odds Ratio",
            value = fmt_metric(if (!is.null(da$null)) da$null$mor else NA, 2),
            showcase = icon("scale-balanced"), theme = "secondary",
            p(class = "mb-0", "Between-stratum heterogeneity on the odds-ratio scale")
          )
        ),
        div(class = "small text-muted mt-2",
            "AUC = 0.5 is chance. A high between-stratum VPC can still translate into only modest individual-level discriminatory accuracy -- the cautionary message at the heart of the 'DA' in MAIHDA.")
      )
    } else {
      NULL
    }

    tagList(
      layout_columns(
        col_widths = c(5, 7),
        class = "maihda-metric-row",
        value_box(
          title = tagList("Variance Partition Coefficient (VPC) / ICC ",
                          tooltip(icon("info-circle"),
                                  "Share of outcome variation that lies between strata. For binary/count outcomes this is on the model's latent scale.")),
          value = sprintf("%.2f%%", vpc$estimate * 100),
          showcase = icon("layer-group"),
          theme = "primary",
          vpc_interval,
          family_line,
          response_vpc_line
        ),
        card(
          card_header("Fit diagnostics & strata overview"),
          diag_ui,
          strata_ui
        )
      ),
      da_card,
      layout_columns(
        card(
          card_header("Variance Components"),
          DTOutput("dt_var_comp")
        ),
        card(
          card_header("Fixed Effects"),
          DTOutput("dt_fix_eff")
        )
      ),
      card(
        card_header("Stratum Estimates (top 10)"),
        DTOutput("dt_stratum")
      )
    )
  })

  output$dt_var_comp <- renderDT({
    req(summary_results())
    datatable(summary_results()$variance_components, options = list(dom = 't', paging = FALSE))
  })

  output$dt_fix_eff <- renderDT({
    req(summary_results())
    datatable(as.data.frame(summary_results()$fixed_effects), options = list(dom = 't', paging = FALSE))
  })

  output$dt_stratum <- renderDT({
    req(summary_results())
    datatable(head(summary_results()$stratum_estimates, 10), options = list(dom = 't', paging = FALSE))
  })

  output$pvc_summary_ui <- renderUI({
    # Cross-classified mode: show the additive/interaction partition read off the
    # single model instead of the two-model PCV.
    if (identical(decomposition_mode(), "cross-classified")) {
      if (is.null(decomposition_results()) || is.null(summary_results())) {
        return(MAIHDA:::maihda_app_empty_state(
          "No decomposition yet",
          "Fit a cross-classified MAIHDA model from the sidebar to see the additive
          (per-dimension) vs. intersectional-interaction split of the between-strata
          variance."))
      }
      d <- decomposition_results()
      fmt_share <- function(est, ci) {
        base <- sprintf("%.2f%%", est * 100)
        if (!is.null(ci) && length(ci) == 2 && all(is.finite(ci))) {
          paste0(base, sprintf(" [%.1f%%, %.1f%%]", ci[1] * 100, ci[2] * 100))
        } else {
          base
        }
      }
      per_dim_df <- data.frame(
        Dimension = names(d$per_dim),
        `Additive variance` = as.numeric(d$per_dim),
        `Share of between-strata` = sprintf("%.1f%%",
                                            as.numeric(d$per_dim) / d$between_var * 100),
        check.names = FALSE
      )
      return(card(
        card_header("Additive vs. Intersectional Decomposition (cross-classified)"),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            class = "maihda-metric-row mb-2",
            value_box(
              title = "Additive variance",
              value = sprintf("%.4f", d$additive_var),
              showcase = icon("layer-group"), theme = "secondary",
              p(class = "mb-0", "Sum of the dimensions' main-effect random variances")
            ),
            value_box(
              title = "Interaction variance",
              value = sprintf("%.4f", d$interaction_var),
              showcase = icon("diagram-project"), theme = "info",
              p(class = "mb-0", "Intersection random effect (interaction beyond additive)")
            ),
            value_box(
              title = tagList("Additive share ",
                              tooltip(shiny::icon("info-circle"),
                                      "The additive (dimension main-effect) variance as a fraction of the total between-strata variance. This is the cross-classified analogue of the PCV, but a different (partial-pooling) estimator; the complement is the intersectional interaction share.")),
              value = fmt_share(d$additive_share, d$additive_share_ci),
              showcase = icon("arrow-down-wide-short"),
              theme = "success",
              p(class = "mb-0",
                sprintf("Interaction share: %s", fmt_share(d$interaction_share, d$interaction_share_ci)))
            )
          ),
          div(class = "small text-muted mb-3",
              "Total between-strata variance: ",
              tags$strong(sprintf("%.4f", d$between_var)),
              HTML("&nbsp;&bull;&nbsp;"),
              "Model: ",
              tags$code(paste(deparse(model_results()$formula), collapse = ""))),
          h6("Per-dimension additive variance"),
          DT::datatable(per_dim_df, options = list(dom = "t", paging = FALSE),
                        rownames = FALSE),
          div(class = "small text-muted mt-3",
              "Dimensions with few levels (e.g. a binary variable) are poorly identified ",
              "and often give a singular fit -- their additive variance can be unstable. ",
              "See the Model Summary tab for fit diagnostics.")
        )
      ))
    }

    if (is.null(pvc_results()) || is.null(model_results())) {
      return(MAIHDA:::maihda_app_empty_state(
        "No PCV results yet",
        "Fit a MAIHDA model from the sidebar to see the proportional change in
        between-stratum variance between the null and adjusted models."))
    }

    pvc <- pvc_results()
    mod <- model_results()

    adjusted_formula <- deparse(mod$formula)
    outcome_var <- all.vars(mod$formula)[1]
    # The null model carries any selected covariates alongside the random intercept
    # (the adjusted model adds the stratum dimensions' main effects), so the displayed
    # null formula must include those covariates to match the fitted model.
    null_covars <- if (!is.null(fit_params())) fit_params()$covariates else character()
    null_formula <- paste(
      deparse(MAIHDA:::maihda_formula_with_stratum(
        outcome_var, if (is.null(null_covars)) character() else null_covars)),
      collapse = ""
    )

    bootstrap_ui <- if (isTRUE(pvc$bootstrap) && !is.null(pvc$ci_lower) && !is.null(pvc$ci_upper)) {
        div(class = "mt-4 text-center text-muted",
            h5("Bootstrap 95% Confidence Interval"),
            tags$p(sprintf("[%.2f%%, %.2f%%]", pvc$ci_lower * 100, pvc$ci_upper * 100))
        )
    } else if (!is.null(pvc$boot_message)) {
        div(class = "mt-4 text-center text-muted",
            tags$p(sprintf("Bootstrap CI unavailable: %s", pvc$boot_message))
        )
    } else {
        NULL
    }

    card(
      card_header("Proportional Change in Variance (PCV)"),
      card_body(
        layout_columns(
          col_widths = c(4, 4, 4),
          class = "maihda-metric-row mb-2",
          value_box(
            title = "Null model variance",
            value = if (!is.null(pvc$var_model1) && is.finite(pvc$var_model1)) sprintf("%.4f", pvc$var_model1) else "N/A",
            showcase = icon("seedling"), theme = "secondary",
            p(class = "mb-0", "Model 1: null (strata + any covariates)")
          ),
          value_box(
            title = "Adjusted model variance",
            value = if (!is.null(pvc$var_model2) && is.finite(pvc$var_model2)) sprintf("%.4f", pvc$var_model2) else "N/A",
            showcase = icon("layer-group"), theme = "info",
            p(class = "mb-0", "Model 2: + strata main effects")
          ),
          value_box(
            title = tagList("Estimated PCV ",
                            tooltip(
                              shiny::icon("info-circle"),
                              "PCV is the proportional change in between-stratum variance from the Null to the Adjusted model. The two models hold any selected covariates fixed and differ only by the strata dimensions' additive main effects, so the PCV is the additive share of those dimensions. A high PCV means the between-stratum variance is much smaller after adding those main effects; a low or negative PCV means little change (or an increase). This is a model-dependent change, not proof that inequality was causally 'explained away' -- it can also reflect suppression, rescaling, sample composition, or uncertainty, not interaction alone.")),
            value = if (is.finite(pvc$pvc)) sprintf("%.2f%%", pvc$pvc * 100) else "N/A",
            showcase = icon("arrow-down-wide-short"),
            theme = if (is.finite(pvc$pvc)) "success" else "warning",
            p(class = "mb-0",
              if (is.finite(pvc$pvc)) "Proportional change vs the null model" else "Undefined for this fit")
          )
        ),
        div(class = "small text-muted mb-3",
            tagList("Null: ", tags$code(null_formula),
                    HTML("&nbsp;&bull;&nbsp;"),
                    "Adjusted: ", tags$code(paste(adjusted_formula, collapse = "")))),
        if (!is.finite(pvc$pvc)) {
          div(class = "alert alert-warning",
              tags$strong("PCV could not be calculated. "),
              if (!is.null(pvc$message)) pvc$message else
                "The baseline between-stratum variance is zero, so the proportional change is undefined. The model fit, VPC and visualizations above remain valid.")
        } else NULL,
        bootstrap_ui
      )
    )
  })

  output$stepwise_pcv_ui <- renderUI({
    if (identical(decomposition_mode(), "cross-classified")) {
      return(MAIHDA:::maihda_app_empty_state(
        "Stepwise PCV not used here",
        "The stepwise PCV decomposition belongs to the **two-model** workflow. In
        cross-classified mode the additive and interaction shares are read directly
        from the single model -- see the **PCV summary** sub-tab."))
    }
    req(stepwise_results())

    card(
      card_header("Stepwise Proportional Change in Variance Decomposition"),
      card_body(
        markdown("
        This table shows the proportional change in between-stratum variance as covariates are added step by step. It is a model-dependent, order-dependent comparison, not a causal decomposition of inequality.

        *   **Step_PCV**: Proportional change in between-stratum variance relative to the *previous* model step.
        *   **Total_PCV**: Proportional change relative to the *null* model (Step 0).
        "),
        shinycssloaders::withSpinner(plotlyOutput("stepwise_pcv_plot", height = "400px")),
        hr(),
        DTOutput("stepwise_pcv_dt")
      )
    )
  })

  output$stepwise_pcv_plot <- renderPlotly({
    req(stepwise_results())
    res <- stepwise_results()

    # Ensure Model column is an ordered factor to maintain step sequence
    res$Model <- factor(res$Model, levels = res$Model)

    # Calculate step variance drop
    res$Step_Variance <- c(0, -diff(res$Variance))

    # Format tooltip text
    hover_text <- paste(
      "<b>Model:</b>", res$Model, "<br>",
      "<b>Added:</b>", ifelse(is.na(res$Added_Variable), "None", res$Added_Variable), "<br>",
      "<b>Step Variance Drop:</b>", round(res$Step_Variance, 4), "<br>",
      "<b>Step PCV:</b>", ifelse(!is.na(res$Step_PCV), paste0(round(res$Step_PCV * 100, 2), "%"), "0%"), "<br>",
      "<b>Total PCV:</b>", ifelse(!is.na(res$Total_PCV), paste0(round(res$Total_PCV * 100, 2), "%"), "0%")
    )

    plot_ly(
      data = res,
      x = ~Model,
      y = ~Total_PCV,
      type = "bar",
      text = hover_text,
      hoverinfo = "text",
      marker = list(color = "#4D9DE0")
    ) |>
      layout(
        title = "Cumulative Change in Between-Stratum Variance",
        xaxis = list(title = "Sequential Model Step", tickangle = -45),
        yaxis = list(title = "Total PCV (Proportional Change in Variance)", tickformat = ".1%")
      )
  })

  output$stepwise_pcv_dt <- renderDT({
    req(stepwise_results())
    res <- stepwise_results()

    # Format the table for the viewer
    df <- res
    df$Variance <- sprintf("%.4f", df$Variance)
    df$Step_PCV <- ifelse(!is.na(df$Step_PCV), sprintf("%.2f%%", df$Step_PCV * 100), "0.00%")
    df$Total_PCV <- ifelse(!is.na(df$Total_PCV), sprintf("%.2f%%", df$Total_PCV * 100), "0.00%")

    datatable(df, options = list(dom = 't', paging = FALSE, ordering = FALSE), rownames = FALSE, escape = FALSE)
  })

  # Visualizations tab: plot picker, ggplot/plotly swap and PNG download.
  MAIHDA:::mod_visualizations_server("viz", model_results = model_results)

  # Interactive Explorer (HUD) tab: key metrics, filterable strata-deviation
  # plot and filtered data export.
  MAIHDA:::mod_explorer_server(
    "explorer",
    model_results = model_results,
    null_summary_results = null_summary_results,
    summary_results = summary_results,
    pvc_results = pvc_results,
    group_vars = reactive(input$group_vars),
    decomposition_results = decomposition_results
  )

  # Model Comparison tab: nested null-vs-adjusted VPC + stratified-by-group MAIHDA.
  MAIHDA:::mod_compare_server(
    "compare",
    comparison_results = comparison_results,
    reactive_data = reactive_data,
    fit_params = fit_params,
    fitted_family = fitted_family
  )

  # --- Reproduce in R: a console script mirroring the last fit ----------------
  # Accessed via a small sidebar link that opens a modal (rather than a full tab).
  # Reads the stored parameters of the last fit (not the live sidebar) and the
  # resolved family, so the emitted script reproduces the model actually fitted.
  observeEvent(input$show_code, {
    showModal(modalDialog(
      title = "Reproduce this analysis in R",
      easyClose = TRUE,
      size = "l",
      markdown(
        "Copy or download the script below to reproduce your **last model fit** from
        the R console -- same models, VPC and PCV (the random seed makes the bootstrap
        intervals reproducible)."),
      div(class = "mb-2",
          downloadButton("download_code", "Download .R Script", class = "btn-secondary",
                         `aria-label` = "Download the reproducible R script")),
      verbatimTextOutput("repro_code"),
      footer = modalButton("Close")
    ))
  })

  repro_code_text <- reactive({
    req(model_results(), fit_params(), fitted_family())
    p <- fit_params()
    MAIHDA:::maihda_app_generate_code(
      outcome_var = p$outcome,
      grouping_vars = p$grouping_vars,
      additional_covars = p$covariates,
      family = fitted_family(),
      autobin = p$autobin,
      use_boot = p$use_boot,
      n_boot = p$n_boot,
      seed = p$seed,
      dataset = p$dataset,
      upload_name = p$upload_name,
      decomposition = if (!is.null(p$decomposition)) p$decomposition else "two-model"
    )
  })

  output$repro_code <- renderText({
    if (is.null(model_results()) || is.null(fit_params())) {
      return("Fit a model to generate a reproducible R script for this analysis.")
    }
    repro_code_text()
  })

  output$download_code <- downloadHandler(
    filename = function() "maihda_analysis.R",
    content = function(file) {
      writeLines(repro_code_text(), file)
    }
  )
}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
