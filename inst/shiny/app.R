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

ui <- page_sidebar(
  shinyjs::useShinyjs(),
  title = "MAIHDA Analysis Dashboard",
  theme = bs_theme(version = 5, primary = "#2C3E50", success = "#6BCF7F", info = "#4D9DE0"),

  sidebar = sidebar(
    title = "Controls",
    selectInput("dataset", "1. Select Dataset:",
                choices = c("Built-in: Simulated Data" = "sim",
                            "Built-in: NHANES Health Data" = "health",
                            "Upload Custom Data" = "upload")),
    conditionalPanel(
      condition = "input.dataset == 'upload'",
      fileInput("upload", "Upload Data (CSV/DTA/SAV)", accept = c(".csv", ".dta", ".sav"))
    ),

    # Model specification
    selectizeInput("outcome", "Outcome Variable", choices = NULL),
    selectizeInput("group_vars", "Strata Grouping Variables", choices = NULL, multiple = TRUE),
    uiOutput("group_var_hint"),
    checkboxInput("autobin", "Auto-bin continuous strata vars (>10 unique values) into 3 groups", value = TRUE),
    selectizeInput("covariates", "Additional Covariates (Fixed Effects)", choices = NULL, multiple = TRUE),

    # Model settings
    selectInput("family", "Family", choices = c("gaussian", "binomial", "poisson"), selected = "gaussian"),
    checkboxInput("use_boot", "Compute Bootstrap CIs (Slower)", value = FALSE),
    conditionalPanel(
      condition = "input.use_boot == true",
      numericInput("n_boot", "Bootstrap Samples", value = 100, min = 10, step = 10),
      numericInput("seed", "Random Seed (reproducible bootstrap)", value = 123, min = 1, step = 1)
    ),

    # Action button to trigger fitting
    actionButton("fit_btn", "Fit MAIHDA Model", class = "btn-primary")
  ),

  navset_card_tab(
    id = "main_tabs",
    nav_panel("Data View",
              shinycssloaders::withSpinner(DTOutput("data_table"))),
    nav_panel("Model Summary",
              uiOutput("model_summary_ui")),
    nav_panel("PCV Results",
              uiOutput("pvc_summary_ui")),
    nav_panel("Stepwise PCV",
              uiOutput("stepwise_pcv_ui")),
    nav_panel("Model Comparison",
              MAIHDA:::mod_compare_ui("compare")),
    nav_panel("Visualizations",
              MAIHDA:::mod_visualizations_ui("viz")),
    nav_panel("Interactive Explorer",
              MAIHDA:::mod_explorer_ui("explorer")),
    nav_panel("Reproduce in R",
              card(
                card_header("Reproducible R script for this analysis"),
                card_body(
                  markdown(
                    "Copy or download the script below to reproduce this dashboard
                    analysis from the R console. It uses the variables, family and
                    settings of your **last model fit**, so re-running it gives the
                    same models, VPC and PCV (the random seed makes the bootstrap
                    intervals reproducible)."
                  ),
                  div(class = "mb-3",
                      downloadButton("download_code", "Download .R Script", class = "btn-secondary")),
                  verbatimTextOutput("repro_code")
                )
              ))
  )
)

server <- function(input, output, session) {

  # Load data: if no file, use maihda_sim_data
  reactive_data <- reactive({
    if (input$dataset == "sim") {
      return(MAIHDA::maihda_sim_data)
    } else if (input$dataset == "health") {
      # Use the new real-world health dataset
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

    new_outcome <- ifelse(!is.null(curr_outcome) && curr_outcome %in% cols, curr_outcome, ifelse("health_outcome" %in% cols, "health_outcome", cols[1]))
    new_group <- if(!is.null(curr_group) && all(curr_group %in% cols)) curr_group else intersect(c("gender", "race"), cols)

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
  fit_params <- reactiveVal(NULL)      # inputs of the last fit (for the "Reproduce in R" tab)
  da_results <- reactiveVal(NULL)      # discriminatory accuracy (binomial only): null + adjusted
  comparison_results <- reactiveVal(NULL)  # nested-model VPC comparison (null vs adjusted)

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
    seed_opt <- if (isTRUE(use_boot)) input$seed else NULL

    # Claim a request token and lock the button so a second click cannot launch a
    # competing fit whose (possibly out-of-order) result would clobber this one.
    this_fit <- fit_id() + 1
    fit_id(this_fit)
    shinyjs::disable("fit_btn")

    # Remember exactly what this fit used so the "Reproduce in R" tab mirrors it,
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
      seed = seed_opt
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
          seed = seed_opt
        )
      }, seed = TRUE) %...>% (function(res) {
        removeNotification(id)
        # A newer fit has superseded this one: drop the stale result and leave the
        # button for the newer fit to re-enable when it finishes.
        if (!identical(this_fit, fit_id())) return(invisible(NULL))
        shinyjs::enable("fit_btn")

        model_results(res$model)
        fitted_family(res$family_used)
        # Dispatch summary() in the main thread, where the MAIHDA S3 method is
        # reliably found, then merge in the VPC/ICC bootstrap intervals computed
        # in the background worker (see maihda_app_bootstrap_vpc_cis()).
        null_summary_results(MAIHDA:::maihda_app_attach_vpc_ci(summary(res$null_model), res$vpc_ci_null))
        summary_results(MAIHDA:::maihda_app_attach_vpc_ci(summary(res$model), res$vpc_ci_adjusted))
        pvc_results(res$pvc)
        stepwise_results(res$stepwise)
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
        if (isTRUE(res$family_autoswitched)) {
          showNotification(
            sprintf("Outcome '%s' is binary -- fitted as 'binomial' (not the selected 'gaussian').",
                    outcome_var),
            type = "warning", duration = 12
          )
        }
        nav_select("main_tabs", "PCV Results")
      }) %...!% (function(err) {
        removeNotification(id)
        if (!identical(this_fit, fit_id())) return(invisible(NULL))
        shinyjs::enable("fit_btn")
        showNotification(paste("Error fitting model:", err$message), type = "error", duration = 15)
      })
  })

  output$model_summary_ui <- renderUI({
    req(summary_results())
    res <- summary_results()

    vpc <- res$vpc
    vpc_interval <- if (MAIHDA:::maihda_vpc_has_interval(vpc)) {
      div(class = "text-muted",
          sprintf("[%.2f%%, %.2f%%] %s",
                  vpc$ci_lower * 100, vpc$ci_upper * 100,
                  MAIHDA:::maihda_vpc_interval_label(vpc)))
    } else {
      NULL
    }

    family_line <- if (!is.null(fitted_family())) {
      div(class = "text-muted small", sprintf("Fitted with family = '%s'", fitted_family()))
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
        div(class = "d-flex justify-content-around text-center",
            div(h5("AUC -- strata only"),
                h3(fmt_metric(if (!is.null(da$null)) da$null$auc else NA)),
                p(class = "text-muted mb-0", "C-statistic of the intersectional strata alone")),
            div(h5("AUC -- adjusted"),
                h3(fmt_metric(if (!is.null(da$adjusted)) da$adjusted$auc else NA)),
                p(class = "text-muted mb-0", "With individual covariates added")),
            div(h5("Median Odds Ratio"),
                h3(fmt_metric(if (!is.null(da$null)) da$null$mor else NA, 2)),
                p(class = "text-muted mb-0", "Between-stratum heterogeneity on the odds-ratio scale"))
        ),
        div(class = "small text-muted mt-2",
            "AUC = 0.5 is chance. A high between-stratum VPC can still translate into only modest individual-level discriminatory accuracy -- the cautionary message at the heart of the 'DA' in MAIHDA.")
      )
    } else {
      NULL
    }

    tagList(
      card(
        card_header("Variance Partition Coefficient (VPC) / ICC"),
        h3(HTML(sprintf("<span class='text-primary'>%.2f%%</span>", vpc$estimate * 100))),
        vpc_interval,
        family_line
      ),
      card(
        card_header("Fit diagnostics & strata overview"),
        diag_ui,
        strata_ui
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
    req(pvc_results())
    req(model_results())

    pvc <- pvc_results()
    mod <- model_results()

    adjusted_formula <- deparse(mod$formula)
    outcome_var <- all.vars(mod$formula)[1]
    null_formula <- paste(
      deparse(MAIHDA:::maihda_formula_with_stratum(outcome_var)),
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
        div(class = "d-flex justify-content-around text-center mb-4",
          div(
            h5("Null Model (Model 1)"),
            tags$code(null_formula),
            br(),br(),
            h5("Variance:"),
            h4(if (!is.null(pvc$var_model1) && is.finite(pvc$var_model1)) sprintf("%.4f", pvc$var_model1) else "N/A")
          ),
          div(
            h5("Adjusted Model (Model 2)"),
            tags$code(paste(adjusted_formula, collapse = "")),
            br(),br(),
            h5("Variance:"),
            h4(if (!is.null(pvc$var_model2) && is.finite(pvc$var_model2)) sprintf("%.4f", pvc$var_model2) else "N/A")
          )
        ),
        hr(),
        div(class = "text-center",
          h3(
            "Estimated PCV ",
            tooltip(
              shiny::icon("info-circle"),
              "PCV is the proportional change in between-stratum variance from the Null to the Adjusted model. A high PCV means the between-stratum variance is much smaller after adding the additive main effects; a low or negative PCV means little change (or an increase). This is a model-dependent change, not proof that inequality was causally 'explained away' -- it can also reflect suppression, rescaling, sample composition, or uncertainty, not interaction alone."
            )
          ),
          if (is.finite(pvc$pvc)) {
            h2(class = "text-success", sprintf("%.2f%%", pvc$pvc * 100))
          } else {
            tagList(
              h2(class = "text-muted", "N/A"),
              div(class = "alert alert-warning text-start",
                  tags$strong("PCV could not be calculated. "),
                  if (!is.null(pvc$message)) pvc$message else
                    "The baseline between-stratum variance is zero, so the proportional change is undefined. The model fit, VPC and visualizations above remain valid.")
            )
          }
        ),
        bootstrap_ui
      )
    )
  })

  output$stepwise_pcv_ui <- renderUI({
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
    group_vars = reactive(input$group_vars)
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
  # Reads the stored parameters of the last fit (not the live sidebar) and the
  # resolved family, so the emitted script reproduces the model actually fitted.
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
      upload_name = p$upload_name
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

shinyApp(ui = ui, server = server)
