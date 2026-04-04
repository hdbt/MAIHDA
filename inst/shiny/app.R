library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(MAIHDA)
library(future)
library(promises)

# Set up multisession for async processing
future::plan(multisession)

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
      fileInput("upload", "Upload Data (CSV/RDS/DTA/SAV)", accept = c(".csv", ".rds", ".dta", ".sav"))
    ),

    # Model specification
    selectizeInput("outcome", "Outcome Variable", choices = NULL),
    selectizeInput("group_vars", "Strata Grouping Variables", choices = NULL, multiple = TRUE),
    selectizeInput("covariates", "Additional Covariates (Fixed Effects)", choices = NULL, multiple = TRUE),

    # Model settings
    selectInput("family", "Family", choices = c("gaussian", "binomial", "poisson"), selected = "gaussian"),
    checkboxInput("use_boot", "Compute Bootstrap CIs (Slower)", value = FALSE),
    conditionalPanel(
      condition = "input.use_boot == true",
      numericInput("n_boot", "Bootstrap Samples", value = 100, min = 10, step = 10)
    ),

    # Action button to trigger fitting
    actionButton("fit_btn", "Fit MAIHDA Model", class = "btn-primary")
  ),

  navset_card_tab(
    nav_panel("Data View",
              DTOutput("data_table")),
    nav_panel("Model Summary",
              uiOutput("model_summary_ui")),
    nav_panel("PVC Results",
              uiOutput("pvc_summary_ui")),
    nav_panel("Stepwise PCV",
              uiOutput("stepwise_pcv_ui")),
    nav_panel("Visualizations",
              div(class = "d-flex justify-content-between align-items-center align-items-md-end mb-3",
                div(class = "flex-grow-1 me-3",
                  selectInput("plot_type", "Select Plot Type:",
                              choices = c("caterpillar", "vpc", "obs_vs_shrunken", "predicted"),
                              width = "100%")
                ),
                div(class = "mb-3",
                  downloadButton("download_plot", "Download Plot", class = "btn-secondary")
                )
              ),
              plotOutput("maihda_plot", height = "500px"))
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
        if (ext == "csv") {
          read.csv(input$upload$datapath)
        } else if (ext == "rds") {
          readRDS(input$upload$datapath)
        } else if (ext == "dta") {
          if (!requireNamespace('haven', quietly = TRUE)) stop("haven package required for DTA files")
          haven::as_factor(haven::read_dta(input$upload$datapath))
        } else if (ext == "sav") {
          if (!requireNamespace('haven', quietly = TRUE)) stop("haven package required for SAV files")
          haven::as_factor(haven::read_sav(input$upload$datapath))
        } else {
          stop("Unsupported format")
        }
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

  output$data_table <- renderDT({
    datatable(reactive_data(), options = list(pageLength = 10, scrollX = TRUE))
  })

  # Reactive values for model and results
  model_results <- reactiveVal(NULL)
  summary_results <- reactiveVal(NULL)
  pvc_results <- reactiveVal(NULL)
  stepwise_results <- reactiveVal(NULL)

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

    # Reset old results
    model_results(NULL)
    summary_results(NULL)
    pvc_results(NULL)
    stepwise_results(NULL)

    id <- showNotification("Creating strata & Fitting Models (May take a moment)...", duration = NULL, type = "message")

    # Formula construction
    if (length(additional_covars) > 0) {
        fmla_str <- paste(outcome_var, "~", paste(c(grouping_vars, additional_covars), collapse = " + "), "+ (1 | stratum)")
    } else {
        fmla_str <- paste(outcome_var, "~", paste(grouping_vars, collapse = " + "), "+ (1 | stratum)")
    }
    fmla <- as.formula(fmla_str)

    # Variables for stepwise PCV
    stepwise_vars <- c(grouping_vars, additional_covars)

    future_promise({
      # Step 1: Handle Missing Data to ensure identical sample sizes across all models.
      # If models have different Ns, variance comparisons (PCV) and parametric bootstraps are invalid.
      all_required_cols <- unique(c(outcome_var, grouping_vars, additional_covars))
      complete_dat <- dat[complete.cases(dat[, all_required_cols, drop = FALSE]), ]

      if (nrow(complete_dat) == 0) {
          stop("Error: No complete cases remaining after omitting missing values (NAs). Please select different variables.")
      }

      # Step 2: Make strata
      strata_dat <- make_strata(complete_dat, vars = grouping_vars)

      # Step 3: Fit model
      fmla_null <- as.formula(paste(outcome_var, "~ 1 + (1 | stratum)"))
      mod1 <- fit_maihda(formula = fmla_null, data = strata_dat$data, engine = eng, family = fam)
      mod2 <- fit_maihda(formula = fmla, data = strata_dat$data, engine = eng, family = fam)

      summ <- summary_maihda(mod2)
      pvc <- calculate_pvc(mod1, mod2, bootstrap = use_boot, n_boot = n_boot)

      stepwise <- stepwise_pcv(strata_dat$data, outcome = outcome_var, vars = stepwise_vars, engine = eng, family = fam)

      list(model = mod2, summary = summ, pvc = pvc, stepwise = stepwise)
    }, seed = TRUE) %...>% (function(res) {
      removeNotification(id)
      model_results(res$model)
      summary_results(res$summary)
      pvc_results(res$pvc)
      stepwise_results(res$stepwise)
    })

  })

  output$model_summary_ui <- renderUI({
    req(summary_results())
    res <- summary_results()

    tagList(
      card(
        card_header("Variance Partition Coefficient (VPC) / ICC"),
        h3(HTML(sprintf("<span class='text-primary'>%.2f%%</span>", res$vpc$estimate * 100)))
      ),
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
    null_formula <- paste(outcome_var, "~ 1 + (1 | stratum)")

    bootstrap_ui <- if (isTRUE(pvc$bootstrap) && !is.null(pvc$ci_lower) && !is.null(pvc$ci_upper)) {
        div(class = "mt-4 text-center text-muted",
            h5("Bootstrap 95% Confidence Interval"),
            tags$p(sprintf("[%.2f%%, %.2f%%]", pvc$ci_lower * 100, pvc$ci_upper * 100))
        )
    } else {
        NULL
    }

    card(
      card_header("Proportional Change in Variance (PVC)"),
      card_body(
        div(class = "d-flex justify-content-around text-center mb-4",
          div(
            h5("Null Model (Model 1)"),
            tags$code(null_formula),
            br(),br(),
            h5("Variance:"),
            h4(if (!is.null(pvc$var_model1)) sprintf("%.4f", pvc$var_model1) else "N/A")
          ),
          div(
            h5("Adjusted Model (Model 2)"),
            tags$code(paste(adjusted_formula, collapse = "")),
            br(),br(),
            h5("Variance:"),
            h4(if (!is.null(pvc$var_model2)) sprintf("%.4f", pvc$var_model2) else "N/A")
          )
        ),
        hr(),
        div(class = "text-center",
          h3(
            "Estimated PVC ",
            tooltip(
              shiny::icon("info-circle"),
              "PVC measures the reduction in between-stratum variance when moving from the Null model to the Adjusted model. High PVC = inequalities explained by additive characteristics. Low PVC = strong intersectional effects."
            )
          ),
          h2(class = "text-success", sprintf("%.2f%%", pvc$pvc * 100))
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
        This table displays how much between-stratum inequality is explained incrementally.

        *   **Step_PCV**: Percentage of variance explained compared to the *previous* model step.
        *   **Total_PCV**: Percentage of variance explained compared to the *null* model (Step 0).
        "),
        hr(),
        DTOutput("stepwise_pcv_dt")
      )
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

  current_plot <- reactive({
    req(model_results())
    req(input$plot_type)

    if (input$plot_type %in% c("caterpillar", "predicted")) {
      plot_maihda(model_results(), type = input$plot_type, n_strata = 20)
    } else {
      plot_maihda(model_results(), type = input$plot_type)
    }
  })

  output$maihda_plot <- renderPlot({
    current_plot()
  })

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("maihda_", input$plot_type, "_plot.png")
    },
    content = function(file) {
      ggsave(file, plot = current_plot(), width = 10, height = 8, dpi = 300)
    }
  )
}

shinyApp(ui = ui, server = server)
