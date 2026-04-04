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
  title = "MAIHDA Analysis Dashboard",
  theme = bs_theme(version = 5, primary = "#2C3E50", success = "#6BCF7F", info = "#4D9DE0"),

  sidebar = sidebar(
    title = "Controls",
    fileInput("upload", "Upload Data (CSV/RDS)", accept = c(".csv", ".rds")),

    # Model specification
    textInput("outcome", "Outcome Variable", value = "outcome"),
    textInput("group_vars", "Grouping Variables (comma-separated)", value = "sex, age_group"),

    # Engine settings
    selectInput("engine", "Engine", choices = c("lme4", "brms"), selected = "lme4"),
    selectInput("family", "Family", choices = c("gaussian", "binomial", "poisson"), selected = "gaussian"),

    # Plot settings
    sliderInput("n_strata", "Strata count for plots", min = 5, max = 50, value = 20),

    # Action button to trigger fitting
    actionButton("fit_btn", "Fit MAIHDA Model", class = "btn-primary")
  ),

  navset_card_tab(
    nav_panel("Data View",
              DTOutput("data_table")),
    nav_panel("Model Summary",
              verbatimTextOutput("model_summary")),
    nav_panel("PVC Results",
              verbatimTextOutput("pvc_summary")),
    nav_panel("Visualizations",
              selectInput("plot_type", "Select Plot Type",
                          choices = c("caterpillar", "vpc", "obs_vs_shrunken", "predicted")),
              plotOutput("maihda_plot", height = "500px"))
  )
)

server <- function(input, output, session) {

  # Load data: if no file, use maihda_sim_data
  reactive_data <- reactive({
    if (is.null(input$upload)) {
      return(MAIHDA::maihda_sim_data)
      } else {
      ext <- tools::file_ext(input$upload$name)
      if (ext == "csv") {
        return(read.csv(input$upload$datapath))
      } else if (ext == "rds") {
        return(readRDS(input$upload$datapath))
      } else {
        stop("Invalid file format")
      }
    }
  })

  output$data_table <- renderDT({
    datatable(reactive_data(), options = list(pageLength = 10, scrollX = TRUE))
  })

  # Reactive values for model and results
  model_results <- reactiveVal(NULL)
  summary_results <- reactiveVal(NULL)
  pvc_results <- reactiveVal(NULL)

  observeEvent(input$fit_btn, {
    dat <- reactive_data()
    req(dat)

    grouping_vars <- trimws(unlist(strsplit(input$group_vars, ",")))
    outcome_var <- input$outcome
    eng <- input$engine
    fam <- input$family

    if (eng == "brms") {
      # Use sequential plan to avoid conflict with brms internal parallelization
      future::plan(sequential)
    } else {
      future::plan(multisession)
    }

    # Reset old results
    model_results(NULL)
    summary_results(NULL)
    pvc_results(NULL)

    id <- showNotification("Creating strata...", duration = NULL, type = "message")

    # Formula construction
    fmla_str <- paste(outcome_var, "~", paste(grouping_vars, collapse = " + "), "+ (1 | stratum)")
    fmla <- as.formula(fmla_str)

    future_promise({
      # Step 1: Make strata
      strata_dat <- make_strata(dat, vars = grouping_vars)

      # Step 2: Fit model
      mod <- fit_maihda(formula = fmla, data = strata_dat, engine = eng, family = fam)

      # Step 3: Summarize and PVC
      summ <- summary(mod)
      pvc <- calculate_pvc(mod)

      list(model = mod, summary = summ, pvc = pvc)
    }, seed = TRUE) %...>% (function(res) {
      removeNotification(id)
      model_results(res$model)
      summary_results(res$summary)
      pvc_results(res$pvc)
    })

  })

  output$model_summary <- renderPrint({
    req(summary_results())
    print(summary_results())
  })

  output$pvc_summary <- renderPrint({
    req(pvc_results())
    print(pvc_results())
  })

  output$maihda_plot <- renderPlot({
    req(model_results())
    req(input$plot_type)

    if (input$plot_type %in% c("caterpillar", "predicted")) {
      plot(model_results(), type = input$plot_type, n_strata = input$n_strata)
    } else {
      plot(model_results(), type = input$plot_type)
    }
  })
}

shinyApp(ui = ui, server = server)
