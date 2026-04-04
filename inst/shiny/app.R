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
    selectizeInput("outcome", "Outcome Variable", choices = NULL),
    selectizeInput("group_vars", "Grouping Variables", choices = NULL, multiple = TRUE),

    # Engine settings
    selectInput("engine", "Engine", choices = c("lme4", "brms"), selected = "lme4"),
    selectInput("family", "Family", choices = c("gaussian", "binomial", "poisson"), selected = "gaussian"),

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
  
  observe({
    req(reactive_data())
    cols <- names(reactive_data())
    updateSelectizeInput(session, "outcome", choices = cols, selected = ifelse("health_outcome" %in% cols, "health_outcome", cols[1]), server = TRUE)
    updateSelectizeInput(session, "group_vars", choices = cols, selected = intersect(c("gender", "race"), cols), server = TRUE)
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

    grouping_vars <- input$group_vars
    req(length(grouping_vars) > 0)
    
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
      fmla_null <- as.formula(paste(outcome_var, "~ 1 + (1 | stratum)"))
      mod1 <- fit_maihda(formula = fmla_null, data = strata_dat$data, engine = eng, family = fam)
      mod2 <- fit_maihda(formula = fmla, data = strata_dat$data, engine = eng, family = fam)

      summ <- summary_maihda(mod2)
      pvc <- calculate_pvc(mod1, mod2)

      list(model = mod2, summary = summ, pvc = pvc)
    }, seed = TRUE) %...>% (function(res) {
      removeNotification(id)
      model_results(res$model)
      summary_results(res$summary)
      pvc_results(res$pvc)
    })

  })

  output$model_summary_ui <- renderUI({
    req(summary_results())
    res <- summary_results()
    
    tagList(
      card(
        card_header("Variance Partition Coefficient (VPC) / ICC"),
        h3(HTML(sprintf("<span class='text-primary'>%.2f%%</span>", res$vpc * 100)))
      ),
      layout_columns(
        card(
          card_header("Variance Components"),
          DT::renderDT(datatable(res$variance_components, options = list(dom = 't', paging = FALSE)))
        ),
        card(
          card_header("Fixed Effects"),
          DT::renderDT(datatable(as.data.frame(res$fixed_effects), options = list(dom = 't', paging = FALSE)))
        )
      ),
      card(
        card_header("Stratum Estimates (top 10)"),
        DT::renderDT(datatable(head(res$stratum_estimates, 10), options = list(dom = 't', paging = FALSE)))
      )
    )
  })

  output$pvc_summary_ui <- renderUI({
    req(pvc_results())
    pvc <- pvc_results()
    
    card(
      card_header("Proportional Change in Variance (PVC)"),
      card_body(
        div(class = "d-flex justify-content-around text-center mb-4",
          div(
            h5("Null Model (Model 1) Variance"),
            h4(sprintf("%.4f", pvc$var_model1))
          ),
          div(
            h5("Adjusted Model (Model 2) Variance"),
            h4(sprintf("%.4f", pvc$var_model2))
          )
        ),
        hr(),
        div(class = "text-center",
          h3("Estimated PVC"),
          h2(class = "text-success", sprintf("%.2f%%", pvc$pvc * 100))
        )
      )
    )
  })

  output$maihda_plot <- renderPlot({
    req(model_results())
    req(input$plot_type)

    if (input$plot_type %in% c("caterpillar", "predicted")) {
      plot_maihda(model_results(), type = input$plot_type, n_strata = 20)
    } else {
      plot_maihda(model_results(), type = input$plot_type)
    }
  })
}

shinyApp(ui = ui, server = server)
