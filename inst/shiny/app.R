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
    checkboxInput("autobin", "Auto-bin continuous strata vars (>10 unique values) into 3 groups", value = TRUE),
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
    id = "main_tabs",
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
                              choices = c(
                                          "Prediction Deviation Panels" = "pred_dev",
                                          "Risk vs. Intersectional Effect (Quadrant)" = "risk_vs_effect",
                                          "Effect Decomposition" = "effect_decomp",
                                          "Effect Decomposition (Ternary)" = "ternary",
                                          "VPC"="vpc", "Observed VS Shrunken" = "obs_vs_shrunken",
                                          "Predicted Values" =  "predicted"),
                              width = "100%")
                ),
                div(class = "mb-3",
                  downloadButton("download_plot", "Download Plot", class = "btn-secondary")
                )
              ),
              uiOutput("maihda_plot_wrapper")),
    nav_panel("Interactive Explorer",
              uiOutput("interactive_explorer_ui"))
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
  null_summary_results <- reactiveVal(NULL)
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
    autobin_opt <- input$autobin

    # Reset old results
    model_results(NULL)
    null_summary_results(NULL)
    summary_results(NULL)
    pvc_results(NULL)
    stepwise_results(NULL)

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
          engine = eng
        )
      }, seed = TRUE) %...>% (function(res) {
        removeNotification(id)
        model_results(res$model)
        # Call S3 dispatch in the main thread where MAIHDA environment is perfectly active
        null_summary_results(summary(res$null_model))
        summary_results(summary(res$model))
      pvc_results(res$pvc)
      stepwise_results(res$stepwise)
      nav_select("main_tabs", "PVC Results")    }) %...!% (function(err) {
      removeNotification(id)
      showNotification(paste("Error fitting model:", err$message), type = "error", duration = 15)    })

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
    null_formula <- paste(
      deparse(MAIHDA:::maihda_formula_with_stratum(outcome_var)),
      collapse = ""
    )

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
        plotlyOutput("stepwise_pcv_plot", height = "400px"),
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
        title = "Cumulative Intersectional Variance Explained",
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

  current_plot <- reactive({
    req(model_results())
    req(input$plot_type)

    if (input$plot_type == "pred_dev") {
      plot_prediction_deviation_panels(model_results(), data = NULL, type = "auto")
    } else if (input$plot_type %in% c("predicted")) {
      plot(model_results(), type = input$plot_type, n_strata = 20)
    } else if (input$plot_type == "ternary") {
      out <- maihda_ternary_plot(model_results())
      out$plot
    } else {
      plot(model_results(), type = input$plot_type)
    }
  })

  output$maihda_plot_wrapper <- renderUI({
    if (input$plot_type == "ternary") {
      plotlyOutput("maihda_plotly", height = "500px")
    } else {
      plotOutput("maihda_plot", height = "500px")
    }
  })

  output$maihda_plot <- renderPlot({
    current_plot()
  })

  output$maihda_plotly <- renderPlotly({
    req(model_results())
    req(input$plot_type == "ternary")

    out <- maihda_ternary_plot(model_results())
    td <- out$data

    marker_sizes <- pmax(sqrt(td$n) * 2, 4)
    marker_colors <- as.numeric(as.factor(td$label))

    plot_ly(
      data = td,
      type = 'scatterternary',
      mode = 'markers',
      a = ~additive_prop,
      b = ~interaction_prop,
      c = ~uncertainty_prop,
      text = ~paste0(
        "<b>Stratum:</b> ", label, "<br>",
        "<b>Size (N):</b> ", n, "<br>",
        "<b>Additive:</b> ", round(additive_prop * 100, 1), "%<br>",
        "<b>Intersection-specific:</b> ", round(interaction_prop * 100, 1), "%<br>",
        "<b>Uncertainty:</b> ", round(uncertainty_prop * 100, 1), "%"
      ),
      hoverinfo = "text",
      marker = list(
        size = marker_sizes,
        color = marker_colors,
        colorscale = "Viridis",
        opacity = 0.8,
        line = list(color = "rgba(0,0,0,0.5)", width = 1)
      )
    ) |>
      layout(
        title = "Interactive MAIHDA Strata Effects Decomposition",
        ternary = list(
          sum = 1,
          aaxis = list(title = "Additive", min = 0.01, linewidth = 2, ticks = "outside", tickvals = seq(0, 1, by = 0.2)),
          baxis = list(title = "Intersection", min = 0.01, linewidth = 2, ticks = "outside", tickvals = seq(0, 1, by = 0.2)),
          caxis = list(title = "Uncertainty", min = 0.01, linewidth = 2, ticks = "outside", tickvals = seq(0, 1, by = 0.2))
        ),
        margin = list(t = 50, b = 50, l = 50, r = 50)
      )
  })

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("maihda_", input$plot_type, "_plot.png")
    },
    content = function(file) {
      ggsave(file, plot = current_plot(), width = 10, height = 8, dpi = 300)
    }
  )

  output$interactive_explorer_ui <- renderUI({
    req(model_results(), null_summary_results(), summary_results(), pvc_results())
    null_res <- null_summary_results()
    res <- summary_results()
    pvc <- pvc_results()

    # Extract metrics for HUD
    vpc_val <- round(null_res$vpc$estimate * 100, 2)
    pvc_val <- round(pvc$pvc * 100, 2)

    layout_columns(
      col_widths = c(12, 12),
      card(
        card_header("HUD: Key MAIHDA Metrics"),
        div(class = "d-flex justify-content-around text-center",
            div(h4("VPC (Null)"), h3(paste0(vpc_val, "%")), p(class="text-muted", "Total Variance b/w Strata")),
            div(h4("PVC (Adjusted)"), h3(paste0(pvc_val, "%")), p(class="text-muted", "Variance Explained by Main Effects")),
            div(h4("Intersectionality"), h3(paste0(100 - pvc_val, "%")), p(class="text-muted", "Unexplained Variance (Interaction Effects)"))
        ),
        markdown("
        **Interpretation Guide**:
        - **VPC** (Variance Partition Coefficient) measures how much of the total outcome variance is due to the strata definitions.
        - **PVC** (Proportional Change in Variance) shows how much of that strata variation is explained by simple additive effects.
        - The remaining percentage represents the true **intersectional effect**, revealing disparities unique to specific strata combinations.
        "),
        uiOutput("dynamic_interpretation")
      ),
      card(
        card_header("Interactive Strata Deviations (Residuals with CIs)"),
        layout_columns(
          col_widths = c(4, 4, 4),
          selectInput("hud_color_var", "Color Points By:",
                      choices = c("Significance (Tolerance)" = "deviant", isolate(input$group_vars))),
          selectInput("hud_sort_var", "Sort Y-Axis By:",
                      choices = c("Effect Size (Magnitude)" = "effect", "Sample Size (N)" = "n", "Alphabetical" = "alpha")),
          sliderInput("hud_top_n", "Show Top Deviant Strata (by Magnitude):",
                      min = 5, max = max(5, nrow(res$stratum_estimates)),
                      value = min(25, nrow(res$stratum_estimates)), step = 1)
        ),
        plotlyOutput("interactive_plot", height = "600px"),
        markdown("
        *Hover over the points to see individual stratum details.*
        - Points far from the zero-line (red) represent **deviant strata**: groups whose outcome significantly departs from what simple additive effects would predict.
        - Error bars represent 95% Confidence Intervals (simulated/bootstrap). If the bar does not cross zero, the intersectional effect is statistically significant.
        - **Point size** represents the total number of individuals (N) within that stratum configuration.
        ")
      ),
      card(
        card_header("Filtered Strata Data Export"),
        div(class = "mb-3", downloadButton("download_hud_data", "Download Highlighted Data (CSV)", class = "btn-secondary")),
        DTOutput("interactive_table")
      )
    )
  })

  # Reactive containing exactly the dataframe filtered for HUD exploring
  hud_plot_data <- reactive({
    req(summary_results(), model_results())

    # Build a simple data frame for plotting
    stratum_df <- as.data.frame(summary_results()$stratum_estimates)

    # Merge with strata_info to get specific variables (N, gender, race, etc.)
    strata_info <- model_results()$strata_info
    if (!is.null(strata_info)) {
      # resolve duplicate column names gracefully
      cols_to_merge <- setdiff(names(strata_info), names(stratum_df))
      stratum_df <- merge(stratum_df, strata_info[, c("stratum", cols_to_merge), drop = FALSE], by = "stratum", all.x = TRUE)
    }

    # Add Absolute Predicted Values via margin average
    mod <- model_results()
    if (!is.null(mod$data)) {
        pred_vals <- tryCatch({
            pred <- predict_maihda(mod)
            agg <- aggregate(pred ~ stratum, data = mod$data, FUN = mean)
            names(agg)[2] <- "abs_pred"
            agg
        }, error = function(e) NULL)
        if (!is.null(pred_vals)) stratum_df <- merge(stratum_df, pred_vals, by = "stratum", all.x = TRUE)
    }

    # Use stratum labels if generated, otherwise default to IDs
    if ("label" %in% names(stratum_df) && !all(is.na(stratum_df$label))) {
      stratum_df$display_label <- paste0(stratum_df$stratum, ": ", stratum_df$label)
    } else {
      stratum_df$display_label <- paste0("Stratum ", stratum_df$stratum)
    }

    # Add a flag for 'deviant' (significant at 95%)
    if (!"lower_95" %in% names(stratum_df)) stratum_df$lower_95 <- stratum_df$random_effect - 1.96 * stratum_df$se
    if (!"upper_95" %in% names(stratum_df)) stratum_df$upper_95 <- stratum_df$random_effect + 1.96 * stratum_df$se

    stratum_df$deviant <- ifelse(stratum_df$lower_95 > 0 | stratum_df$upper_95 < 0, "Significant", "Not Significant")

    # Filter the Top N Deviant strata (by highest absolute effect, retaining original signs)
    if (!is.null(input$hud_top_n)) {
      stratum_df <- stratum_df[order(abs(stratum_df$random_effect), decreasing = TRUE), ]
      stratum_df <- head(stratum_df, input$hud_top_n)
    }

    stratum_df
  })

  output$interactive_plot <- renderPlotly({
    req(hud_plot_data())
    stratum_df <- hud_plot_data()

    # Y-axis Sorting Control
    sort_by <- if (!is.null(input$hud_sort_var)) input$hud_sort_var else "effect"
    if (sort_by == "n" && "n" %in% names(stratum_df)) {
        stratum_df$display_label <- factor(stratum_df$display_label, levels = stratum_df$display_label[order(stratum_df$n)])
    } else if (sort_by == "alpha") {
        stratum_df$display_label <- factor(stratum_df$display_label, levels = rev(stratum_df$display_label[order(as.character(stratum_df$display_label))]))
    } else {
        stratum_df$display_label <- factor(stratum_df$display_label, levels = stratum_df$display_label[order(stratum_df$random_effect)])
    }

    # Create tooltip format
    n_text <- if ("n" %in% names(stratum_df)) paste("<br>Sample Size (N):", stratum_df$n) else ""
    abs_text <- if ("abs_pred" %in% names(stratum_df)) paste("<br>Absolute Pred. Outcome:", round(stratum_df$abs_pred, 3)) else ""
    stratum_df$tooltip <- paste0("<b>", stratum_df$display_label, "</b>",
                                  n_text,
                                  abs_text,
                                  "<br>Effect:", round(stratum_df$random_effect, 3),
                                  "<br>95% CI:", round(stratum_df$lower_95, 3), " to ", round(stratum_df$upper_95, 3))

    # Choose mapping variables
    color_var <- if(!is.null(input$hud_color_var)) input$hud_color_var else "deviant"
    size_mapped <- "n" %in% names(stratum_df)

    p <- ggplot(stratum_df, aes(x = random_effect, y = display_label,
                                color = .data[[color_var]],
                                text = tooltip)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

    if (size_mapped) {
       p <- p + geom_point(aes(size = n), alpha = 0.8) + scale_size_continuous(range = c(2, 6))
    } else {
       p <- p + geom_point(size = 3)
    }

    p <- p + geom_errorbar(aes(xmin = lower_95, xmax = upper_95), width = 0.2, orientation = "y") +
      theme_minimal() +
      labs(x = "Intersectional Intercept / Effect (Deviation)",
           y = "Stratum", color = tools::toTitleCase(color_var), size = "Sample Size (N)") +
      theme(axis.text.y = element_text(size = 8))

    # If using standard deviant coloring, retain manual scale
    if (color_var == "deviant") {
        p <- p + scale_color_manual(values = c("Significant" = "#E74C3C", "Not Significant" = "#34495E"))
    }

    # Disable tooltip for size parameter so it doesn't double-up and break Plotly's rendering gracefully
    ggplotly(p, tooltip = "text")
  })

  output$dynamic_interpretation <- renderUI({
    req(null_summary_results(), summary_results(), pvc_results())

    # Grab data
    null_res <- null_summary_results()
    res <- summary_results()
    pvc <- pvc_results()
    df <- as.data.frame(res$stratum_estimates)
    strata_info <- model_results()$strata_info

    # Merge N into df if available
    if (!is.null(strata_info) && "n" %in% names(strata_info)) {
      df <- merge(df, strata_info[, c("stratum", "n")], by = "stratum", all.x = TRUE)
    } else {
      df$n <- "Unknown"
    }

    # Find the most deviant stratum
    df <- df[order(abs(df$random_effect), decreasing = TRUE), ]
    most_deviant <- df[1, ]

    # Stratum Label
    if ("label" %in% names(most_deviant) && !is.na(most_deviant$label)) {
      dev_label <- paste0(most_deviant$stratum, " (", most_deviant$label, ")")
    } else {
      dev_label <- paste0("Stratum ", most_deviant$stratum)
    }

    dev_effect <- round(most_deviant$random_effect, 3)
    dev_n <- most_deviant$n

    # Extract metrics
    vpc_val <- round(null_res$vpc$estimate * 100, 2)
    pvc_val <- round(pvc$pvc * 100, 2)
    interaction_val <- 100 - pvc_val

    # Construct the summary paragraph dynamically
    tags$div(class = "alert alert-info mt-3",
      tags$strong("Automated Research Summary: "),
      "In this analysis, ", tags$strong(paste0(vpc_val, "%")),
      " of the total variance in the outcome is attributable to the defined intersecting demographic or social strata. ",
      "When considering simple additive (main) effects, ", tags$strong(paste0(pvc_val, "%")),
      " of this between-strata disparity is explained away, meaning that ", tags$strong(paste0(interaction_val, "%")),
      " of the disparity represents unique intersectional interaction effects not captured by standard main-effect modeling. ",
      "Exploring these residuals reveals that the most prominent intersectional disparity occurs in ",
      tags$strong(dev_label), " (N = ", dev_n, "), which shows an intersectional deviation score of ",
      tags$strong(dev_effect), " from what simple additive effects would predict."
    )
  })

  output$interactive_table <- renderDT({
    req(hud_plot_data())
    df <- hud_plot_data()

    # Drop tooltip and internal parsing columns before showing table
    cols_to_drop <- c("tooltip", "display_label")
    df <- df[, !names(df) %in% cols_to_drop]

    # Round numerics
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], round, 3)

    datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$download_hud_data <- downloadHandler(
    filename = function() {
      paste0("maihda_highlighted_strata_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- hud_plot_data()
      cols_to_drop <- c("tooltip", "display_label")
      df <- df[, !names(df) %in% cols_to_drop]
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
