# Shiny modules for run_maihda_app().
#
# These live in R/ (rather than inline in inst/shiny/app.R) so they are unit
# testable with shiny::testServer() and covered by R CMD check. They run in the
# MAIHDA namespace, which does not attach shiny/bslib/DT/plotly, so every UI and
# server helper is fully qualified; the suggested packages are guaranteed present
# because run_maihda_app() gates on maihda_app_required_packages(). IDs created
# inside server-side renderUI() are namespaced with session$ns(), as Shiny does
# not auto-namespace dynamically generated UI.

# --- Visualizations tab ------------------------------------------------------
# Static plot-type picker + download button + a wrapper that swaps a static
# ggplot output for an interactive plotly one (the ternary view). Depends only on
# the fitted model.

#' @noRd
mod_visualizations_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "d-flex justify-content-between align-items-center align-items-md-end mb-3",
      shiny::div(
        class = "flex-grow-1 me-3",
        shiny::selectInput(ns("plot_type"), "Select Plot Type:",
          choices = c(
            "Prediction Deviation Panels" = "pred_dev",
            "Mean Prediction vs. Stratum Effect (Quadrant)" = "risk_vs_effect",
            "Effect Decomposition" = "effect_decomp",
            "Effect Decomposition (Ternary)" = "ternary",
            "VPC" = "vpc", "Observed VS Shrunken" = "obs_vs_shrunken",
            "Predicted Values" = "predicted"),
          width = "100%")
      ),
      shiny::div(
        class = "mb-3",
        shiny::downloadButton(ns("download_plot"), "Download Plot", class = "btn-secondary")
      )
    ),
    shiny::uiOutput(ns("maihda_plot_wrapper"))
  )
}

#' @noRd
mod_visualizations_server <- function(id, model_results) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_plot <- shiny::reactive({
      shiny::req(model_results())
      shiny::req(input$plot_type)

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

    output$maihda_plot_wrapper <- shiny::renderUI({
      if (input$plot_type == "ternary") {
        shinycssloaders::withSpinner(plotly::plotlyOutput(ns("maihda_plotly"), height = "500px"))
      } else {
        shinycssloaders::withSpinner(shiny::plotOutput(ns("maihda_plot"), height = "500px"))
      }
    })

    output$maihda_plot <- shiny::renderPlot({
      current_plot()
    })

    output$maihda_plotly <- plotly::renderPlotly({
      shiny::req(model_results())
      shiny::req(input$plot_type == "ternary")

      out <- maihda_ternary_plot(model_results())
      maihda_app_ternary_plotly(out$data)
    })

    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        paste0("maihda_", input$plot_type, "_plot.png")
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = current_plot(), width = 10, height = 8, dpi = 300)
      }
    )

    current_plot
  })
}

# --- Interactive Explorer (HUD) tab ------------------------------------------
# Key-metrics header + interactive strata-deviation plot with filters + filtered
# data export. Depends on the fitted model and the null/adjusted summaries and
# PCV. `group_vars` is a reactive giving the grouping-variable names (for the
# "colour by" choices). Returns the filtered hud_plot_data reactive for testing.

#' @noRd
mod_explorer_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("interactive_explorer_ui"))
}

#' @noRd
mod_explorer_server <- function(id, model_results, null_summary_results,
                                summary_results, pvc_results, group_vars) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$interactive_explorer_ui <- shiny::renderUI({
      shiny::req(model_results(), null_summary_results(), summary_results(), pvc_results())
      null_res <- null_summary_results()
      res <- summary_results()
      pvc <- pvc_results()

      # Extract metrics for HUD
      null_vpc <- null_res$vpc
      vpc_val <- round(null_vpc$estimate * 100, 2)
      vpc_ci_text <- if (maihda_vpc_has_interval(null_vpc)) {
        sprintf("95%% CI [%.2f%%, %.2f%%]", null_vpc$ci_lower * 100, null_vpc$ci_upper * 100)
      } else {
        NULL
      }
      pvc_val <- round(pvc$pvc * 100, 2)
      pvc_val_display <- if (is.finite(pvc_val)) paste0(pvc_val, "%") else "N/A"
      pvc_display <- maihda_app_pvc_display(pvc_val)

      bslib::layout_columns(
        col_widths = c(12, 12),
        bslib::card(
          bslib::card_header("HUD: Key MAIHDA Metrics"),
          shiny::div(class = "d-flex justify-content-around text-center",
              shiny::div(shiny::h4("VPC (Null)"), shiny::h3(paste0(vpc_val, "%")),
                  shiny::p(class = "text-muted mb-0", "Total Variance b/w Strata"),
                  if (!is.null(vpc_ci_text)) shiny::p(class = "text-muted small mb-0", vpc_ci_text) else NULL),
              shiny::div(shiny::h4("PCV (Adjusted)"), shiny::h3(pvc_val_display), shiny::p(class = "text-muted", "Between-Stratum Variance Change with Main Effects")),
              shiny::div(shiny::h4(pvc_display$label), shiny::h3(pvc_display$value), shiny::p(class = "text-muted", pvc_display$description))
          ),
          shiny::markdown("
          **Interpretation Guide**:
          - **VPC** (Variance Partition Coefficient) measures the share of the unexplained outcome variance that lies between strata.
          - **PCV** (Proportional Change in Variance) is the proportional change in between-stratum variance when the additive main effects are added. It is a model-dependent comparison, not a causal measure of variance 'explained'.
          - The remaining between-stratum variation is often read as the intersectional component, but it is model-dependent and should be interpreted cautiously (a negative PCV does not by itself prove hidden structural inequality).
          "),
          shiny::uiOutput(ns("dynamic_interpretation"))
        ),
        bslib::card(
          bslib::card_header("Interactive Strata Deviations (Residuals with CIs)"),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            shiny::selectInput(ns("hud_color_var"), "Color Points By:",
                        choices = c("Conditional Interval Status" = "deviant", shiny::isolate(group_vars()))),
            shiny::selectInput(ns("hud_sort_var"), "Sort Y-Axis By:",
                        choices = c("Effect Size (Magnitude)" = "effect", "Sample Size (N)" = "n", "Alphabetical" = "alpha")),
            shiny::sliderInput(ns("hud_top_n"), "Show Top Strata (by Effect Magnitude):",
                        min = 5, max = max(5, nrow(res$stratum_estimates)),
                        value = min(25, nrow(res$stratum_estimates)), step = 1)
          ),
          shinycssloaders::withSpinner(plotly::plotlyOutput(ns("interactive_plot"), height = "600px")),
          shiny::markdown("
          *Hover over the points to see individual stratum details.*
          - Points far from the zero-line (red) are the **most extreme strata** for exploration: groups whose outcome departs most from what the additive main effects alone would predict (a descriptive flag, not a model-misfit diagnosis).
          - Error bars are approximate conditional intervals for stratum random effects. If the bar does not cross zero, treat it as a screening signal, not a formal bootstrap or posterior significance test.
          - **Point size** represents the total number of individuals (N) within that stratum configuration.
          ")
        ),
        bslib::card(
          bslib::card_header("Filtered Strata Data Export"),
          shiny::div(class = "mb-3", shiny::downloadButton(ns("download_hud_data"), "Download Highlighted Data (CSV)", class = "btn-secondary")),
          DT::DTOutput(ns("interactive_table"))
        )
      )
    })

    # Reactive containing exactly the dataframe filtered for HUD exploring
    hud_plot_data <- shiny::reactive({
      shiny::req(summary_results(), model_results())

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

      # Add an exploratory flag for conditional intervals that exclude zero.
      if (!"lower_95" %in% names(stratum_df)) stratum_df$lower_95 <- stratum_df$random_effect - 1.96 * stratum_df$se
      if (!"upper_95" %in% names(stratum_df)) stratum_df$upper_95 <- stratum_df$random_effect + 1.96 * stratum_df$se

      stratum_df$deviant <- ifelse(stratum_df$lower_95 > 0 | stratum_df$upper_95 < 0, "Excludes zero", "Includes zero")

      # Filter the Top N Deviant strata (by highest absolute effect, retaining original signs)
      if (!is.null(input$hud_top_n)) {
        stratum_df <- stratum_df[order(abs(stratum_df$random_effect), decreasing = TRUE), ]
        stratum_df <- head(stratum_df, input$hud_top_n)
      }

      stratum_df
    })

    output$interactive_plot <- plotly::renderPlotly({
      shiny::req(hud_plot_data())
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
                                    "<br>Approx. conditional interval:", round(stratum_df$lower_95, 3), " to ", round(stratum_df$upper_95, 3))

      # Choose mapping variables
      color_var <- if (!is.null(input$hud_color_var)) input$hud_color_var else "deviant"
      size_mapped <- "n" %in% names(stratum_df)

      p <- ggplot2::ggplot(stratum_df, ggplot2::aes(x = .data$random_effect, y = .data$display_label,
                                  color = .data[[color_var]],
                                  text = .data$tooltip)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

      if (size_mapped) {
        p <- p + ggplot2::geom_point(ggplot2::aes(size = .data$n), alpha = 0.8) + ggplot2::scale_size_continuous(range = c(2, 6))
      } else {
        p <- p + ggplot2::geom_point(size = 3)
      }

      p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95), width = 0.2, orientation = "y") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Intersectional Intercept / Effect (Deviation)",
             y = "Stratum", color = tools::toTitleCase(color_var), size = "Sample Size (N)") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

      # If using standard deviant coloring, retain manual scale
      if (color_var == "deviant") {
        p <- p + ggplot2::scale_color_manual(values = c("Excludes zero" = "#E74C3C", "Includes zero" = "#34495E"))
      }

      # Disable tooltip for size parameter so it doesn't double-up and break Plotly's rendering gracefully
      plotly::ggplotly(p, tooltip = "text")
    })

    output$dynamic_interpretation <- shiny::renderUI({
      shiny::req(null_summary_results(), summary_results(), pvc_results())

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
      pvc_display <- maihda_app_pvc_display(pvc_val)
      pvc_interpretation <- if (identical(pvc_display$status, "negative")) {
        shiny::tagList(
          "After adding the selected main effects, between-strata variance increases by ",
          shiny::tags$strong(pvc_display$value),
          ", a suppression or unmasking pattern. The adjusted model therefore has ",
          shiny::tags$strong(pvc_display$remaining_value),
          " of the null between-strata variance, rather than an explained-away share. "
        )
      } else if (identical(pvc_display$status, "unknown")) {
        shiny::tagList(
          "The proportional change in variance could not be summarized for this fit, ",
          "so the adjusted-model share of between-strata variance is not available. "
        )
      } else {
        shiny::tagList(
          "After adding the simple additive (main) effects, the between-strata variance is ",
          shiny::tags$strong(paste0(pvc_val, "%")),
          " smaller, leaving ",
          shiny::tags$strong(pvc_display$remaining_value),
          " of the original between-strata variance in the adjusted model. This is a ",
          "model-dependent change, not necessarily variance causally explained by those effects. "
        )
      }

      # Construct the summary paragraph dynamically
      shiny::tags$div(class = "alert alert-info mt-3",
        shiny::tags$strong("Automated Research Summary: "),
        "In this analysis, ", shiny::tags$strong(paste0(vpc_val, "%")),
        " of the (null-model) variance in the outcome lies between the defined intersecting demographic or social strata",
        " -- a between-stratum share of variance (on the model's latent scale for binary or count outcomes), not variance causally attributable to those strata. ",
        pvc_interpretation,
        "Among these residuals, the largest between-stratum departure from what the additive main effects alone would predict is in ",
        shiny::tags$strong(dev_label), " (N = ", dev_n, "), with an intersectional deviation score of ",
        shiny::tags$strong(dev_effect), ". Treat this as a descriptive screening flag rather than a confirmed disparity, especially for non-representative data, latent-scale models, or conditional intervals."
      )
    })

    output$interactive_table <- DT::renderDT({
      shiny::req(hud_plot_data())
      df <- hud_plot_data()

      # Drop tooltip and internal parsing columns before showing table
      cols_to_drop <- c("tooltip", "display_label")
      df <- df[, !names(df) %in% cols_to_drop]

      # Round numerics
      num_cols <- vapply(df, is.numeric, logical(1))
      df[num_cols] <- lapply(df[num_cols], round, 3)

      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
    })

    output$download_hud_data <- shiny::downloadHandler(
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

    hud_plot_data
  })
}
