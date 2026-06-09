# Shiny modules for run_maihda_app().
#
# These live in R/ (rather than inline in inst/shiny/app.R) so they are unit
# testable with shiny::testServer() and covered by R CMD check. They run in the
# MAIHDA namespace, which does not attach shiny/bslib/DT/plotly, so every UI and
# server helper is fully qualified; the suggested packages are guaranteed present
# because run_maihda_app() gates on maihda_app_required_packages(). IDs created
# inside server-side renderUI() are namespaced with session$ns(), as Shiny does
# not auto-namespace dynamically generated UI.

# --- Shared empty-state card -------------------------------------------------
# A friendly placeholder shown on a results surface *before* a model is fitted,
# instead of a blank panel. Centralised here so every tab/module shows the same
# styling (.maihda-empty in inst/shiny/www/custom.css) and a consistent nudge
# back to the sidebar workflow.

#' @noRd
maihda_app_empty_state <- function(title, message,
                                   icon_name = "wand-magic-sparkles") {
  bslib::card(
    class = "maihda-empty text-center",
    bslib::card_body(
      shiny::icon(icon_name, class = "maihda-empty-icon"),
      shiny::h5(title, class = "mt-2"),
      shiny::div(class = "text-muted", shiny::markdown(message))
    )
  )
}

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
        shiny::downloadButton(ns("download_plot"), "Download Plot", class = "btn-secondary",
                              `aria-label` = "Download the current plot as a PNG image")
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
      if (is.null(model_results())) {
        return(maihda_app_empty_state(
          "No plots yet",
          "Fit a MAIHDA model from the sidebar, then pick a plot type above to
           visualise the strata effects."))
      }
      if (input$plot_type == "ternary") {
        shinycssloaders::withSpinner(plotly::plotlyOutput(ns("maihda_plotly"), height = "500px"))
      } else {
        shinycssloaders::withSpinner(shiny::plotOutput(ns("maihda_plot"), height = "500px"))
      }
    })

    output$maihda_plot <- shiny::renderPlot(
      current_plot(),
      alt = function() {
        paste0("MAIHDA ", if (!is.null(input$plot_type)) input$plot_type else "",
               " plot for the fitted model; see the accompanying tables for exact values.")
      }
    )

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
      if (is.null(model_results()) || is.null(null_summary_results()) ||
          is.null(summary_results()) || is.null(pvc_results())) {
        return(maihda_app_empty_state(
          "Nothing to explore yet",
          "Fit a MAIHDA model from the sidebar to unlock the interactive strata
           explorer, key metrics and filtered data export."))
      }
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
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            class = "maihda-metric-row",
            bslib::value_box(
              title = "VPC (Null)",
              value = paste0(vpc_val, "%"),
              showcase = shiny::icon("layer-group"),
              theme = "primary",
              shiny::p(class = "mb-0", "Total variance between strata"),
              if (!is.null(vpc_ci_text)) shiny::p(class = "small mb-0", vpc_ci_text) else NULL
            ),
            bslib::value_box(
              title = "PCV (Adjusted)",
              value = pvc_val_display,
              showcase = shiny::icon("arrow-down-wide-short"),
              theme = "info",
              shiny::p(class = "mb-0", "Between-stratum variance change with main effects")
            ),
            bslib::value_box(
              title = pvc_display$label,
              value = pvc_display$value,
              showcase = shiny::icon("chart-pie"),
              theme = switch(pvc_display$status,
                             negative = "warning", unknown = "secondary", "success"),
              shiny::p(class = "mb-0", pvc_display$description)
            )
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
          shiny::div(class = "mb-3", shiny::downloadButton(ns("download_hud_data"), "Download Highlighted Data (CSV)", class = "btn-secondary",
                              `aria-label` = "Download the highlighted strata as a CSV file")),
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

      # When colouring by conditional-interval status, encode it by shape as well
      # as colour so the signal is not conveyed by colour alone (colourblind-safe).
      deviant_mode <- identical(color_var, "deviant")
      point_aes <- if (deviant_mode) {
        ggplot2::aes(x = .data$random_effect, y = .data$display_label,
                     color = .data[[color_var]], shape = .data[[color_var]],
                     text = .data$tooltip)
      } else {
        ggplot2::aes(x = .data$random_effect, y = .data$display_label,
                     color = .data[[color_var]], text = .data$tooltip)
      }

      p <- ggplot2::ggplot(stratum_df, point_aes) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

      if (size_mapped) {
        p <- p + ggplot2::geom_point(ggplot2::aes(size = .data$n), alpha = 0.8) + ggplot2::scale_size_continuous(range = c(2, 6))
      } else {
        p <- p + ggplot2::geom_point(size = 3)
      }

      p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower_95, xmax = .data$upper_95), width = 0.2, orientation = "y") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Intersectional Intercept / Effect (Deviation)",
             y = "Stratum", color = tools::toTitleCase(color_var),
             shape = tools::toTitleCase(color_var), size = "Sample Size (N)") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

      # Colourblind-safe (Okabe-Ito) colours + distinct shapes for the deviant flag.
      if (deviant_mode) {
        p <- p +
          ggplot2::scale_color_manual(values = c("Excludes zero" = "#D55E00", "Includes zero" = "#0072B2")) +
          ggplot2::scale_shape_manual(values = c("Excludes zero" = 17, "Includes zero" = 16))
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

# --- Model Comparison tab ----------------------------------------------------
# Wires in two exported, tested comparison functions the app otherwise never
# surfaces: compare_maihda() (nested-model VPC forest plot, computed from the
# null + adjusted models already fitted) and compare_maihda_groups() (a full
# MAIHDA refitted within each level of a higher-level group). The latter is slow,
# so it has its own button + async future. `comparison_results` is a reactive
# carrying the precomputed maihda_comparison; `reactive_data` / `fit_params` /
# `fitted_family` describe the last fit so the stratified comparison can be built.

#' @noRd
mod_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    bslib::card(
      bslib::card_header("Nested models: VPC comparison (null vs adjusted)"),
      bslib::card_body(
        shiny::markdown(
          "Compares the between-stratum VPC of the **null** (strata-only) model with
          the **adjusted** model (main effects added) -- computed automatically from
          your last fit."),
        shiny::uiOutput(ns("nested_ui"))
      )
    ),
    bslib::card(
      bslib::card_header("Stratified MAIHDA: VPC by a higher-level group"),
      bslib::card_body(
        shiny::markdown(
          "Fits a *separate* MAIHDA within each level of a higher-level group
          (e.g. country, region) and compares the between-stratum VPC across groups.
          This refits a model per group, so it can be slow."),
        bslib::layout_columns(
          col_widths = c(5, 4, 3),
          shiny::selectInput(ns("group_var"), "Group by:", choices = NULL),
          shiny::selectInput(ns("plot_type"), "Plot:",
                             choices = c("VPC by group" = "vpc",
                                         "Variance components" = "components",
                                         "Between-stratum variance" = "between_variance")),
          shiny::div(class = "mt-3",
                     shiny::actionButton(ns("run_group"), "Run group comparison",
                                         class = "btn-primary"))
        ),
        shinycssloaders::withSpinner(shiny::plotOutput(ns("group_plot"), height = "420px")),
        DT::DTOutput(ns("group_table"))
      )
    )
  )
}

#' @noRd
mod_compare_server <- function(id, comparison_results, reactive_data, fit_params,
                               fitted_family) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Don't let a restored bookmark auto-launch the (slow) group comparison.
    shiny::setBookmarkExclude("run_group")

    round_numeric <- function(df, digits = 4) {
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], round, digits)
      df
    }

    # ---- Nested comparison (auto, from the last fit) ----
    # Swap a friendly empty state for the table/plot until a fit exists. The
    # nested_table / nested_plot renderers below stay defined unconditionally
    # (testServer accesses output$nested_table directly), so this only governs
    # what is *mounted* in the UI.
    output$nested_ui <- shiny::renderUI({
      if (is.null(comparison_results())) {
        return(maihda_app_empty_state(
          "No comparison yet",
          "Fit a MAIHDA model from the sidebar; the null-vs-adjusted VPC
           comparison is computed automatically from that fit."))
      }
      shiny::tagList(
        DT::DTOutput(ns("nested_table")),
        shinycssloaders::withSpinner(shiny::plotOutput(ns("nested_plot"), height = "320px"))
      )
    })

    output$nested_table <- DT::renderDT({
      shiny::req(comparison_results())
      DT::datatable(round_numeric(as.data.frame(comparison_results())),
                    options = list(dom = "t", paging = FALSE), rownames = FALSE)
    })

    output$nested_plot <- shiny::renderPlot(
      {
        shiny::req(comparison_results())
        plot(comparison_results())
      },
      alt = "Forest plot comparing the between-stratum VPC of the null and adjusted models."
    )

    # ---- Group-variable choices: data columns minus outcome + strata + covars ----
    shiny::observe({
      dat <- reactive_data()
      p <- fit_params()
      shiny::req(dat, p)
      used <- c(p$outcome, p$grouping_vars, p$covariates)
      shiny::updateSelectInput(session, "group_var", choices = setdiff(names(dat), used))
    })

    # ---- Stratified group comparison (own async fit) ----
    group_cmp <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run_group, {
      dat <- reactive_data()
      p <- fit_params()
      grp <- input$group_var
      shiny::req(dat, p, grp, length(p$grouping_vars) > 0)

      fam <- if (!is.null(fitted_family())) fitted_family() else "gaussian"
      autobin_opt <- isTRUE(p$autobin)

      # Build outcome ~ 1 + (1 | g1:g2:...); compare_maihda_groups() forms strata
      # itself. Detach the formula environment so the future does not serialise the
      # whole reactive context.
      rand <- paste(vapply(p$grouping_vars, maihda_quote_name, character(1)), collapse = ":")
      fml <- stats::as.formula(
        paste0(maihda_quote_name(p$outcome), " ~ 1 + (1 | ", rand, ")")
      )
      environment(fml) <- globalenv()

      group_cmp(NULL)
      shinyjs::disable(ns("run_group"))
      note <- shiny::showNotification(
        sprintf("Fitting a MAIHDA per level of '%s' (may take a while)...", grp),
        duration = NULL, type = "message")

      finish <- function() {
        shiny::removeNotification(note)
        shinyjs::enable(ns("run_group"))
      }

      promises::then(
        promises::future_promise({
          compare_maihda_groups(formula = fml, data = dat, group = grp,
                                family = fam, autobin = autobin_opt)
        }, seed = TRUE),
        onFulfilled = function(res) {
          finish()
          group_cmp(res)
        },
        onRejected = function(err) {
          finish()
          shiny::showNotification(
            paste("Group comparison failed:", conditionMessage(err)),
            type = "error", duration = 12)
        }
      )
    })

    output$group_table <- DT::renderDT({
      shiny::req(group_cmp())
      DT::datatable(round_numeric(as.data.frame(group_cmp())),
                    options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
    })

    output$group_plot <- shiny::renderPlot(
      {
        shiny::req(group_cmp())
        plot(group_cmp(), type = input$plot_type)
      },
      alt = function() {
        paste0("MAIHDA group-comparison plot (", input$plot_type,
               ") across levels of the selected grouping variable.")
      }
    )
  })
}
