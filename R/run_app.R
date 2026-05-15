maihda_app_required_packages <- function() {
  c("shiny", "bslib", "DT", "future", "promises", "haven", "shinyjs", "plotly", "ggtern")
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
