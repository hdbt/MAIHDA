# MAIHDA brand theme -----------------------------------------------------------
#
# This file carries the package's *brand identity* into plot chrome (titles,
# text, gridlines) so the figures match the hex logo, the pkgdown site, and the
# Shiny app. It deliberately does NOT touch any data-encoding colours: every
# scale_*_manual() value in the plot functions keeps the colourblind-safe
# Okabe-Ito / Paul Tol palettes. Brand colours and data colours are two
# different jobs and are kept strictly separate.

#' MAIHDA brand colours (internal)
#'
#' Brand-identity colours used only for plot *chrome* (titles, text,
#' gridlines) and never for encoding data. The navy matches the hex logo, the
#' pkgdown theme, and the Shiny app; the violet is the logo's "emergent
#' intersection" accent and is reserved for chrome.
#'
#' @keywords internal
#' @noRd
maihda_brand <- list(
  navy   = "#2C3E50", # primary brand colour (logo background, pkgdown, Shiny)
  slate  = "#5D6D7E", # secondary text
  grid   = "#D5DBE0", # subtle gridlines
  violet = "#7E33A8"  # logo "emergent intersection" accent (chrome only)
)

#' MAIHDA plot theme
#'
#' A \pkg{ggplot2} theme that applies the MAIHDA brand identity to plot
#' \emph{chrome} only -- navy titles and axis labels and soft slate gridlines.
#' It encodes no data: the colourblind-safe data palettes used by the package's
#' plotting functions are left untouched.
#'
#' Built on [ggplot2::theme_minimal()], so it composes with additional
#' `+ theme()` calls in the usual way (later settings win).
#'
#' @section Fonts:
#' By default the theme uses the graphics device's default font, which is safe
#' on every device (including the PostScript/PDF devices used by `R CMD check`
#' and many rendering back-ends). To render figures in the brand font
#' (Montserrat) -- matching the hex logo -- set `options(maihda.font =
#' "Montserrat")` or pass `base_family = "Montserrat"`, and use a graphics
#' device that can resolve that family (e.g. \pkg{ragg} or \pkg{showtext}).
#' Forcing a font that the active device cannot find produces "invalid font
#' type" errors, which is why it is opt-in rather than automatic.
#'
#' @param base_size Base font size in points. Default `11`.
#' @param base_family Base font family. Defaults to
#'   `getOption("maihda.font", "")`, i.e. the graphics device default unless you
#'   have opted into a brand font globally. Pass a string to force a family for
#'   a single plot.
#'
#' @return A \pkg{ggplot2} theme object that can be added to a plot with `+`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point(colour = "#0072B2") +
#'   labs(title = "MAIHDA brand theme", x = "MPG", y = "Weight") +
#'   theme_maihda()
theme_maihda <- function(base_size = 11, base_family = NULL) {
  if (is.null(base_family)) {
    base_family <- getOption("maihda.font", "")
  }

  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(colour = maihda_brand$navy, face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = maihda_brand$slate),
      plot.caption     = ggplot2::element_text(colour = maihda_brand$slate),
      axis.title       = ggplot2::element_text(colour = maihda_brand$navy),
      axis.text        = ggplot2::element_text(colour = maihda_brand$slate),
      legend.title     = ggplot2::element_text(colour = maihda_brand$navy),
      legend.text      = ggplot2::element_text(colour = maihda_brand$slate),
      strip.text       = ggplot2::element_text(colour = maihda_brand$navy, face = "bold"),
      panel.grid.major = ggplot2::element_line(colour = maihda_brand$grid),
      panel.grid.minor = ggplot2::element_line(colour = maihda_brand$grid, linewidth = 0.25),
      plot.title.position = "plot"
    )
}
