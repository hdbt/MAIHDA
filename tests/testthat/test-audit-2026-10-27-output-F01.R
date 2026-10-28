# Audit 2026-10-27 (batch 2026-10-27-output, F01): the subtitle of the longitudinal
# plot(type = "trajectories") ran past the right edge of a 7 x 4.5 in figure, the size the
# longitudinal vignette draws it at. theme_maihda() aligns titles with the whole plot, so the
# subtitle's row is the figure's width less the plot margins, 6.85 in; drawn as one line at the
# theme's 11 pt it was 7.5 in wide and was cut at "referenc". The sentence now breaks after its
# semicolon, and the longer line is 4.8 in.
#
# Nothing here draws random numbers.

# Width of a plot's drawn subtitle as a share of the row it is drawn in, read from the plot's
# gtable on a width x height in pdf device, whose Helvetica metrics are the same on every
# platform.
subtitle_share <- function(p, width = 7, height = 4.5) {
  old <- grDevices::dev.cur()
  grDevices::pdf(NULL, width = width, height = height)
  on.exit({
    grDevices::dev.off()
    if (old > 1) grDevices::dev.set(old)
  })
  g <- ggplot2::ggplotGrob(p)
  i <- which(g$layout$name == "subtitle")
  ws <- grid::convertWidth(g$widths, "in", valueOnly = TRUE)
  cols <- seq_along(ws)
  outside <- sum(ws[cols < g$layout$l[i] | cols > g$layout$r[i]])
  grid::convertWidth(grid::grobWidth(g$grobs[[i]]), "in", valueOnly = TRUE) / (width - outside)
}

test_that("the trajectories subtitle fits a 7 x 4.5 in figure", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA", envir = environment())
  m <- suppressWarnings(suppressMessages(
    fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
               data = maihda_long_data, id = "id", time = "wave")))
  p <- plot(m, type = "trajectories")
  # The vignette's twelve strata keep the legend on the right.
  expect_identical(p$theme$legend.position, "right")
  # Under 0.85 in Helvetica, so it also fits in a wider system font on a png device: Verdana
  # sets this subtitle 16% wider than Helvetica.
  expect_lt(subtitle_share(p), 0.85)
  # The measure sees a subtitle that is too long.
  expect_gt(subtitle_share(p + ggplot2::labs(subtitle = strrep("x", 150))), 1)
})
