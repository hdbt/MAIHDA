# Audit 2026-09-15: plot(type = "context_vpc") refused a crossed-dimensions fit that
# carries a context -- CONFIRMED.
#
# plot_context_vpc() accepted only a variance table tagged kind = "contextual", which
# the two-model contextual summary builds. A crossed-dimensions fit with context =
# builds its table as kind = "cross_classified" -- the same "Context: <var>" rows beside
# the additive dimension and interaction rows -- and carries the same $context
# partition, yet plot(type = "context_vpc") stopped with "No contextual partition is
# available. Fit the model with fit_maihda(context = )", naming the argument that had
# been given, and plot(type = "all") on the model and on the maihda() analysis warned
# that the panel "could not be computed and was omitted" for every such fit, lme4 and
# brms alike. plot(type = "vpc") drew the same table without trouble.
#
# FIX: plot_context_vpc() also takes a cross_classified table with at least one context
# row, drawing the between-stratum variance as its additive dimension and interaction
# bars in the colours plot_vpc() gives them (both now read
# maihda_cc_component_colors()). A fit without a context, crossed-dimensions or not, is
# refused with the same message as before. plot_vpc() and the two-model context_vpc
# view were identical() to the old ones on the fits checked.

d15_data <- function() {
  set.seed(762)
  dg <- expand.grid(a = factor(1:4), b = factor(1:4))
  d <- dg[rep(seq_len(nrow(dg)), each = 60), ]
  d$site <- factor(sample(1:16, nrow(d), TRUE))
  d$y <- stats::rnorm(16)[d$site] + stats::rnorm(nrow(d))
  d
}

d15_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# Fill colour of each bar, named by its component.
d15_fills <- function(p) {
  stats::setNames(ggplot2::layer_data(p, 1)$fill, as.character(p$data$component))
}

# The warnings an expression raises, muffled. plot(type = "all") prints every panel,
# so the assertions below look only for the warning that names the context_vpc
# panel, not for silence from panels this fix does not touch.
d15_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  w
}

test_that("context_vpc draws a crossed-dimensions fit that carries a context", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  a <- d15_quiet(maihda(y ~ 1 + (1 | a:b), data = d15_data(),
                        decomposition = "crossed-dimensions", context = "site"))
  vc <- a$summary$variance_components
  expect_identical(attr(vc, "kind"), "cross_classified")
  vc <- vc[vc$component != "Total", ]

  p <- plot(a, type = "context_vpc")
  expect_s3_class(p, "ggplot")
  expect_identical(levels(p$data$component), vc$component)
  expect_equal(ggplot2::layer_data(p, 1)$y, vc$variance, tolerance = 1e-12)
  expect_identical(ggplot2::layer_data(p, 2)$label, sprintf("%.1f%%", vc$proportion * 100))
  expect_identical(p$labels$title, sprintf("Stratum vs. Context Variance (VPC/ICC = %.3f)",
                                           a$summary$vpc$estimate))
  # Each component keeps the colour the stacked VPC bar gives it.
  fills <- d15_fills(p)
  expect_identical(fills, d15_fills(plot(a, type = "vpc"))[names(fills)])

  expect_s3_class(plot(a$model, type = "context_vpc", summary_obj = a$summary), "ggplot")
})

test_that("plot(type = 'all') draws context_vpc for a crossed-dimensions fit without dropping it", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  a <- d15_quiet(maihda(y ~ 1 + (1 | a:b), data = d15_data(),
                        decomposition = "crossed-dimensions", context = "site"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  w_model <- d15_warnings(m_plots <- plot(a$model, type = "all", summary_obj = a$summary))
  expect_identical(grep("context_vpc", w_model, fixed = TRUE, value = TRUE), character(0))
  expect_s3_class(m_plots$context_vpc, "ggplot")
  w_analysis <- d15_warnings(a_plots <- plot(a, type = "all"))
  expect_identical(grep("context_vpc", w_analysis, fixed = TRUE, value = TRUE), character(0))
  expect_s3_class(a_plots$context_vpc, "ggplot")
})

test_that("context_vpc still refuses a fit without a context", {
  # Negative control: passes before and after the fix.
  skip_on_cran()
  cc <- d15_quiet(maihda(y ~ 1 + (1 | a:b), data = d15_data(),
                         decomposition = "crossed-dimensions"))
  expect_error(plot(cc, type = "context_vpc"), "No contextual partition")
  expect_error(plot(cc$model, type = "context_vpc", summary_obj = cc$summary),
               "No contextual partition")
})

test_that("the two-model context_vpc view keeps its bars, colours and caption", {
  # Negative control: passes before and after the fix.
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  b <- d15_quiet(maihda(y ~ 1 + (1 | a:b), data = d15_data(), context = "site"))
  p <- plot(b, type = "context_vpc")
  expect_identical(d15_fills(p),
                   c("Between-stratum (random)" = "#E69F00", "Context: site" = "#117733",
                     "Within-stratum (residual)" = "#56B4E9"))
  expect_match(p$labels$caption, "^Contextual cross-classified MAIHDA:")
})

test_that("a crossed-dimensions table with several dimensions and contexts is drawn in the VPC colours", {
  skip_if_not_installed("ggplot2")
  # The table shape both engines build (maihda_cc_summary_lme4() / _brms()).
  s <- list(
    variance_components = maihda_cc_components_table(
      per_dim = c(gender = 0.3, race = 0.2, education = 0.1), interaction_var = 0.05,
      within_var = 1, per_context = c(school = 0.4, region = 0.2)),
    context = list(context_vars = c("school", "region")),
    vpc = list(estimate = 0.1)
  )
  p <- plot_context_vpc(s)
  expect_identical(d15_fills(p),
                   c("Additive: gender" = "#CC79A7", "Additive: race" = "#009E73",
                     "Additive: education" = "#0072B2",
                     "Intersectional interaction" = "#E69F00",
                     "Context: school" = "#117733", "Context: region" = "#44AA99",
                     "Within-stratum (residual)" = "#56B4E9"))
  expect_match(p$labels$caption, "^Crossed-dimensions contextual MAIHDA:")
  expect_identical(d15_fills(plot_vpc(s)), d15_fills(p))
})
