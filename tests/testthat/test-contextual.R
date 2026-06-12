# Contextual cross-classified MAIHDA: individuals cross-classified by their
# intersectional stratum AND a higher-level context (school, hospital, region),
# outcome ~ covars + (1 | stratum) + (1 | context). This is the literature's
# "cross-classified MAIHDA"; the summary partitions the unexplained variance into
# between-stratum vs. between-context vs. residual. Distinct from the
# crossed-DIMENSIONS decomposition tested in test-cross-classified.R.

# Simulated data with a genuine stratum signal AND a genuine context signal. A
# many-level context (30 sites) keeps the context variance well identified, so the
# partition assertions are stable (the bundled 6-country data is too few-level for
# tight tests).
make_context_data <- function(seed = 8101, n = 1800, n_sites = 30) {
  set.seed(seed)
  d <- data.frame(
    g1 = sample(c("m", "f"), n, replace = TRUE),
    g2 = sample(c("low", "mid", "high"), n, replace = TRUE),
    site = sample(paste0("S", seq_len(n_sites)), n, replace = TRUE),
    region = sample(paste0("R", seq_len(10)), n, replace = TRUE),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$g1, d$g2, drop = TRUE)
  u_stratum <- rnorm(nlevels(stratum), sd = 1.0)[stratum]
  u_site <- stats::setNames(rnorm(n_sites, sd = 0.8), paste0("S", seq_len(n_sites)))
  d$y <- 1 + 0.4 * d$x + u_stratum + u_site[d$site] + rnorm(n, sd = 1.2)
  d
}

# ---- pure partition arithmetic ----------------------------------------------

test_that("maihda_context_partition computes the stratum/context/residual split", {
  part <- maihda_context_partition(2, list(site = 1, region = 1), 4)
  expect_equal(part$context_total, 2)
  expect_equal(part$total, 8)
  expect_equal(part$vpc_stratum, 0.25)
  expect_equal(part$vpc_context_total, 0.25)
  expect_equal(part$vpc_context$site, 1 / 8)

  # var_other enters the denominator but not the context total.
  part2 <- maihda_context_partition(2, c(site = 2), 4, var_other = 2)
  expect_equal(part2$total, 10)
  expect_equal(part2$vpc_stratum, 0.2)
  expect_equal(part2$vpc_context_total, 0.2)

  # Elementwise over draw vectors (the brms path).
  part3 <- maihda_context_partition(c(1, 2), list(site = c(1, 2)), c(2, 4))
  expect_equal(part3$vpc_stratum, c(0.25, 0.25))
})

test_that("maihda_context_components_table rows are labelled and sum to 1", {
  tab <- maihda_context_components_table(2, c(site = 1), 0, 4)
  expect_identical(attr(tab, "kind"), "contextual")
  expect_true("Between-stratum (random)" %in% tab$component)
  expect_true("Context: site" %in% tab$component)
  expect_true("Within-stratum (residual)" %in% tab$component)
  expect_false("Other random effects" %in% tab$component)
  non_total <- tab[tab$component != "Total", ]
  expect_equal(sum(non_total$proportion), 1, tolerance = 1e-8)

  tab2 <- maihda_context_components_table(2, c(site = 1), 0.5, 4)
  expect_true("Other random effects" %in% tab2$component)
})

test_that("maihda_validate_context rejects bad input", {
  d <- data.frame(y = 1, site = "a", stratum = "s")
  expect_null(maihda_validate_context(NULL, d))
  expect_identical(maihda_validate_context("site", d), "site")
  expect_error(maihda_validate_context(1, d), "character")
  expect_error(maihda_validate_context("nope", d), "not found")
  expect_error(maihda_validate_context("stratum", d), "stratum")
})

# ---- fit_maihda(context = ) --------------------------------------------------

test_that("fit_maihda(context = ) appends the context RE and tags the model", {
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  rhs <- paste(deparse(m$formula), collapse = " ")
  expect_true(grepl("1 | stratum", rhs, fixed = TRUE))
  expect_true(grepl("1 | site", rhs, fixed = TRUE))
  expect_identical(m$context_vars, "site")
  expect_identical(m$context_info$context_vars, "site")
})

test_that("fit_maihda(context = ) works with a pre-built stratum column", {
  d <- make_context_data()
  s <- make_strata(d, vars = c("g1", "g2"))
  m <- fit_maihda(y ~ x + (1 | stratum), data = s$data, context = "site")
  rhs <- paste(deparse(m$formula), collapse = " ")
  expect_true(grepl("1 | site", rhs, fixed = TRUE))
  expect_identical(m$context_vars, "site")
})

test_that("fit_maihda(context = ) is idempotent when the RE is already present", {
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  # Refit the derived formula (already carries (1 | site)) with context again --
  # the maihda() refit path. The RE must not be duplicated.
  m2 <- fit_maihda(m$formula, m$original_data, context = "site")
  rhs <- paste(deparse(m2$formula), collapse = " ")
  expect_identical(lengths(regmatches(rhs, gregexpr("1 | site", rhs, fixed = TRUE))), 1L)
  expect_identical(m2$context_vars, "site")
})

test_that("fit_maihda(context = ) validation errors", {
  d <- make_context_data()
  f <- y ~ x + (1 | g1:g2)
  expect_error(fit_maihda(f, d, context = "nope"), "not found")
  expect_error(fit_maihda(f, d, context = "g1"), "also define the intersectional strata")
  expect_error(fit_maihda(y ~ x + (1 | g1:g2), d, context = "x"),
               "fixed part")
  expect_error(fit_maihda(y ~ x, d, context = "site"), "no stratum random effect")
})

# ---- summary: the contextual partition ----------------------------------------

test_that("summary() partitions stratum vs. context vs. residual coherently", {
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  s <- summary(m)

  expect_false(is.null(s$context))
  expect_identical(s$context$context_vars, "site")
  expect_identical(attr(s$variance_components, "kind"), "contextual")
  expect_true("Context: site" %in% s$variance_components$component)

  # Proportions sum to 1 and the headline VPC is the stratum share.
  vc <- s$variance_components
  non_total <- vc[vc$component != "Total", ]
  expect_equal(sum(non_total$proportion), 1, tolerance = 1e-8)
  expect_equal(s$vpc$estimate, s$context$vpc_stratum, tolerance = 1e-10)
  expect_equal(s$context$vpc_stratum + s$context$vpc_context_total +
                 non_total$proportion[non_total$component == "Within-stratum (residual)"],
               1, tolerance = 1e-8)

  # With a genuine site signal the context share is comfortably positive.
  expect_gt(s$context$vpc_context_total, 0.05)
})

test_that("the contextual headline VPC equals the generic multi-RE VPC", {
  # The generic single-stratum path already puts extra RE variance in the VPC
  # denominator ("Other random effects"); the contextual path must agree on the
  # number and only improve the labelling. The raw multi-RE formula needs the
  # pre-built stratum column (auto-strata supports a single RE term only).
  d <- make_context_data()
  s <- make_strata(d, vars = c("g1", "g2"))
  m_ctx <- fit_maihda(y ~ x + (1 | stratum), data = s$data, context = "site")
  m_raw <- fit_maihda(y ~ x + (1 | stratum) + (1 | site), data = s$data)
  expect_null(m_raw$context_info)
  s_ctx <- summary(m_ctx)
  s_raw <- summary(m_raw)
  expect_equal(s_ctx$vpc$estimate, s_raw$vpc$estimate, tolerance = 1e-10)
  expect_true("Other random effects" %in% s_raw$variance_components$component)
  expect_null(s_raw$context)
})

test_that("the stratum VPC shrinks when the context is partitioned out", {
  d <- make_context_data()
  m0 <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  m1 <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  expect_lt(summary(m1)$vpc$estimate, summary(m0)$vpc$estimate)
})

test_that("two context variables each get their own component and share", {
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = c("site", "region"))
  s <- summary(m)
  expect_setequal(s$context$context_vars, c("site", "region"))
  expect_true(all(c("Context: site", "Context: region") %in%
                    s$variance_components$component))
  expect_equal(sum(s$context$vpc_context), s$context$vpc_context_total,
               tolerance = 1e-10)
  vc <- s$variance_components
  expect_equal(sum(vc$proportion[vc$component != "Total"]), 1, tolerance = 1e-8)
})

test_that("contextual bootstrap returns intervals for both shares", {
  d <- make_context_data(n = 900, n_sites = 20)
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  sb <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 25))
  expect_true(is.finite(sb$vpc$ci_lower) && is.finite(sb$vpc$ci_upper))
  expect_length(sb$context$vpc_context_total_ci, 2)
  expect_true(all(is.finite(sb$context$vpc_context_total_ci)))
  expect_true(sb$context$bootstrap)
})

test_that("print methods surface the contextual partition", {
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  expect_output(print(m), "Context: site")
  expect_output(print(summary(m)), "Contextual Cross-Classified Partition")
  expect_output(print(summary(m)), "Context 'site'")
})

# ---- maihda(context = ) --------------------------------------------------------

test_that("maihda(context = ) carries the context through null and adjusted", {
  d <- make_context_data()
  a <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site"))
  expect_identical(a$context_vars, "site")
  expect_identical(a$mode, "two-model")
  for (f in list(a$formula, a$adjusted_formula)) {
    expect_true(grepl("1 | site", paste(deparse(f), collapse = " "), fixed = TRUE))
  }
  expect_false(is.null(a$summary$context))
  expect_false(is.null(a$pcv))
  expect_output(print(a), "Context: site")
  expect_output(print(a), "Context share \\(null\\)")
})

test_that("the PCV is essentially unchanged by adding an orthogonal context", {
  # The context RE sits in both the null and adjusted models, and site membership
  # is independent of the strata by construction, so the PCV (a ratio of
  # between-stratum variances) should match the no-context PCV closely.
  d <- make_context_data()
  a0 <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d))
  a1 <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site"))
  expect_equal(a1$pcv$pvc, a0$pcv$pvc, tolerance = 0.05)
})

test_that("maihda(group, context) errors: they are different designs", {
  d <- make_context_data()
  expect_error(
    maihda(y ~ x + (1 | g1:g2), data = d, group = "region", context = "site"),
    "either 'group'.*or 'context'")
})

test_that("maihda(decomposition = 'crossed-dimensions', context = ) composes", {
  d <- make_context_data()
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | g1:g2), data = d,
           decomposition = "crossed-dimensions", context = "site")))
  expect_identical(cc$mode, "crossed-dimensions")
  rhs <- paste(deparse(cc$formula), collapse = " ")
  expect_true(grepl("1 | site", rhs, fixed = TRUE))
  expect_true(grepl("1 | stratum", rhs, fixed = TRUE))

  # The components table carries dimensions, interaction, AND the context, and the
  # VPC denominator includes the context variance.
  vc <- cc$summary$variance_components
  expect_true("Context: site" %in% vc$component)
  expect_equal(sum(vc$proportion[vc$component != "Total"]), 1, tolerance = 1e-8)
  expect_false(is.null(cc$summary$context))
  dcmp <- cc$decomposition
  expect_equal(cc$summary$vpc$estimate,
               dcmp$between_var / (dcmp$between_var + dcmp$within_var +
                                     cc$summary$context$context_var_total),
               tolerance = 1e-8)
})

# ---- plots ---------------------------------------------------------------------

test_that("contextual plots render", {
  skip_if_not_installed("ggplot2")
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d, context = "site")
  s <- summary(m)
  expect_s3_class(plot(m, type = "vpc", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "context_vpc", summary_obj = s), "ggplot")

  a <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site"))
  expect_s3_class(plot(a, type = "context_vpc"), "ggplot")
})

test_that("context_vpc errors without a contextual fit", {
  d <- make_context_data()
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  expect_error(plot(m, type = "context_vpc"), "No contextual partition")
  a <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d))
  expect_error(plot(a, type = "context_vpc"), "No contextual partition")
})

# ---- brms parity -----------------------------------------------------------------

test_that("contextual brms summary returns posterior shares with intervals", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). The draws-based
  # partition arithmetic is covered Stan-free above.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- make_context_data(seed = 8200, n = 900, n_sites = 20)
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, engine = "brms", context = "site",
               chains = 2, iter = 500, refresh = 0, seed = 1)))
  s <- summary(m)
  expect_identical(s$vpc$method, "posterior")
  expect_true(s$vpc$ci_lower < s$vpc$ci_upper)
  expect_identical(s$context$method, "posterior")
  expect_length(s$context$vpc_context_total_ci, 2)
  expect_true(all(is.finite(s$context$vpc_context_total_ci)))
  vc <- s$variance_components
  expect_true("Context: site" %in% vc$component)
  expect_equal(sum(vc$proportion[vc$component != "Total"]), 1, tolerance = 1e-6)
})
