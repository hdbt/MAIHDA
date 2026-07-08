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

  # A degenerate all-zero partition has no defined proportions (0/0), not zeros.
  tab0 <- maihda_context_components_table(0, c(site = 0), 0, 0)
  expect_true(all(is.na(tab0$proportion[tab0$component != "Total"])))
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

test_that("summary() of a binomial contextual fit attaches no discriminatory accuracy or response VPC", {
  # The DA's AUC is built from full predictions that INCLUDE the context random
  # effects, while the response-scale VPC simulates only the stratum variance -- so
  # neither matches the stratum-vs-context partition the contextual summary reports.
  # Both are therefore skipped for a contextual fit (as they are for crossed-dimensions
  # and longitudinal fits), rather than pinning a mismatched estimand to the partition.
  set.seed(8131)
  d <- make_context_data(n = 2000, n_sites = 25)
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(d$y - mean(d$y)))

  m <- suppressWarnings(suppressMessages(
    fit_maihda(yb ~ x + (1 | g1:g2), data = d, context = "site", family = "binomial")
  ))
  expect_false(is.null(m$context_info))

  s <- suppressWarnings(suppressMessages(summary(m, response_vpc = TRUE, seed = 1)))
  # The contextual partition is still produced...
  expect_false(is.null(s$context))
  # ...but the binomial companions that would carry a mismatched estimand are not.
  expect_null(s$discriminatory_accuracy)
  expect_null(s$vpc_response)

  # A single-stratum binomial fit on the same data DOES surface the DA -- confirming
  # the skip is specific to the contextual structure, not the family.
  m_plain <- suppressWarnings(suppressMessages(
    fit_maihda(yb ~ x + (1 | g1:g2), data = d, family = "binomial")
  ))
  expect_s3_class(summary(m_plain)$discriminatory_accuracy, "maihda_da")
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

test_that("maihda_context_summary_lme4 errors when the context RE is absent", {
  # A model tagged for a context whose random effect is not in the fit (only
  # reachable through internal misuse, but the guard must name the missing RE).
  d <- make_context_data(n = 600, n_sites = 10)
  m <- fit_maihda(y ~ x + (1 | g1:g2), data = d)
  expect_error(
    maihda_context_summary_lme4(m, list(context_vars = "site"),
                                lme4::VarCorr(m$model), FALSE, 10, 0.95),
    "missing the random effect")
})

test_that("an explicit family is threaded through maihda(context = )", {
  d <- make_context_data(n = 900, n_sites = 15)
  a <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site",
                               family = "gaussian"))
  expect_identical(a$context_vars, "site")
  expect_false(is.null(a$summary$context))
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

test_that("plot(type = 'vpc', model = ) works for a contextual analysis", {
  d <- make_context_data()
  a <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site"))

  # model = "both" -> one change plot carrying BOTH models, the broken-out context
  # slice, and the PCV in the subtitle.
  p_both <- plot(a, type = "vpc", model = "both")
  expect_s3_class(p_both, "ggplot")
  expect_false(inherits(p_both, "patchwork"))
  expect_true("Context: site" %in% as.character(p_both$data$component))
  expect_setequal(as.character(unique(p_both$data$model)),
                  c("Null model", "Adjusted model"))
  expect_match(p_both$labels$subtitle, "PCV")

  # The null and adjusted single views keep the contextual bar and their labels.
  expect_identical(plot(a, type = "vpc")$labels$subtitle, "Null model")
  expect_identical(plot(a, type = "vpc", model = "adjusted")$labels$subtitle,
                   "Adjusted model")
})

test_that("the PCV is essentially unchanged by adding an orthogonal context", {
  # The context RE sits in both the null and adjusted models, and site membership
  # is independent of the strata by construction, so the PCV (a ratio of
  # between-stratum variances) should match the no-context PCV closely.
  d <- make_context_data()
  a0 <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d))
  a1 <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site"))
  expect_equal(a1$pcv$pcv, a0$pcv$pcv, tolerance = 0.05)
})

# ---- group + context compose (stratified x contextual) -----------------------

test_that("maihda(group, context) composes: per-group contextual fits", {
  d <- make_context_data()
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | g1:g2), data = d, group = "region", context = "site")))
  expect_s3_class(a, "maihda_analysis")
  expect_false(is.null(a$groups))
  # The per-group table carries the contextual partition columns and attribute.
  expect_true(all(c("var_context", "vpc_context") %in% names(a$groups)))
  expect_identical(attr(a$groups, "context_var"), "site")
  expect_false(is.null(attr(a$groups, "context_per")))
  # Net-of-context share is a valid proportion for every fitted group.
  ok <- a$groups[a$groups$status == "ok", , drop = FALSE]
  expect_gt(nrow(ok), 0)
  expect_true(all(is.finite(ok$vpc_context)))
  expect_true(all(ok$vpc_context >= 0 & ok$vpc_context <= 1))
  # The headline analysis still prints and reports the group comparison.
  expect_output(print(a), "Group comparison")
})

test_that("compare_maihda_groups(context =) reports the per-group context partition", {
  d <- make_context_data()
  cmp <- suppressWarnings(compare_maihda_groups(
    y ~ x + (1 | g1:g2), data = d, group = "region", context = "site"))
  expect_true(all(c("var_context", "vpc_context") %in% names(cmp)))
  per <- attr(cmp, "context_per")
  expect_false(is.null(per))
  ok <- cmp[cmp$status == "ok", , drop = FALSE]
  g1 <- as.character(ok$group[1])
  # The summed per-context split equals the reported (summed) var_context.
  expect_equal(sum(per[[g1]]), ok$var_context[as.character(ok$group) == g1],
               tolerance = 1e-6)
  expect_output(print(cmp), "Context: site")
})

test_that("per-group PCV is essentially unchanged by an orthogonal context", {
  # site membership is independent of the strata by construction, so adding the
  # (1 | site) context RE to each per-group fit should barely move the between-
  # stratum PCV (which reads only the stratum variance component).
  d <- make_context_data()
  base <- suppressWarnings(as.data.frame(compare_maihda_groups(
    y ~ x + (1 | g1:g2), data = d, group = "region")))
  ctx <- suppressWarnings(as.data.frame(compare_maihda_groups(
    y ~ x + (1 | g1:g2), data = d, group = "region", context = "site")))
  m <- merge(base[, c("group", "pcv")], ctx[, c("group", "pcv")],
             by = "group", suffixes = c("_base", "_ctx"))
  ok <- is.finite(m$pcv_base) & is.finite(m$pcv_ctx)
  expect_gt(sum(ok), 0)
  expect_equal(mean(m$pcv_ctx[ok]), mean(m$pcv_base[ok]), tolerance = 0.1)
})

test_that("group + context is rejected up front for wemix / ordinal", {
  d <- make_context_data()
  expect_error(
    maihda(y ~ x + (1 | g1:g2), data = d, group = "region", context = "site",
           engine = "wemix"),
    "does not support 'context'")
  expect_error(
    maihda(y ~ x + (1 | g1:g2), data = d, group = "region", context = "site",
           engine = "ordinal"),
    "does not support 'context'")
  expect_error(
    compare_maihda_groups(y ~ x + (1 | g1:g2), data = d, group = "region",
                          context = "site", engine = "wemix"),
    "does not support 'context'")
})

test_that("context may not name the group variable", {
  d <- make_context_data()
  expect_error(
    compare_maihda_groups(y ~ x + (1 | g1:g2), data = d, group = "region",
                          context = "region"),
    "cannot include the group variable")
})

test_that("group + context warns when a group has too few context levels", {
  set.seed(11)
  n <- 600
  region <- sample(c("A", "B"), n, replace = TRUE)
  # Region A sees only 3 sites (weakly identified); region B sees 20.
  site <- ifelse(region == "A",
                 sample(paste0("a", 1:3), n, replace = TRUE),
                 sample(paste0("b", 1:20), n, replace = TRUE))
  d <- data.frame(
    g1 = sample(c("m", "f"), n, replace = TRUE),
    g2 = sample(c("lo", "hi"), n, replace = TRUE),
    region = region, site = site, x = rnorm(n),
    stringsAsFactors = FALSE
  )
  d$y <- 0.3 * d$x + rnorm(n)
  # Other (singular-fit) warnings may fire too; capture all and assert the
  # weak-identification one names region A x site.
  w <- capture_warnings(
    compare_maihda_groups(y ~ x + (1 | g1:g2), data = d, group = "region",
                          context = "site"))
  expect_true(any(grepl("few context levels", w)))
  expect_true(any(grepl("A x site", w)))
})

test_that("plot(type = 'components') gains a context slice that still sums to 1", {
  d <- make_context_data()
  cmp <- suppressWarnings(compare_maihda_groups(
    y ~ x + (1 | g1:g2), data = d, group = "region", context = "site"))
  p <- plot(cmp, type = "components")
  expect_s3_class(p, "ggplot")
  expect_true("Context (higher level)" %in% as.character(p$data$component))
  agg <- tapply(p$data$proportion, p$data$group, sum)
  expect_true(all(abs(agg - 1) < 1e-6))
})

test_that("a non-contextual group comparison keeps no context columns", {
  d <- make_context_data()
  cmp <- suppressWarnings(compare_maihda_groups(
    y ~ x + (1 | g1:g2), data = d, group = "region"))
  expect_false(any(c("var_context", "vpc_context") %in% names(cmp)))
  expect_null(attr(cmp, "context_var"))
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

test_that("crossed-dimensions + context: print, vpc plot, and bootstrap intervals", {
  d <- make_context_data(n = 900, n_sites = 15)
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | g1:g2), data = d,
           decomposition = "crossed-dimensions", context = "site")))

  expect_output(print(cc), "crossed-dimensions \\(single model\\)")
  expect_output(print(cc), "Context:")
  expect_output(print(cc), "Contextual Cross-Classified Partition")
  expect_output(print(cc$summary), "Contextual Cross-Classified Partition")

  skip_if_not_installed("ggplot2")
  # The cc VPC bar gains green context slice(s).
  p <- plot(cc$model, type = "vpc", summary_obj = cc$summary)
  expect_s3_class(p, "ggplot")

  # The cc parametric bootstrap also intervals the context share.
  sb <- suppressWarnings(summary(cc$model, bootstrap = TRUE, n_boot = 15))
  expect_true(sb$context$bootstrap)
  expect_length(sb$context$vpc_context_total_ci, 2)
})

# ---- plots ---------------------------------------------------------------------

test_that("plot(type = 'all') includes the context_vpc view for a contextual fit", {
  skip_if_not_installed("ggplot2")
  d <- make_context_data(n = 900, n_sites = 15)
  a <- suppressMessages(maihda(y ~ x + (1 | g1:g2), data = d, context = "site"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  plots <- suppressWarnings(plot(a, type = "all"))
  expect_s3_class(plots$context_vpc, "ggplot")
  m_plots <- suppressWarnings(plot(a$model, summary_obj = a$summary))
  expect_s3_class(m_plots$context_vpc, "ggplot")
})

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

# ---- brms contextual summary (Stan-free) -----------------------------------------
# maihda_context_summary_brms() is a pure function of the posterior draws once
# maihda_posterior_draws_brms() is mocked: everything downstream
# (maihda_group_variance_draws_brms, maihda_residual_variance_draws_brms on a
# gaussian stub, the partition and table builders) needs no fitted Stan model.
# Mirrors the brmsfit-stub pattern in test-summary_variance.R.

make_context_brms_stub <- function() {
  structure(
    list(formula = y ~ x + (1 | stratum) + (1 | site),
         family = list(family = "gaussian", link = "identity")),
    class = "brmsfit"
  )
}

test_that("maihda_context_summary_brms partitions per draw (Stan-free)", {
  set.seed(4711)
  n <- 500
  draws <- data.frame(
    sd_stratum__Intercept = sqrt(1.0 + 0.2 * runif(n)),
    sd_site__Intercept    = sqrt(0.8 + 0.2 * runif(n)),
    sigma                 = sqrt(1.4 + 0.2 * runif(n))
  )
  object <- list(model = make_context_brms_stub(), engine = "brms")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)

  res <- maihda_context_summary_brms(object, list(context_vars = "site"),
                                     conf_level = 0.9)

  v_s <- draws$sd_stratum__Intercept^2
  v_c <- draws$sd_site__Intercept^2
  v_e <- draws$sigma^2
  vpc_draws <- v_s / (v_s + v_c + v_e)
  expect_equal(res$vpc_result$estimate, stats::median(vpc_draws))
  expect_equal(res$vpc_result$ci_lower,
               stats::quantile(vpc_draws, 0.05, names = FALSE))
  expect_equal(res$vpc_result$ci_upper,
               stats::quantile(vpc_draws, 0.95, names = FALSE))
  expect_identical(res$vpc_result$method, "posterior")

  ctx_draws <- v_c / (v_s + v_c + v_e)
  expect_equal(res$context$vpc_context_total, stats::median(ctx_draws))
  expect_equal(res$context$vpc_context_total_ci,
               stats::quantile(ctx_draws, c(0.05, 0.95), names = FALSE))
  expect_identical(res$context$method, "posterior")
  expect_equal(res$context$var_stratum, mean(v_s))
  expect_equal(unname(res$context$per_context["site"]), mean(v_c))
  expect_equal(res$context$within_var, mean(v_e))

  tab <- res$variance_components
  expect_identical(attr(tab, "kind"), "contextual")
  expect_true("Context: site" %in% tab$component)
  expect_false("Other random effects" %in% tab$component)
  expect_equal(sum(tab$proportion[tab$component != "Total"]), 1, tolerance = 1e-8)
})

test_that("maihda_context_summary_brms keeps extra REs in the denominator; point = 'mean'", {
  set.seed(4712)
  n <- 300
  draws <- data.frame(
    sd_stratum__Intercept = sqrt(1.0 + 0.2 * runif(n)),
    sd_site__Intercept    = sqrt(0.8 + 0.2 * runif(n)),
    sd_extra__Intercept   = sqrt(0.3 + 0.2 * runif(n)),
    sigma                 = sqrt(1.4 + 0.2 * runif(n))
  )
  object <- list(model = make_context_brms_stub(), engine = "brms")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)

  res <- maihda_context_summary_brms(object, list(context_vars = "site"),
                                     conf_level = 0.95, point = "mean")

  v_s <- draws$sd_stratum__Intercept^2
  v_c <- draws$sd_site__Intercept^2
  v_x <- draws$sd_extra__Intercept^2
  v_e <- draws$sigma^2
  expect_equal(res$vpc_result$estimate, mean(v_s / (v_s + v_c + v_x + v_e)))
  expect_equal(res$context$other_var, mean(v_x))
  expect_true("Other random effects" %in% res$variance_components$component)
  tab <- res$variance_components
  expect_equal(sum(tab$proportion[tab$component != "Total"]), 1, tolerance = 1e-8)
})

test_that("maihda_context_summary_brms reports NA when no draw is finite", {
  # All-zero variances make every per-draw share 0/0 = NaN; the summariser must
  # come back NA rather than erroring.
  draws <- data.frame(
    sd_stratum__Intercept = c(0, 0),
    sd_site__Intercept = c(0, 0),
    sigma = c(0, 0)
  )
  object <- list(model = make_context_brms_stub(), engine = "brms")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)
  res <- maihda_context_summary_brms(object, list(context_vars = "site"), 0.95)
  expect_true(is.na(res$vpc_result$estimate))
  expect_true(is.na(res$context$vpc_context_total))
})

test_that("maihda_context_summary_brms errors when an RE is missing from the draws", {
  draws <- data.frame(sd_stratum__Intercept = c(1, 1.1), sigma = c(1.2, 1.3))
  object <- list(model = make_context_brms_stub(), engine = "brms")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)
  expect_error(
    maihda_context_summary_brms(object, list(context_vars = "site"), 0.95),
    "missing the random effect")
})

# ---- print branches of the contextual partition ----------------------------------

test_that("maihda_print_context_partition prints multi-context and interval branches", {
  shared <- list(
    var_stratum = 1, within_var = 2, other_var = 0, vpc_stratum = 1 / 4.2,
    bootstrap = FALSE
  )
  ctx_multi <- c(shared, list(
    context_vars = c("site", "region"),
    per_context = c(site = 0.8, region = 0.4),
    context_var_total = 1.2,
    vpc_context = c(site = 0.8 / 4.2, region = 0.4 / 4.2),
    vpc_context_total = 1.2 / 4.2,
    vpc_context_total_ci = c(0.2, 0.35)
  ))
  expect_output(maihda_print_context_partition(ctx_multi), "All contexts combined")
  expect_output(maihda_print_context_partition(ctx_multi), "Context 'region'")

  ctx_single_ci <- c(shared, list(
    context_vars = "site",
    per_context = c(site = 0.8),
    context_var_total = 0.8,
    vpc_context = c(site = 0.8 / 3.8),
    vpc_context_total = 0.8 / 3.8,
    vpc_context_total_ci = c(0.15, 0.3)
  ))
  expect_output(maihda_print_context_partition(ctx_single_ci),
                "Context share interval")
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

# ---- stepwise_pcv(context = ) --------------------------------------------------
# The stepwise between-stratum PCV as a contextual cross-classified model: the
# (1 | context) intercept is held in the null model and every step, so Step_PCV /
# Total_PCV are the between-stratum PCV NET OF the context. Every step is fit through
# fit_maihda(context = ), and extract_between_variance() already reads only the
# stratum component, so the isolation is automatic; these tests pin the forwarding,
# the shared complete-case sample, the Context_Variance column, and the engine guards.

test_that("stepwise_pcv(context = ) isolates the between-stratum PCV net of context (lme4 gaussian)", {
  d <- make_context_data()
  s <- make_strata(d, vars = c("g1", "g2"))

  out <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "y", c("g1", "g2", "x"), context = "site")))
  expect_s3_class(out, "maihda_stepwise")
  expect_true(all(is.finite(out$Total_PCV)))
  expect_true("Context_Variance" %in% names(out))

  # The final step's between-stratum Variance equals the direct contextual fit's
  # (ML-refit) stratum variance -- proving the (1 | site) effect is held every step
  # and the stratum variance is read net of it.
  m_ctx <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ g1 + g2 + x + (1 | stratum), data = s$data, context = "site")))
  v_ctx <- extract_between_variance(maihda_pcv_refit_ml(m_ctx))
  expect_equal(out$Variance[nrow(out)], v_ctx, tolerance = 1e-6)

  # Context genuinely matters for this data: the same fit WITHOUT the context has a
  # materially different stratum variance, so the match above is discriminating.
  m_noctx <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ g1 + g2 + x + (1 | stratum), data = s$data)))
  v_noctx <- extract_between_variance(maihda_pcv_refit_ml(m_noctx))
  expect_false(isTRUE(all.equal(v_ctx, v_noctx, tolerance = 1e-4)))
})

test_that("stepwise_pcv Context_Variance reports the between-context variance and leaves other tables unchanged", {
  d <- make_context_data()
  s <- make_strata(d, vars = c("g1", "g2"))

  out    <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "y", c("g1", "g2", "x"), context = "site")))
  out_no <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "y", c("g1", "g2", "x"))))

  # Present only with context =; the non-contextual gaussian table is byte-for-byte
  # the historical six columns.
  expect_true("Context_Variance" %in% names(out))
  expect_identical(names(out_no),
                   c("Step", "Model", "Added_Variable", "Variance",
                     "Step_PCV", "Total_PCV"))
  # ...and Context_Variance sits directly after the stratum Variance.
  expect_equal(which(names(out) == "Context_Variance"),
               which(names(out) == "Variance") + 1L)

  # Positive at every step and equal to the fitted site variance at the final step.
  expect_true(all(out$Context_Variance > 0))
  m_ctx <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ g1 + g2 + x + (1 | stratum), data = s$data, context = "site")))
  site_var <- unname(maihda_random_variances_lme4(
    maihda_pcv_refit_ml(m_ctx)$model)["site"])
  expect_equal(out$Context_Variance[nrow(out)], site_var, tolerance = 1e-6)
})

test_that("stepwise_pcv(context = ) filters context to one shared complete-case sample", {
  d <- make_context_data()
  d$site[1:50] <- NA
  s <- make_strata(d, vars = c("g1", "g2"))

  out <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "y", c("g1", "g2", "x"), context = "site")))

  # The 50 NA-site rows are dropped once, up front, so every step (here the null)
  # uses the same reduced analytic sample as a manual complete-site fit.
  s_cc <- s$data[!is.na(s$data$site), , drop = FALSE]
  m_null <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ (1 | stratum), data = s_cc, context = "site")))
  expect_equal(out$Variance[1],
               extract_between_variance(maihda_pcv_refit_ml(m_null)),
               tolerance = 1e-6)
})

test_that("stepwise_pcv omits the AUC/MOR trajectory for a binary contextual fit", {
  set.seed(4242)
  d <- make_context_data(n = 1500, n_sites = 20)
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(d$y - mean(d$y)))
  s <- make_strata(d, vars = c("g1", "g2"))

  # No discriminatory-accuracy columns (the AUC would include the context effect),
  # but the net-of-context PCV trajectory and the context variance are reported.
  outb <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "yb", c("g1", "g2"), family = "binomial",
                 context = "site")))
  expect_false(any(c("AUC", "Step_AUC", "Total_AUC", "MOR") %in% names(outb)))
  expect_true("Context_Variance" %in% names(outb))
  expect_true(all(is.finite(outb$Total_PCV)))

  # The same binary stepwise WITHOUT context still carries the AUC trajectory --
  # confirming the omission is specific to the contextual structure, not the family.
  outb0 <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "yb", c("g1", "g2"), family = "binomial")))
  expect_true(all(c("AUC", "MOR") %in% names(outb0)))
})

test_that("stepwise_pcv(context = ) rejects the wemix and ordinal engines", {
  d <- make_context_data()
  d$w <- stats::runif(nrow(d), 0.5, 2)
  s <- make_strata(d, vars = c("g1", "g2"))

  # sampling_weights routes to wemix, which fits no crossed random effects.
  expect_error(
    suppressMessages(stepwise_pcv(s$data, "y", c("g1", "g2"), context = "site",
                                  sampling_weights = "w")),
    "does not support 'context'")
  # An ordinal family routes to the clmm engine, likewise no crossed REs.
  expect_error(
    suppressMessages(stepwise_pcv(s$data, "y", c("g1", "g2"), context = "site",
                                  family = "ordinal")),
    "does not support 'context'")
})

test_that("stepwise_pcv(context = ) works on the brms engine", {
  # Compiles a Stan model per step, so OPT-IN (set MAIHDA_TEST_BRMS=true).
  # stepwise_pcv() takes no `...`, so the per-step brms fits use default sampling;
  # a two-fit trajectory (null + one step) keeps the opt-in run bounded.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- make_context_data(seed = 8300, n = 400, n_sites = 12)
  s <- make_strata(d, vars = c("g1", "g2"))
  out <- suppressWarnings(suppressMessages(
    stepwise_pcv(s$data, "y", c("g1"), engine = "brms", context = "site")))
  expect_s3_class(out, "maihda_stepwise")
  expect_true("Context_Variance" %in% names(out))
  expect_true(all(is.finite(out$Context_Variance)))
})
