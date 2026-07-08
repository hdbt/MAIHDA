# Cross-classified MAIHDA: the single-model additive/interaction decomposition
# (each dimension's additive main effect as a random intercept, the intersection RE
# as the interaction). Data with genuine additive main effects AND an interaction so
# the additive share is strictly inside (0, 1).
make_cc_data <- function(seed = 7001, n = 1500) {
  set.seed(seed)
  d <- data.frame(
    a = sample(c("a1", "a2", "a3"), n, replace = TRUE),
    b = sample(c("b1", "b2", "b3"), n, replace = TRUE),
    cc = sample(c("c1", "c2", "c3"), n, replace = TRUE),
    grp = sample(c("G1", "G2"), n, replace = TRUE),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )
  ua <- stats::setNames(rnorm(3, sd = 0.9), c("a1", "a2", "a3"))
  ub <- stats::setNames(rnorm(3, sd = 0.7), c("b1", "b2", "b3"))
  uc <- stats::setNames(rnorm(3, sd = 0.5), c("c1", "c2", "c3"))
  stratum <- interaction(d$a, d$b, d$cc, drop = TRUE)
  uint <- rnorm(nlevels(stratum), sd = 0.8)[stratum]
  d$y <- 2 + 0.3 * d$x + ua[d$a] + ub[d$b] + uc[d$cc] + uint + rnorm(n, sd = 1)
  d
}

# ---- pure partition arithmetic ---------------------------------------------

test_that("maihda_cc_variance_split and maihda_cc_partition compute the partition", {
  v <- c(a = 1, b = 2, stratum = 3, Residual = 4)
  sp <- maihda_cc_variance_split(v, c(A = "a", B = "b"), "stratum")
  expect_equal(sp$additive, 3)
  expect_equal(sp$interaction, 3)
  expect_equal(unname(sp$per_dim), c(1, 2))
  expect_equal(names(sp$per_dim), c("A", "B"))

  part <- maihda_cc_partition(sp$additive, sp$interaction, 4)
  expect_equal(part$between, 6)
  expect_equal(part$vpc, 6 / 10)
  expect_equal(part$additive_share, 0.5)
  expect_equal(part$interaction_share, 0.5)
})

test_that("maihda_cc_partition returns NA (not NaN) shares with zero between variance", {
  # A degenerate fit with no between-strata variance: the additive / interaction
  # shares split that (zero) variance, so they are undefined -- 0 / 0 = NaN. The
  # partition must report NA, not leak NaN into the public summary.
  part <- maihda_cc_partition(0, 0, 1)
  expect_equal(part$between, 0)
  expect_true(is.na(part$additive_share))
  expect_true(is.na(part$interaction_share))
  expect_false(is.nan(part$additive_share))    # NA, specifically not NaN
  expect_false(is.nan(part$interaction_share))
  expect_equal(part$vpc, 0)                     # between / total = 0 / 1 is well defined

  # When the TOTAL variance is also zero the VPC is undefined too -> NA, not NaN.
  part0 <- maihda_cc_partition(0, 0, 0)
  expect_true(is.na(part0$vpc))
  expect_false(is.nan(part0$vpc))

  # A non-zero between variance is unchanged (regression guard for the ifelse).
  part2 <- maihda_cc_partition(2, 1, 1)
  expect_equal(part2$additive_share, 2 / 3)
  expect_equal(part2$interaction_share, 1 / 3)
  expect_equal(part2$vpc, 3 / 4)
})

test_that("maihda_cc_variance_split errors when a random effect is missing", {
  v <- c(a = 1, stratum = 3, Residual = 4)
  expect_error(maihda_cc_variance_split(v, c(A = "a", B = "b"), "stratum"),
               "missing")
})

test_that("maihda_group_variance_draws_brms splits per group, underscored names ok", {
  draws <- data.frame(
    sd_Gender__Intercept = c(0.2, 0.4),
    sd_.maihda_dim_age__Intercept = c(0.1, 0.3),
    sd_stratum__Intercept = c(0.5, 0.5),
    sigma = c(1, 1),
    check.names = FALSE
  )
  gv <- maihda_group_variance_draws_brms(draws)
  expect_setequal(names(gv), c("Gender", ".maihda_dim_age", "stratum"))
  expect_equal(gv$Gender, c(0.04, 0.16))
  expect_equal(gv[[".maihda_dim_age"]], c(0.01, 0.09))
  expect_equal(gv$stratum, c(0.25, 0.25))
})

# ---- formula builder --------------------------------------------------------

test_that("maihda_cross_classified_formula builds dimension REs + intersection RE", {
  d <- make_cc_data()
  s <- make_strata(d, vars = c("a", "b", "cc"))
  cc <- maihda_cross_classified_formula(y ~ x + (1 | stratum), c("a", "b", "cc"),
                                        list(), s$data)
  expect_equal(cc$interaction_group, "stratum")
  expect_setequal(unname(cc$dim_groups), c("a", "b", "cc"))
  rhs <- paste(deparse(cc$formula), collapse = " ")
  expect_true(grepl("1 | a", rhs, fixed = TRUE))
  expect_true(grepl("1 | b", rhs, fixed = TRUE))
  expect_true(grepl("1 | stratum", rhs, fixed = TRUE))
  # Fewer than two dimensions -> no decomposition.
  expect_null(maihda_cross_classified_formula(y ~ x + (1 | stratum), c("a"),
                                              list(), s$data))
})

# ---- maihda() crossed-dimensions mode -----------------------------------------

test_that("maihda(decomposition = 'crossed-dimensions') returns a coherent partition", {
  d <- make_cc_data()
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, decomposition = "crossed-dimensions")))

  expect_s3_class(cc, "maihda_analysis")
  expect_identical(cc$mode, "crossed-dimensions")
  expect_null(cc$pcv)
  expect_null(cc$model_adjusted)
  expect_s3_class(cc$model, "maihda_model")
  expect_false(is.null(cc$model$cc_info))

  dcmp <- cc$decomposition
  expect_false(is.null(dcmp))
  # Invariants.
  expect_equal(dcmp$additive_var + dcmp$interaction_var, dcmp$between_var,
               tolerance = 1e-8)
  expect_equal(dcmp$additive_share + dcmp$interaction_share, 1, tolerance = 1e-8)
  vpc <- cc$summary$vpc$estimate
  expect_equal(vpc, dcmp$between_var / (dcmp$between_var + dcmp$within_var),
               tolerance = 1e-8)
  # Per-dimension additive variances are labelled by the dimension names and sum to
  # the additive variance.
  expect_setequal(names(dcmp$per_dim), c("a", "b", "cc"))
  expect_equal(sum(dcmp$per_dim), dcmp$additive_var, tolerance = 1e-8)
  # With a real interaction the shares are strictly inside (0, 1).
  expect_gt(dcmp$interaction_share, 0)
  expect_lt(dcmp$interaction_share, 1)
})

test_that("the crossed-dimensions formula carries dimension REs + the intersection RE", {
  d <- make_cc_data()
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, decomposition = "crossed-dimensions")))
  rhs <- paste(deparse(cc$formula), collapse = " ")
  expect_true(grepl("1 | a", rhs, fixed = TRUE))
  expect_true(grepl("1 | b", rhs, fixed = TRUE))
  expect_true(grepl("1 | cc", rhs, fixed = TRUE))
  expect_true(grepl("1 | stratum", rhs, fixed = TRUE))
})

test_that("the crossed-dimensions summary table has one row per dimension + interaction", {
  d <- make_cc_data()
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, decomposition = "crossed-dimensions")))
  vc <- cc$summary$variance_components
  expect_identical(attr(vc, "kind"), "cross_classified")
  expect_equal(sum(grepl("^Additive: ", vc$component)), 3L)
  expect_true("Intersectional interaction" %in% vc$component)
  expect_true("Within-stratum (residual)" %in% vc$component)
  # Non-total proportions sum to 1.
  non_total <- vc[vc$component != "Total", ]
  expect_equal(sum(non_total$proportion), 1, tolerance = 1e-8)
})

test_that("crossed-dimensions mode needs at least two dimensions", {
  d <- make_cc_data()
  expect_error(
    suppressWarnings(suppressMessages(
      maihda(y ~ x + (1 | a), data = d, decomposition = "crossed-dimensions"))),
    "two stratum dimensions|single dimension")
})

# ---- bootstrap --------------------------------------------------------------

test_that("crossed-dimensions bootstrap returns finite VPC and share intervals", {
  d <- make_cc_data()
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, decomposition = "crossed-dimensions")))
  sb <- suppressWarnings(summary(cc$model, bootstrap = TRUE, n_boot = 25))
  expect_true(is.finite(sb$vpc$ci_lower) && is.finite(sb$vpc$ci_upper))
  # A percentile bootstrap interval is well-ordered but need not contain the point
  # estimate (the bootstrap VPC distribution can be skewed at small n_boot).
  expect_true(sb$vpc$ci_lower <= sb$vpc$ci_upper)
  asci <- sb$decomposition$additive_share_ci
  expect_length(asci, 2)
  expect_true(all(is.finite(asci)))
})

# ---- plots ------------------------------------------------------------------

test_that("crossed-dimensions plots render", {
  skip_if_not_installed("ggplot2")
  d <- make_cc_data()
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, decomposition = "crossed-dimensions")))
  expect_s3_class(plot(cc, type = "vpc"), "ggplot")
  expect_s3_class(suppressWarnings(plot(cc, type = "effect_decomp")), "ggplot")
  expect_s3_class(suppressWarnings(plot(cc, type = "predicted")), "ggplot")
})

# ---- group comparison -------------------------------------------------------

test_that("compare_maihda_groups(decomposition = 'crossed-dimensions') reports shares", {
  d <- make_cc_data()
  g <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ x + (1 | a:b:cc), data = d, group = "grp",
                          decomposition = "crossed-dimensions")))
  expect_s3_class(g, "maihda_group_comparison")
  expect_true(all(c("var_additive", "var_interaction", "additive_share",
                    "interaction_share") %in% names(g)))
  expect_false(any(c("pcv", "var_between_adjusted", "var_between_adjusted_ml") %in%
                     names(g)))
  ok <- g[g$status == "ok", ]
  expect_true(all(ok$additive_share >= -1e-8 & ok$additive_share <= 1 + 1e-8))
  expect_identical(attr(g, "decomposition"), "crossed-dimensions")
})

test_that("plot(type = 'all') includes group_additive_share for a crossed-dimensions group analysis", {
  skip_if_not_installed("ggplot2")
  d <- make_cc_data()
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, group = "grp",
           decomposition = "crossed-dimensions")))

  # The group comparison carries the crossed-dimensions additive share, not the
  # two-model PCV -- the two are mutually exclusive by decomposition mode.
  expect_true("additive_share" %in% names(a$groups))
  expect_false("pcv" %in% names(a$groups))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  plots <- suppressWarnings(plot(a))          # type = "all"

  # Regression: the additive-share-by-group view (the one estimable for this mode)
  # must appear in the default plot(); the inapplicable group_pcv drops out.
  expect_s3_class(plots$group_additive_share, "ggplot")
  expect_null(plots$group_pcv)
})

test_that("plot(type = 'all') keeps group_pcv (not additive_share) for a two-model group analysis", {
  skip_if_not_installed("ggplot2")
  d <- make_cc_data()
  a <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b), data = d, group = "grp")))   # default: two-model

  expect_true("pcv" %in% names(a$groups))
  expect_false("additive_share" %in% names(a$groups))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  plots <- suppressWarnings(plot(a))          # type = "all"

  # The mutual exclusion holds in the other direction: two-model shows the PCV
  # view and drops the crossed-dimensions additive share.
  expect_s3_class(plots$group_pcv, "ggplot")
  expect_null(plots$group_additive_share)
})

test_that("crossed-dimensions group comparison requires shared_strata = TRUE", {
  d <- make_cc_data()
  expect_error(
    compare_maihda_groups(y ~ x + (1 | a:b:cc), data = d, group = "grp",
                          shared_strata = FALSE, decomposition = "crossed-dimensions"),
    "shared_strata")
})

# ---- brms parity ------------------------------------------------------------

test_that("crossed-dimensions brms summary returns posterior shares with intervals", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). The draws-based
  # partition logic is covered Stan-free by the helper test above.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- make_cc_data(seed = 7100, n = 900)
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b:cc), data = d, engine = "brms",
           decomposition = "crossed-dimensions",
           chains = 2, iter = 500, refresh = 0, seed = 1)))

  dcmp <- cc$decomposition
  expect_identical(dcmp$method, "posterior")
  expect_equal(dcmp$additive_share + dcmp$interaction_share, 1, tolerance = 1e-6)
  expect_length(dcmp$additive_share_ci, 2)
  expect_true(all(is.finite(dcmp$additive_share_ci)))
  expect_identical(cc$summary$vpc$method, "posterior")
  expect_true(cc$summary$vpc$ci_lower < cc$summary$vpc$ci_upper)
})

# ---- regression: the default two-model path is unchanged --------------------

test_that("default maihda() stays two-model with a PCV and no decomposition", {
  d <- make_cc_data()
  a <- suppressWarnings(suppressMessages(maihda(y ~ x + (1 | a:b:cc), data = d)))
  expect_identical(a$mode, "two-model")
  expect_false(is.null(a$pcv))
  expect_null(a$decomposition)
  expect_false(is.null(a$model_adjusted))
})

test_that("default compare_maihda_groups() keeps the two-model PCV schema", {
  d <- make_cc_data()
  g <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ x + (1 | a:b:cc), data = d, group = "grp")))
  expect_true("pcv" %in% names(g))
  expect_false("additive_share" %in% names(g))
})

# ---- deprecated alias --------------------------------------------------------

test_that("decomposition = 'cross-classified' warns and maps to 'crossed-dimensions'", {
  d <- make_cc_data()
  # capture_warnings() so the deprecation warning can be asserted while incidental
  # lme4 fit warnings (singular dimension REs) are tolerated.
  w1 <- capture_warnings(
    a <- suppressMessages(
      maihda(y ~ x + (1 | a:b:cc), data = d, decomposition = "cross-classified")))
  expect_true(any(grepl("renamed.*crossed-dimensions", w1)))
  expect_identical(a$mode, "crossed-dimensions")
  expect_false(is.null(a$decomposition))

  w2 <- capture_warnings(
    g <- suppressMessages(
      compare_maihda_groups(y ~ x + (1 | a:b:cc), data = d, group = "grp",
                            decomposition = "cross-classified")))
  expect_true(any(grepl("renamed.*crossed-dimensions", w2)))
  expect_identical(attr(g, "decomposition"), "crossed-dimensions")
  expect_true("additive_share" %in% names(g))
})

# ---- brms crossed-dimensions summary with a context (Stan-free) ---------------
# maihda_cc_summary_brms() is a pure function of the posterior draws once
# maihda_posterior_draws_brms() is mocked (the same pattern as the contextual
# brms tests in test-contextual.R), so the context-aware partition is testable
# without compiling a Stan model.

test_that("maihda_cc_summary_brms folds a context into the partition (Stan-free)", {
  set.seed(7311)
  n <- 400
  draws <- data.frame(
    sd_a__Intercept = sqrt(0.9 + 0.2 * runif(n)),
    sd_b__Intercept = sqrt(0.6 + 0.2 * runif(n)),
    sd_stratum__Intercept = sqrt(0.7 + 0.2 * runif(n)),
    sd_site__Intercept = sqrt(0.5 + 0.2 * runif(n)),
    sigma = sqrt(1.1 + 0.2 * runif(n))
  )
  stub <- structure(
    list(formula = y ~ x + (1 | a) + (1 | b) + (1 | stratum) + (1 | site),
         family = list(family = "gaussian", link = "identity")),
    class = "brmsfit")
  object <- list(model = stub, engine = "brms",
                 context_info = list(context_vars = "site"))
  cc_info <- list(dim_groups = c(a = "a", b = "b"), interaction_group = "stratum")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)

  res <- maihda_cc_summary_brms(object, cc_info, conf_level = 0.9)

  v_a <- draws$sd_a__Intercept^2
  v_b <- draws$sd_b__Intercept^2
  v_s <- draws$sd_stratum__Intercept^2
  v_c <- draws$sd_site__Intercept^2
  v_e <- draws$sigma^2
  between <- v_a + v_b + v_s
  total <- between + v_c + v_e
  expect_equal(res$vpc_result$estimate, stats::median(between / total))
  expect_equal(res$decomposition$additive_share,
               stats::median((v_a + v_b) / between))

  # The context element mirrors the contextual summary.
  expect_identical(res$context$context_vars, "site")
  expect_equal(res$context$vpc_context_total, stats::median(v_c / total))
  expect_length(res$context$vpc_context_total_ci, 2)
  expect_equal(unname(res$context$per_context["site"]), mean(v_c))

  tab <- res$variance_components
  expect_true("Context: site" %in% tab$component)
  expect_equal(sum(tab$proportion[tab$component != "Total"]), 1, tolerance = 1e-8)
})

test_that("maihda_cc_summary_brms names a context RE missing from the draws", {
  draws <- data.frame(
    sd_a__Intercept = c(1, 1.1),
    sd_b__Intercept = c(0.8, 0.9),
    sd_stratum__Intercept = c(0.7, 0.8),
    sigma = c(1, 1.1)
  )
  stub <- structure(
    list(formula = y ~ 1, family = list(family = "gaussian", link = "identity")),
    class = "brmsfit")
  object <- list(model = stub, engine = "brms",
                 context_info = list(context_vars = "site"))
  cc_info <- list(dim_groups = c(a = "a", b = "b"), interaction_group = "stratum")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)
  expect_error(maihda_cc_summary_brms(object, cc_info, 0.9),
               "missing the random effect")
})

# ---- decomposition print paths -------------------------------------------------

test_that("maihda_print_cc_decomposition prints shares, per-dim variances, and CIs", {
  d <- list(
    additive_var = 1.5, interaction_var = 0.5, between_var = 2, within_var = 2,
    additive_share = 0.75, interaction_share = 0.25,
    per_dim = c(a = 1, b = 0.5),
    additive_share_ci = NULL, interaction_share_ci = NULL
  )
  expect_output(maihda_print_cc_decomposition(d), "crossed-dimensions")
  expect_output(maihda_print_cc_decomposition(d), "different estimator")
  expect_output(maihda_print_cc_decomposition(d), "b: 0.5000")

  d$additive_share_ci <- c(0.6, 0.9)
  d$interaction_share_ci <- c(0.1, 0.4)
  expect_output(maihda_print_cc_decomposition(d), "[60.0%, 90.0%]", fixed = TRUE)
})

# ---- per-group crossed-dimensions via maihda(group = ) --------------------------

test_that("maihda(crossed-dimensions, group = ) attaches the per-group decomposition", {
  d <- make_cc_data(n = 1200)
  cc <- suppressWarnings(suppressMessages(
    maihda(y ~ x + (1 | a:b), data = d, decomposition = "crossed-dimensions",
           group = "grp")))
  expect_s3_class(cc$groups, "maihda_group_comparison")
  expect_true(all(c("additive_share", "interaction_share") %in% names(cc$groups)))

  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(cc$groups, type = "additive_share"), "ggplot")
  expect_s3_class(plot(cc$groups, type = "components"), "ggplot")
})

test_that("plot(type = 'additive_share') errors on a two-model group comparison", {
  d <- make_cc_data(n = 800)
  g <- suppressWarnings(suppressMessages(
    compare_maihda_groups(y ~ x + (1 | a:b), data = d, group = "grp")))
  expect_error(plot(g, type = "additive_share"), "crossed-dimensions")
})

test_that("compare_maihda_groups(crossed-dimensions) needs two stratum dimensions", {
  d <- make_cc_data(n = 800)
  expect_error(
    suppressWarnings(compare_maihda_groups(
      y ~ x + (1 | a), data = d, group = "grp",
      decomposition = "crossed-dimensions")),
    "at least two stratum dimensions")
})

test_that("maihda_cc_summary_brms without a context returns context = NULL", {
  set.seed(7312)
  n <- 200
  draws <- data.frame(
    sd_a__Intercept = sqrt(0.9 + 0.2 * runif(n)),
    sd_b__Intercept = sqrt(0.6 + 0.2 * runif(n)),
    sd_stratum__Intercept = sqrt(0.7 + 0.2 * runif(n)),
    sigma = sqrt(1.1 + 0.2 * runif(n))
  )
  stub <- structure(
    list(formula = y ~ 1, family = list(family = "gaussian", link = "identity")),
    class = "brmsfit")
  object <- list(model = stub, engine = "brms")
  cc_info <- list(dim_groups = c(a = "a", b = "b"), interaction_group = "stratum")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) draws)

  res <- maihda_cc_summary_brms(object, cc_info, conf_level = 0.95)
  expect_null(res$context)
  expect_false(any(grepl("^Context: ", res$variance_components$component)))
})
