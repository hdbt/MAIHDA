# Audit 2026-10-12 (brms-mor-variance-table): an engine = "brms" analysis reports ONE
# between-stratum variance, the posterior median of its draws.
#
# The PCV compares posterior medians. The Median Odds Ratio and the variance-components
# table of summary() -- with what is read from it: maihda_table()'s rows,
# compare_maihda_groups()'s variance columns, stepwise_pcv()'s Context_Variance, the VPC
# plots -- hold the same medians, with the table's Total and proportions computed from
# its rows. A crossed-dimensions fit follows the same rule. Measured against a known
# truth at 8 strata, the posterior mean of the null model's variance averages 0.45 and
# its median 0.34 for a true 0.26 (MOR 1.89 and 1.73 for 1.62).
#
# The blocks below pin each summary on hand-built draws (no Stan) and repeat the checks
# on real fits when MAIHDA_TEST_BRMS=true.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
mvt_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

mvt_mor <- function(v) exp(sqrt(2 * v) * stats::qnorm(0.75))

# Skewed SD draws, built without random numbers. One far draw each puts the mean of
# every variance well above its median. An odd number of draws, so each median is one
# of the draws and (median SD)^2 is the median variance exactly. The fixed shuffles keep
# the components from rising together, so that a ratio or a sum of medians is not the
# median of the per-draw ratio or sum.
mvt_n <- 201L
mvt_shuffle <- order(sin(seq_len(mvt_n) * 12.9898))
mvt_sd_null <- c(seq(0.35, 0.75, length.out = 200), 4)
mvt_sd_adj  <- c(seq(0.05, 0.45, length.out = 200), 5)
mvt_sd_ctx  <- c(seq(0.20, 0.60, length.out = 200), 4.5)[mvt_shuffle]
mvt_sd_ctx2 <- c(seq(0.10, 0.90, length.out = 200), 3)[rev(mvt_shuffle)]
mvt_sd_dim  <- c(seq(0.30, 1.10, length.out = 200), 6)[order(cos(seq_len(mvt_n) * 78.233))]
mvt_sigma   <- seq(0.80, 1.20, length.out = mvt_n)[rev(order(cos(seq_len(mvt_n) * 78.233)))]

mvt_logit <- list(family = "bernoulli", link = "logit")
mvt_gauss <- list(family = "gaussian", link = "identity")

# A maihda_model shell around a brmsfit stub that carries its own posterior draws.
mvt_stub <- function(draws, formula = y ~ 1 + (1 | stratum), family = mvt_logit, ...) {
  structure(
    c(list(model = structure(list(formula = formula, family = family, draws = draws),
                             class = "brmsfit"),
           engine = "brms", formula = formula, family = family),
      list(...)),
    class = "maihda_model"
  )
}

# The draws and brms's own SD summary, both read off the stub. VarCorr(robust = TRUE)
# reports the median SD as its Estimate, as brms does.
mvt_mock_brms <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) model$draws,
    .env = env
  )
  local_mocked_bindings(
    VarCorr = function(x, ...) {
      robust <- isTRUE(list(...)$robust)
      sd_cols <- grep("^sd_", names(x$draws), value = TRUE)
      out <- lapply(sd_cols, function(cn) {
        est <- if (robust) stats::median(x$draws[[cn]]) else mean(x$draws[[cn]])
        list(sd = matrix(
          c(est, 0, 0, 0), nrow = 1,
          dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5"))))
      })
      stats::setNames(out, sub("__Intercept$", "", sub("^sd_", "", sd_cols)))
    },
    .package = "brms", .env = env
  )
}

test_that("the table of a logistic brms fit holds the posterior-median variance", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  m <- mvt_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = mvt_sd_null))
  v <- mvt_sd_null^2
  r <- pi^2 / 3

  res <- maihda_vpc_draws_brms(m$model)
  expect_identical(res$var_stratum, stats::median(v))
  expect_identical(res$var_other_random, 0)
  expect_identical(res$var_residual, r)
  # Fixture premise: the mean of the variance lies well above its median.
  expect_gt(mean(v), 1.2 * stats::median(v))

  tab <- maihda_variance_components_table(res$var_stratum, res$var_other_random,
                                          res$var_residual)
  between <- tab$component == "Between-stratum (random)"
  expect_identical(tab$variance[between], stats::median(v))
  # The SD row is the posterior median SD, not the root of a mean variance.
  expect_equal(tab$sd[between], stats::median(mvt_sd_null), tolerance = 1e-12)
  expect_gt(sqrt(mean(v)) - tab$sd[between], 0.05)
  expect_identical(tab$variance[tab$component == "Total"], stats::median(v) + r)
  # With the level-1 variance fixed by the link and one random effect, the VPC is a
  # monotone function of the variance: the table's stratum proportion IS the headline
  # VPC, the posterior median of the per-draw VPC.
  expect_equal(res$vpc$estimate, stats::median(v / (v + r)), tolerance = 1e-12)
  expect_equal(tab$proportion[between], res$vpc$estimate, tolerance = 1e-12)
  # A table built from posterior means would show a larger share.
  expect_gt(mean(v) / (mean(v) + r) - tab$proportion[between], 0.02)

  # `point` summarises the VPC; the variance components stay posterior medians.
  res_mean <- maihda_vpc_draws_brms(m$model, point = "mean")
  expect_equal(res_mean$vpc$estimate, mean(v / (v + r)), tolerance = 1e-12)
  expect_identical(res_mean$var_stratum, stats::median(v))
})

test_that("each row of a brms table is its own posterior median; totals come from the rows", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  # A Gaussian fit with a second random effect: three components that all vary.
  m <- mvt_stub(
    data.frame(b_Intercept = 0, sd_stratum__Intercept = mvt_sd_null,
               sd_site__Intercept = mvt_sd_ctx, sigma = mvt_sigma),
    formula = y ~ 1 + (1 | stratum) + (1 | site), family = mvt_gauss)
  v_s <- mvt_sd_null^2
  v_o <- mvt_sd_ctx^2
  v_e <- mvt_sigma^2

  res <- maihda_vpc_draws_brms(m$model)
  expect_identical(res$var_stratum, stats::median(v_s))
  expect_equal(res$var_other_random, stats::median(v_o), tolerance = 1e-12)
  expect_identical(res$var_residual, stats::median(v_e))

  tab <- maihda_variance_components_table(res$var_stratum, res$var_other_random,
                                          res$var_residual)
  rows <- tab$component != "Total"
  total <- stats::median(v_s) + stats::median(v_o) + stats::median(v_e)
  expect_equal(tab$variance[!rows], total, tolerance = 1e-12)
  expect_equal(tab$variance[!rows], sum(tab$variance[rows]), tolerance = 1e-12)
  expect_equal(tab$proportion[rows], tab$variance[rows] / total, tolerance = 1e-12)
  expect_equal(sum(tab$proportion[rows]), 1, tolerance = 1e-12)

  # The headline VPC stays the posterior median of the per-draw share. The table's
  # stratum proportion is built from the rows, so here the two are not the same number,
  # and the Total is not the median of the per-draw total.
  share <- v_s / (v_s + v_o + v_e)
  expect_equal(res$vpc$estimate, stats::median(share), tolerance = 1e-12)
  expect_gt(abs(tab$proportion[1] - stats::median(share)), 1e-4)
  expect_gt(abs(total - stats::median(v_s + v_o + v_e)), 1e-4)
})

test_that("a contextual brms summary tabulates posterior medians that add up", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  m <- mvt_stub(
    data.frame(b_Intercept = 0, sd_stratum__Intercept = mvt_sd_null,
               sd_site__Intercept = mvt_sd_ctx, sd_region__Intercept = mvt_sd_ctx2,
               sd_extra__Intercept = mvt_sd_adj),
    formula = y ~ 1 + (1 | stratum) + (1 | site) + (1 | region) + (1 | extra))
  v_s <- mvt_sd_null^2
  v_c1 <- mvt_sd_ctx^2
  v_c2 <- mvt_sd_ctx2^2
  v_x <- mvt_sd_adj^2
  r <- pi^2 / 3

  res <- maihda_context_summary_brms(m, list(context_vars = c("site", "region")), 0.95)
  cx <- res$context
  expect_identical(cx$var_stratum, stats::median(v_s))
  expect_identical(unname(cx$per_context), c(stats::median(v_c1), stats::median(v_c2)))
  expect_identical(cx$other_var, stats::median(v_x))
  expect_identical(cx$within_var, r)
  # The contexts' total is the sum of the context rows, which is not the median of the
  # per-draw sum on these draws.
  expect_equal(cx$context_var_total, stats::median(v_c1) + stats::median(v_c2),
               tolerance = 1e-12)
  expect_gt(abs(stats::median(v_c1 + v_c2) - cx$context_var_total), 0.005)

  # The variance a PCV of this fit compares is the stratum row.
  expect_identical(maihda_pcv_between_variance(m), cx$var_stratum)

  tab <- res$variance_components
  rows <- tab$component != "Total"
  expect_identical(tab$variance[tab$component == "Between-stratum (random)"],
                   stats::median(v_s))
  expect_identical(tab$variance[tab$component == "Context: site"], stats::median(v_c1))
  expect_identical(tab$variance[tab$component == "Context: region"], stats::median(v_c2))
  expect_identical(tab$variance[tab$component == "Other random effects"],
                   stats::median(v_x))
  expect_equal(tab$variance[!rows], sum(tab$variance[rows]), tolerance = 1e-12)
  expect_equal(sum(tab$proportion[rows]), 1, tolerance = 1e-12)

  # The shares keep their per-draw summary.
  total <- v_s + v_c1 + v_c2 + v_x + r
  expect_equal(res$vpc_result$estimate, stats::median(v_s / total), tolerance = 1e-12)
  expect_equal(cx$vpc_context_total, stats::median((v_c1 + v_c2) / total),
               tolerance = 1e-12)
})

test_that("a crossed-dimensions brms summary tabulates posterior medians that add up", {
  skip_on_cran()
  skip_if_not_installed("brms")
  mvt_mock_brms()
  m <- mvt_stub(
    data.frame(b_Intercept = 0, sd_ga__Intercept = mvt_sd_dim,
               sd_gb__Intercept = mvt_sd_ctx2, sd_stratum__Intercept = mvt_sd_adj,
               sd_site__Intercept = mvt_sd_ctx),
    formula = y ~ 1 + (1 | ga) + (1 | gb) + (1 | stratum) + (1 | site),
    context_info = list(context_vars = "site"))
  cc <- list(dim_groups = c(a = "ga", b = "gb"), interaction_group = "stratum")
  v_a <- mvt_sd_dim^2
  v_b <- mvt_sd_ctx2^2
  v_i <- mvt_sd_adj^2
  v_c <- mvt_sd_ctx^2
  r <- pi^2 / 3

  res <- maihda_cc_summary_brms(m, cc, conf_level = 0.95)
  dc <- res$decomposition
  expect_identical(dc$per_dim, c(a = stats::median(v_a), b = stats::median(v_b)))
  expect_identical(dc$interaction_var, stats::median(v_i))
  expect_identical(dc$within_var, r)
  # Additive and between-strata totals are sums of the rows shown.
  expect_equal(dc$additive_var, stats::median(v_a) + stats::median(v_b), tolerance = 1e-12)
  expect_equal(dc$between_var, dc$additive_var + dc$interaction_var, tolerance = 1e-12)
  expect_gt(abs(stats::median(v_a + v_b + v_i) - dc$between_var), 0.005)
  # Well below the posterior mean the far draw carries.
  expect_gt(mean(v_a), 1.3 * dc$per_dim[["a"]])

  expect_identical(res$context$var_stratum, dc$between_var)
  expect_identical(unname(res$context$per_context), stats::median(v_c))
  expect_equal(res$context$context_var_total, stats::median(v_c), tolerance = 1e-12)

  tab <- res$variance_components
  rows <- tab$component != "Total"
  expect_equal(tab$variance[!rows], sum(tab$variance[rows]), tolerance = 1e-12)
  expect_equal(sum(tab$variance[rows]),
               dc$between_var + stats::median(v_c) + r, tolerance = 1e-12)

  # The VPC and the shares keep their per-draw summary, so the table's between-strata
  # proportion is not the headline VPC on these draws.
  between <- v_a + v_b + v_i
  total <- between + v_c + r
  expect_equal(res$vpc_result$estimate, stats::median(between / total), tolerance = 1e-12)
  expect_equal(dc$additive_share, stats::median((v_a + v_b) / between), tolerance = 1e-12)
  expect_gt(abs(dc$between_var / sum(tab$variance[rows]) - res$vpc_result$estimate), 1e-3)

  # A Gaussian outcome and two contexts: the residual row is the median of sigma^2 (a
  # logit residual is a constant and could not tell), and the contexts' total is the
  # sum of the two context rows, not the median of the per-draw sum.
  m2 <- mvt_stub(
    data.frame(b_Intercept = 0, sd_ga__Intercept = mvt_sd_dim,
               sd_gb__Intercept = mvt_sd_ctx2, sd_stratum__Intercept = mvt_sd_adj,
               sd_site__Intercept = mvt_sd_ctx, sd_region__Intercept = mvt_sd_null,
               sigma = mvt_sigma),
    formula = y ~ 1 + (1 | ga) + (1 | gb) + (1 | stratum) + (1 | site) + (1 | region),
    family = mvt_gauss, context_info = list(context_vars = c("site", "region")))
  v_r <- mvt_sd_null^2
  v_e <- mvt_sigma^2
  res2 <- maihda_cc_summary_brms(m2, cc, conf_level = 0.95)
  expect_identical(res2$decomposition$within_var, stats::median(v_e))
  expect_gt(abs(mean(v_e) - stats::median(v_e)), 0.005)
  tab2 <- res2$variance_components
  expect_identical(tab2$variance[tab2$component == "Within-stratum (residual)"],
                   stats::median(v_e))
  expect_identical(unname(res2$context$per_context),
                   c(stats::median(v_c), stats::median(v_r)))
  expect_equal(res2$context$context_var_total, stats::median(v_c) + stats::median(v_r),
               tolerance = 1e-12)
  expect_gt(abs(stats::median(v_c + v_r) - res2$context$context_var_total), 0.005)
  rows2 <- tab2$component != "Total"
  expect_equal(tab2$variance[!rows2], sum(tab2$variance[rows2]), tolerance = 1e-12)
})

test_that("the brms MOR is the posterior median of the per-draw MOR", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  v <- mvt_sd_null^2
  for (fam in list(mvt_logit, list(family = "cumulative", link = "logit"))) {
    m <- mvt_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = mvt_sd_null),
                  family = fam)
    expect_identical(maihda_mor_between_variance(m), stats::median(v))
    expect_equal(maihda_mor(m), mvt_mor(stats::median(v)), tolerance = 1e-12)
    # The closed form is monotone in the variance, so that is the median of the MOR.
    expect_equal(maihda_mor(m), stats::median(mvt_mor(v)), tolerance = 1e-12)
    # Not the MOR at the posterior-mean variance, nor the posterior mean of the MOR.
    expect_gt(mvt_mor(mean(v)) - maihda_mor(m), 0.05)
    expect_gt(mean(mvt_mor(v)) - maihda_mor(m), 0.05)
  }
})

test_that("the crossed-dimensions brms MOR is evaluated at the posterior-median variances", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  cells <- expand.grid(ga = c("p", "q"), gb = c("x", "y", "z"), stringsAsFactors = FALSE)
  m <- mvt_stub(
    data.frame(b_Intercept = 0, sd_ga__Intercept = mvt_sd_dim,
               sd_gb__Intercept = mvt_sd_ctx2, sd_stratum__Intercept = mvt_sd_adj),
    formula = y ~ 1 + (1 | ga) + (1 | gb) + (1 | stratum),
    cc_info = list(dim_groups = c(a = "ga", b = "gb"), interaction_group = "stratum"),
    data = cells[rep(seq_len(nrow(cells)), each = 3), ])
  med <- c(stats::median(mvt_sd_dim^2), stats::median(mvt_sd_ctx2^2))
  avg <- c(mean(mvt_sd_dim^2), mean(mvt_sd_ctx2^2))

  cv <- maihda_cc_variances(m)
  expect_identical(unname(cv$dims), med)
  expect_identical(cv$interaction, stats::median(mvt_sd_adj^2))
  expect_equal(maihda_mor_between_variance(m), sum(med) + stats::median(mvt_sd_adj^2),
               tolerance = 1e-12)

  parts <- maihda_mor_crossed_parts(m)
  expect_identical(parts$dims, med)
  expect_identical(maihda_mor(m), maihda_mor_crossed_from_parts(parts))
  # The posterior means would give a markedly larger MOR.
  parts_mean <- parts
  parts_mean$dims <- avg
  parts_mean$interaction <- mean(mvt_sd_adj^2)
  expect_gt(maihda_mor_crossed_from_parts(parts_mean) - maihda_mor(m), 0.2)
})

test_that("the PCV, the table and the group comparison share one pair of variances", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  # The stubs are not fitted models, so the same-sample checks have nothing to read.
  local_mocked_bindings(validate_pcv_models = function(model1, model2) invisible(TRUE))
  m0 <- mvt_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = mvt_sd_null))
  m1 <- mvt_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = mvt_sd_adj),
                 formula = y ~ a + b + (1 | stratum))

  row0 <- maihda_vpc_draws_brms(m0$model)$var_stratum
  row1 <- maihda_vpc_draws_brms(m1$model)$var_stratum
  pcv <- calculate_pcv(m0, m1)
  expect_identical(pcv$var_model1, row0)
  expect_identical(pcv$var_model2, row1)
  # The PCV is the reduction the table's two variance rows show.
  expect_equal(pcv$pcv, 1 - row1 / row0, tolerance = 1e-12)
  # compare_maihda_groups() derives var_between_adjusted as var_between * (1 - pcv):
  # with var_between the table's row, that is the adjusted fit's own variance, the
  # value it also returns as var_between_adjusted_ml.
  expect_equal(row0 * (1 - pcv$pcv), pcv$var_model2, tolerance = 1e-12)
  # So maihda_table() has no second basis to flag.
  shown <- list(null = list("Between-stratum variance" = c(row0, NA_real_, NA_real_)),
                adjusted = list("Between-stratum variance" = c(row1, NA_real_, NA_real_)))
  expect_null(maihda_table_pcv_note(shown, pcv, c("null", "adjusted"), "brms"))
})

test_that("stepwise_pcv() holds a brms context variance as the sum of posterior medians", {
  skip_if_not_installed("brms")
  mvt_mock_brms()
  sd_by_step <- list(mvt_sd_null, mvt_sd_adj)
  local_mocked_bindings(fit_maihda = function(formula, data, ...) {
    k <- length(attr(stats::terms(reformulas::nobars(formula)), "term.labels"))
    mvt_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = sd_by_step[[k + 1L]],
                        sd_site__Intercept = mvt_sd_ctx, sd_region__Intercept = mvt_sd_ctx2),
             formula = formula, family = mvt_gauss)
  })
  d <- data.frame(stratum = factor(rep(seq_len(6), each = 10)),
                  site = factor(rep(seq_len(5), times = 12)),
                  region = factor(rep(seq_len(3), times = 20)),
                  a = rep(c(-1, 0, 1), times = 20), y = rep(seq(-2, 2, length.out = 12), 5))

  out <- suppressMessages(suppressWarnings(
    stepwise_pcv(d, "y", "a", engine = "brms", context = c("site", "region"))))
  expect_identical(out$Variance,
                   c(stats::median(mvt_sd_null^2), stats::median(mvt_sd_adj^2)))
  ctx_total <- stats::median(mvt_sd_ctx^2) + stats::median(mvt_sd_ctx2^2)
  expect_equal(out$Context_Variance, c(ctx_total, ctx_total), tolerance = 1e-12)
  # The posterior means would give a larger context variance.
  expect_gt(mean(mvt_sd_ctx^2) + mean(mvt_sd_ctx2^2) - ctx_total, 0.05)

  printed <- capture.output(print(out))
  expect_true(any(grepl("Variance basis: posterior medians", printed, fixed = TRUE)))
  expect_false(any(grepl("posterior means", printed, fixed = TRUE)))
})

test_that("the help pages say which summary a brms variance and MOR are", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  read_rd <- function(f) {
    gsub("[[:space:]]+", " ", paste(readLines(file.path(man, f), warn = FALSE), collapse = " "))
  }
  smry <- read_rd("summary.maihda_model.Rd")
  expect_true(grepl("posterior median of each variance component", smry, fixed = TRUE))
  expect_true(grepl("computed from those rows", smry, fixed = TRUE))
  expect_true(grepl("median overstates the variance", smry, fixed = TRUE))
  expect_false(grepl("posterior-mean variance components", smry, fixed = TRUE))

  mor <- read_rd("maihda_mor.Rd")
  expect_true(grepl("\\strong{brms fits.}", mor, fixed = TRUE))
  expect_true(grepl("posterior median of the Median Odds Ratio", mor, fixed = TRUE))

  # No page describes a brms variance, MOR or context variance as a posterior mean.
  for (f in c("calculate_pcv.Rd", "stepwise_pcv.Rd", "compare_maihda_groups.Rd",
              "maihda_table.Rd")) {
    expect_false(grepl("posterior mean", read_rd(f), fixed = TRUE), info = f)
  }
})

test_that("real brms fits: one between-stratum variance in the table, the MOR and the PCV", {
  # Compiles three Stan models, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # Eight strata, 100 Bernoulli rows each: additive main effects plus a small stratum
  # effect. brms draws its sampler seed from R's stream, so the fits sit inside the
  # seeded block too.
  fits <- mvt_with_seed(20261012, {
    g <- expand.grid(d1 = c("a", "b"), d2 = c("c", "d"), d3 = c("e", "f"),
                     stringsAsFactors = FALSE)
    eta <- -0.9 + 0.5 * (g$d1 == "b") + 0.45 * (g$d2 == "d") + 0.55 * (g$d3 == "f") +
      stats::rnorm(8, 0, 0.2)
    d <- g[rep(seq_len(8), each = 100), ]
    d$y <- stats::rbinom(800, 1, stats::plogis(eta[rep(seq_len(8), each = 100)]))
    rownames(d) <- NULL
    list(
      two = suppressWarnings(suppressMessages(
        maihda(y ~ 1 + (1 | d1:d2:d3), data = d, engine = "brms", family = "binomial",
               chains = 2, iter = 1000, refresh = 0))),
      crossed = suppressWarnings(suppressMessages(
        maihda(y ~ 1 + (1 | d1:d2:d3), data = d, engine = "brms", family = "binomial",
               decomposition = "crossed-dimensions", chains = 2, iter = 1000,
               refresh = 0))))
  })
  a <- fits$two
  v0 <- as.data.frame(a$model$model)$sd_stratum__Intercept^2
  v1 <- as.data.frame(a$model_adjusted$model)$sd_stratum__Intercept^2
  between <- function(s) {
    s$variance_components$variance[s$variance_components$component ==
                                     "Between-stratum (random)"]
  }
  share <- function(s) {
    s$variance_components$proportion[s$variance_components$component ==
                                       "Between-stratum (random)"]
  }

  # The table's rows are the medians the PCV compares, not the means.
  expect_identical(between(a$summary), stats::median(v0))
  expect_identical(between(a$summary_adjusted), stats::median(v1))
  expect_identical(a$pcv$var_model1, between(a$summary))
  expect_identical(a$pcv$var_model2, between(a$summary_adjusted))
  expect_gt(mean(v0), stats::median(v0))
  # The stratum proportion is the headline VPC. The number of draws is even here, so a
  # median of a transform and the transform of the median differ by the interpolation
  # between the two middle draws; hence the tolerances below.
  expect_equal(share(a$summary), a$summary$vpc$estimate, tolerance = 1e-5)
  expect_equal(share(a$summary_adjusted), a$summary_adjusted$vpc$estimate,
               tolerance = 1e-5)

  # The MOR is the posterior median of the per-draw MOR.
  expect_equal(maihda_mor(a$model), mvt_mor(stats::median(v0)), tolerance = 1e-12)
  expect_equal(maihda_mor(a$model), stats::median(mvt_mor(v0)), tolerance = 1e-5)
  expect_identical(a$summary$discriminatory_accuracy$mor, maihda_mor(a$model))
  expect_lt(maihda_mor(a$model), mvt_mor(mean(v0)))

  # The report table: its PCV row is the reduction its variance rows show, its MOR row
  # the MOR at its variance row, and it carries no note.
  tab <- suppressWarnings(suppressMessages(maihda_table(a)))
  rw <- function(stat) tab$models[tab$models$statistic == stat, ]
  bv <- rw("Between-stratum variance")
  expect_identical(c(bv$null, bv$adjusted), c(stats::median(v0), stats::median(v1)))
  expect_equal(rw("Between-stratum SD")$null, stats::median(sqrt(v0)), tolerance = 1e-5)
  expect_equal(rw("PCV (null -> adjusted)")$adjusted, 1 - bv$adjusted / bv$null,
               tolerance = 1e-12)
  expect_equal(rw("MOR")$null, mvt_mor(bv$null), tolerance = 1e-12)
  expect_equal(rw("VPC/ICC")$null, bv$null / (bv$null + pi^2 / 3), tolerance = 1e-5)
  expect_null(tab$models_note)

  # print() of the analysis shows the same pair of variances as the table.
  out <- capture.output(print(a))
  expect_true(any(grepl(
    sprintf("Between-stratum variance: %.4f (null) -> %.4f (adjusted)",
            between(a$summary), between(a$summary_adjusted)), out, fixed = TRUE)))

  # The VPC bar's stratum slice is the VPC its title quotes.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- plot_vpc(a$summary)
    expect_equal(p$data$proportion[p$data$component == "Between-stratum (random)"],
                 a$summary$vpc$estimate, tolerance = 1e-5)
  }

  # Crossed-dimensions fit: rows, totals and the MOR from the posterior medians.
  cc <- fits$crossed
  dr <- as.data.frame(cc$model$model)
  gv <- maihda_group_variance_draws_brms(dr)
  info <- cc$model$cc_info
  med_dim <- vapply(unname(info$dim_groups), function(grp) stats::median(gv[[grp]]),
                    numeric(1))
  dc <- cc$summary$decomposition
  expect_identical(unname(dc$per_dim), unname(med_dim))
  expect_identical(dc$interaction_var, stats::median(gv[[info$interaction_group]]))
  expect_equal(dc$between_var, sum(med_dim) + dc$interaction_var, tolerance = 1e-12)
  parts <- maihda_mor_crossed_parts(cc$model)
  expect_identical(unname(parts$dims), unname(med_dim))
  expect_identical(maihda_mor(cc$model), maihda_mor_crossed_from_parts(parts))
})
