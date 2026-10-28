# Audit 2026-10-07 (brms-pcv-posterior-mean): the PCV of engine = "brms" fits compares
# the posterior MEDIANS of the between-stratum variances.
#
# The posterior of a between-stratum variance is right-skewed when strata are few, the
# adjusted model's more than the null's, so one minus a ratio of two posterior MEANS
# understates the PCV; and a posterior mean is set by an upper tail the sampler visits
# rarely, so that ratio also moves from one sampler seed to the next. Measured against
# a known truth of 0.84 at 8 strata: 0.58 from posterior means, 0.74 from medians.
#
# The blocks below pin the summary on hand-built draws (no Stan), check that the MOR
# reads the same variance and that the other engines keep their fitted one, that the
# printed tables carry one basis line, and repeat the check on real fits when
# MAIHDA_TEST_BRMS=true. The variance-components table and the MOR have their own
# file, test-audit-2026-10-12-brms-mor-variance-table.R.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
pcv_with_seed <- function(seed, expr) {
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

# Skewed SD draws: one far draw each, so the mean of the variance sits well above its
# median, and further above it for the "adjusted" set. An odd number of draws, so the
# median is one of them and (median SD)^2 is the median variance exactly.
pcv_sd_null <- c(seq(0.35, 0.75, length.out = 200), 4)
pcv_sd_mid  <- c(seq(0.20, 0.60, length.out = 200), 4.5)
pcv_sd_adj  <- c(seq(0.05, 0.45, length.out = 200), 5)

# A maihda_model shell around a brmsfit stub that carries its own SD draws.
pcv_brms_stub <- function(sd_draws, formula = y ~ 1 + (1 | stratum),
                          family = list(family = "gaussian", link = "identity")) {
  structure(
    list(model = structure(list(formula = formula, family = family, sd_draws = sd_draws),
                           class = "brmsfit"),
         engine = "brms", formula = formula, family = family),
    class = "maihda_model"
  )
}

# The draws and brms's own SD summary, both read off the stub. VarCorr(robust = TRUE)
# reports the median SD as its Estimate, as brms does.
pcv_mock_brms <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) {
      data.frame(b_Intercept = 0, sd_stratum__Intercept = model$sd_draws)
    },
    .env = env
  )
  local_mocked_bindings(
    VarCorr = function(x, ...) {
      est <- if (isTRUE(list(...)$robust)) stats::median(x$sd_draws) else mean(x$sd_draws)
      list(stratum = list(sd = matrix(
        c(est, 0, 0, 0), nrow = 1,
        dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5")))))
    },
    .package = "brms", .env = env
  )
}

test_that("the brms stratum variance is the posterior median of its draws", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  fit <- pcv_brms_stub(pcv_sd_adj)$model
  v <- pcv_sd_adj^2

  expect_identical(maihda_stratum_variance_brms(fit), stats::median(v))
  # The median does not depend on whether it is taken on the SD or on the variance.
  expect_equal(maihda_stratum_variance_brms(fit), stats::median(pcv_sd_adj)^2,
               tolerance = 1e-12)
  # Fixture premise: skewed enough that the median cannot be taken for the mean.
  expect_gt(mean(v), 3 * stats::median(v))
})

test_that("without draws the brms stratum variance falls back to brms's median SD", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) stop("no draws"))
  fit <- pcv_brms_stub(pcv_sd_adj)$model

  expect_equal(maihda_stratum_variance_brms(fit), stats::median(pcv_sd_adj)^2,
               tolerance = 1e-12)
  # Not the squared mean SD, which is what VarCorr() reports without robust = TRUE.
  expect_gt(abs(maihda_stratum_variance_brms(fit) - mean(pcv_sd_adj)^2), 0.01)
})

test_that("calculate_pcv() compares posterior medians for two brms fits", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  # The stubs are not fitted models, so the same-sample checks have nothing to read.
  local_mocked_bindings(validate_pcv_models = function(model1, model2) invisible(TRUE))
  m0 <- pcv_brms_stub(pcv_sd_null)
  m1 <- pcv_brms_stub(pcv_sd_adj, y ~ a + b + (1 | stratum))
  v0 <- pcv_sd_null^2
  v1 <- pcv_sd_adj^2

  res <- calculate_pcv(m0, m1)
  expect_identical(res$var_model1, stats::median(v0))
  expect_identical(res$var_model2, stats::median(v1))
  expect_equal(res$pcv, 1 - stats::median(v1) / stats::median(v0), tolerance = 1e-12)
  expect_s3_class(tryCatch(res$pvc, warning = function(w) w), "deprecatedWarning")
  expect_identical(unclass(res)$pvc, res$pcv)
  expect_identical(res$estimation_used, "posterior")
  expect_false(res$adjusted_at_boundary)

  # Not one minus the ratio of posterior means, which is far lower on these draws.
  pcv_means <- 1 - mean(v1) / mean(v0)
  expect_lt(pcv_means, 0.55)
  expect_gt(res$pcv, 0.75)

  # estimation = "ML" refits nothing on a brms fit.
  res_ml <- calculate_pcv(m0, m1, estimation = "ML")
  expect_identical(res_ml$pcv, res$pcv)
  expect_identical(res_ml$estimation_used, "posterior")

  out <- capture.output(print(res))
  expect_true(any(grepl("Variance basis: posterior medians", out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Model 1: %.6f", stats::median(v0)), out, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Model 2: %.6f", stats::median(v1)), out, fixed = TRUE)))
})

test_that("the PCV and the MOR read one between-stratum variance, the posterior median", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  m0 <- pcv_brms_stub(pcv_sd_null, family = list(family = "bernoulli", link = "logit"))
  v0 <- pcv_sd_null^2

  expect_identical(maihda_pcv_between_variance(m0), stats::median(v0))
  expect_identical(extract_between_variance(m0), stats::median(v0))
  expect_identical(maihda_mor_between_variance(m0), stats::median(v0))
  expect_equal(maihda_mor(m0), exp(sqrt(2 * stats::median(v0)) * stats::qnorm(0.75)),
               tolerance = 1e-12)
})

test_that("the PCV variance of an lme4 or ordinal fit is the fitted one, exactly", {
  skip_on_cran()
  fits <- pcv_with_seed(20261007, {
    d <- data.frame(g = sample(c("a", "b"), 400, replace = TRUE),
                    r = sample(c("x", "y", "z"), 400, replace = TRUE),
                    age = rnorm(400))
    d <- make_strata(d, vars = c("g", "r"))$data
    u <- stats::rnorm(nlevels(factor(d$stratum)), sd = 0.6)[as.integer(factor(d$stratum))]
    d$y <- 1 + 0.3 * d$age + u + rnorm(400)
    d$yb <- stats::rbinom(400, 1, stats::plogis(-0.2 + 1.5 * u))
    d$yo <- cut(d$y, stats::quantile(d$y, c(0, 0.3, 0.7, 1)), include.lowest = TRUE,
                ordered_result = TRUE)
    list(
      gaussian = fit_maihda(y ~ age + (1 | stratum), data = d),
      binomial = suppressMessages(
        fit_maihda(yb ~ age + (1 | stratum), data = d, family = "binomial")),
      ordinal = if (requireNamespace("ordinal", quietly = TRUE)) {
        suppressMessages(suppressWarnings(
          fit_maihda(yo ~ age + (1 | stratum), data = d, engine = "ordinal")))
      }
    )
  })
  for (m in fits[c("gaussian", "binomial")]) {
    expect_identical(maihda_pcv_between_variance(m), extract_between_variance(m))
    expect_identical(maihda_pcv_between_variance(m), maihda_stratum_variance_lme4(m$model))
  }

  skip_if_not_installed("ordinal")
  m_o <- fits$ordinal
  expect_identical(maihda_pcv_between_variance(m_o), extract_between_variance(m_o))
  expect_identical(maihda_pcv_between_variance(m_o), maihda_clmm_variances(m_o)$stratum)
})

test_that("stepwise_pcv() and pcv_importance() tabulate posterior medians for brms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  pcv_mock_brms()
  # Stan-free stand-in for the fits: a stub whose draws depend on how many predictors
  # the step's formula carries.
  fits <- list(pcv_sd_null, pcv_sd_mid, pcv_sd_adj)
  local_mocked_bindings(fit_maihda = function(formula, data, ...) {
    k <- length(attr(stats::terms(reformulas::nobars(formula)), "term.labels"))
    pcv_brms_stub(fits[[k + 1L]], formula)
  })
  d <- pcv_with_seed(4, data.frame(stratum = factor(rep(seq_len(6), each = 10)),
                                   a = rnorm(60), b = rnorm(60), y = rnorm(60)))
  med <- vapply(fits, function(s) stats::median(s^2), numeric(1))
  avg <- vapply(fits, function(s) mean(s^2), numeric(1))

  out <- suppressMessages(stepwise_pcv(d, "y", c("a", "b"), engine = "brms"))
  expect_identical(out$Variance, med)
  expect_equal(out$Total_PCV, 1 - med / med[1], tolerance = 1e-12)
  expect_equal(out$Step_PCV[-1], 1 - med[-1] / med[-3], tolerance = 1e-12)
  expect_identical(attr(out, "estimation_used"), "posterior")
  # The posterior means would give a different table.
  expect_gt(max(abs(out$Variance - avg)), 0.05)
  printed <- capture.output(print(out))
  expect_true(any(grepl("Variance basis: posterior medians", printed, fixed = TRUE)))
  # Nothing in this table is a posterior mean, so nothing is said about one.
  expect_false(any(grepl("posterior means", printed, fixed = TRUE)))

  imp <- suppressMessages(pcv_importance(d, "y", c("a", "b"), engine = "brms"))
  expect_equal(imp$total_pcv, 1 - med[3] / med[1], tolerance = 1e-12)
})

test_that("printed tables of a brms fit carry one basis line", {
  # stepwise_pcv(): Variance, Context_Variance and the MOR all rest on posterior
  # medians, so the basis line covers the table and no column is singled out.
  sw <- structure(
    data.frame(Step = 0:1, Model = c("Null Model", "Model 1"),
               Added_Variable = c("None (Intercept only)", "a"),
               Variance = c(0.30, 0.06), Context_Variance = c(0.20, 0.21),
               Step_PCV = c(0, 0.8), Total_PCV = c(0, 0.8), MOR = c(1.7, 1.3),
               stringsAsFactors = FALSE),
    estimation = "fitted", estimation_used = "posterior",
    boundary_steps = character(0), undefined_step_pcv = character(0),
    class = c("maihda_stepwise", "data.frame"))
  out <- capture.output(print(sw))
  expect_true(any(grepl("Variance basis: posterior medians", out, fixed = TRUE)))
  expect_false(any(grepl("posterior means", out, fixed = TRUE)))
  no_context <- sw
  no_context$Context_Variance <- NULL      # `$<-` keeps the table's attributes
  out_nc <- capture.output(print(no_context))
  expect_true(any(grepl("Variance basis: posterior medians", out_nc, fixed = TRUE)))
  expect_false(any(grepl("posterior means", out_nc, fixed = TRUE)))
  attr(sw, "estimation_used") <- "fitted"
  expect_false(any(grepl("posterior", capture.output(print(sw)), fixed = TRUE)))

  # compare_maihda_groups(): the variance columns are each group's summary() variances,
  # posterior medians on brms as the PCV's are, so the PCV's label serves the table.
  basis_line <- paste0("Variance basis: ", maihda_pcv_basis_label("posterior"))

  cmp <- structure(
    data.frame(group = c("A", "B"), n = c(100L, 120L), n_strata = c(8L, 8L),
               vpc = c(0.10, 0.12), var_between = c(0.40, 0.45), var_other = c(0, 0),
               var_residual = c(3.29, 3.29), status = c("ok", "ok"), pcv = c(0.7, 0.6),
               var_between_adjusted = c(0.12, 0.18),
               var_between_adjusted_ml = c(0.12, 0.18), pcv_status = c("ok", "ok"),
               stringsAsFactors = FALSE),
    group_var = "grp", engine = "brms", family = "binomial", shared_strata = TRUE,
    estimation = "fitted", estimation_used = "posterior",
    class = c("maihda_group_comparison", "data.frame"))
  out <- capture.output(print(cmp))
  expect_true(any(grepl(basis_line, out, fixed = TRUE)))
  expect_false(any(grepl("posterior mean", out, fixed = TRUE)))
  # A comparison without the PCV columns (strata from one dimension) reads the same:
  # its var_between is a posterior median too.
  one_dim <- cmp
  one_dim$pcv <- one_dim$var_between_adjusted <- one_dim$var_between_adjusted_ml <-
    one_dim$pcv_status <- NULL
  out1 <- capture.output(print(one_dim))
  expect_true(any(grepl(basis_line, out1, fixed = TRUE)))
  expect_false(any(grepl("posterior mean", out1, fixed = TRUE)))
})

test_that("the printed basis and the bootstrap refusal name the posterior medians", {
  lbl <- maihda_pcv_basis_label("posterior")
  expect_true(grepl("posterior medians", lbl, fixed = TRUE))
  expect_false(grepl("ML-refit", lbl, fixed = TRUE))
  # The other bases keep their wording.
  expect_identical(maihda_pcv_basis_label("fitted"),
                   "as fitted (REML for Gaussian lmer, matching summary())")
  expect_true(grepl("ML-refit", maihda_pcv_basis_label("ML"), fixed = TRUE))

  fake <- structure(list(engine = "brms"), class = "maihda_model")
  msg <- tryCatch(bootstrap_pcv(fake, fake, 10, 0.95), error = conditionMessage)
  expect_true(grepl("point estimate", msg, fixed = TRUE))
  expect_true(grepl("posterior-median", msg, fixed = TRUE))
  expect_false(grepl("posterior-mean", msg, fixed = TRUE))
})

test_that("maihda_table() attaches no PCV note to a brms analysis", {
  # The table's variance rows are summary()'s posterior medians, the values the PCV
  # compares, so there is no second basis to explain.
  rows <- function(bv) list("Between-stratum variance" = c(bv, NA_real_, NA_real_))
  shown <- list(null = rows(0.40), adjusted = rows(0.20))
  agree <- list(pcv = 0.5, var_model1 = 0.40, var_model2 = 0.20,
                estimation_used = "posterior")
  keys <- c("null", "adjusted")
  expect_null(maihda_table_pcv_note(shown, agree, keys, "brms"))

  # The lme4 note, for a PCV on another basis than the rows, keeps its wording, and it
  # is a note about lme4's estimation bases: no other engine is ever handed it, even
  # with variances that differ from the rows.
  other <- list(pcv = 0.79, var_model1 = 0.30, var_model2 = 0.063, estimation_used = "ML")
  ml <- maihda_table_pcv_note(shown, other, keys, "lme4")
  expect_true(grepl("maximum-likelihood", ml, fixed = TRUE))
  expect_false(grepl("posterior", ml, fixed = TRUE))
  for (engine in c("brms", "wemix", "ordinal")) {
    expect_null(maihda_table_pcv_note(shown, other, keys, engine), info = engine)
  }
  expect_null(maihda_table_pcv_note(shown, other, "estimate", "lme4"))
  expect_null(maihda_table_pcv_note(shown, NULL, keys, "lme4"))
})

test_that("?calculate_pcv says what a brms PCV compares and which way it leans", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  rd <- gsub("[[:space:]]+", " ",
             paste(readLines(file.path(man, "calculate_pcv.Rd"), warn = FALSE),
                   collapse = " "))
  expect_true(grepl("\\strong{brms fits.}", rd, fixed = TRUE))
  sec <- sub(".*\\\\strong\\{brms fits\\.\\}", "", rd)
  sec <- sub("\\\\examples.*", "", sec)
  expect_lt(nchar(sec), nchar(rd))               # the paragraph was actually isolated
  expect_true(grepl("posterior medians", sec, fixed = TRUE))
  expect_true(grepl("\\emph{understated} with few strata", sec, fixed = TRUE))
  expect_false(grepl("posterior-mean", rd, fixed = TRUE))
})

test_that("real brms fits: the PCV is built from the draws' medians", {
  # Compiles two Stan models, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # Eight strata, 100 Bernoulli rows each: additive main effects plus a small
  # stratum effect. brms draws its sampler seed from R's stream, so the fits sit
  # inside the seeded block too.
  a <- pcv_with_seed(20261007, {
    g <- expand.grid(d1 = c("a", "b"), d2 = c("c", "d"), d3 = c("e", "f"),
                     stringsAsFactors = FALSE)
    eta <- -0.9 + 0.5 * (g$d1 == "b") + 0.45 * (g$d2 == "d") + 0.55 * (g$d3 == "f") +
      stats::rnorm(8, 0, 0.2)
    d <- g[rep(seq_len(8), each = 100), ]
    d$y <- stats::rbinom(800, 1, stats::plogis(eta[rep(seq_len(8), each = 100)]))
    rownames(d) <- NULL
    suppressWarnings(suppressMessages(
      maihda(y ~ 1 + (1 | d1:d2:d3), data = d, engine = "brms", family = "binomial",
             chains = 2, iter = 1000, refresh = 0)))
  })
  v0 <- as.data.frame(a$model$model)$sd_stratum__Intercept^2
  v1 <- as.data.frame(a$model_adjusted$model)$sd_stratum__Intercept^2

  expect_identical(a$pcv$var_model1, stats::median(v0))
  expect_identical(a$pcv$var_model2, stats::median(v1))
  expect_equal(a$pcv$pcv, 1 - stats::median(v1) / stats::median(v0), tolerance = 1e-12)
  expect_identical(a$pcv$estimation_used, "posterior")
  expect_identical(calculate_pcv(a$model, a$model_adjusted)$pcv, a$pcv$pcv)
  # The skew is real on a fit like this one: the mean of each variance lies above its
  # median, by more in the adjusted model.
  expect_gt(mean(v0), stats::median(v0))
  expect_gt(mean(v1) / stats::median(v1), mean(v0) / stats::median(v0))

  # summary()'s table and the MOR hold the same posterior medians.
  vb <- a$summary$variance_components
  expect_identical(vb$variance[vb$component == "Between-stratum (random)"],
                   stats::median(v0))
  expect_identical(extract_between_variance(a$model), stats::median(v0))
  expect_equal(maihda_mor(a$model), exp(sqrt(2 * stats::median(v0)) * stats::qnorm(0.75)),
               tolerance = 1e-12)

  # The report table shows the variances the PCV compares, so its PCV row is the
  # reduction its two variance rows imply and it needs no note; the PCV has no
  # interval to report.
  tab <- suppressWarnings(suppressMessages(maihda_table(a)))
  bv <- tab$models[tab$models$statistic == "Between-stratum variance", ]
  expect_identical(c(bv$null, bv$adjusted), c(stats::median(v0), stats::median(v1)))
  expect_equal(tab$models$adjusted[tab$models$statistic == "PCV (null -> adjusted)"],
               1 - bv$adjusted / bv$null, tolerance = 1e-12)
  expect_null(tab$models_note)
  g <- suppressWarnings(glance(a))
  expect_true(is.na(g$pcv.conf.low) && is.na(g$pcv.conf.high))

  # Without draws the variance comes from brms's own median SD (robust = TRUE).
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) stop("no draws"))
  expect_equal(maihda_pcv_between_variance(a$model), stats::median(sqrt(v0))^2,
               tolerance = 1e-10)
  expect_equal(extract_between_variance(a$model), stats::median(sqrt(v0))^2,
               tolerance = 1e-10)
})
