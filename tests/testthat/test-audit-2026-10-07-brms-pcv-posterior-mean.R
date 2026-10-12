# Audit 2026-10-07 (brms-pcv-posterior-mean): the PCV of engine = "brms" fits compares
# the posterior MEDIANS of the between-stratum variances.
#
# The posterior of a between-stratum variance is right-skewed when strata are few, the
# adjusted model's more than the null's, so one minus a ratio of two posterior MEANS
# understates the PCV; and a posterior mean is set by an upper tail the sampler visits
# rarely, so that ratio also moves from one sampler seed to the next. Measured against
# a known truth of 0.84 at 8 strata: 0.58 from posterior means, 0.74 from medians.
#
# The blocks below pin the summary on hand-built draws (no Stan), check that only the
# PCV's variances moved -- the MOR's plug-in, summary()'s table and the other engines
# did not -- that the printed tables say which summary each column holds, and repeat
# the check on real fits when MAIHDA_TEST_BRMS=true.

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

test_that("the brms stratum variance is its posterior mean or median, as asked", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  fit <- pcv_brms_stub(pcv_sd_adj)$model
  v <- pcv_sd_adj^2

  expect_identical(maihda_stratum_variance_brms(fit), mean(v))
  expect_identical(maihda_stratum_variance_brms(fit, point = "mean"), mean(v))
  expect_identical(maihda_stratum_variance_brms(fit, point = "median"), stats::median(v))
  # The median does not depend on whether it is taken on the SD or on the variance.
  expect_equal(maihda_stratum_variance_brms(fit, point = "median"),
               stats::median(pcv_sd_adj)^2, tolerance = 1e-12)
  # Fixture premise: skewed enough that the two summaries cannot be confused.
  expect_gt(mean(v), 3 * stats::median(v))
  # No pattern: match.arg()'s message is translated by R.
  expect_error(maihda_stratum_variance_brms(fit, point = "mode"))
})

test_that("without draws the brms stratum variance falls back to brms's SD summary", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) stop("no draws"))
  fit <- pcv_brms_stub(pcv_sd_adj)$model

  expect_equal(maihda_stratum_variance_brms(fit), mean(pcv_sd_adj)^2, tolerance = 1e-12)
  expect_equal(maihda_stratum_variance_brms(fit, point = "median"),
               stats::median(pcv_sd_adj)^2, tolerance = 1e-12)
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
  expect_identical(res$pvc, res$pcv)
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

test_that("only the PCV reads the median: the MOR's variance keeps the posterior mean", {
  skip_if_not_installed("brms")
  pcv_mock_brms()
  m0 <- pcv_brms_stub(pcv_sd_null, family = list(family = "bernoulli", link = "logit"))
  v0 <- pcv_sd_null^2

  expect_identical(maihda_pcv_between_variance(m0), stats::median(v0))
  expect_identical(extract_between_variance(m0), mean(v0))
  expect_identical(extract_between_variance(m0, brms_point = "mean"), mean(v0))
  expect_identical(extract_between_variance(m0, brms_point = "median"), stats::median(v0))
  expect_identical(maihda_mor_between_variance(m0), mean(v0))
  expect_equal(maihda_mor(m0), exp(sqrt(2 * mean(v0)) * stats::qnorm(0.75)),
               tolerance = 1e-12)
})

test_that("the PCV variance of an lme4 or ordinal fit is the fitted one, exactly", {
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
    expect_identical(extract_between_variance(m, brms_point = "median"),
                     extract_between_variance(m))
    expect_identical(maihda_pcv_between_variance(m), maihda_stratum_variance_lme4(m$model))
  }

  skip_if_not_installed("ordinal")
  m_o <- fits$ordinal
  expect_identical(maihda_pcv_between_variance(m_o), extract_between_variance(m_o))
  expect_identical(maihda_pcv_between_variance(m_o), maihda_clmm_variances(m_o)$stratum)
})

test_that("stepwise_pcv() and pcv_importance() tabulate posterior medians for brms", {
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

  imp <- suppressMessages(suppressWarnings(
    pcv_importance(d, "y", c("a", "b"), engine = "brms", method = "sequential")))
  expect_equal(imp$total_pcv, 1 - med[3] / med[1], tolerance = 1e-12)
})

test_that("printed tables do not pass a posterior-mean variance off as a median", {
  # stepwise_pcv(): the basis line is about the Variance column; on a brms fit the
  # context variance and the MOR are read from posterior means and are named as such.
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
  expect_true(any(grepl("Context_Variance and MOR: from posterior means.", out, fixed = TRUE)))
  no_context <- sw
  no_context$Context_Variance <- NULL      # `$<-` keeps the table's attributes
  expect_true(any(grepl("^MOR: from posterior means\\.", capture.output(print(no_context)))))
  attr(sw, "estimation_used") <- "fitted"
  expect_false(any(grepl("posterior", capture.output(print(sw)), fixed = TRUE)))

  # compare_maihda_groups(): var_between comes from each group's summary(), a posterior
  # mean on brms, so the header does not call the table's variances medians.
  with_pcv <- maihda_group_basis_label("posterior", has_pcv = TRUE)
  expect_true(grepl("var_between is a posterior mean", with_pcv, fixed = TRUE))
  expect_true(grepl("the PCV compares posterior medians", with_pcv, fixed = TRUE))
  no_pcv <- maihda_group_basis_label("posterior", has_pcv = FALSE)
  expect_identical(no_pcv, "as fitted (brms posterior summaries, not maximum likelihood)")
  for (basis in c("fitted", "ML", "mixed")) {
    expect_identical(maihda_group_basis_label(basis, TRUE), maihda_pcv_basis_label(basis))
    expect_identical(maihda_group_basis_label(basis, FALSE), maihda_pcv_basis_label(basis))
  }

  cmp <- structure(
    data.frame(group = c("A", "B"), n = c(100L, 120L), n_strata = c(8L, 8L),
               vpc = c(0.10, 0.12), var_between = c(0.40, 0.45), var_other = c(0, 0),
               var_residual = c(3.29, 3.29), status = c("ok", "ok"), pcv = c(0.7, 0.6),
               var_between_adjusted = c(0.12, 0.18),
               var_between_adjusted_ml = c(0.06, 0.09), pcv_status = c("ok", "ok"),
               stringsAsFactors = FALSE),
    group_var = "grp", engine = "brms", family = "binomial", shared_strata = TRUE,
    estimation = "fitted", estimation_used = "posterior",
    class = c("maihda_group_comparison", "data.frame"))
  out <- capture.output(print(cmp))
  expect_true(any(grepl(paste0("Variance basis: ", with_pcv), out, fixed = TRUE)))
  # A comparison without the PCV columns (strata from one dimension) mentions no median.
  one_dim <- cmp
  one_dim$pcv <- one_dim$var_between_adjusted <- one_dim$var_between_adjusted_ml <-
    one_dim$pcv_status <- NULL
  out1 <- capture.output(print(one_dim))
  expect_true(any(grepl(paste0("Variance basis: ", no_pcv), out1, fixed = TRUE)))
  expect_false(any(grepl("median", out1, fixed = TRUE)))
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

test_that("maihda_table() says a brms PCV is not the ratio of its variance rows", {
  # The table's variance rows are summary()'s posterior means; the PCV compares medians.
  rows <- function(bv) list("Between-stratum variance" = c(bv, NA_real_, NA_real_))
  shown <- list(null = rows(0.40), adjusted = rows(0.20))
  pcv <- list(pcv = 0.79, var_model1 = 0.30, var_model2 = 0.063,
              estimation_used = "posterior")
  keys <- c("null", "adjusted")

  note <- maihda_table_pcv_note(shown, pcv, keys, "brms")
  expect_type(note, "character")
  expect_true(grepl("posterior medians", note, fixed = TRUE))
  # The variance row is a posterior mean; the SD row is its square root, not a mean SD.
  expect_true(grepl("variance row is that variance's posterior mean", note, fixed = TRUE))
  expect_true(grepl("SD row its square root", note, fixed = TRUE))
  expect_false(grepl("maximum-likelihood", note, fixed = TRUE))

  # No note when the PCV's variances are the rows shown, without a null/adjusted pair,
  # without a PCV, or on the engines whose PCV reads the fitted variances.
  agree <- list(pcv = 0.5, var_model1 = 0.40, var_model2 = 0.20)
  expect_null(maihda_table_pcv_note(shown, agree, keys, "brms"))
  expect_null(maihda_table_pcv_note(shown, pcv, "estimate", "brms"))
  expect_null(maihda_table_pcv_note(shown, NULL, keys, "brms"))
  expect_null(maihda_table_pcv_note(shown, pcv, keys, "wemix"))
  expect_null(maihda_table_pcv_note(shown, pcv, keys, "ordinal"))

  # The lme4 note keeps its maximum-likelihood wording.
  ml <- maihda_table_pcv_note(shown, c(pcv[1:3], list(estimation_used = "ML")), keys, "lme4")
  expect_true(grepl("maximum-likelihood", ml, fixed = TRUE))
  expect_false(grepl("posterior", ml, fixed = TRUE))
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

  # Nothing else moved: summary()'s table and the MOR still read the posterior mean.
  vb <- a$summary$variance_components
  expect_equal(vb$variance[vb$component == "Between-stratum (random)"], mean(v0),
               tolerance = 1e-8)
  expect_identical(extract_between_variance(a$model), mean(v0))
  expect_equal(maihda_mor(a$model), exp(sqrt(2 * mean(v0)) * stats::qnorm(0.75)),
               tolerance = 1e-12)

  # The report table shows those posterior-mean variances, so it carries the note; and
  # the PCV has no interval to report.
  tab <- suppressWarnings(suppressMessages(maihda_table(a)))
  bv <- tab$models[tab$models$statistic == "Between-stratum variance", ]
  expect_equal(c(bv$null, bv$adjusted), c(mean(v0), mean(v1)), tolerance = 1e-8)
  expect_true(grepl("posterior medians", tab$models_note, fixed = TRUE))
  g <- suppressWarnings(glance(a))
  expect_true(is.na(g$pcv.conf.low) && is.na(g$pcv.conf.high))

  # Without draws the two summaries come from brms's own SD summary: its median SD
  # (robust = TRUE) for the PCV, its mean SD otherwise.
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) stop("no draws"))
  expect_equal(maihda_pcv_between_variance(a$model), stats::median(sqrt(v0))^2,
               tolerance = 1e-10)
  expect_equal(extract_between_variance(a$model), mean(sqrt(v0))^2, tolerance = 1e-10)
})
