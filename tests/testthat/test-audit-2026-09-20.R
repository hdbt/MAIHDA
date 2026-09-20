# Audit 2026-09-20 -- na.exclude predictions combined with unpadded analytic rows.
#
# FINDING (external audit, item 6, P2). A fit made with na.action = na.exclude keeps
# an "exclude"-class na.action on its model frame, so base R PADS every fitted-row
# accessor -- stats::predict(), fitted(), residuals(), weights(type = "prior") -- back
# out to the ORIGINAL input rows, reinserting an NA at each dropped position. But
# object$data IS the analytic model frame. Consumers pairing the two were misaligned.
#
# Reproduced exactly as reported (720 input rows, 717 analytic):
#   predict_maihda()                720 values against 717 rows of object$data
#   stratum predictions             "arguments imply differing number of rows: 717, 720"
#   maihda_discriminatory_accuracy  "'prob' and 'y' must have the same length (720 vs 717)"
#   maihda_table()                  ranked-strata table dropped by its error handler
#
# ...and WIDER than reported, in ways that are SILENT rather than loud:
#   * when the analytic count DIVIDES the padded count (720 input / 360 analytic) R
#     recycles instead of erroring: stratum predictions came back with no error and no
#     warning, off by up to 0.0157 on the probability scale, with every stratum's n and
#     w_sum exactly DOUBLED (114 vs 57, 118 vs 59, ...);
#   * the zero-inflation adequacy check vanished entirely (getME(model, "y") is not
#     padded but fitted() is, so every expected-zero count went NA);
#   * a growth count fit's summary()$count_vpc was biased -- lambda 7.052 against the
#     correct 7.826 (9.9%), level1_variance 0.1326 against 0.1203 -- because its
#     per-row latent variance v_i is row-varying there and got recycled;
#   * plot(type = "all") silently dropped four panels (obs_vs_shrunken, predicted,
#     effect_decomp, prediction_deviation).
#
# NOT a defect, measured rather than assumed: the INTERCEPT-ONLY count VPC is
# unaffected (lambda 2.313946712603 either way). mu and the prior weights are BOTH
# padded there and v_i is constant, so maihda_weighted_obs_mean()'s is.finite() filter
# drops exactly the same positions. The auditor was right not to count their
# weighted-Poisson probe: those recycling warnings are real but numerically inert
# until a random slope makes v_i vary by row.
#
# Engines: lme4 only. WeMix refuses na.action outright; the ordinal engine accepts it
# but is immune (it builds object$data as complete cases and predicts from object$data).
# brms is untested here -- no Stan compiler available.
#
# FIX: one internal fitted-row accessor family -- maihda_fit_rows_predict/fitted/
# residuals/weights(), all built on maihda_unpad_fit_rows() -- returning the ANALYTIC
# rows, routed through all 23 internal call sites. A length that is neither the
# analytic nor the padded one is an error naming both counts -- but only for values
# that REACH the accessor; a site still calling stats::predict() directly never does,
# which is what the na.exclude-equals-na.omit equivalence test below is for. Per the
# maintainer the
# public predict_maihda() also returns the analytic rows, so every predict_maihda()
# return has one length and can be bound to object$data.

# -- shared fixtures ----------------------------------------------------------

# A binomial MAIHDA over `n` rows with `miss` of them missing a covariate. `scatter`
# picks whether the missing rows are contiguous or spread through the frame: with
# `miss` = n/2 the padded length is an exact multiple of the analytic length, which is
# the case R recycles silently, and only SCATTERED missingness then mispairs the rows.
audit0920_binomial <- function(n = 720, miss = 3, scatter = TRUE, seed = 1) {
  set.seed(seed)
  d <- data.frame(
    sex = factor(sample(c("F", "M"), n, TRUE)),
    edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
    x   = rnorm(n))
  st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.7)
  d$y <- rbinom(n, 1, stats::plogis(-0.2 + 0.5 * d$x + u[as.integer(st)]))
  if (miss > 0) {
    set.seed(seed + 98)
    d$x[if (scatter) sample(n, miss) else seq_len(miss)] <- NA
  }
  d
}

audit0920_pair <- function(d, ...) {
  list(
    omit = fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial(),
                      na.action = stats::na.omit, ...),
    exclude = fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial(),
                         na.action = stats::na.exclude, ...))
}

# -- the accessor and its guard ----------------------------------------------

test_that("maihda_unpad_fit_rows() strips na.exclude padding and is a no-op otherwise", {
  skip_on_cran()
  f <- audit0920_pair(audit0920_binomial())

  # Premise of the whole finding: the fit frame carries an "exclude" na.action, is the
  # ANALYTIC rows, and base R pads the accessors back out to the input rows.
  expect_s3_class(attr(f$exclude$model@frame, "na.action"), "exclude")
  expect_identical(nrow(f$exclude$data), 717L)
  expect_identical(length(stats::predict(f$exclude$model)), 720L)
  expect_identical(sum(is.na(stats::predict(f$exclude$model))), 3L)
  # ... and that na.omit records the SAME dropped rows but never pads, which is why
  # it was always correct and why the accessor keys off the class, not the presence.
  expect_s3_class(attr(f$omit$model@frame, "na.action"), "omit")
  expect_false(inherits(attr(f$omit$model@frame, "na.action"), "exclude"))
  expect_identical(length(stats::predict(f$omit$model)), 717L)

  ex <- MAIHDA:::maihda_na_exclude_rows(f$exclude$model)
  expect_identical(ex$n_fit, 717L)
  expect_identical(length(ex$omit), 3L)
  expect_null(MAIHDA:::maihda_na_exclude_rows(f$omit$model))

  padded <- stats::predict(f$exclude$model)
  stripped <- MAIHDA:::maihda_unpad_fit_rows(padded, f$exclude$model)
  expect_identical(length(stripped), 717L)
  expect_false(anyNA(stripped))
  expect_equal(unname(stripped), unname(padded[-ex$omit]))
  # Positional correctness, not just the count: the survivors keep the frame's rows.
  expect_identical(names(stripped), rownames(f$exclude$model@frame))

  # Already-analytic input is returned untouched (getME(model, "y") never pads).
  y <- as.numeric(lme4::getME(f$exclude$model, "y"))
  expect_identical(MAIHDA:::maihda_unpad_fit_rows(y, f$exclude$model), y)

  # A matrix is stripped by rows...
  m <- matrix(seq_len(720 * 2), nrow = 720)
  expect_identical(nrow(MAIHDA:::maihda_unpad_fit_rows(m, f$exclude$model)), 717L)
  # ... and so is a data frame. length() of one is its COLUMN count, which would slip
  # past both the analytic and the padded test and return it unstripped.
  df <- data.frame(a = seq_len(720), b = seq_len(720))
  expect_identical(nrow(MAIHDA:::maihda_unpad_fit_rows(df, f$exclude$model)), 717L)

  # A length that is neither analytic nor padded cannot be aligned by any rule, so it
  # errors naming both counts instead of recycling into a wrong answer. This guards
  # values passed THROUGH the accessor, not call sites that bypass it.
  expect_error(
    MAIHDA:::maihda_unpad_fit_rows(rep(1, 42), f$exclude$model, "prediction"),
    "42 row\\(s\\).*717 analytic row\\(s\\).*720")
  # Under na.omit the accessor never intervenes, so nothing can be perturbed.
  expect_identical(MAIHDA:::maihda_unpad_fit_rows(rep(1, 42), f$omit$model), rep(1, 42))
  expect_identical(MAIHDA:::maihda_unpad_fit_rows(NULL, f$exclude$model), NULL)
})

# -- the three sites the finding named ---------------------------------------

test_that("predict_maihda() returns the analytic rows under na.exclude", {
  skip_on_cran()
  f <- audit0920_pair(audit0920_binomial())

  for (ty in c("response", "link")) {
    p <- predict_maihda(f$exclude, type = ty)
    # The bug: 720 values (base R's na.exclude padding) against 717 rows of $data.
    expect_identical(length(p), nrow(f$exclude$data))
    expect_false(anyNA(p))
    expect_equal(unname(p), unname(predict_maihda(f$omit, type = ty)))
  }
  # The point of aligning them: this is what a caller actually does.
  expect_identical(
    nrow(cbind(f$exclude$data, pred = predict_maihda(f$exclude))),
    717L)
})

test_that("stratum predictions no longer error under na.exclude", {
  skip_on_cran()
  f <- audit0920_pair(audit0920_binomial())
  # The bug: "arguments imply differing number of rows: 717, 720".
  got <- expect_silent(
    MAIHDA:::maihda_stratum_predictions_lme4(f$exclude, summary(f$exclude)))
  want <- MAIHDA:::maihda_stratum_predictions_lme4(f$omit, summary(f$omit))
  expect_equal(got, want)
})

test_that("maihda_table() keeps its ranked-strata table under na.exclude", {
  skip_on_cran()
  f <- audit0920_pair(audit0920_binomial())
  tb <- maihda_table(f$exclude)
  # The bug: the ranking error was swallowed by tryCatch(error = function(e) NULL),
  # so the table came back with its model rows and no strata at all.
  expect_false(is.null(tb$strata))
  expect_equal(tb$strata, maihda_table(f$omit)$strata)
  expect_identical(tb$n_obs, maihda_table(f$omit)$n_obs)
})

test_that("discriminatory accuracy no longer errors under na.exclude", {
  skip_on_cran()
  f <- audit0920_pair(audit0920_binomial())
  # The bug: "'prob' and 'y' must have the same length (720 vs 717)", which summary()
  # turned into a warning and an omitted AUC/MOR.
  da <- maihda_discriminatory_accuracy(f$exclude)
  expect_equal(da$auc, maihda_discriminatory_accuracy(f$omit)$auc)
  expect_true(is.finite(da$auc))
  s <- expect_silent(summary(f$exclude))
  expect_false(is.null(s$discriminatory_accuracy))
  expect_equal(s$discriminatory_accuracy$auc, summary(f$omit)$discriminatory_accuracy$auc)
})

# -- the SILENT cases the finding did not reach ------------------------------

test_that("a padded length that divides the analytic one no longer recycles silently", {
  skip_on_cran()
  # 720 input rows, 360 analytic: 720 = 2 * 360, so R recycled WITHOUT a warning and
  # the stratum table came back looking perfectly healthy.
  d <- audit0920_binomial(n = 720, miss = 360, scatter = TRUE)
  f <- audit0920_pair(d)
  expect_identical(nrow(f$exclude$data), 360L)

  got <- expect_silent(
    MAIHDA:::maihda_stratum_predictions_lme4(f$exclude, summary(f$exclude)))
  want <- MAIHDA:::maihda_stratum_predictions_lme4(f$omit, summary(f$omit))

  # The bug's two fingerprints, both silent: displaced predictions (up to 0.0157 on
  # the probability scale) and stratum sizes that were exactly DOUBLED.
  expect_equal(got$predicted_row, want$predicted_row)
  expect_equal(got$n, want$n)
  expect_equal(got$n, c(57, 59, 66, 41, 61, 76))
  expect_equal(sum(got$n), nrow(f$exclude$data))
  expect_equal(got$w_sum, want$w_sum)

  # Contiguous missingness recycled into the CORRECT pairing by luck (the NA half sat
  # ahead of the data half), so only `n` was wrong there. Pin that it stays correct.
  f2 <- audit0920_pair(audit0920_binomial(n = 720, miss = 360, scatter = FALSE))
  g2 <- MAIHDA:::maihda_stratum_predictions_lme4(f2$exclude, summary(f2$exclude))
  w2 <- MAIHDA:::maihda_stratum_predictions_lme4(f2$omit, summary(f2$omit))
  expect_equal(g2, w2)
  expect_equal(sum(g2$n), nrow(f2$exclude$data))
})

test_that("the zero-inflation adequacy check survives na.exclude", {
  skip_on_cran()
  set.seed(3)
  n <- 720
  d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                  x   = rnorm(n))
  st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.7)
  d$y <- rpois(n, exp(1.2 + 0.5 * d$x + u[as.integer(st)]))
  d$x[c(5, 100, 500)] <- NA

  p <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::poisson(),
               na.action = na))
  ad <- lapply(p, function(f) MAIHDA:::maihda_adequacy_checks(f$model, NULL))

  # The bug: getME(model, "y") is NOT padded but fitted() IS, so every dpois(0, mu)
  # term went NA, exp0 was NA, and the whole check was dropped without a word.
  expect_false(is.null(ad$exclude$zeroinflation))
  expect_equal(ad$exclude$zeroinflation, ad$omit$zeroinflation)
  expect_identical(ad$exclude$zeroinflation$n, 717L)
  # Overdispersion was never wrong (its statistic drops NAs); pin that it stays put.
  expect_equal(ad$exclude$overdispersion, ad$omit$overdispersion)
})

test_that("a growth count fit's count_vpc is not biased by na.exclude", {
  skip_on_cran()
  # A random SLOPE is what makes the per-row latent variance v_i vary, and only then
  # does recycling eta against v move a number: lambda was 7.052 against 7.826.
  set.seed(21)
  N <- 300; TT <- 4
  ld <- data.frame(
    id   = rep(seq_len(N), each = TT),
    time = rep(0:(TT - 1), times = N),
    sex  = factor(rep(sample(c("F", "M"), N, TRUE), each = TT)),
    edu  = factor(rep(sample(c("low", "mid", "high"), N, TRUE), each = TT)))
  st <- interaction(ld$sex, ld$edu, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.4); b <- rnorm(nlevels(st), 0, 0.25)
  ri <- rnorm(N, 0, 0.3)
  ld$cnt <- rpois(nrow(ld), exp(1.1 + 0.15 * ld$time + u[as.integer(st)] +
                                  b[as.integer(st)] * ld$time + ri[ld$id]))
  ld$cov <- rnorm(nrow(ld)); ld$cov[c(7, 300, 900)] <- NA

  f <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    suppressWarnings(suppressMessages(
      fit_maihda(cnt ~ cov + (1 | sex:edu), data = ld, family = stats::poisson(),
                 id = "id", time = "time", na.action = na))))
  s <- lapply(f, function(x) suppressWarnings(summary(x)))

  # v_i really is row-varying here -- otherwise this test would pass vacuously.
  bars <- reformulas::findbars(stats::formula(f$omit$model))
  vc <- lme4::VarCorr(f$omit$model)
  v <- MAIHDA:::maihda_latent_re_variance_rows(
    bars, lapply(bars, function(bb) vc[[MAIHDA:::maihda_bar_group_name(bb)]]),
    MAIHDA:::maihda_model_frame(f$omit$model))
  expect_gt(stats::sd(v), 0.1)

  expect_equal(s$exclude$count_vpc$lambda, s$omit$count_vpc$lambda)
  expect_equal(s$omit$count_vpc$lambda, 7.825722551708, tolerance = 1e-8)
  expect_equal(s$exclude$count_vpc$level1_variance, s$omit$count_vpc$level1_variance)
  expect_equal(s$omit$count_vpc$level1_variance, 0.120254404, tolerance = 1e-7)
  expect_equal(s$exclude$count_vpc$alternatives, s$omit$count_vpc$alternatives)
})

test_that("the intercept-only count VPC was never affected", {
  skip_on_cran()
  # The auditor's own weighted-Poisson probe, settled: mu and the prior weights are
  # BOTH padded and v_i is constant, so the is.finite() filter drops the same
  # positions and lambda is identical. Recorded so a future change cannot quietly
  # turn this inert case into a live one.
  set.seed(3)
  n <- 720
  d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                  x   = rnorm(n))
  st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.7)
  eta <- -0.2 + 0.5 * d$x + u[as.integer(st)]
  d$y <- rpois(n, exp(eta + 1.2))
  set.seed(11); d$w <- round(runif(n, 0.5, 4), 3)
  d$x[c(5, 100, 500)] <- NA

  f <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::poisson(),
               weights = w, na.action = na))
  lam <- vapply(f, function(x) MAIHDA:::maihda_weighted_obs_mean(
    MAIHDA:::maihda_count_marginal_mu_lme4(x$model),
    MAIHDA:::maihda_fit_prior_weights(x$model)), numeric(1))
  expect_equal(unname(lam[["exclude"]]), unname(lam[["omit"]]))
  expect_equal(unname(lam[["omit"]]), 2.313946712603, tolerance = 1e-9)
})

test_that("no plot panel is dropped under na.exclude", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  set.seed(5)
  n <- 720
  d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                  x   = rnorm(n))
  st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.7)
  d$y <- 1 + 0.5 * d$x + u[as.integer(st)] + rnorm(n)
  d$x[c(5, 100, 500)] <- NA
  f <- list(
    omit = fit_maihda(y ~ x + (1 | sex:edu), data = d, na.action = stats::na.omit),
    exclude = fit_maihda(y ~ x + (1 | sex:edu), data = d,
                         na.action = stats::na.exclude))

  # All four errored before the fix; plot(type = "all") warned and left them out.
  for (ty in c("obs_vs_shrunken", "predicted", "effect_decomp", "prediction_deviation")) {
    p <- expect_silent(plot(f$exclude, type = ty))
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  expect_silent(plot(f$exclude, type = "all"))
  expect_equal(summary(f$exclude)$vpc$estimate, summary(f$omit)$vpc$estimate)
})

# -- the routes into na.exclude, and the end-to-end invariant ----------------

test_that("na.exclude reaches the fit by argument, abbreviation and global option", {
  skip_on_cran()
  d <- audit0920_binomial(n = 400, miss = 3)
  ref <- fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial(),
                    na.action = stats::na.omit)

  # (a) spelled out, (b) abbreviated -- R's own partial matching binds na.act to
  # glmer's na.action, so the padding arrives without the exact name appearing.
  by_arg <- fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial(),
                       na.action = stats::na.exclude)
  by_abbrev <- fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial(),
                          na.act = stats::na.exclude)
  # (c) never named at all: model.frame() defaults to getOption("na.action").
  old <- options(na.action = "na.exclude")
  on.exit(options(old), add = TRUE)
  by_option <- fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial())
  options(old)

  for (f in list(by_arg, by_abbrev, by_option)) {
    expect_s3_class(attr(f$model@frame, "na.action"), "exclude")
    expect_identical(length(predict_maihda(f)), nrow(f$data))
    expect_equal(unname(predict_maihda(f)), unname(predict_maihda(ref)))
    expect_equal(maihda_discriminatory_accuracy(f)$auc,
                 maihda_discriminatory_accuracy(ref)$auc)
  }
})

test_that("na.exclude is observationally equivalent to na.omit across the surface", {
  skip_on_cran()
  # The invariant the accessor family exists to hold, and the test that catches a
  # NEW unrouted accessor: nothing the package reports may depend on which of the two
  # NA actions was used, because they fit exactly the same rows.
  f <- audit0920_pair(audit0920_binomial(n = 400, miss = 5))
  expect_equal(lme4::fixef(f$exclude$model), lme4::fixef(f$omit$model))
  expect_identical(nrow(f$exclude$data), nrow(f$omit$data))

  se <- summary(f$exclude); so <- summary(f$omit)
  for (k in c("vpc", "variance_components", "decomposition", "stratum_estimates",
              "fixed_effects", "discriminatory_accuracy", "vpc_response")) {
    expect_equal(se[[k]], so[[k]], info = k)
  }
  expect_equal(maihda_mor(f$exclude), maihda_mor(f$omit))
  expect_equal(maihda_vpc_response(f$exclude, seed = 1),
               maihda_vpc_response(f$omit, seed = 1))
  expect_equal(maihda_auc(predict_maihda(f$exclude), f$exclude$data$y),
               maihda_auc(predict_maihda(f$omit), f$omit$data$y))
  expect_equal(maihda_table(f$exclude)$strata, maihda_table(f$omit)$strata)
  expect_equal(maihda_table(f$exclude)$models, maihda_table(f$omit)$models)
  expect_equal(MAIHDA:::maihda_weight_fingerprint(f$exclude$model),
               MAIHDA:::maihda_weight_fingerprint(f$omit$model))
  expect_equal(MAIHDA:::maihda_fit_prior_weights(f$exclude$model),
               MAIHDA:::maihda_fit_prior_weights(f$omit$model))
})

test_that("the padding is stripped correctly alongside subset and NA weights", {
  skip_on_cran()
  # An na.action records positions in the POST-subset, post-weight-drop frame, and
  # napredict() pads back only to that frame's length -- not to the original input.
  # Stripping has to use the very indices R padded with, so pin the awkward cases.
  mk <- function(seed) {
    set.seed(seed)
    n <- 720
    d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                    edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                    x   = rnorm(n))
    st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
    d$y <- 1 + 0.5 * d$x + rnorm(nlevels(st), 0, 0.7)[as.integer(st)] + rnorm(n)
    d
  }

  # (a) subset = : 720 rows in, 480 kept, 477 analytic. The padded length is 480.
  d <- mk(1); d$x[c(5, 100, 500)] <- NA
  d$keep <- rep(c(TRUE, TRUE, FALSE), length.out = 720)
  f <- list(
    omit = fit_maihda(y ~ x + (1 | sex:edu), data = d, subset = keep,
                      na.action = stats::na.omit),
    exclude = fit_maihda(y ~ x + (1 | sex:edu), data = d, subset = keep,
                         na.action = stats::na.exclude))
  expect_identical(nrow(f$exclude$data), 477L)
  expect_identical(length(stats::predict(f$exclude$model)), 480L)  # padded to the SUBSET
  expect_identical(length(predict_maihda(f$exclude)), 477L)
  expect_equal(unname(predict_maihda(f$exclude)), unname(predict_maihda(f$omit)))
  expect_equal(MAIHDA:::maihda_stratum_predictions_lme4(f$exclude, summary(f$exclude)),
               MAIHDA:::maihda_stratum_predictions_lme4(f$omit, summary(f$omit)))

  # (b) NA precision weights: lme4 drops those rows too, so the na.action covers both
  # kinds and the padding runs to the full 720.
  d2 <- mk(2); d2$x[c(7, 200)] <- NA
  set.seed(31); d2$w <- round(runif(720, 0.5, 4), 3); d2$w[c(11, 400, 650)] <- NA
  g <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    suppressWarnings(fit_maihda(y ~ x + (1 | sex:edu), data = d2, weights = w,
                                na.action = na)))
  expect_identical(nrow(g$exclude$data), 715L)
  expect_identical(length(stats::predict(g$exclude$model)), 720L)
  expect_equal(MAIHDA:::maihda_fit_prior_weights(g$exclude$model),
               MAIHDA:::maihda_fit_prior_weights(g$omit$model))
  expect_equal(unname(predict_maihda(g$exclude)), unname(predict_maihda(g$omit)))
  expect_equal(summary(g$exclude)$vpc, summary(g$omit)$vpc)

  # (c) no missing data at all: model.frame() records no na.action, so na.exclude is
  # a total no-op and cannot perturb the overwhelmingly common case.
  d3 <- mk(5)
  k <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    fit_maihda(y ~ x + (1 | sex:edu), data = d3, na.action = na))
  expect_null(attr(k$exclude$model@frame, "na.action"))
  expect_null(MAIHDA:::maihda_na_exclude_rows(k$exclude$model))
  expect_equal(unname(predict_maihda(k$exclude)), unname(predict_maihda(k$omit)))
  expect_equal(summary(k$exclude)[c("vpc", "variance_components", "stratum_estimates")],
               summary(k$omit)[c("vpc", "variance_components", "stratum_estimates")])
})

test_that("a contextual fit's intersectional-scope AUC survives na.exclude", {
  skip_on_cran()
  # maihda_da_scope_scores() is a separate predict() site, reached only by a fit that
  # carries a non-intersectional random effect -- so a plain (1 | sex:edu) fixture
  # never touches it.
  set.seed(4)
  n <- 900
  d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                  site = factor(sample(paste0("s", 1:12), n, TRUE)),
                  x = rnorm(n))
  st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
  u <- rnorm(nlevels(st), 0, 0.7); v <- rnorm(nlevels(d$site), 0, 0.5)
  d$y <- rbinom(n, 1, stats::plogis(-0.2 + 0.5 * d$x + u[as.integer(st)] +
                                      v[as.integer(d$site)]))
  d$x[c(5, 100, 500)] <- NA
  f <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    fit_maihda(y ~ x + (1 | sex:edu), data = d, family = stats::binomial(),
               context = "site", na.action = na))
  da <- lapply(f, maihda_discriminatory_accuracy)
  expect_identical(da$exclude$auc_scope, "intersectional")
  expect_equal(da$exclude$auc, da$omit$auc)
  expect_equal(da$exclude$auc_full, da$omit$auc_full)
  expect_equal(da$exclude$mor, da$omit$mor)
})

test_that("the parametric bootstrap runs and agrees under na.exclude", {
  skip_on_cran()
  # simulate.merMod() passes its draws through napredict(), so under na.exclude they
  # come back on the ORIGINAL rows. maihda_refit_draw() then tagged such a draw as
  # "already on the model frame", suppressing the very truncation it needed, and
  # EVERY refit failed: summary(bootstrap = TRUE) stopped with "All VPC bootstrap
  # refits failed; no interval can be computed." na.omit was unaffected throughout.
  set.seed(7)
  n <- 600
  d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                  age = factor(sample(c("y", "o"), n, TRUE)),
                  x   = rnorm(n))
  st <- interaction(d$sex, d$edu, d$age, sep = "_", drop = TRUE)
  d$y <- 1 + 0.5 * d$x + rnorm(nlevels(st), 0, 0.9)[as.integer(st)] + rnorm(n)
  d$x[c(5, 100, 300)] <- NA
  f <- lapply(list(omit = stats::na.omit, exclude = stats::na.exclude), function(na)
    fit_maihda(y ~ x + (1 | sex:edu:age), data = d, na.action = na))
  expect_identical(nrow(f$exclude$data), nrow(f$omit$data))

  # The premise: lme4 draws on the model frame and then pads.
  set.seed(99); raw_om <- stats::simulate(f$omit$model, 3)
  set.seed(99); raw_ex <- stats::simulate(f$exclude$model, 3)
  expect_identical(nrow(raw_om), 597L)
  expect_identical(nrow(raw_ex), 600L)
  omit_idx <- as.integer(attr(f$exclude$model@frame, "na.action"))
  expect_equal(unname(as.matrix(raw_ex[-omit_idx, , drop = FALSE])),
               unname(as.matrix(raw_om)))

  # The package's own simulator now returns the fitted rows, and the SAME draws.
  set.seed(99); sim_om <- MAIHDA:::maihda_simulate_lme4(f$omit$model, 3)
  set.seed(99); sim_ex <- MAIHDA:::maihda_simulate_lme4(f$exclude$model, 3)
  expect_identical(nrow(sim_ex), 597L)
  expect_equal(unname(as.matrix(sim_ex)), unname(as.matrix(sim_om)))

  # End to end. The bootstrap draws from the ambient RNG (summary()'s `seed` is
  # documented for the response-scale VPC only), so seed it the conventional way;
  # done so, the two NA actions agree exactly rather than merely to within
  # Monte-Carlo error.
  set.seed(42); b_om <- suppressWarnings(summary(f$omit, bootstrap = TRUE, n_boot = 12))
  set.seed(42); b_ex <- suppressWarnings(summary(f$exclude, bootstrap = TRUE, n_boot = 12))
  expect_true(is.finite(b_ex$vpc$ci_lower) && is.finite(b_ex$vpc$ci_upper))
  expect_equal(b_ex$vpc, b_om$vpc)
  expect_equal(b_ex$decomposition, b_om$decomposition)
  expect_equal(b_ex$variance_components, b_om$variance_components)

  set.seed(5); p_om <- suppressWarnings(summary(f$omit, df_method = "bootstrap",
                                                n_boot = 10))
  set.seed(5); p_ex <- suppressWarnings(summary(f$exclude, df_method = "bootstrap",
                                                n_boot = 10))
  expect_equal(p_ex$fixed_effects, p_om$fixed_effects)
})

test_that("the ordinal engine was immune and stays immune", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # It accepts na.action but builds object$data as complete cases itself and predicts
  # from object$data, so it never carried the padding. WeMix refuses na.action.
  set.seed(7)
  m <- 400
  d <- data.frame(sex = factor(sample(c("F", "M"), m, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), m, TRUE)),
                  x   = rnorm(m))
  st <- interaction(d$sex, d$edu, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.7)
  d$yo <- factor(cut(0.5 * d$x + u[as.integer(st)] + rnorm(m),
                     breaks = c(-Inf, -0.5, 0.5, Inf), labels = c("a", "b", "c")),
                 ordered = TRUE)
  d$x[1:3] <- NA
  f <- suppressWarnings(fit_maihda(yo ~ x + (1 | sex:edu), data = d,
                                   engine = "ordinal", na.action = stats::na.exclude))
  expect_identical(nrow(f$data), 397L)
  expect_identical(length(predict_maihda(f)), 397L)
})
