# Audit 2026-09-06 (second pass), finding A2: "prediction helpers discard the fitted
# transformation basis".
#
# maihda_lme4_fixed_link(), maihda_wemix_linpred() and maihda_clmm_linpred() rebuilt
# terms() from a BARE formula, which carries no "predvars" attribute. A data-dependent
# term -- scale(x), poly(x, 2), splines::ns(x, 3) -- was therefore RE-EVALUATED on the
# prediction batch, so:
#   * predicting all rows and subsetting disagreed with predicting the subset,
#   * a single-row (or otherwise constant) grid divided by a zero SD and gave NaN,
#   * the longitudinal count VPC(t) grid, which fixes every row to one time, returned
#     VPC = NaN at every time for a fit containing scale(<time>).
# The fix takes the basis from the FIT: terms(model, fixed.only = TRUE) for lme4, and
# maihda_fitted_predict_terms() -- the model frame rebuilt on the stored analytic rows
# -- for the WeMix and ordinal engines.

test_that("maihda_fitted_predict_terms recovers the fitted scale() basis", {
  set.seed(101)
  d <- data.frame(x = stats::rnorm(60, mean = 5, sd = 3),
                  g = factor(rep(c("a", "b"), 30)))
  d$y <- stats::rnorm(60)

  basis <- MAIHDA:::maihda_fitted_predict_terms(y ~ scale(x) + g, d)
  pv <- attr(basis$terms, "predvars")
  expect_false(is.null(pv))
  # The centre and scale are the FITTED data's, not any prediction batch's.
  expect_equal(as.numeric(pv[[2]]$center), mean(d$x), tolerance = 1e-12)
  expect_equal(as.numeric(pv[[2]]$scale), stats::sd(d$x), tolerance = 1e-12)
  expect_equal(basis$xlev, list(g = c("a", "b")))
  # Response removed, so the terms apply to newdata that lacks the response.
  expect_equal(as.integer(attr(basis$terms, "response")), 0L)

  # THE TEETH: the same expression evaluated WITHOUT predvars (the old path) picks up
  # the prediction batch's centre instead of the fitted one.
  bare <- stats::delete.response(stats::terms(y ~ scale(x) + g))
  expect_null(attr(bare, "predvars"))
  sub <- d[d$x > 5, , drop = FALSE]
  old_sub <- stats::model.matrix(bare, stats::model.frame(bare, sub))
  new_sub <- stats::model.matrix(basis$terms,
                                 stats::model.frame(basis$terms, sub,
                                                    xlev = basis$xlev))
  expect_gt(max(abs(old_sub[, "scale(x)"] - new_sub[, "scale(x)"])), 0.5)
})

test_that("maihda_lme4_fixed_link is batch-invariant and matches predict(re.form=NA)", {
  skip_on_cran()
  set.seed(102)
  n <- 300
  d <- data.frame(x = stats::rnorm(n, mean = 3, sd = 2),
                  g = factor(sample(c("p", "q"), n, TRUE)),
                  st = factor(sample(paste0("s", 1:6), n, TRUE)))
  d$y <- 1 + 0.7 * as.numeric(scale(d$x)) + stats::rnorm(n, sd = 0.6)
  sub <- d$x > 3.2
  expect_gt(sum(sub), 10)

  for (rhs in c("scale(x)", "poly(x, 2)", "splines::ns(x, 3)")) {
    f <- stats::as.formula(paste("y ~", rhs, "+ g + (1 | st)"))
    m <- suppressWarnings(suppressMessages(lme4::lmer(f, data = d)))

    full <- MAIHDA:::maihda_lme4_fixed_link(m, d)
    part <- MAIHDA:::maihda_lme4_fixed_link(m, d[sub, , drop = FALSE])
    native <- as.numeric(stats::predict(m, newdata = d[sub, , drop = FALSE],
                                        re.form = NA))
    # Subsetting after predicting == predicting the subset.
    expect_equal(full[sub], part, tolerance = 1e-10, ignore_attr = TRUE)
    # And both agree with the engine's own safe-prediction route.
    expect_equal(part, native, tolerance = 1e-10, ignore_attr = TRUE)
    # A single-row grid is finite (the old path gave NaN: sd of one value is 0).
    one <- MAIHDA:::maihda_lme4_fixed_link(m, d[1, , drop = FALSE])
    expect_length(one, 1L)
    expect_true(is.finite(one))
    expect_equal(one, full[1], tolerance = 1e-10, ignore_attr = TRUE)
    # A grid holding the transformed variable CONSTANT (the VPC(t) shape) is finite.
    dc <- d
    dc$x <- 2
    expect_true(all(is.finite(MAIHDA:::maihda_lme4_fixed_link(m, dc))))
  }
})

test_that("the fitted-basis fix leaves untransformed fixed parts bit-identical", {
  skip_on_cran()
  set.seed(103)
  n <- 240
  d <- data.frame(x = stats::rnorm(n, 5, 3), z = stats::runif(n, 1, 9),
                  g = factor(sample(c("p", "q", "r"), n, TRUE)),
                  st = factor(sample(paste0("s", 1:6), n, TRUE)))
  d$y <- 1 + 0.4 * d$x + stats::rnorm(n, sd = 0.6)

  # The PRE-FIX code path, reproduced inline: terms() of the bare formula.
  old_path <- function(model, newdata) {
    beta <- lme4::fixef(model)
    tt <- stats::delete.response(
      stats::terms(MAIHDA:::maihda_nobars(stats::formula(model))))
    mfr <- MAIHDA:::maihda_model_frame(model)
    term_vars <- all.vars(tt)
    fac <- names(mfr)[vapply(mfr, is.factor, logical(1))]
    xlev <- lapply(mfr[intersect(fac, term_vars)], levels)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = xlev)
    X <- stats::model.matrix(tt, mf,
                             contrasts.arg = attr(lme4::getME(model, "X"),
                                                  "contrasts"))
    as.numeric(X[, names(beta), drop = FALSE] %*% beta)
  }

  # Terms with NO data-dependent transformation must be untouched by the fix --
  # including poly(raw = TRUE) and I(), which do not depend on the data.
  forms <- list(y ~ x + g + (1 | st),
                y ~ x * g + (1 | st),
                y ~ log(z) + I(x^2) + (1 | st),
                y ~ poly(x, 2, raw = TRUE) + (1 | st))
  for (f in forms) {
    m <- suppressWarnings(suppressMessages(lme4::lmer(f, data = d)))
    expect_equal(MAIHDA:::maihda_lme4_fixed_link(m, d), old_path(m, d),
                 tolerance = 0)
  }
})

test_that("maihda_lme4_fixed_link still carries the offset it is handed", {
  skip_on_cran()
  set.seed(104)
  n <- 200
  d <- data.frame(x = stats::rnorm(n, 5, 3), e = stats::runif(n, 1, 4),
                  st = factor(sample(paste0("s", 1:5), n, TRUE)))
  d$cnt <- stats::rpois(n, exp(0.2 + 0.3 * as.numeric(scale(d$x)) + log(d$e)))

  # Formula offset: the helper reproduces predict(newdata =, re.form = NA), which
  # keeps a formula offset. scale(x) must not perturb that agreement.
  mf_off <- suppressWarnings(lme4::glmer(
    cnt ~ scale(x) + offset(log(e)) + (1 | st), data = d, family = stats::poisson))
  off <- MAIHDA:::maihda_lme4_formula_offset_at(mf_off, d)
  expect_equal(MAIHDA:::maihda_lme4_fixed_link(mf_off, d, offset = off),
               as.numeric(stats::predict(mf_off, newdata = d, re.form = NA)),
               tolerance = 1e-10)

  # External offset=: predict(newdata =) DROPS it, so the helper must differ from
  # that by exactly the offset -- and match the no-newdata prediction, which keeps it.
  me <- suppressWarnings(lme4::glmer(cnt ~ scale(x) + (1 | st), data = d,
                                     family = stats::poisson, offset = log(d$e)))
  ext <- MAIHDA:::maihda_fitted_offset_external(me)
  got <- MAIHDA:::maihda_lme4_fixed_link(me, d, offset = ext)
  expect_equal(got - as.numeric(stats::predict(me, newdata = d, re.form = NA)),
               as.numeric(ext), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(got, as.numeric(stats::predict(me, re.form = NA)),
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("maihda_clmm_linpred uses the fitted transformation basis", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  set.seed(105)
  n <- 400
  d <- data.frame(x = stats::rnorm(n, 5, 3),
                  gender = factor(sample(c("m", "f"), n, TRUE)),
                  race = factor(sample(c("a", "b"), n, TRUE)))
  cell <- interaction(d$gender, d$race, drop = TRUE)
  eta <- 0.9 * as.numeric(scale(d$x)) + stats::rnorm(4, sd = 0.5)[as.integer(cell)]
  d$y <- factor(cut(eta + stats::rlogis(n), breaks = c(-Inf, -1, 0.5, Inf),
                    labels = 1:3), ordered = TRUE)

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ scale(x) + (1 | gender:race), data = d, engine = "ordinal")))
  sub <- fit$data$x > 5.5
  expect_gt(sum(sub), 10)

  full <- maihda_clmm_linpred(fit, newdata = fit$data, include_re = FALSE)
  part <- maihda_clmm_linpred(fit, newdata = fit$data[sub, , drop = FALSE],
                              include_re = FALSE)
  expect_equal(full[sub], part, tolerance = 1e-10, ignore_attr = TRUE)

  one <- maihda_clmm_linpred(fit, newdata = fit$data[1, , drop = FALSE],
                             include_re = FALSE)
  expect_true(is.finite(one))
  expect_equal(one, full[1], tolerance = 1e-10, ignore_attr = TRUE)

  # The rebuilt basis equals the basis clmm itself stored (its terms carry predvars;
  # element 1 there is the response, which the rebuilt terms have dropped).
  basis <- MAIHDA:::maihda_fitted_predict_terms(fit$formula, fit$data)
  expect_equal(as.numeric(attr(basis$terms, "predvars")[[2]]$center),
               as.numeric(attr(fit$model$terms, "predvars")[[3]]$center),
               tolerance = 1e-12)
})

test_that("maihda_wemix_linpred uses the fitted transformation basis", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  set.seed(106)
  n <- 400
  d <- data.frame(x = stats::rnorm(n, 5, 3),
                  gender = factor(sample(c("m", "f"), n, TRUE)),
                  race = factor(sample(c("a", "b"), n, TRUE)))
  cell <- interaction(d$gender, d$race, drop = TRUE)
  d$y <- 2 + 0.8 * as.numeric(scale(d$x)) +
    stats::rnorm(4, sd = 0.6)[as.integer(cell)] + stats::rnorm(n, sd = 0.7)
  d$w <- stats::runif(n, 0.5, 3)

  fit <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ scale(x) + (1 | gender:race), data = d, engine = "wemix",
               sampling_weights = "w")))
  sub <- fit$data$x > 5.5
  expect_gt(sum(sub), 10)

  full <- maihda_wemix_linpred(fit, newdata = fit$data, include_re = FALSE)
  part <- maihda_wemix_linpred(fit, newdata = fit$data[sub, , drop = FALSE],
                               include_re = FALSE)
  expect_equal(full[sub], part, tolerance = 1e-10, ignore_attr = TRUE)

  one <- maihda_wemix_linpred(fit, newdata = fit$data[1, , drop = FALSE],
                              include_re = FALSE)
  expect_true(is.finite(one))
  expect_equal(one, full[1], tolerance = 1e-10, ignore_attr = TRUE)

  # WeMix stores no terms object at all, so the basis MUST come from the stored
  # analytic rows -- which are exactly the rows WeMix::mix() fitted.
  expect_null(fit$model$terms)
  basis <- MAIHDA:::maihda_fitted_predict_terms(fit$formula, fit$data)
  expect_equal(as.numeric(attr(basis$terms, "predvars")[[2]]$center),
               mean(fit$data$x), tolerance = 1e-12)
})

test_that("a longitudinal count VPC(t) is finite when the fit scales time", {
  skip_on_cran()
  set.seed(107)
  nid <- 150
  d <- expand.grid(wave = 0:4, id = seq_len(nid))
  d$gender <- factor(c("m", "f")[1 + (d$id %% 2)])
  d$ethnicity <- factor(c("a", "b", "c")[1 + (d$id %% 3)])
  cellf <- interaction(d$gender, d$ethnicity, drop = TRUE)
  eta <- 0.6 + 0.25 * d$wave +
    stats::rnorm(nlevels(cellf), sd = 0.35)[as.integer(cellf)] +
    stats::rnorm(nid, sd = 0.30)[d$id]
  d$count <- stats::rpois(nrow(d), exp(eta))
  d$id <- factor(d$id)

  fit <- suppressWarnings(suppressMessages(
    maihda(count ~ scale(wave) + (1 | gender:ethnicity), data = d,
           id = "id", time = "wave", family = "poisson",
           decomposition = "longitudinal")))
  vt <- suppressWarnings(summary(fit))$longitudinal$vpc_t

  # Every time point is a finite proportion. The old path fixed all rows to one
  # time, so scale(wave) divided by a zero SD and every estimate was NaN.
  expect_equal(nrow(vt), 5L)
  expect_true(all(is.finite(vt$estimate)))
  expect_true(all(vt$estimate > 0 & vt$estimate < 1))
})

# ---------------------------------------------------------------------------
# Follow-up in the same pass: a longitudinal prediction grid is built FROM the
# stored model frame, which holds a transformed term only under its derived name
# ("scale(z)") and never as the raw z. A covariate appearing only inside a
# transformation therefore had no column to evaluate against and the call died with
# "object 'z' not found". maihda_lme4_fixed_link() now falls back to the values the
# FIT stored for such a term -- per fitted row (fallback = "rows", the VPC(t) grid)
# or at a representative value (fallback = "mean", the trajectory) -- exactly as
# maihda_lme4_formula_offset_at() already did for an offset() term.
# ---------------------------------------------------------------------------

test_that("the stored-column fallback reproduces the re-evaluated term exactly", {
  skip_on_cran()
  set.seed(108)
  n <- 300
  d <- data.frame(z = stats::runif(n, 2, 20), q = stats::rnorm(n, 2, 1),
                  st = factor(sample(paste0("s", 1:6), n, TRUE)))
  d$y <- 1 + 0.6 * as.numeric(scale(d$z)) + 0.3 * d$q + stats::rnorm(n, sd = 0.6)

  for (rhs in c("scale(z)", "poly(z, 2)", "splines::ns(z, 3)", "log(z)",
                "scale(z) + q")) {
    f <- stats::as.formula(paste("y ~", rhs, "+ (1 | st)"))
    m <- suppressWarnings(suppressMessages(lme4::lmer(f, data = d)))
    fr <- MAIHDA:::maihda_model_frame(m)
    # The premise: the model frame does NOT carry the raw z.
    expect_false("z" %in% names(fr))

    got <- MAIHDA:::maihda_lme4_fixed_link(m, fr)          # stored-column fallback
    fr2 <- fr
    fr2$z <- d$z[as.integer(rownames(fr))]                 # raw z restored
    want <- MAIHDA:::maihda_lme4_fixed_link(m, fr2)        # re-evaluated via predvars
    expect_equal(got, want, tolerance = 1e-10)
    expect_equal(got, as.numeric(stats::predict(m, re.form = NA)), tolerance = 1e-10)
  }
})

test_that("the mean fallback holds the term at one representative value", {
  skip_on_cran()
  set.seed(109)
  n <- 300
  d <- data.frame(z = stats::runif(n, 2, 20), q = stats::rnorm(n, 2, 1),
                  st = factor(sample(paste0("s", 1:6), n, TRUE)))
  d$y <- 1 + 0.6 * as.numeric(scale(d$z)) + 0.3 * d$q + stats::rnorm(n, sd = 0.6)
  m <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ scale(z) + q + (1 | st), data = d)))
  fr <- MAIHDA:::maihda_model_frame(m)
  nd <- fr[rep(1L, 20), , drop = FALSE]
  nd$q <- mean(d$q)

  got <- MAIHDA:::maihda_lme4_fixed_link(m, nd, fallback = "mean")
  expect_length(got, 20L)
  expect_true(all(is.finite(got)))
  # Held constant down the grid, and -- scale() being linear -- EXACTLY the value the
  # raw representative covariate would have produced.
  expect_lt(diff(range(got)), 1e-12)
  nd2 <- nd
  nd2$z <- mean(d$z)
  expect_equal(got, MAIHDA:::maihda_lme4_fixed_link(m, nd2), tolerance = 1e-10)

  # The "rows" fallback refuses a grid that does not align with the fitted rows,
  # rather than silently recycling the stored column.
  expect_error(MAIHDA:::maihda_lme4_fixed_link(m, fr[1:10, , drop = FALSE]),
               "prediction grid")

  # THE MISALIGNMENT TRAP: the stored values are placed POSITIONALLY, so a grid with
  # the same NUMBER of rows but different (here permuted) rows must be refused, not
  # silently mis-paired. Row count alone does not settle alignment.
  set.seed(112)
  shuffled <- fr[sample(nrow(fr)), , drop = FALSE]
  expect_false(identical(rownames(shuffled), rownames(fr)))
  expect_error(MAIHDA:::maihda_lme4_fixed_link(m, shuffled),
               "fitted rows in the fitted order")
  # ... and the same grid predicts fine once it carries the raw variable, so the guard
  # is about the FALLBACK, not about permuted grids in general.
  s2 <- shuffled
  s2$z <- d$z[as.integer(rownames(shuffled))]
  expect_true(all(is.finite(MAIHDA:::maihda_lme4_fixed_link(m, s2))))

  # The grid the package itself builds keeps the fitted row names, so the guard never
  # fires on the real VPC(t) path.
  expect_identical(rownames(MAIHDA:::maihda_longitudinal_set_time(fr, "q", 1)),
                   rownames(fr))
})

test_that("a grid carrying every raw variable is untouched by the fallback", {
  skip_on_cran()
  set.seed(110)
  n <- 240
  d <- data.frame(z = stats::runif(n, 2, 20), q = stats::rnorm(n, 2, 1),
                  g = factor(sample(c("p", "r"), n, TRUE)),
                  st = factor(sample(paste0("s", 1:6), n, TRUE)))
  d$y <- 1 + 0.4 * d$z + stats::rnorm(n, sd = 0.6)
  # No stale term can arise here, so the result must equal what the fit's own
  # safe-prediction route gives.
  for (f in list(y ~ scale(z) + g + (1 | st), y ~ poly(z, 2) + (1 | st),
                 y ~ z * g + (1 | st), y ~ log(z) + I(q^2) + (1 | st))) {
    m <- suppressWarnings(suppressMessages(lme4::lmer(f, data = d)))
    expect_equal(MAIHDA:::maihda_lme4_fixed_link(m, d),
                 as.numeric(stats::predict(m, newdata = d, re.form = NA)),
                 tolerance = 1e-10)
  }
})

test_that("a longitudinal fit with a transformed covariate no longer errors", {
  skip_on_cran()
  set.seed(111)
  nid <- 120
  d <- expand.grid(wave = 0:4, id = seq_len(nid))
  d$gender <- factor(c("m", "f")[1 + (d$id %% 2)])
  d$ethnicity <- factor(c("a", "b", "c")[1 + (d$id %% 3)])
  # z appears ONLY inside scale(), so the model frame carries "scale(z)" and no z.
  d$z <- stats::rnorm(nrow(d), 10, 4)
  cf <- interaction(d$gender, d$ethnicity, drop = TRUE)
  lin <- 0.6 + 0.25 * d$wave + 0.1 * as.numeric(scale(d$z)) +
    stats::rnorm(nlevels(cf), sd = 0.35)[as.integer(cf)] +
    stats::rnorm(nid, sd = 0.30)[d$id]
  d$count <- stats::rpois(nrow(d), exp(lin))
  d$score <- lin + stats::rnorm(nrow(d), sd = 0.5)
  d$id <- factor(d$id)

  # Poisson: the VPC(t) grid fixes every row to one time (fallback = "rows").
  fp <- suppressWarnings(suppressMessages(
    maihda(count ~ scale(wave) + scale(z) + (1 | gender:ethnicity), data = d,
           id = "id", time = "wave", family = "poisson",
           decomposition = "longitudinal")))
  vt <- suppressWarnings(summary(fp))$longitudinal$vpc_t
  expect_equal(nrow(vt), 5L)
  expect_true(all(is.finite(vt$estimate)))
  expect_true(all(vt$estimate > 0 & vt$estimate < 1))

  # Gaussian: the residual grid returns early for a non-count family, so it is the
  # TRAJECTORY plot that builds a grid here (fallback = "mean").
  fg <- suppressWarnings(suppressMessages(
    maihda(score ~ scale(wave) + scale(z) + (1 | gender:ethnicity), data = d,
           id = "id", time = "wave", decomposition = "longitudinal")))
  expect_s3_class(suppressWarnings(plot(fg, type = "trajectories")), "ggplot")
  expect_no_error(ggplot2::ggplot_build(
    suppressWarnings(plot(fg, type = "trajectories"))))
})

test_that("staleness is judged before any dummy is written", {
  skip_on_cran()
  set.seed(113)
  n <- 300
  d <- data.frame(z = stats::runif(n, 2, 20), q = stats::rnorm(n, 2, 1),
                  g = factor(sample(c("p", "r"), n, TRUE)),
                  st = factor(sample(paste0("s", 1:6), n, TRUE)))
  d$cnt <- stats::rpois(n, exp(0.2 + 0.3 * as.numeric(scale(d$z)) + log(d$z)))
  d$y <- 1 + 0.6 * as.numeric(scale(d$z)) + stats::rnorm(n, sd = 0.6)
  raw_z <- function(fr) d$z[as.integer(rownames(fr))]

  # A variable SHARED between an offset() and a model-matrix term. The offset's raw
  # variable is dummy-filled first; judging staleness against the PATCHED newdata would
  # make scale(z) look present and silently evaluate it on the dummy, giving wrong
  # numbers with no error.
  m <- suppressWarnings(lme4::glmer(cnt ~ scale(z) + offset(log(z)) + (1 | st),
                                    data = d, family = stats::poisson))
  fr <- MAIHDA:::maihda_model_frame(m)
  expect_false("z" %in% names(fr))
  got <- MAIHDA:::maihda_lme4_fixed_link(
    m, fr, offset = MAIHDA:::maihda_lme4_formula_offset_at(m, fr))
  expect_equal(got, as.numeric(stats::predict(m, re.form = NA)), tolerance = 1e-10)

  # A stale term inside an interaction, and two stale terms with different bases.
  mi <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ scale(z) * g + (1 | st), data = d)))
  fri <- MAIHDA:::maihda_model_frame(mi)
  expect_equal(MAIHDA:::maihda_lme4_fixed_link(mi, fri),
               as.numeric(stats::predict(mi, re.form = NA)), tolerance = 1e-10)

  m2 <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ scale(z) + poly(q, 2) + g + (1 | st), data = d)))
  fr3 <- MAIHDA:::maihda_model_frame(m2)
  expect_false("z" %in% names(fr3))
  expect_false("q" %in% names(fr3))
  expect_equal(MAIHDA:::maihda_lme4_fixed_link(m2, fr3),
               as.numeric(stats::predict(m2, re.form = NA)), tolerance = 1e-10)

  # A stale term alongside an EXTERNAL offset=, which has no expression to re-evaluate.
  me <- suppressWarnings(lme4::glmer(cnt ~ scale(z) + (1 | st), data = d,
                                     family = stats::poisson, offset = log(d$z)))
  fre <- MAIHDA:::maihda_model_frame(me)
  expect_equal(MAIHDA:::maihda_lme4_fixed_link(
                 me, fre, offset = MAIHDA:::maihda_fitted_offset_external(me)),
               as.numeric(stats::predict(me, re.form = NA)), tolerance = 1e-10)
})
