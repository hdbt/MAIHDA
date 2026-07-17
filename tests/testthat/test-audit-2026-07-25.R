# Regression tests for the 2026-07-25 audit findings.
#
#   1 [P1] calculate_pcv(bootstrap=TRUE) simulated responses from model2 (its row
#          order) and fed them UNCHANGED into model1's refit. When the two fits share
#          the same observations in a different row order (which validate_pcv_models
#          accepts) each draw was paired with the wrong row for model1, corrupting
#          var1 and the whole interval. Draws are now permuted into model1's row
#          order (by persistent row id) before refitting model1.
#   2 [P1] a caller-supplied 'stratum' that CONTRADICTS the dimension columns was
#          accepted, pairing one intersection's fixed effects with another's random
#          effect. It is now cross-checked against the dimensions.
#   3 [P1] maihda_fit_diagnostics() read only conv$lme4$messages, so an optimizer
#          that stopped early (conv$opt != 0) but landed at a small-gradient point
#          (empty lme4 messages) was reported converged. It now also consults
#          conv$opt / optinfo$message.
#   4 [P2] brms individual predictions zeroed only the unseen STRATUM effect; an
#          unseen context/longitudinal level stayed in re_formula and was SAMPLED
#          (brms default sample_new_levels="uncertainty"), diverging from lme4's
#          zero. Every unseen grouping level is now dropped per row.
#   5 [P2] maihda_longitudinal_pcv() stored only ml_refit; a boundary skip or failed
#          refitML() left a mixed REML/ML comparison indistinguishable from a clean
#          fitted one. It now records estimation_used ("fitted"/"ML"/"mixed").

# ---- Finding 1: bootstrap response alignment ---------------------------------

test_that("maihda_reorder_response permutes a vector and a matrix by row", {
  f <- MAIHDA:::maihda_reorder_response
  expect_identical(f(c(10, 20, 30), c(3L, 1L, 2L)), c(30, 10, 20))
  m <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)          # 3 rows x 2 cols
  out <- f(m, c(3L, 1L, 2L))
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(3L, 2L))
  expect_identical(out, matrix(c(3, 1, 2, 6, 4, 5), ncol = 2))
})

maihda_reorder_fixture <- function() {
  set.seed(42); N <- 200
  d1 <- factor(sample(c("m", "f"), N, TRUE))
  d2 <- factor(sample(c("lo", "mid", "hi"), N, TRUE))
  sid <- droplevels(interaction(d1, d2))
  u  <- rnorm(nlevels(sid), 0, 1.2)[as.integer(sid)]
  y  <- 1 + 0.5 * (d1 == "m") + 0.3 * (d2 == "hi") + u + rnorm(N)
  dat <- data.frame(y = y, d1 = d1, d2 = d2)
  rownames(dat) <- paste0("obs", seq_len(N)); dat
}

test_that("PCV bootstrap aligns a reordered model2 draw to model1's rows", {
  skip_on_cran()
  dat <- maihda_reorder_fixture()
  null_fit <- suppressMessages(suppressWarnings(fit_maihda(y ~ 1 + (1 | d1:d2), data = dat)))
  adj_fit  <- suppressMessages(suppressWarnings(fit_maihda(y ~ d1 + d2 + (1 | d1:d2), data = dat)))
  set.seed(7); shuf <- dat[sample(nrow(dat)), ]
  adj_shuf <- suppressMessages(suppressWarnings(fit_maihda(y ~ d1 + d2 + (1 | d1:d2), data = shuf)))

  # Same row order -> identity alignment (NULL); reordered -> a real permutation.
  expect_null(MAIHDA:::maihda_pcv_boot_align(null_fit$model, adj_fit$model))
  perm <- MAIHDA:::maihda_pcv_boot_align(null_fit$model, adj_shuf$model)
  expect_false(is.null(perm))
  expect_identical(sort(perm), seq_len(nrow(dat)))

  set.seed(101)
  base  <- suppressWarnings(calculate_pcv(null_fit, adj_fit,  bootstrap = TRUE, n_boot = 40))
  set.seed(101)
  reord <- suppressWarnings(calculate_pcv(null_fit, adj_shuf, bootstrap = TRUE, n_boot = 40))
  # Point PCV is order-invariant; the reordered CI is now sane (pre-fix it was a
  # huge-magnitude interval with many spurious zero-null-variance boundary draws).
  expect_equal(reord$pcv, base$pcv, tolerance = 1e-3)
  expect_true(is.finite(reord$ci_lower) && reord$ci_lower > -5)
  expect_true(reord$n_boot_boundary <= base$n_boot_boundary + 2L)
})

# ---- Finding 2: stratum vs dimensions consistency ----------------------------

test_that("prediction rejects a supplied stratum that contradicts the dimensions", {
  skip_on_cran()
  dat <- maihda_reorder_fixture()
  adj <- suppressMessages(suppressWarnings(fit_maihda(y ~ d1 + d2 + (1 | d1:d2), data = dat)))

  nd  <- adj$data[1, , drop = FALSE]
  bad <- nd
  bad$stratum <- setdiff(as.character(unique(adj$data$stratum)),
                         as.character(nd$stratum))[1]
  expect_error(predict_maihda(adj, newdata = bad, type = "link"),
               "must match the dimension columns")
  # A consistent supplied stratum, and predicting on the fitted data, still work.
  expect_true(is.finite(as.numeric(predict_maihda(adj, newdata = nd, type = "link"))))
  expect_length(predict_maihda(adj, newdata = adj$data, type = "link"), nrow(adj$data))
})

test_that("stratum/dimension check is a no-op when the dimensions are absent", {
  # A helper-level check that needs no fit: without the dimension columns there is
  # nothing to cross-check, so a supplied stratum is trusted (returns invisibly).
  obj <- list(
    strata_vars = c("d1", "d2"),
    strata_info = data.frame(stratum = c("1", "2"), label = c("m x lo", "f x hi"),
                             stringsAsFactors = FALSE),
    strata_sep = " x ", strata_autobin_info = NULL)
  nd_no_dims <- data.frame(stratum = "1")
  expect_null(MAIHDA:::maihda_check_stratum_matches_dims(obj, nd_no_dims))
  # With dims present and matching -> no error; contradicting -> error.
  nd_ok  <- data.frame(stratum = "1", d1 = "m", d2 = "lo")
  nd_bad <- data.frame(stratum = "2", d1 = "m", d2 = "lo")   # dims say stratum 1
  expect_null(MAIHDA:::maihda_check_stratum_matches_dims(obj, nd_ok))
  expect_error(MAIHDA:::maihda_check_stratum_matches_dims(obj, nd_bad),
               "must match the dimension columns")
})

# ---- Finding 3: optimizer convergence code -----------------------------------

test_that("maihda_fit_diagnostics consults the optimizer return code", {
  skip_on_cran()
  dat <- maihda_reorder_fixture()
  sid <- droplevels(interaction(dat$d1, dat$d2)); d3 <- cbind(dat, sid = sid)

  # A clean fit is converged.
  good <- lme4::lmer(y ~ d1 + d2 + (1 | sid), data = d3)
  expect_true(MAIHDA:::maihda_fit_diagnostics(good)$converged)

  # An optimizer that stops early (conv$opt != 0) at a small-gradient point (empty
  # lme4 messages) must NOT be reported converged, and the optimizer message shows.
  bad <- good
  bad@optinfo$conv$opt <- 1L
  bad@optinfo$message  <- "bobyqa -- maximum number of function evaluations exceeded"
  bad@optinfo$conv$lme4 <- list(code = 0L, messages = character(0))
  diag_bad <- MAIHDA:::maihda_fit_diagnostics(bad)
  expect_false(diag_bad$converged)
  expect_true(any(grepl("function evaluations", diag_bad$messages)))

  # A natural early stop (bobyqa maxfun) is likewise flagged.
  m8 <- suppressWarnings(lme4::lmer(y ~ d1 + d2 + (1 | sid), data = d3,
          control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 8))))
  expect_false(MAIHDA:::maihda_fit_diagnostics(m8)$converged)
})

# ---- Finding 4: brms unseen non-stratum levels (Stan-free) --------------------

test_that("brms individual prediction zeroes every unseen grouping level", {
  grp <- function(re) if (length(re) == 1 && is.na(re)) "NA" else
    paste(sort(unlist(lapply(reformulas::findbars(re), function(b) all.vars(b[[3]])))),
          collapse = "+")
  mk <- function(formula, data = NULL)
    structure(list(model = structure(list(), class = "brmsfit"),
                   engine = "brms", formula = formula, data = data),
              class = "maihda_model")
  train <- data.frame(stratum = c("1", "2", "1", "2"), school = c("A", "B", "A", "B"),
                      x = 0, y = 0, stringsAsFactors = FALSE)
  obj <- mk(y ~ x + (1 | stratum) + (1 | school), train)

  run <- function(object, nd, allow) testthat::with_mocked_bindings(
    MAIHDA:::maihda_brms_individual_prediction(object, nd, "link", allow,
                                               list(allow_new_levels = allow)),
    maihda_brms_predict_rows = function(object, nd, scale, dots)
      rep(if ("re_formula" %in% names(dots)) grp(dots$re_formula) else "FULL", nrow(nd)),
    .package = "MAIHDA")

  nd <- data.frame(x = 0,
    stratum = c("1", "2", "ZZ", "ZZ"),   # seen, seen, unseen, unseen
    school  = c("A", "ZZ", "A", "ZZ"),   # seen, unseen, seen, unseen
    stringsAsFactors = FALSE)
  expect_identical(run(obj, nd, TRUE), c("FULL", "stratum", "school", "NA"))

  # Stratum-only model reduces to the previous behaviour exactly.
  obj2 <- mk(y ~ x + (1 | stratum), train)
  nd2 <- data.frame(x = 0, stratum = c("1", "ZZ"), stringsAsFactors = FALSE)
  expect_identical(run(obj2, nd2, TRUE), c("FULL", "NA"))

  # allow_new_levels = FALSE -> plain full prediction, no per-row re_formula.
  expect_true(all(run(obj, nd, FALSE) == "FULL"))
})

# ---- Finding 5: longitudinal PCV estimation basis ----------------------------

test_that("longitudinal PCV print surfaces a mixed REML/ML basis (no fit)", {
  # Build with no duplicate keys (x$field returns the FIRST match, so an appended
  # ml_refit must not shadow a base one).
  base_obj <- list(
    pcv_intercept = 0.4, pcv_slope = 0.3,
    var_baseline_null = 0.5, var_baseline_adjusted = 0.3,
    var_slope_null = 0.2, var_slope_adjusted = 0.14,
    ref_time = 0, time = "wave",
    Sigma_stratum_null = matrix(c(0.5, 0, 0, 0.2), 2),
    Sigma_stratum_adjusted = matrix(c(0.3, 0, 0, 0.14), 2),
    null_at_boundary = FALSE)
  mk_lp <- function(...) structure(c(base_obj, list(...)), class = "maihda_long_pcv")
  mixed  <- mk_lp(estimation = "ML",     estimation_used = "mixed",  ml_refit = FALSE)
  fitted <- mk_lp(estimation = "fitted", estimation_used = "fitted", ml_refit = FALSE)
  ml     <- mk_lp(estimation = "ML",     estimation_used = "ML",     ml_refit = TRUE)
  expect_output(print(mixed), "kept its REML fit")
  expect_false(any(grepl("kept its REML fit", capture.output(print(fitted)))))
  expect_false(any(grepl("refitted with maximum likelihood", capture.output(print(fitted)))))
  expect_output(print(ml), "refitted with maximum likelihood")
})

test_that("maihda_longitudinal_pcv records fitted / ML / mixed bases", {
  skip_on_cran()
  a <- suppressMessages(suppressWarnings(maihda(
    wellbeing ~ wave + (1 | gender:ethnicity:education),
    data = maihda_long_data, id = "id", time = "wave", decomposition = "longitudinal")))
  nullm <- a$model; adjm <- a$model_adjusted

  fit_pcv <- MAIHDA:::maihda_longitudinal_pcv(nullm, adjm, estimation = "fitted")
  expect_identical(fit_pcv$estimation_used, "fitted")
  expect_false(isTRUE(fit_pcv$ml_refit))

  ml_pcv <- MAIHDA:::maihda_longitudinal_pcv(nullm, adjm, estimation = "ML")
  expect_identical(ml_pcv$estimation_used, "ML")
  expect_true(isTRUE(ml_pcv$ml_refit))

  # An ML request whose refit does not apply (mocked here as a no-op, standing in
  # for a boundary skip / failed refitML) leaves a mixed comparison, recorded as
  # such and NOT as ml_refit.
  mixed <- testthat::with_mocked_bindings(
    MAIHDA:::maihda_longitudinal_pcv(nullm, adjm, estimation = "ML"),
    maihda_longitudinal_refit_ml = function(model) model,
    .package = "MAIHDA")
  expect_identical(mixed$estimation_used, "mixed")
  expect_false(isTRUE(mixed$ml_refit))
})
