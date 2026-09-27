# Audit 2026-09-26 (a third pass of that date, run after the F06 observed-vs-shrunken
# pass on R/plot_maihda.R; this is its sibling in the discriminatory-accuracy path).
#
# FINDING: maihda_da_observed_response() read the observed response for every
# non-lme4 engine as
#
#     resp <- all.vars(model$formula)[1]
#     y    <- model$data[[resp]]
#
# -- the first variable NAMED in the formula, taken raw out of the frame, rather than
# the response the engine actually FITTED. When the outcome is an expression the two
# are different columns. A brms model frame is model.frame.brmsfit() = the brmsfit's
# own `data` slot, which keeps every raw input column alongside the evaluated
# response, so the wrong column was silently there to be read.
#
# MEASURED on ONE real 240-row Stan fit of
# round(raw) | trials(ntr) ~ x + (1 | stratum), read with TWO spellings of the decoy
# column. The two spellings give the two failure modes, and only the second has an
# AUC at all -- do not quote these numbers against one fixture:
#
#   * raw = successes + 0.4 -- exceeds the trial count on each all-success row (3 of
#     240), and the call HARD STOPS with "Internal error: negative case/control mass
#     in the weighted AUC";
#   * raw = max(successes - 0.4, 0) -- stays inside [0, trials], so nothing complains:
#     AUC 0.6587315786 against a true 0.6562245743, with fractional "counts" reported
#     as the totals (909.6 cases / 1456.4 controls for data holding 1002 / 1364).
#
# The fixture below is the FIRST spelling, shrunk to 4 rows.
#
# The value reaches maihda_discriminatory_accuracy() (the AUC) and the
# aggregated-binomial success counts in maihda_da_brms_aggregated_counts().
#
# SCOPE (re-verified, not assumed): lme4 never reached the branch (it returns
# getME(, "y")); wemix refuses a non-symbol response outright ("binary (Bernoulli)
# 0/1 outcome only"), as does the ordinal engine ("no cbind()/addition terms"); brms
# itself refuses a logical expression response and a bare numeric 0/1 expression. The
# live route is the aggregated binomial with a transformed success count -- including
# its DESIGN-WEIGHTED variant, whose response then carries two addition terms
# (`y | trials(n) + weights(.maihda_sw)`) and was equally affected. A bare-symbol
# response is a no-op on every engine, two addition terms included.
#
# Fixed by reading the response through the new maihda_response_from_frame(), which
# asks the frame's own terms (stats::model.response()) before falling back to
# evaluating the response expression.

# A FAITHFUL brms model frame. Measured on the real Stan fit above, the frame's terms
# deparse to `round(raw) ~ raw + ntr + x + stratum` with response = 1 -- brms has
# already stripped the `| trials(ntr)` addition term, and the untransformed `raw`
# column survives beside the fitted `round(raw)`. stats::model.frame() reproduces that
# structure exactly, so the bug path runs with no Stan compile.
audit_0926c_frame <- function(d) {
  stats::model.frame(round(raw) ~ raw + ntr + x + stratum, data = d,
                     na.action = stats::na.pass)
}

audit_0926c_data <- function() {
  d <- data.frame(
    x       = c(-1, -0.5, 0.5, 1),
    stratum = factor(c("s1", "s1", "s2", "s2")),
    ntr     = c(8, 13, 11, 7)
  )
  # The last row is ALL successes, so the decoy column exceeds its trial count there
  # -- the row that turned the silently-wrong AUC into a hard stop.
  d$succ <- c(2, 6, 4, 7)
  d$raw  <- d$succ + 0.4
  d
}

audit_0926c_model <- function(data, formula, engine = "brms") {
  structure(
    list(
      model = structure(list(data = data), class = "brmsfit"),
      engine = engine,
      formula = formula,
      data = data,
      original_data = data,
      family = list(family = "binomial", link = "logit"),
      strata_vars = "stratum",
      strata_sep = " x ",
      strata_info = NULL,
      strata_autobin_info = NULL,
      context_vars = NULL,
      sampling_weights = NULL,
      longitudinal_info = NULL,
      response_recoding = NULL
    ),
    class = "maihda_model"
  )
}

test_that("the DA response is the FITTED outcome, not the first formula variable", {
  d <- audit_0926c_data()
  f <- round(raw) | trials(ntr) ~ x + (1 | stratum)
  m <- audit_0926c_model(audit_0926c_frame(d), f)

  # Teeth: the old spelling's target really is a different column, and it really is
  # present in the frame -- so the wrong read succeeded rather than erroring.
  expect_identical(all.vars(f)[1], "raw")
  expect_true("raw" %in% names(m$data))
  expect_equal(as.numeric(m$data[["raw"]]), d$raw)

  got <- as.numeric(MAIHDA:::maihda_da_observed_response(m))
  expect_equal(got, d$succ)                               # the response that was fitted
  expect_false(isTRUE(all.equal(got, d$raw)))             # was: the raw decoy column
})

test_that("the aggregated-binomial success counts stay inside [0, trials]", {
  d <- audit_0926c_data()
  m <- audit_0926c_model(audit_0926c_frame(d),
                         round(raw) | trials(ntr) ~ x + (1 | stratum))

  counts <- MAIHDA:::maihda_da_brms_aggregated_counts(m)
  expect_equal(counts$successes, d$succ)
  expect_equal(counts$trials, d$ntr)
  # The contract the weighted AUC depends on. The decoy column broke it on the
  # all-success row (7.4 successes out of 7 trials).
  expect_true(all(counts$successes >= 0 & counts$successes <= counts$trials))
  expect_true(any(d$raw > d$ntr))
})

test_that("the weighted AUC no longer hard-stops on an out-of-range response", {
  d <- audit_0926c_data()
  prob <- c(0.2, 0.5, 0.4, 0.8)

  # What the old read handed maihda_auc_weighted(): a negative control mass.
  expect_error(MAIHDA:::maihda_auc_weighted(prob, d$raw, d$ntr),
               "negative case/control mass")

  auc <- MAIHDA:::maihda_auc_weighted(prob, d$succ, d$ntr)
  expect_true(is.finite(auc) && auc >= 0 && auc <= 1)
})

test_that("a terms-free frame with a bare-symbol response is read exactly as before", {
  # The wemix / ordinal shape: fit_maihda() stores the pre-built analytic data, a
  # plain data frame with no terms attribute, so model.response() is NULL. Both
  # engines refuse a non-symbol response at fit time, so this is the whole of their
  # surface -- and it must return the same column the old all.vars() read did.
  d <- data.frame(y = c(0, 1, 1, 0), x = c(-1, -0.5, 0.5, 1),
                  stratum = factor(c("s1", "s1", "s2", "s2")))
  m <- audit_0926c_model(d, y ~ x + (1 | stratum), engine = "wemix")

  expect_null(tryCatch(stats::model.response(m$data), error = function(e) NULL))
  expect_equal(as.numeric(MAIHDA:::maihda_da_observed_response(m)), d$y)

  # The 0/1 coercions are preserved.
  df <- d; df$y <- factor(c("no", "yes", "yes", "no"), levels = c("no", "yes"))
  mf <- audit_0926c_model(df, y ~ x + (1 | stratum), engine = "wemix")
  expect_equal(MAIHDA:::maihda_da_observed_response(mf), c(0L, 1L, 1L, 0L))

  dl <- d; dl$y <- c(FALSE, TRUE, TRUE, FALSE)
  ml <- audit_0926c_model(dl, y ~ x + (1 | stratum), engine = "wemix")
  expect_equal(MAIHDA:::maihda_da_observed_response(ml), c(0L, 1L, 1L, 0L))
})

test_that("maihda_response_from_frame resolves each shape, and NULLs the unresolvable", {
  d <- audit_0926c_data()

  # 1. A frame with terms: the engine's own evaluated response.
  expect_equal(as.numeric(MAIHDA:::maihda_response_from_frame(
    audit_0926c_frame(d), round(raw) | trials(ntr) ~ x + (1 | stratum))), d$succ)

  # 2. A plain frame, bare symbol: the column.
  expect_equal(MAIHDA:::maihda_response_from_frame(d, succ ~ x + (1 | stratum)), d$succ)

  # 3. A plain frame, expression response: rebuilt.
  expect_equal(as.numeric(MAIHDA:::maihda_response_from_frame(
    d, round(raw) ~ x + (1 | stratum))), d$succ)

  # 4. na.pass keeps the rebuild aligned to the frame instead of silently shortening
  #    it -- `data` is already the analytic sample.
  dn <- d; dn$raw[2] <- NA_real_
  out <- MAIHDA:::maihda_response_from_frame(dn, round(raw) ~ x + (1 | stratum))
  expect_length(out, nrow(dn))
  expect_true(is.na(out[2]))

  # 5. Unresolvable: NULL, so a caller keeps its own fallback rather than a wrong
  #    answer. The DA reader turns that into a named error.
  expect_null(MAIHDA:::maihda_response_from_frame(d, absent ~ x + (1 | stratum)))
  m <- audit_0926c_model(d, absent ~ x + (1 | stratum), engine = "wemix")
  expect_error(MAIHDA:::maihda_da_observed_response(m), "Could not recover the response")
})

test_that("the rebuild never answers from the formula's enclosure", {
  # Caught in this pass's own self-check, in the fix rather than in the finding:
  # model.frame() resolves a variable it cannot find in `data` against the enclosure,
  # and environment(formula) chains out to globalenv. An intermediate version of
  # maihda_response_from_frame() therefore answered a rebuild from a same-named object
  # that had nothing to do with the fit -- the failure mode the 2026-09-26 ordinal pass
  # names as the worst of its kind (a published number from the wrong data).
  dcoy <- c(100, 200, 300, 400)          # the decoy, visible to the formula below
  d <- data.frame(x = c(-1, -0.5, 0.5, 1))   # and deliberately NOT a column

  f <- round(dcoy) ~ x                   # environment(): this test's frame
  expect_true("dcoy" %in% all.vars(f))
  expect_false("dcoy" %in% names(d))
  expect_true(exists("dcoy", envir = environment(f), inherits = TRUE))
  expect_null(MAIHDA:::maihda_response_from_frame(d, f))

  # The bare-symbol branch was always guarded by its names(data) check.
  expect_null(MAIHDA:::maihda_response_from_frame(d, dcoy ~ x))

  # The legitimate rebuild still works: the VARIABLE comes from the frame, only the
  # FUNCTION from the environment.
  d2 <- data.frame(x = c(-1, -0.5, 0.5, 1), dcoy = c(2.4, 6.4, 4.4, 7.4))
  expect_equal(as.numeric(MAIHDA:::maihda_response_from_frame(d2, round(dcoy) ~ x)),
               c(2, 6, 4, 7))

  # A partially-resolvable expression is refused rather than half-read.
  d3 <- data.frame(x = c(-1, -0.5, 0.5, 1), a = c(1, 2, 3, 4))
  bvar <- c(9, 9, 9, 9)
  expect_null(MAIHDA:::maihda_response_from_frame(d3, I(a + bvar) ~ x))
})

test_that("lme4 still reads its response from getME(, 'y')", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  set.seed(926)
  n <- 240
  d <- data.frame(x = stats::rnorm(n),
                  stratum = factor(sample(sprintf("s%d", 1:8), n, TRUE)))
  d$y  <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * d$x))
  d$yr <- d$y + 0.4                       # the same decoy, on the lme4 path
  m <- suppressWarnings(suppressMessages(
    fit_maihda(round(yr) ~ x + (1 | stratum), data = d, family = "binomial",
               engine = "lme4")))

  # lme4's frame does not even carry the decoy column, so the branch was never
  # reachable here -- but pin the engine-native read anyway.
  expect_false("yr" %in% names(m$data))
  expect_equal(as.numeric(MAIHDA:::maihda_da_observed_response(m)),
               as.numeric(lme4::getME(m$model, "y")))
})

test_that("a real brms fit with a transformed response gets the fitted successes", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  set.seed(1)
  K <- 24
  d <- data.frame(stratum = factor(rep(sprintf("s%02d", seq_len(K)), each = 10)),
                  x = stats::rnorm(240))
  d$ntr  <- sample(5:15, 240, TRUE)
  d$succ <- stats::rbinom(240, d$ntr, stats::plogis(-0.3 + 0.6 * d$x))
  d$raw  <- d$succ + 0.4

  m <- suppressWarnings(suppressMessages(
    fit_maihda(round(raw) | trials(ntr) ~ x + (1 | stratum), data = d,
               engine = "brms", family = "binomial",
               chains = 1, iter = 400, warmup = 200, refresh = 0, seed = 1)))

  expect_identical(all.vars(m$formula)[1], "raw")
  got <- as.numeric(MAIHDA:::maihda_da_observed_response(m))
  expect_equal(got, as.numeric(d$succ))
  expect_false(isTRUE(all.equal(got, as.numeric(d$raw))))

  # The AUC that previously aborted, and its unweighted-count totals.
  prob <- predict_maihda(m, type = "individual", scale = "response")
  da <- maihda_discriminatory_accuracy(m)
  expect_equal(da$auc, MAIHDA:::maihda_auc_weighted(prob, d$succ, d$ntr))
  expect_equal(da$n_case, sum(d$succ))
  expect_equal(da$n_control, sum(d$ntr) - sum(d$succ))
})

test_that("a design-weighted brms fit with a transformed response is also correct", {
  # The composite the finding did not name and the suite never runs: brms aggregated
  # binomial + sampling_weights puts TWO addition terms on the response
  # (`round(rawsucc) | trials(total) + weights(.maihda_sw)`). Measured: the frame's
  # terms are `round(rawsucc) ~ rawsucc + total + .maihda_sw + stratum` with
  # response = 1, so model.response() still resolves, and all.vars()[1] still did not.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  set.seed(909)
  K <- 12
  agg <- data.frame(stratum = factor(sprintf("s%02d", seq_len(K))),
                    total = sample(20:40, K, replace = TRUE))
  agg$success <- stats::rbinom(K, agg$total,
                               stats::plogis(seq(-1.6, 1.6, length.out = K)))
  agg$w <- rev(seq(0.2, 3, length.out = K))
  agg$rawsucc <- agg$success + 0.4

  m <- suppressWarnings(suppressMessages(
    fit_maihda(round(rawsucc) | trials(total) ~ (1 | stratum), data = agg,
               engine = "brms", family = "binomial", sampling_weights = "w",
               chains = 1, iter = 300, refresh = 0)))

  expect_identical(all.vars(m$formula)[1], "rawsucc")
  expect_equal(as.numeric(MAIHDA:::maihda_da_observed_response(m)),
               as.numeric(agg$success))
  expect_equal(MAIHDA:::maihda_da_brms_aggregated_counts(m)$successes,
               as.numeric(agg$success))

  da <- suppressWarnings(maihda_discriminatory_accuracy(m))
  expect_true(is.finite(da$auc))
  # Reported totals stay unweighted observation counts, as the 2026-08 weighted-AUC
  # test pins for the bare-symbol spelling.
  expect_equal(da$n_case, sum(agg$success))
  expect_equal(da$n_control, sum(agg$total) - sum(agg$success))
})

test_that("the obs_vs_shrunken reader is guarded against the enclosure too", {
  # The same escape, found in R/plot_maihda.R's maihda_observed_response_from_model_frame()
  # after this pass fixed it in its own helper. That function's BARE-SYMBOL branch was
  # already guarded by a names(data) check; its model.frame() rebuild was not, and its
  # one-value-per-row check catches a namesake of the WRONG length only. Measured before
  # the guard: a frame with no `dcoy` column plus a four-element `dcoy` in scope returned
  # 100 200 300 400 as the observed response -- a plotted axis from the wrong data.
  #
  # Only reachable at the unit level: the expression branch needs a terms-free frame
  # (wemix / ordinal) whose response variable is not one of its columns, and
  # fit_maihda() refuses a response that lives outside `data` ("undefined columns
  # selected"). Hardening, not a reproduced field defect.
  skip_if_not(exists("maihda_observed_response_from_model_frame",
                     envir = asNamespace("MAIHDA"), inherits = FALSE),
              "obs_vs_shrunken reader not present in this tree")
  h <- get("maihda_observed_response_from_model_frame", envir = asNamespace("MAIHDA"))

  d <- data.frame(x = c(-1, -0.5, 0.5, 1),
                  stratum = factor(c("s1", "s1", "s2", "s2")))
  dcoy <- c(100, 200, 300, 400)          # same length as nrow(d): the dangerous case
  expect_false("dcoy" %in% names(d))
  expect_error(h(d, round(dcoy) ~ x), "not a column of the model frame")

  # A bare symbol was already refused, and stays refused.
  expect_error(h(d, dcoy ~ x), "Outcome variable not found in data")

  # The legitimate rebuild is untouched: the VARIABLE comes from the frame, only the
  # FUNCTION from the environment. This is F06's own reproducer shape.
  d2 <- data.frame(x = c(-1, -0.5, 0.5, 1), dcoy = c(2.4, 6.4, 4.4, 7.4),
                   stratum = factor(c("s1", "s1", "s2", "s2")))
  expect_equal(as.numeric(h(d2, round(dcoy) ~ x)), c(2, 6, 4, 7))
  d3 <- data.frame(x = c(-1, -0.5, 0.5, 1), yv = c(1, 10, 100, 1000),
                   stratum = factor(c("s1", "s1", "s2", "s2")))
  expect_equal(as.numeric(h(d3, log(yv) ~ x)), log(c(1, 10, 100, 1000)))
})

test_that("both readers are one implementation with two contracts", {
  # Two concurrent passes on 2026-09-26 fixed the SAME defect in two copies of this
  # logic, and the copies drifted -- one had the brmsformula unwrapping and the
  # one-value-per-row check, the other the enclosure guard. They are now a single
  # maihda_response_from_frame(); the plot reader is a strict = TRUE wrapper. This
  # block is what makes a future re-divergence fail rather than pass quietly.
  skip_if_not(exists("maihda_observed_response_from_model_frame",
                     envir = asNamespace("MAIHDA"), inherits = FALSE),
              "obs_vs_shrunken reader not present in this tree")
  strict <- get("maihda_observed_response_from_model_frame",
                envir = asNamespace("MAIHDA"))
  loose <- MAIHDA:::maihda_response_from_frame

  d <- data.frame(x = c(-1, -0.5, 0.5, 1),
                  y = c(1, 10, 100, 1000),
                  n = c(8, 13, 11, 7),
                  f = c(7, 7, 7, 7),
                  stratum = factor(c("s1", "s1", "s2", "s2")))

  # SAME VALUE on every shape either one resolves -- that is what "one
  # implementation" has to mean in practice.
  shapes <- list(y ~ x, log(y) ~ x, I(y > 3) ~ x, y | trials(n) ~ x,
                 stratum ~ x, cbind(y, f) ~ x)
  for (fo in shapes) {
    expect_equal(loose(d, fo), strict(d, fo),
                 info = paste(deparse(fo), collapse = " "))
  }
  # ... including the shared extras each copy used to have alone.
  expect_true(is.matrix(strict(d, cbind(y, f) ~ x)))          # was F06's only
  expect_identical(strict(d, stratum ~ x), d$stratum)         # class preserved

  # DIFFERENT CONTRACT, same decision: where one returns NULL the other throws.
  unresolvable <- list(
    list(f = absent ~ x,            pattern = "Outcome variable not found"),
    list(f = mean(y) ~ x,           pattern = "one value per row"))
  for (u in unresolvable) {
    expect_null(loose(d, u$f), info = paste(deparse(u$f), collapse = " "))
    expect_error(strict(d, u$f), u$pattern)
  }

  # The enclosure guard is in the shared core, so BOTH readers have it now -- this
  # is the half F06's copy was missing and the DA copy had.
  dcoy <- c(100, 200, 300, 400)
  dg <- data.frame(x = c(-1, -0.5, 0.5, 1))
  expect_null(loose(dg, round(dcoy) ~ x))
  expect_error(strict(dg, round(dcoy) ~ x), "not a column of the model frame")
})

test_that("the DA reader inherited the brmsformula unwrapping from the fold", {
  # The half the DA copy was missing. Before the fold, a brmsformula reached
  # maihda_describe_response_expr() unwrapped and the read failed; it now resolves
  # on both readers, from one place, so the outcome and the trials() denominator
  # cannot come from two different objects.
  skip_if_not_installed("brms")
  d <- data.frame(x = c(-1, -0.5, 0.5, 1), y = c(1, 10, 100, 1000),
                  n = c(8, 13, 11, 7), stratum = factor(c("s1", "s1", "s2", "s2")))
  bf <- brms::bf(log(y) | trials(n) ~ x + (1 | stratum))

  expect_equal(unname(as.numeric(MAIHDA:::maihda_response_from_frame(d, bf))),
               log(d$y))
  expect_equal(MAIHDA:::maihda_trials_from_formula(bf, d), d$n)
})
