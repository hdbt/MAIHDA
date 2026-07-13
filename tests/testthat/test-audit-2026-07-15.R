# Regression tests for the 2026-07-15 audit findings.
#
#   P1 #1  brms silently ignored lme4-style subset/weights/offset dots
#   P1 #2  a design-weighted brms fit's derived null/adjusted refit collided
#          with the reserved .maihda_sw column and aborted
#   P2 #3  zero/negative Gaussian precision weights produced a degenerate lmer
#          fit (logLik -Inf) yet a finite VPC
#   P2 #4  WeMix/ordinal individual predictions silently mapped a missing
#          stratum to a zero random effect (population-average prediction)
#
# The brms cases are exercised WITHOUT a Stan compile: #1 is rejected during
# argument handling (before any brms call), and #2 drives the package's pure
# weight-preparation helpers directly.

# ---- P1 #1: brms rejects lme4-style data-masked fitting arguments ------------

test_that("engine = 'brms' rejects subset/weights/offset dots (audit P1 #1)", {
  d <- data.frame(y = rnorm(20), x = rnorm(20),
                  stratum = factor(rep(letters[1:4], 5)))
  # The guard fires while arguments are handled, before brms/Stan is touched --
  # so this needs neither brms installed nor a sampler.
  expect_error(
    fit_maihda(y ~ x + (1 | stratum), data = d, engine = "brms", subset = 1:7),
    "not supported by engine = \"brms\"", fixed = TRUE)
  expect_error(
    fit_maihda(y ~ x + (1 | stratum), data = d, engine = "brms",
               weights = rep(1, 20)),
    "not supported by engine = \"brms\"", fixed = TRUE)
  expect_error(
    fit_maihda(y ~ x + (1 | stratum), data = d, engine = "brms",
               offset = rep(0, 20)),
    "not supported by engine = \"brms\"", fixed = TRUE)
})

# ---- P1 #2: a derived weighted brms refit re-prepares cleanly ---------------

test_that("maihda_strip_brms_weights_term removes only the internal weights term", {
  # Sole addition term -> the `|` is dropped, leaving the bare response.
  f1 <- MAIHDA:::maihda_strip_brms_weights_term(y | weights(.maihda_sw) ~ x + (1 | s))
  expect_equal(paste(deparse(f1[[2]]), collapse = " "), "y")

  # trials(n) + weights(.maihda_sw) -> the trials() addition term is preserved.
  f2 <- MAIHDA:::maihda_strip_brms_weights_term(
    y | trials(n) + weights(.maihda_sw) ~ x + (1 | s))
  expect_equal(paste(deparse(f2[[2]]), collapse = " "), "y | trials(n)")

  # A user's own weights(w0) term is a genuine conflict and is left in place.
  f3 <- MAIHDA:::maihda_strip_brms_weights_term(y | weights(w0) ~ x)
  expect_equal(paste(deparse(f3[[2]]), collapse = " "), "y | weights(w0)")

  # A formula with no addition term is returned unchanged.
  f4 <- MAIHDA:::maihda_strip_brms_weights_term(y ~ x + (1 | s))
  expect_equal(paste(deparse(f4), collapse = " "), "y ~ x + (1 | s)")
})

test_that("maihda_prepare_brms_sampling_weights re-prepares a derived fit (audit P1 #2)", {
  d <- data.frame(y = rnorm(6), x = rnorm(6),
                  stratum = rep(c("a", "b"), 3),
                  w = c(1, 2, 1, 3, 1, 2))

  # First preparation injects weights(.maihda_sw) and the .maihda_sw column.
  prep1 <- MAIHDA:::maihda_prepare_brms_sampling_weights(
    d, y ~ x + (1 | stratum), "w")
  expect_true(".maihda_sw" %in% names(prep1$data))
  expect_match(paste(deparse(prep1$formula), collapse = " "),
               "weights(.maihda_sw)", fixed = TRUE)

  # A maihda()-/compare_maihda_groups()-derived refit re-enters with BOTH the
  # injected formula term and the .maihda_sw column already present (copied from
  # the prior fit). Previously the reserved-column guard aborted; now it strips
  # and re-prepares cleanly, re-normalizing from the ORIGINAL 'w' column.
  prep2 <- MAIHDA:::maihda_prepare_brms_sampling_weights(
    prep1$data, prep1$formula, "w")
  expect_true(".maihda_sw" %in% names(prep2$data))
  expect_equal(prep2$data$.maihda_sw, prep1$data$.maihda_sw)

  # Exactly one weights() term survives -- no duplication.
  fstr <- paste(deparse(prep2$formula), collapse = " ")
  n_weights <- length(regmatches(fstr, gregexpr("weights(", fstr, fixed = TRUE))[[1]])
  expect_equal(n_weights, 1L)
})

# ---- P2 #3: non-positive precision weights are dropped before the fit --------

test_that("zero/negative precision weights are dropped before the lmer fit (audit P2 #3)", {
  set.seed(11)
  n <- 120
  g <- factor(rep(1:12, each = 10))
  d <- data.frame(x = rnorm(n), stratum = g)
  d$y <- 0.5 * d$x + rnorm(12, sd = 1)[g] + rnorm(n, sd = 0.5)
  w <- rep(1, n); w[3] <- 0; w[7] <- -2          # one zero, one negative

  wmsgs <- testthat::capture_warnings(
    m <- fit_maihda(y ~ x + (1 | stratum), data = d, weights = w))
  expect_true(any(grepl("non-positive or non-finite precision weight", wmsgs)))

  # The two offending rows are dropped, giving a fit IDENTICAL to removing them
  # by hand -- not the degenerate zero-weight lmer fit (logLik -Inf).
  keep <- w > 0
  m_ref <- fit_maihda(y ~ x + (1 | stratum), data = d[keep, ], weights = w[keep])
  expect_equal(as.integer(MAIHDA:::maihda_nobs(m$model)), sum(keep))
  expect_true(is.finite(as.numeric(stats::logLik(m$model))))
  expect_equal(as.numeric(stats::logLik(m$model)),
               as.numeric(stats::logLik(m_ref$model)))

  # The reported VPC is the clean-fit VPC, not a finite value read off a
  # degenerate fit.
  expect_equal(summary(m)$vpc$estimate, summary(m_ref)$vpc$estimate)
})

# ---- P2 #4: missing strata are rejected for individual predictions -----------

test_that("maihda_check_known_strata rejects NA strata only for individual predictions (audit P2 #4)", {
  known <- c("a", "b", "c")
  # Individual: an NA stratum has no random effect and must be rejected (it would
  # otherwise become a silent zero-RE, population-average prediction).
  expect_error(
    MAIHDA:::maihda_check_known_strata(c("a", NA), known, type = "individual"),
    "missing", ignore.case = TRUE)
  # Stratum-level: an NA stratum simply yields no row (empty result), matching
  # lme4, so it is tolerated.
  expect_silent(
    MAIHDA:::maihda_check_known_strata(c("a", NA), known, type = "strata"))
  # A genuinely unknown (non-NA) stratum is still rejected.
  expect_error(
    MAIHDA:::maihda_check_known_strata(c("a", "zzz"), known, type = "individual"),
    "not present in the fitted model")
})

test_that("missing strata are rejected for lme4 individual predictions (audit P2 #4)", {
  set.seed(4)
  n <- 200
  d <- data.frame(
    y = rnorm(n), x = rnorm(n),
    gender = rep(c("F", "M"), n / 2),
    race = rep(c("A", "B"), each = n / 2))
  # The strata shorthand records the defining dimensions, so predict can rebuild
  # 'stratum' from newdata (exercising the missing-dimension path, Case B).
  m <- fit_maihda(y ~ x + (1 | gender:race), data = d)

  # Case A: a supplied NA stratum column is rejected for individual predictions.
  nd_na <- data.frame(x = 0, stratum = NA_character_)
  expect_error(predict_maihda(m, newdata = nd_na, type = "individual"),
               "missing", ignore.case = TRUE)

  # allow_new_levels = TRUE opts into the population-average (fixed-only) value
  # for an unseen (non-NA) stratum -- the opt-in path stays functional. (An
  # all-NA grouping column is a separate lme4 predict limitation, out of scope
  # for this finding, which concerns the silent DEFAULT.)
  nd_new <- data.frame(x = 0, stratum = "ZZ-unseen")
  pa <- predict_maihda(m, newdata = nd_new, type = "individual",
                       allow_new_levels = TRUE)
  expect_length(pa, 1L)
  expect_true(is.finite(as.numeric(pa)))

  # Case B: a missing stratum-defining dimension is rejected.
  nd_dim <- data.frame(x = 0, gender = NA_character_, race = "A")
  expect_error(predict_maihda(m, newdata = nd_dim, type = "individual"),
               "missing", ignore.case = TRUE)

  # type = "strata" with an NA stratum still returns an empty result (unchanged).
  strata_na <- predict_maihda(
    m, newdata = data.frame(x = 0, stratum = NA_character_), type = "strata")
  expect_equal(nrow(strata_na), 0L)
})

test_that("WeMix individual predictions reject a missing stratum (audit P2 #4)", {
  skip_if_not_installed("WeMix")
  set.seed(41)
  n <- 300
  d <- data.frame(
    y = rnorm(n), x = rnorm(n),
    gender = rep(c("F", "M"), n / 2),
    race = rep(c("A", "B"), each = n / 2),
    w = runif(n, 0.5, 2))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, sampling_weights = "w")))
  expect_equal(m$engine, "wemix")

  # Previously silently returned a fixed-only prediction; now rejected like lme4.
  nd <- data.frame(x = 0, stratum = NA_character_)
  expect_error(predict_maihda(m, newdata = nd, type = "individual"),
               "missing", ignore.case = TRUE)
})

test_that("ordinal individual predictions reject a missing stratum (audit P2 #4)", {
  skip_if_not_installed("ordinal")
  set.seed(42)
  n <- 240
  d <- data.frame(
    y = factor(sample(1:3, n, replace = TRUE), ordered = TRUE),
    x = rnorm(n),
    gender = rep(c("F", "M"), n / 2),
    race = rep(c("A", "B"), each = n / 2))
  strata <- make_strata(d, vars = c("gender", "race"))
  d$stratum <- strata$data$stratum
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | stratum), data = d, engine = "ordinal")))

  nd <- data.frame(x = 0, stratum = NA_character_)
  expect_error(predict_maihda(m, newdata = nd, type = "individual"),
               "missing", ignore.case = TRUE)
})
