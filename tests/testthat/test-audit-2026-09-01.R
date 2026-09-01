# Regression tests for the 2026-09-01 statistical audit pass.
#
# FINDING: the ordinal proportional-odds diagnostic tested the wrong model. It
# dropped the stratum random intercept, compared two fixed-only ordinal::clm()
# fits, and referred the resulting nominal-effects LRT to an ordinary chi-squared
# distribution -- then used that p-value to raise an automatic caveat on the
# fitted clmm.
#
# A conditional cumulative model with a normal random intercept does NOT in general
# remain an ordinary proportional-odds model once the random intercept is
# marginalised away: the implied marginal cumulative-logit slopes differ across
# thresholds for any non-zero stratum variance. The chi-squared null was therefore
# false under the correctly specified model, and the flag's rejection rate grew
# without bound in n (measured: 5% at tau = 0, but 24% at a realistic 7% stratum
# VPC with n = 96,000, and 58% at tau = 3 with n = 4,800).
#
# CAVEAT worth knowing when reading these fixtures: an EXACTLY SYMMETRIC threshold
# configuration is the exception -- the marginal slopes then coincide and the
# fixed-only statistic is valid. Three categories cut at -c and +c is that case,
# and it is what make_po_data() in test-audit-2026-08-03.R uses, which is why the
# old screen behaved there. The fixtures below deliberately use four categories at
# (-1.5, 0, 1.5) so the confounding is present.
#
# FIX: the automatic flag and the chi-squared p-value are gone. The statistic is
# still computed and stored, under a name that says what it is
# ($adequacy$marginal_po_proxy, carrying lrt/df/n_terms only), and
# maihda_proportional_odds_test() calibrates it by parametric bootstrap under the
# fitted clmm for callers who want an actual test.

# Data generated EXACTLY from the conditional proportional-odds mixed model that
# clmm fits: proportional odds holds conditionally, by construction. Any flag
# raised on this data is a false positive. `violate = TRUE` instead gives genuinely
# threshold-specific slopes, where a rejection is correct.
make_po_audit_data <- function(seed, sigma, n_per = 100, beta = 0.8,
                               thr = c(-1.5, 0, 1.5), violate = FALSE) {
  set.seed(seed)
  strata <- expand.grid(a = factor(1:4), b = factor(1:3))
  d <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i)
    data.frame(a = strata$a[i], b = strata$b[i], x = stats::rnorm(n_per))))
  idx <- interaction(d$a, d$b, drop = TRUE)
  u <- stats::rnorm(nlevels(idx), 0, sigma)[as.integer(idx)]
  if (violate) {
    lat <- u + stats::rlogis(nrow(d))
    y <- ifelse(lat < thr[1] - 2.0 * d$x, 1L,
         ifelse(lat < thr[2] - 0.8 * d$x, 2L,
         ifelse(lat < thr[3] + 1.2 * d$x, 3L, 4L)))
  } else {
    lat <- beta * d$x + u + stats::rlogis(nrow(d))
    y <- cut(lat, c(-Inf, thr, Inf), labels = FALSE)
  }
  d$y <- factor(y, levels = 1:4, ordered = TRUE)
  d
}

# ---- the removed flag: deterministic, no model fitting ----------------------

test_that("the fixed-only proportional-odds statistic carries no p-value or flag", {
  skip_if_not_installed("ordinal")
  d <- make_po_audit_data(1004, sigma = 3)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))

  po <- maihda_ordinal_po_stat(m$model)
  expect_false(is.null(po))
  expect_setequal(names(po), c("lrt", "df", "n_terms"))
  # a chi-squared p-value must not be reconstructible from the stored object:
  # its presence is what let the confounded statistic be read as evidence
  expect_null(po$min_p)
  expect_null(po$p)
  expect_null(po$flag)
})

test_that("the proportional-odds proxy is stored descriptively and never flagged", {
  skip_if_not_installed("ordinal")
  # tau = 3 with n = 1200: the exact regime in which the old chi-squared screen
  # false-flagged correctly specified models (3 of 10 seeds; this is one of them)
  d <- make_po_audit_data(1004, sigma = 3)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))
  ad <- m$diagnostics$adequacy

  expect_true("marginal_po_proxy" %in% names(ad))
  expect_null(ad$proportional_odds)          # the old, flag-bearing key is gone
  expect_null(ad$marginal_po_proxy$flag)
  expect_gt(ad$marginal_po_proxy$lrt, 0)

  # and no caveat reaches the user, on the object or in print()
  expect_false(any(grepl("Proportional odds", maihda_format_adequacy(ad),
                         fixed = TRUE)))
  expect_false(any(grepl("Proportional odds",
                         utils::capture.output(print(m)), fixed = TRUE)))
})

test_that("the dead chi-squared cut-off is no longer part of the thresholds", {
  th <- maihda_adequacy_thresholds()
  expect_null(th$ordinal_po_p)
  # the checks that DO still flag keep their cut-offs untouched
  expect_equal(th$overdispersion_ratio, 1.5)
  expect_equal(th$autocorr_min, 0.3)
})

# ---- the calibrated replacement ---------------------------------------------

test_that("maihda_proportional_odds_test calibrates away the random-effect confounding", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # correctly specified conditional PO model with a LARGE stratum variance --
  # precisely where the old screen broke. The chi-squared reference calls this
  # fit a violation; the bootstrap, simulating under the fitted clmm, does not.
  d <- make_po_audit_data(1004, sigma = 3)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))
  tt <- maihda_proportional_odds_test(m, n_sim = 199, seed = 1)

  expect_s3_class(tt, "maihda_po_test")
  expect_lt(tt$p_chisq, 0.01)               # what the old screen saw: "violated"
  # The claim under test is that CALIBRATION MOVES THE VERDICT, not that this one
  # seed lands above any particular cut-off. Under a correct model a calibrated
  # p is ~Uniform(0, 1), so pinning "p >= 0.05" here would pin a coincidence --
  # and it sat exactly on the boundary (9 of 199 null draws at or above the
  # observed statistic, where 9 is the smallest count that clears 0.05). Assert
  # the order-of-magnitude shift, and that the fit is no longer condemned at the
  # level the chi-squared screen condemned it at, both with real margin.
  expect_gt(tt$p_value / tt$p_chisq, 10)
  expect_gt(tt$p_value, 0.01)
  expect_equal(tt$n_sim + tt$n_failed, 199L)
  expect_length(tt$null_lrt, tt$n_sim)

  # p_chisq is retained on the object -- the comparison above needs it -- but the
  # print method must show ONE p-value. Printing the chi-squared one beside the
  # bootstrap one offers two answers to a single question, and it is an inference
  # against a null this pass established is false.
  txt <- paste(utils::capture.output(print(tt)), collapse = " ")
  expect_match(txt, "Bootstrap p-value")
  expect_match(txt, "Nominal-effects LRT")
  expect_false(grepl("chisq", txt, ignore.case = TRUE))
  expect_false(grepl(sprintf("%.4f", tt$p_chisq), txt, fixed = TRUE))
})

test_that("maihda_proportional_odds_test still rejects a genuine violation", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_audit_data(311, sigma = 0.5, violate = TRUE)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))
  tt <- maihda_proportional_odds_test(m, n_sim = 99, seed = 311)
  expect_lt(tt$p_value, 0.05)
  expect_match(paste(utils::capture.output(print(tt)), collapse = " "),
               "Bootstrap p-value")
})

test_that("maihda_proportional_odds_test is reproducible and validates its inputs", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_po_audit_data(312, sigma = 0.5)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))

  expect_equal(maihda_proportional_odds_test(m, n_sim = 49, seed = 5)$p_value,
               maihda_proportional_odds_test(m, n_sim = 49, seed = 5)$p_value)
  expect_error(maihda_proportional_odds_test(m, n_sim = 0), "whole number")
  expect_error(maihda_proportional_odds_test(m, n_sim = 2.5), "whole number")

  # a null (covariate-free) cumulative fit has no covariate slopes to test
  m0 <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ (1 | a:b), data = d, family = "ordinal")))
  expect_null(maihda_ordinal_po_stat(m0$model))
  expect_error(maihda_proportional_odds_test(m0, n_sim = 5), "no covariate slopes")

  # and a non-cumulative fit is refused rather than silently mis-analysed
  expect_error(maihda_proportional_odds_test(structure(list(), class = "maihda_model")),
               "engine = 'ordinal'")
})

# ---- the shared statistic, so both consumers stay in step -------------------

test_that("maihda_po_lrt is the single source of the nominal-effects statistic", {
  skip_if_not_installed("ordinal")
  d <- make_po_audit_data(1004, sigma = 3)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))

  direct <- maihda_po_lrt(stats::model.frame(m$model), "y", "x")
  via_stat <- maihda_ordinal_po_stat(m$model)
  expect_equal(direct$lrt, via_stat$lrt)
  expect_equal(direct$df, via_stat$df)

  # degrades to NULL rather than erroring when the refit cannot be done
  expect_null(maihda_po_lrt(data.frame(y = 1, x = 1), "y", "nonexistent_term"))
})

test_that("the bootstrap's frame/data row alignment holds under NAs and reordering", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  # maihda_proportional_odds_test() pairs eta_fixed and the stratum index (built
  # from $data) with the response it overwrites in the clmm model frame. If those
  # two ever fell out of row order the null distribution would be silently wrong,
  # so pin the assumption rather than trusting it.
  d <- make_po_audit_data(7, sigma = 0.6, n_per = 80)
  set.seed(7)
  d$x[sample(nrow(d), 40)] <- NA          # rows the fit must drop
  d <- d[sample(nrow(d)), ]               # and an order the frame must preserve

  m <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ x + (1 | a:b), data = d, family = "ordinal")))
  dat <- stats::model.frame(m$model)

  expect_equal(nrow(dat), nrow(m$data))
  expect_identical(rownames(dat), rownames(m$data))
  expect_identical(as.character(dat$y), as.character(m$data$y))
  expect_equal(as.numeric(dat$x), as.numeric(m$data$x))

  expect_s3_class(maihda_proportional_odds_test(m, n_sim = 19, seed = 1),
                  "maihda_po_test")
})
