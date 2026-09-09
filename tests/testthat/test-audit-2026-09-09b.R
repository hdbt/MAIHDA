# Audit 2026-09-09 (b) -- the brms half of the structured-threshold problem.
#
# The clmm half is in test-audit-2026-09-09.R: ordinal::clmm stores FREE threshold
# parameters, so the package had to learn to expand them. brms fails the other way
# round. Its generated quantities already write the full b_Intercept[1..K-1] vector
# whatever the threshold structure, so READING is fine -- but fit_maihda() rebuilt
# the family as brms::cumulative(link = family$link), which threw away every other
# option the user had set on it. brms::cumulative(threshold = "equidistant") was
# therefore fitted with brms's flexible default, silently: not wrong arithmetic,
# but a different model from the one asked for, with nothing in the output saying so.
#
# maihda_brms_cumulative_family() now carries `threshold` and `link_disc` across.
# The marker lists that the family-string and maihda_cumulative() paths produce
# carry neither field, so those paths still get brms's defaults -- asserted below,
# because that is the route almost every ordinal brms fit takes.
#
# The rebuild is deliberately a small pure function so its contract is testable
# without a Stan toolchain; the end-to-end fit is opt-in like every other brms test.

test_that("the cumulative family rebuild preserves user-set options", {
  skip_on_cran()
  skip_if_not_installed("brms")

  # A real brmsfamily keeps what the user asked for.
  eq <- maihda_brms_cumulative_family(brms::cumulative(threshold = "equidistant"))
  expect_identical(eq$family, "cumulative")
  expect_identical(eq$threshold, "equidistant")

  stz <- maihda_brms_cumulative_family(brms::cumulative(threshold = "sum_to_zero"))
  expect_identical(stz$threshold, "sum_to_zero")

  # The link is carried across independently of the threshold.
  pr <- maihda_brms_cumulative_family(brms::cumulative(link = "probit"))
  expect_identical(pr$link, "probit")
  expect_identical(pr$threshold, "flexible")

  # ...and a non-default link_disc survives too.
  ld <- maihda_brms_cumulative_family(
    brms::cumulative(link = "logit", link_disc = "identity"))
  expect_identical(ld$link_disc, "identity")
})

test_that("the plain cumulative marker lists still get brms defaults", {
  skip_on_cran()
  skip_if_not_installed("brms")

  # maihda_cumulative() and the family-string path produce a bare
  # list(family, link) with no threshold field at all; the rebuild must turn that
  # into a usable brmsfamily rather than propagating a NULL into cumulative().
  for (fam in list(maihda_cumulative("logit"), maihda_cumulative("probit"),
                   list(family = "cumulative", link = "logit"))) {
    out <- maihda_brms_cumulative_family(fam)
    expect_s3_class(out, "brmsfamily")
    expect_identical(out$family, "cumulative")
    expect_identical(out$link, fam$link)
    expect_identical(out$threshold, "flexible")
  }
})

test_that("the rebuild leaves the default path's model untouched", {
  skip_on_cran()
  skip_if_not_installed("brms")

  # Confinement for the route almost every ordinal brms fit takes. Compare the
  # MODEL, not the family object: brms::cumulative() builds fresh linkfun/linkinv
  # closures on every call, so it is not even identical() to itself and that test
  # would fail for a reason unrelated to the rebuild. make_stancode()/
  # make_standata() need no Stan toolchain and settle it exactly -- same code and
  # same data means the same model.
  set.seed(11)
  n <- 120
  d <- data.frame(
    x = stats::rnorm(n),
    g = factor(sample(letters[1:4], n, replace = TRUE)),
    y = factor(sample(1:4, n, replace = TRUE), levels = 1:4, ordered = TRUE)
  )

  for (lk in c("logit", "probit")) {
    before <- brms::cumulative(link = lk)                       # the pre-fix expression
    after  <- maihda_brms_cumulative_family(brms::cumulative(link = lk))

    # Every field except the two closures is untouched...
    keep <- setdiff(names(before), c("linkfun", "linkinv"))
    expect_identical(before[keep], after[keep])
    # ...and the closures agree numerically.
    z <- c(-2, -0.5, 0, 0.5, 2)
    expect_equal(before$linkinv(z), after$linkinv(z), tolerance = 0)

    # The model brms would actually build is byte-identical.
    expect_identical(
      as.character(brms::make_stancode(y ~ x + (1 | g), data = d, family = before)),
      as.character(brms::make_stancode(y ~ x + (1 | g), data = d, family = after))
    )
    expect_equal(
      brms::make_standata(y ~ x + (1 | g), data = d, family = before),
      brms::make_standata(y ~ x + (1 | g), data = d, family = after)
    )
  }

  # A structured request, by contrast, must genuinely change the generated model
  # rather than only relabelling the family -- otherwise the fix would be inert.
  sc_flex <- as.character(brms::make_stancode(y ~ x + (1 | g), data = d,
    family = maihda_brms_cumulative_family(brms::cumulative())))
  sc_eq <- as.character(brms::make_stancode(y ~ x + (1 | g), data = d,
    family = maihda_brms_cumulative_family(brms::cumulative(threshold = "equidistant"))))
  expect_false(identical(sc_flex, sc_eq))
  expect_true(grepl("real<lower=0> delta", sc_eq, fixed = TRUE))
  expect_false(grepl("real<lower=0> delta", sc_flex, fixed = TRUE))

  # Both still expand to the full threshold vector in generated quantities, which
  # is why maihda_brms_ordinal_thresholds() needs no clmm-style expansion.
  expect_true(grepl("vector[nthres] b_Intercept", sc_eq, fixed = TRUE))
  expect_true(grepl("vector[nthres] b_Intercept", sc_flex, fixed = TRUE))
})

test_that("an unusable threshold fails loudly rather than being dropped", {
  skip_on_cran()
  skip_if_not_installed("brms")

  # The whole point of the fix is that the request is not silently discarded, so
  # a request brms cannot honour must error rather than fall back to flexible.
  expect_error(
    maihda_brms_cumulative_family(
      list(family = "cumulative", link = "logit", threshold = "nonsense")),
    "flexible"
  )

  # A multi-value threshold is malformed, not a request to use the default, so it
  # must reach brms and be rejected there rather than quietly becoming flexible.
  expect_error(
    maihda_brms_cumulative_family(
      list(family = "cumulative", link = "logit",
           threshold = c("flexible", "equidistant"))))

  # Absent / NA / empty genuinely mean "not set" and fall back to the default.
  for (unset in list(NULL, NA_character_, "")) {
    out <- maihda_brms_cumulative_family(
      list(family = "cumulative", link = "logit", threshold = unset))
    expect_identical(out$threshold, "flexible")
  }
})

test_that("a brms cumulative fit honours a structured threshold end to end", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")

  set.seed(11)
  n <- 300
  d <- data.frame(
    x  = stats::rnorm(n),
    d1 = sample(c("a", "b"), n, replace = TRUE),
    d2 = sample(c("p", "q"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  st  <- interaction(d$d1, d$d2, drop = TRUE)
  eta <- 0.7 * d$x + stats::rnorm(nlevels(st), 0, 0.6)[st]
  cum <- vapply(c(-0.8, 0.15, 1.10, 2.05),
                function(a) stats::plogis(a - eta), numeric(n))
  cum <- cbind(cum, 1)
  cum <- t(apply(cum - cbind(0, cum[, -5, drop = FALSE]), 1, cumsum))
  cum[, 5] <- 1
  d$y <- factor(max.col(stats::runif(n) <= cum, ties.method = "first"),
                levels = 1:5, ordered = TRUE)

  fit <- function(fam) suppressMessages(suppressWarnings(fit_maihda(
    y ~ x + (1 | d1:d2), data = d, family = fam, engine = "brms",
    chains = 1, iter = 600, warmup = 300, seed = 1, refresh = 0, silent = 2)))

  m <- fit(brms::cumulative(threshold = "equidistant"))
  expect_identical(m$model$family$threshold, "equidistant")

  # brms expands the constrained thresholds itself, so the reader needs no
  # counterpart to maihda_clmm_cutpoints(): fixef() carries all K-1 Intercept[k]
  # rows (never the raw `delta`) under every threshold structure.
  thr <- maihda_brms_ordinal_thresholds(m$model)
  expect_length(thr, nlevels(d$y) - 1L)
  expect_false(is.unsorted(thr, strictly = TRUE))

  # The constraint was really imposed: equal spacing to numerical tolerance.
  gaps <- diff(thr)
  expect_equal(max(abs(gaps - mean(gaps))), 0, tolerance = 1e-8)

  # And the response scale reaches the top categories, as for the clmm engine.
  sc <- predict_maihda(m, type = "individual", scale = "response")
  expect_true(all(sc >= 1 & sc <= nlevels(d$y)))
  expect_gt(max(sc), 3.5)

  # Control: the family-string path is untouched and its thresholds stay free.
  mf <- fit("ordinal")
  expect_identical(mf$model$family$threshold, "flexible")
  gaps_f <- diff(maihda_brms_ordinal_thresholds(mf$model))
  expect_gt(max(abs(gaps_f - mean(gaps_f))), 1e-3)
})
