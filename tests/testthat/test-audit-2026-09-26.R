# Audit 2026-09-26 -- the observed axis has to be on the FITTED response scale.
#
# FINDING (external audit, F06). plot_obs_vs_shrunken() reads the observed outcome
# through maihda_observed_response_from_model_frame(), which tries
# stats::model.response(object$data) first and, when that fails, fell back to
#
#   outcome_var <- all.vars(formula_obj)[1]
#   data[[outcome_var]]
#
# -- the first VARIABLE named in the formula, taken raw out of the data. For a
# response that is an EXPRESSION that is the wrong column: log(y) ~ x has
# all.vars()[1] == "y", so the x-axis carried raw y while the y-axis (fitted and
# shrunken) stayed on the log scale. The two axes were not comparable and the
# y = x diagonal, which the panel draws as its only reference, was meaningless.
#
# WHICH ENGINES REACH THE FALLBACK (measured, not inferred):
#
#   lme4     object$data is a real model frame  -> model.response() works. SAFE.
#   brms     model.frame.brmsfit() returns $data, which DOES carry a "terms"
#            attribute and the evaluated response as its first column
#            (names: "log(y)", "y", "x", "stratum") -> model.response() works. SAFE.
#   ordinal  object$data is the analytic data frame, no terms -> fallback taken,
#            but maihda_ordinal_prepare_response() refuses any non-symbol response
#            up front ("single outcome column"), so it never carries an expression.
#   wemix    object$data is the analytic data frame, no terms -> fallback taken.
#            The binomial path is also symbol-only (maihda_response_is_binary()
#            returns FALSE for a call, and fit_maihda() then refuses the fit), but
#            the GAUSSIAN path accepts log(y) ~ x and is where this bites.
#
# So the reachable defect is wemix + a transformed Gaussian response. The tell is
# that the SAME model written two ways disagreed: fitting log(y) ~ x and fitting
# ly ~ x after ly <- log(y) gave identical coefficients and identical shrunken
# estimates, but observed values differing by up to 2.992642 on the fixture below
# (stratum 1: 3.5385491 against a correct 1.1000563, where the shrunken values span
# only 0.4957 .. 1.2452).
#
# maihda_describe() already got its own outcome right, via
# maihda_describe_response_expr() + eval(); the plot fallback now reuses that same
# expression extractor, so a bare symbol still takes the column verbatim and brms
# addition terms (y | trials(n)) still resolve to the leftmost leaf.
#
# It rebuilds the value with model.frame(resp ~ 1, na.action = na.pass) rather than
# a bare eval(), CHOSEN BY MEASUREMENT: eval() gives the right values but the wrong
# object -- I(y > 3) keeps the "AsIs" wrapper that model.frame() drops, and carries
# no row names -- where the model.frame route matches what lme4's model.response()
# returns in class and names alike, for log(y), I(y > 3) and a bare symbol.

audit_0926_frame <- function(seed = 1, n = 600) {
  set.seed(seed)
  d <- data.frame(
    g1 = factor(sample(c("a", "b"), n, TRUE)),
    g2 = factor(sample(c("p", "q", "r"), n, TRUE)),
    x  = rnorm(n)
  )
  d <- make_strata(d, c("g1", "g2"))$data
  st <- as.integer(factor(d$stratum))
  d$y  <- exp(0.4 * d$x + rnorm(length(unique(st)), 0, 0.5)[st] +
                rnorm(n, 0, 0.4) + 1)
  d$ly <- log(d$y)
  d$w  <- runif(n, 0.5, 2)
  d
}

# A frame shaped like the one wemix/ordinal store on object$data: a plain
# data.frame with no "terms" attribute, so model.response() cannot find the
# response and the formula fallback is what runs.
audit_0926_plain <- function() {
  data.frame(
    y = c(3, 7, 11, 2),
    n = c(10, 20, 12, 20),
    f = c(7, 13, 1, 18),
    x = c(-1, 0, 1, 2),
    stratum = factor(c("A", "B", "C", "D"))
  )
}


test_that("the observed-response fallback evaluates a transformed response", {
  d <- audit_0926_plain()
  expect_null(attr(d, "terms"))
  expect_error(stats::model.response(d))          # the fallback really is reached

  # model.frame() names the response by row, exactly as the lme4 route does (see the
  # class/names parity check below), so compare values here.
  got <- function(fo) unname(MAIHDA:::maihda_observed_response_from_model_frame(d, fo))

  # THE DEFECT: log(y) came back as raw y.
  expect_equal(got(log(y) ~ x + (1 | stratum)), log(d$y))
  expect_false(isTRUE(all.equal(got(log(y) ~ x + (1 | stratum)), d$y)))

  # Other expression spellings of a response.
  expect_equal(got(sqrt(y) ~ x + (1 | stratum)), sqrt(d$y))
  # Rebuilt through model.frame(), so I() does NOT leak its "AsIs" wrapper -- a bare
  # eval() would keep it, where the lme4 route returns plain logical (measured).
  expect_equal(got(I(y > 3) ~ x + (1 | stratum)), d$y > 3)
  expect_equal(got(I(y / n) ~ x + (1 | stratum)), d$y / d$n)

  # A two-column response stays a matrix, which the binomial branch of
  # maihda_observed_outcome_for_plot() needs; all.vars()[1] gave it the successes
  # alone, i.e. a mean success COUNT where the y-axis is a probability.
  cb <- MAIHDA:::maihda_observed_response_from_model_frame(
    d, cbind(y, f) ~ x + (1 | stratum))
  expect_true(is.matrix(cb))
  expect_equal(dim(cb), c(4L, 2L))
  expect_equal(unname(cb[, 1]), d$y)
  od <- MAIHDA:::maihda_observed_outcome_for_plot(
    cb, list(family = "binomial", link = "logit"))
  expect_equal(od$numerator / od$denominator, d$y / (d$y + d$f))

  # The rebuild runs in the FORMULA's environment, so a transformation the user
  # defined beside the model resolves the same way it did at fit time.
  local({
    halve <- function(v) v / 2
    expect_equal(got(halve(y) ~ x + (1 | stratum)), d$y / 2)
  })
})


test_that("the fallback leaves bare-symbol and trials() responses untouched", {
  d <- audit_0926_plain()

  # A bare symbol is still the column verbatim -- including its class, so a factor
  # outcome does not get coerced on its way to the extractor.
  expect_identical(
    MAIHDA:::maihda_observed_response_from_model_frame(d, y ~ x + (1 | stratum)),
    d$y)
  expect_identical(
    MAIHDA:::maihda_observed_response_from_model_frame(d, stratum ~ x + (1 | stratum)),
    d$stratum)

  # brms addition terms hang off the response with `|`; the outcome is the leftmost
  # leaf and the trials() denominator is recovered separately at the call site.
  expect_identical(
    MAIHDA:::maihda_observed_response_from_model_frame(d, y | trials(n) ~ (1 | stratum)),
    d$y)
  expect_equal(MAIHDA:::maihda_trials_from_formula(y | trials(n) ~ (1 | stratum), d), d$n)

  # An absent outcome column is still refused rather than resolved against the
  # formula environment, where a same-named object could be sitting -- `zz` is
  # defined right here, so this frame IS the formula's environment.
  zz <- rep(99, 4)
  expect_error(
    MAIHDA:::maihda_observed_response_from_model_frame(d, zz ~ x + (1 | stratum)),
    "Outcome variable not found")

  # The caller pairs the result with object$data row by row, so a response that is
  # not one value per row is refused instead of being recycled into a wrong plot.
  expect_error(
    MAIHDA:::maihda_observed_response_from_model_frame(d, mean(y) ~ x + (1 | stratum)),
    "one value per row")
})


test_that("a brmsformula is unwrapped the way the trials() helper unwraps it", {
  # Both helpers read the same object$formula at the same call site, so both must
  # unwrap a brmsformula -- otherwise the outcome and its denominator could come
  # from two different places. Its own block, so the skip cannot swallow the
  # assertions above.
  skip_if_not_installed("brms")
  d <- audit_0926_plain()
  bf <- brms::bf(log(y) | trials(n) ~ x + (1 | stratum))
  expect_equal(
    unname(MAIHDA:::maihda_observed_response_from_model_frame(d, bf)), log(d$y))
  expect_equal(MAIHDA:::maihda_trials_from_formula(bf, d), d$n)
})


test_that("a wemix fit plots the same observed values however the response is spelled", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926_frame()

  q <- function(e) suppressMessages(suppressWarnings(e))
  ft <- q(fit_maihda(log(y) ~ x + (1 | stratum), data = d, sampling_weights = "w"))
  fp <- q(fit_maihda(ly ~ x + (1 | stratum), data = d, sampling_weights = "w"))

  # The two fits ARE the same model: same analytic rows, same coefficients.
  expect_equal(unname(ft$model$coef), unname(fp$model$coef))

  # The fallback is what runs here -- object$data is the analytic data frame.
  expect_null(attr(ft$data, "terms"))

  pt <- q(MAIHDA:::plot_obs_vs_shrunken(ft, summary(ft)))
  pp <- q(MAIHDA:::plot_obs_vs_shrunken(fp, summary(fp)))
  m <- merge(pt$data[, c("stratum", "observed", "n", "shrunken")],
             pp$data[, c("stratum", "observed", "n", "shrunken")],
             by = "stratum", suffixes = c(".t", ".p"))
  expect_gt(nrow(m), 1)

  # The y-axis never moved -- the defect was entirely on the observed axis.
  expect_equal(m$shrunken.t, m$shrunken.p)
  expect_equal(m$n.t, m$n.p)
  # ... and the x-axis now agrees too.
  expect_equal(m$observed.t, m$observed.p)

  # Against an outside computation, not just against the sibling spelling: the
  # observed value is the prior-weighted mean of log(y) within the stratum.
  hand <- vapply(split(seq_len(nrow(ft$data)), ft$data$stratum), function(i)
    sum(ft$data$w[i] * log(ft$data$y[i])) / sum(ft$data$w[i]), numeric(1))
  expect_equal(pt$data$observed, unname(hand[as.character(pt$data$stratum)]))

  # Both axes are on the log scale, so the panel's y = x diagonal means something:
  # the observed values sit in the same range as the shrunken ones, not e^that.
  expect_lt(max(abs(pt$data$observed - pt$data$shrunken)), 0.5)

  # The same through the PUBLIC route users actually call, not just the internal
  # panel builder.
  pub <- q(plot(ft, type = "obs_vs_shrunken", summary_obj = summary(ft)))
  expect_s3_class(pub, "ggplot")
  expect_equal(pub$data$observed[order(as.character(pub$data$stratum))],
               pt$data$observed[order(as.character(pt$data$stratum))])
})


test_that("the lme4 path keeps taking model.response() for a transformed outcome", {
  skip_on_cran()
  d <- audit_0926_frame()
  q <- function(e) suppressMessages(suppressWarnings(e))

  fl <- q(fit_maihda(log(y) ~ x + (1 | stratum), data = d, engine = "lme4"))
  # A real model frame: the early return, not the formula fallback, is what runs.
  expect_false(is.null(attr(fl$data, "terms")))
  expect_equal(
    as.numeric(MAIHDA:::maihda_observed_response_from_model_frame(fl$data, fl$formula)),
    log(d$y))

  pl <- q(MAIHDA:::plot_obs_vs_shrunken(fl, summary(fl)))
  hand <- vapply(split(seq_len(nrow(d)), d$stratum), function(i)
    mean(log(d$y[i])), numeric(1))
  expect_equal(pl$data$observed, unname(hand[as.character(pl$data$stratum)]))

  # The two routes -- model.response() on a real model frame, and the formula
  # fallback on a frame without one -- must agree in CLASS as well as value, or a
  # downstream branch of maihda_observed_outcome_for_plot() could pick differently
  # for the same model depending only on which engine stored the data.
  bare <- as.data.frame(d[, c("y", "x", "stratum")])
  for (fo in list(log(y) ~ x + (1 | stratum),
                  I(y > 3) ~ x + (1 | stratum),
                  y ~ x + (1 | stratum))) {
    fit <- q(fit_maihda(fo, data = d, engine = "lme4",
                        family = if (identical(deparse(fo[[2]]), "I(y > 3)"))
                          "binomial" else "gaussian"))
    via_frame <- MAIHDA:::maihda_observed_response_from_model_frame(fit$data, fit$formula)
    via_fallback <- MAIHDA:::maihda_observed_response_from_model_frame(bare, fo)
    expect_identical(class(via_frame), class(via_fallback))
    expect_equal(unname(via_frame), unname(via_fallback))
    # A bare symbol short-circuits to the column itself on the fallback side (HEAD
    # behaviour, kept); an EXPRESSION goes through model.frame() and is named by row
    # just as the lme4 frame is.
    if (!is.symbol(fo[[2]])) {
      expect_identical(names(via_frame), names(via_fallback))
    }
  }
})


test_that("the ordinal engine still refuses an expression response", {
  # The ordinal guard is unchanged: clmm needs a single outcome COLUMN, and
  # maihda_ordinal_prepare_response() refuses any non-symbol response up front, so
  # the observed-response fallback above never carries one on that path.
  skip_on_cran()
  d <- audit_0926_frame()
  q <- function(e) suppressMessages(suppressWarnings(e))

  labs <- c("low", "mid", "high")
  set.seed(3)
  d$yo <- factor(labs[1 + (d$x > -0.4) + (d$x > 0.7)], levels = labs, ordered = TRUE)

  skip_if_not_installed("ordinal")
  expect_error(
    q(fit_maihda(factor(yo, levels = rev(labs), ordered = TRUE) ~ x + (1 | stratum),
                 data = d, engine = "ordinal", family = "ordinal")),
    "single outcome column")
})


test_that("a wemix binomial fit now ACCEPTS an expression response", {
  # CHANGED (audit 2026-09-27). This block used to pin the opposite: the wemix
  # engine refused I(y > 3) with a message about AGGREGATED binomial responses,
  # because maihda_analytic_response() returned NULL for any non-symbol response
  # and maihda_response_is_binary() was therefore FALSE. That was a false negative,
  # not a real guard -- the response is Bernoulli -- and the fallback above turns
  # out to handle it correctly, which is what this block now checks rather than
  # assumes.
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0926_frame()
  q <- function(e) suppressMessages(suppressWarnings(e))

  fe <- q(fit_maihda(I(y > 3) ~ x + (1 | stratum), data = d,
                     family = "binomial", sampling_weights = "w"))
  expect_s3_class(fe, "maihda_model")

  # It is the SAME model as the precomputed column, to the last bit.
  d$b3 <- as.integer(d$y > 3)
  fs <- q(fit_maihda(b3 ~ x + (1 | stratum), data = d,
                     family = "binomial", sampling_weights = "w"))
  expect_equal(unname(fe$model$coef), unname(fs$model$coef), tolerance = 1e-12)

  # ... and the observed-response fallback -- the helper this file exists for --
  # carries the expression correctly on the binomial path: the evaluated logical,
  # equal row for row to the precomputed 0/1 column.
  expect_null(attr(fe$data, "terms"))
  oe <- MAIHDA:::maihda_observed_response_from_model_frame(fe$data, fe$formula)
  expect_true(is.logical(oe))
  expect_equal(as.numeric(oe),
               as.numeric(MAIHDA:::maihda_observed_response_from_model_frame(
                 fs$data, fs$formula)))

  # An AGGREGATED response is still refused -- now because of what it evaluates to
  # (a matrix, and a non-integral proportion) rather than how it is spelled.
  d$ns <- rbinom(nrow(d), 10, 0.5)
  d$nf <- 10L - d$ns
  expect_error(
    q(fit_maihda(cbind(ns, nf) ~ x + (1 | stratum), data = d,
                 family = "binomial", sampling_weights = "w")),
    "binary \\(Bernoulli\\) 0/1 outcome only")
  expect_error(
    q(fit_maihda(I(ns / 10) ~ x + (1 | stratum), data = d,
                 family = "binomial", sampling_weights = "w")),
    "binary \\(Bernoulli\\) 0/1 outcome only")
})
