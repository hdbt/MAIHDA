# Audit 2026-09-27 -- a two-level outcome must be classified by what it IS, not by
# how it is spelled.
#
# FINDING. fit_maihda(I(ly > 0.9) ~ x + (1 | stratum), sampling_weights = "w",
# family = binomial()) was refused with
#
#   engine = "wemix" supports a binary (Bernoulli) 0/1 outcome only; aggregated
#   binomial responses (cbind(success, failure) or trials) are not supported.
#
# I(ly > 0.9) is a Bernoulli outcome, not an aggregated binomial. The same outcome
# precomputed into a column fitted, and engine = "lme4" fitted the expression form.
#
# ROOT CAUSE. NOT maihda_is_binary_vector() rejecting a logical -- that was the
# obvious hypothesis and it is REFUTED below: it returns TRUE for a logical vector.
# maihda_analytic_response() short-circuited on !is.symbol(formula[[2]]), so EVERY
# expression response evaluated to NULL and maihda_response_is_binary() was FALSE
# for all of them alike -- I(ly > 0.9), factor(g1) and a genuinely aggregated
# cbind(s, f). The wemix gate reads that flag to mean "aggregated", so a false
# negative surfaced as a message about the wrong thing.
#
# WIDER THAN THE REPORT, measured: the same NULL also drove
#   * family auto-detection -- fit_maihda(I(ly > 0.9) ~ x) with no family= chose
#     gaussian, where the identical outcome as a column chose binomial; and
#   * maihda_describe(), which reported the outcome as continuous/gaussian with no
#     levels and summarised a 0/1 variable by its mean and SD.
# Both are fixed by the same change and are pinned here.
#
# FIX. maihda_analytic_response() evaluates the response instead of refusing it.
# Exclusion of aggregated responses is preserved but now follows from the VALUE:
# cbind(s, f) is a matrix, which maihda_is_binary_vector() rejects on dim(), and
# I(s / n) is non-integral, which it rejects on its values.
#
# THE ONE ASYMMETRY THAT REMAINS, deliberately. maihda_prepare_binomial_response()
# recodes data[[outcome]] BY NAME, so it reaches bare columns only; an expression
# arrives at the engine as it evaluates. A logical, a 0/1 value and a two-level
# factor are all fine, but a character or a 1/2 coding is not -- glmer and WeMix
# stop with "response must be numeric or factor" and "y values must be 0 <= y <= 1",
# neither of which mentions the spelling. fit_maihda() now refuses those two cases
# itself, naming the cause and the remedy. Making them work instead would mean
# materialising the recoded value under a generated name and rewriting the stored
# formula, which changes what print(), predict() and the observed-response fallback
# all see -- a redesign, not a fix, so it is out of scope here.
#
# THE ORDINAL ENGINE DOES NOT MOVE: maihda_ordinal_prepare_response() has its own
# "single outcome column" guard, measured unchanged (pinned in
# test-audit-2026-09-26.R).

audit_0927_frame <- function() {
  set.seed(1); n <- 600
  d <- data.frame(g1 = factor(sample(c("a", "b"), n, TRUE)),
                  g2 = factor(sample(c("p", "q", "r"), n, TRUE)),
                  x  = rnorm(n))
  d <- make_strata(d, c("g1", "g2"))$data
  st <- as.integer(factor(d$stratum))
  d$y  <- exp(0.4 * d$x + rnorm(length(unique(st)), 0, 0.5)[st] +
                rnorm(n, 0, 0.4) + 1)
  d$w  <- runif(n, 0.5, 2)
  d$ly <- log(d$y)
  d$bb <- as.integer(d$ly > 0.9)                 # the precomputed twin
  d$chr  <- ifelse(d$ly > 0.9, "case", "control")
  d$num2 <- ifelse(d$ly > 0.9, 2, 1)             # binary, but coded 1/2
  d$ns <- rbinom(n, 10, 0.5)
  d$nf <- 10L - d$ns
  d
}

audit_0927_q <- function(e) suppressMessages(suppressWarnings(e))


test_that("the premise: the vector test was never the problem", {
  d <- audit_0927_frame()
  # The reported hypothesis, refuted: a logical vector IS binary to the helper.
  expect_true(maihda_is_binary_vector(d$ly > 0.9))
  expect_true(maihda_is_binary_vector(d$chr))
  expect_true(maihda_is_binary_vector(d$num2))
  # ... and a matrix or a non-integral proportion is not, which is what keeps
  # aggregated responses out now that the spelling no longer does.
  expect_false(maihda_is_binary_vector(cbind(d$ns, d$nf)))
  expect_false(maihda_is_binary_vector(d$ns / 10))
})


test_that("an expression response is evaluated, not refused", {
  d <- audit_0927_frame()

  # THE DEFECT: each of these returned NULL / FALSE.
  expect_false(is.null(maihda_analytic_response(I(ly > 0.9) ~ x + (1 | stratum), d)))
  expect_true(maihda_response_is_binary(I(ly > 0.9) ~ x + (1 | stratum), d))
  expect_true(maihda_response_is_binary((ly > 0.9) ~ x + (1 | stratum), d))
  expect_true(maihda_response_is_binary(factor(g1) ~ x + (1 | stratum), d))

  # It evaluates to what it should, over the analytic sample.
  r <- maihda_analytic_response(I(ly > 0.9) ~ x + (1 | stratum), d)
  expect_true(is.logical(r))
  expect_identical(unname(table(r)[["TRUE"]]), 265L)
  expect_identical(unname(table(r)[["FALSE"]]), 335L)

  # A bare symbol is unchanged, and an aggregated response is STILL not binary --
  # now because of its value, not its spelling.
  expect_true(maihda_response_is_binary(bb ~ x + (1 | stratum), d))
  expect_false(maihda_response_is_binary(cbind(ns, nf) ~ x + (1 | stratum), d))
  expect_false(maihda_response_is_binary(I(ns / 10) ~ x + (1 | stratum), d))
  expect_false(maihda_response_is_binary(ly ~ x + (1 | stratum), d))
  expect_false(maihda_response_is_binary(log(y) ~ x + (1 | stratum), d))

  # An ordered-factor EXPRESSION is now seen by the ordinal detector too.
  labs <- c("low", "mid", "high")
  d$yo <- factor(labs[1 + (d$x > -0.4) + (d$x > 0.7)], levels = labs, ordered = TRUE)
  expect_true(maihda_response_is_ordinal(
    factor(yo, levels = rev(labs), ordered = TRUE) ~ x + (1 | stratum), d))
})


test_that("the wemix fit is the same model however the outcome is spelled", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- audit_0927_frame()
  q <- audit_0927_q

  fe <- q(fit_maihda(I(ly > 0.9) ~ x + (1 | stratum), data = d,
                     sampling_weights = "w", family = binomial()))
  fs <- q(fit_maihda(bb ~ x + (1 | stratum), data = d,
                     sampling_weights = "w", family = binomial()))
  expect_s3_class(fe, "maihda_model")

  # Bit-identical coefficients: one model, two spellings. Measured, not invented.
  expect_equal(unname(fe$model$coef), unname(fs$model$coef), tolerance = 1e-12)
  expect_equal(unname(fe$model$coef), c(-0.375252321, 1.322369647), tolerance = 1e-8)
  expect_identical(nrow(fe$data), nrow(fs$data))

  # A parenthesised spelling is the same fit again.
  fp <- q(fit_maihda((ly > 0.9) ~ x + (1 | stratum), data = d,
                     sampling_weights = "w", family = binomial()))
  expect_equal(unname(fp$model$coef), unname(fs$model$coef), tolerance = 1e-12)

  # Downstream reporting agrees across the two spellings.
  expect_equal(unname(q(predict_maihda(fe))), unname(q(predict_maihda(fs))),
               tolerance = 1e-10)
  expect_equal(as.numeric(unlist(q(summary(fe))$vpc)),
               as.numeric(unlist(q(summary(fs))$vpc)), tolerance = 1e-10)
})


test_that("family auto-detection no longer depends on the spelling", {
  skip_on_cran()
  d <- audit_0927_frame()
  q <- audit_0927_q
  fam <- function(f) {
    m <- q(fit_maihda(f, data = d, engine = "lme4"))
    paste0(maihda_model_family(m)$family, "(", maihda_model_family(m)$link, ")")
  }
  # THE DEFECT: the expression form auto-detected gaussian(identity) -- a linear
  # probability model -- where the identical outcome as a column chose binomial.
  expect_identical(fam(bb ~ x + (1 | stratum)), "binomial(logit)")
  expect_identical(fam(I(ly > 0.9) ~ x + (1 | stratum)), "binomial(logit)")
  expect_identical(fam((ly > 0.9) ~ x + (1 | stratum)), "binomial(logit)")

  # A continuous response is still gaussian, transformed or not.
  expect_identical(fam(ly ~ x + (1 | stratum)), "gaussian(identity)")
  expect_identical(fam(log(y) ~ x + (1 | stratum)), "gaussian(identity)")
  # A proportion stays gaussian: it is not a Bernoulli outcome.
  expect_identical(fam(I(ns / 10) ~ x + (1 | stratum)), "gaussian(identity)")
})


test_that("maihda_describe reports an expression outcome as binary", {
  skip_on_cran()
  d <- audit_0927_frame()
  q <- audit_0927_q
  de <- q(maihda_describe(I(ly > 0.9) ~ x + (1 | stratum), data = d))
  ds <- q(maihda_describe(bb ~ x + (1 | stratum), data = d))

  # THE DEFECT: the expression form reported family "gaussian", outcome_kind
  # "continuous", family_detected FALSE and no outcome levels at all.
  expect_identical(de$family, "binomial")
  expect_identical(de$family_link, "logit")
  expect_true(de$family_detected)
  expect_identical(de$outcome_kind, "binomial")
  expect_identical(de$family, ds$family)
  expect_identical(de$outcome_kind, ds$outcome_kind)

  # The event count matches the column spelling exactly; only the LEVEL NAMES
  # differ, and honestly so -- a logical is FALSE/TRUE where the column is 0/1.
  expect_equal(as.numeric(unlist(de$outcome_overall)),
               as.numeric(unlist(ds$outcome_overall)), tolerance = 1e-12)
  expect_identical(as.character(de$outcome_levels$level), c("FALSE", "TRUE"))
  expect_identical(as.character(ds$outcome_levels$level), c("0", "1"))
  expect_identical(de$outcome_levels$n, ds$outcome_levels$n)
})


test_that("an expression the engine cannot consume is refused by name", {
  skip_on_cran()
  d <- audit_0927_frame()
  q <- audit_0927_q

  # The recoding that turns "case"/"control" or 1/2 into 0/1 applies to a bare
  # COLUMN only. Those two spellings would otherwise reach the engine as they are
  # and draw an error that never mentions the spelling, so they are refused here
  # with the cause and the remedy named. Both engines, same message.
  for (eng in list(list(engine = "lme4"), list(sampling_weights = "w"))) {
    for (f in list(I(chr) ~ x + (1 | stratum), I(num2) ~ x + (1 | stratum))) {
      expect_error(
        do.call(fit_maihda, c(list(f, data = d, family = binomial()), eng)),
        "written as an expression", info = deparse(f[[2]]))
      expect_error(
        do.call(fit_maihda, c(list(f, data = d, family = binomial()), eng)),
        "recoding applies to a bare outcome column only", info = deparse(f[[2]]))
    }
  }
  # The message names the offending value, not just the shape.
  expect_error(q(fit_maihda(I(chr) ~ x + (1 | stratum), data = d,
                            engine = "lme4", family = binomial())),
               "character values")
  expect_error(q(fit_maihda(I(num2) ~ x + (1 | stratum), data = d,
                            engine = "lme4", family = binomial())),
               "the values 1/2")

  # ... and the SAME outcomes stored as columns still fit, because those ARE
  # recoded. That contrast is the whole point of the message.
  expect_s3_class(q(fit_maihda(chr ~ x + (1 | stratum), data = d,
                               engine = "lme4", family = binomial())),
                  "maihda_model")
  expect_s3_class(q(fit_maihda(num2 ~ x + (1 | stratum), data = d,
                               engine = "lme4", family = binomial())),
                  "maihda_model")

  # A SYMBOL response never reaches this branch, so nothing that fitted before is
  # newly refused.
  expect_s3_class(q(fit_maihda(bb ~ x + (1 | stratum), data = d,
                               engine = "lme4", family = binomial())),
                  "maihda_model")
})


test_that("the engine-readiness predicate accepts exactly what the engines take", {
  # Direct pins on the helper, so its contract does not drift from the message above.
  expect_true(maihda_binomial_expression_is_engine_ready(c(TRUE, FALSE, NA)))
  expect_true(maihda_binomial_expression_is_engine_ready(c(0L, 1L)))
  expect_true(maihda_binomial_expression_is_engine_ready(c(0, 1)))
  expect_true(maihda_binomial_expression_is_engine_ready(c(0, 0)))
  expect_true(maihda_binomial_expression_is_engine_ready(factor(c("a", "b"))))
  expect_false(maihda_binomial_expression_is_engine_ready(c("a", "b")))
  expect_false(maihda_binomial_expression_is_engine_ready(c(1, 2)))
  expect_false(maihda_binomial_expression_is_engine_ready(cbind(1:2, 2:1)))
  expect_false(maihda_binomial_expression_is_engine_ready(NULL))
})
