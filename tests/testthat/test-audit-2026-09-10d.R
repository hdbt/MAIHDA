# Audit 2026-09-10 (fourth pass): a zero-random-effect prediction was documented
# and reported to users as a "population average".
#
# predict_maihda(allow_new_levels = TRUE) deliberately sets an unseen stratum's
# random effect to zero. The help page called the result "the population-average
# (fixed-effects-only) prediction", and the same phrase appeared in the three
# error messages that offer the escape. On a non-identity link that name is wrong.
# The response-scale marginal mean integrates over the random-effect distribution,
# so for a log link with stratum variance tau^2 it is exp(eta + tau^2/2) --
# exp(tau^2/2) times the exp(eta) the package returns. On the Poisson fit below
# the returned value is 1.5956 while the marginal mean is 2.2115 (confirmed by
# quadrature and by the package's own maihda_latent_re_variance_rows() path),
# 38.6% larger. A logit fit misses the other way, attenuated towards 0.5
# (0.6895 returned vs 0.6708 marginal). Only the Gaussian identity link, and the
# link scale generally, make the two coincide.
#
# The computation was correct and is left EXACTLY as it was; this pass renamed it.
# The assertions therefore run in BOTH directions: the wording must no longer
# claim a population average, AND the numbers must stay the conditional u = 0
# values, so that nobody later swaps in a marginal mean to make the old name true.

# Shared fixture: 12 strata from a 2 x 6 intersection, with a Poisson, a Gaussian
# and a Bernoulli outcome driven by the same stratum effects.
maihda_a9_data <- function(seed = 1, n_g = 12, n_i = 60) {
  set.seed(seed)
  u <- stats::rnorm(n_g, 0, 1)
  d <- data.frame(g = rep(seq_len(n_g), each = n_i),
                  x = stats::rnorm(n_g * n_i))
  d$gender <- factor(ifelse(d$g %% 2 == 0, "F", "M"))
  d$race <- factor(c("A", "B", "C", "D", "E", "F")[((d$g - 1) %/% 2) + 1])
  d$cnt <- stats::rpois(nrow(d), exp(0.2 + 0.3 * d$x + u[d$g]))
  d$y <- 1 + 0.3 * d$x + u[d$g] + stats::rnorm(nrow(d))
  d$b <- stats::rbinom(nrow(d), 1, stats::plogis(0.4 + 0.3 * d$x + u[d$g]))
  d
}

# One row in an intersection the model never saw (race "Z"), covariate profile
# x = 0, with no 'stratum' column so the package rebuilds the unseen label itself.
maihda_a9_unseen_row <- function(train) {
  nd <- train[1, , drop = FALSE]
  nd$stratum <- NULL
  nd$gender <- factor("F", levels = levels(train$gender))
  nd$race <- factor("Z", levels = c(levels(train$race), "Z"))
  nd$x <- 0
  nd
}

# E_u[ g^-1(eta + u) ] for u ~ N(0, tau2), by quadrature -- the estimand the old
# name promised. Finite limits because exp() overflows integrate()'s infinite
# transformation.
maihda_a9_marginal <- function(linkinv, eta, tau2) {
  s <- sqrt(tau2)
  stats::integrate(function(z) linkinv(eta + z) * stats::dnorm(z, 0, s),
                   -20 * s, 20 * s,
                   rel.tol = 1e-10, subdivisions = 2000L)$value
}

maihda_a9_tau2 <- function(fit) {
  as.numeric(lme4::VarCorr(fit$model)$stratum[1, 1])
}


test_that("an unseen-stratum count prediction is the conditional value at u = 0", {
  skip_on_cran()
  d <- maihda_a9_data()
  st <- make_strata(d, vars = c("gender", "race"))
  fit <- fit_maihda(cnt ~ x + (1 | stratum), data = st$data, family = "poisson")

  nd <- maihda_a9_unseen_row(st$data)
  eta <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "link", allow_new_levels = TRUE))
  got <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "response", allow_new_levels = TRUE))
  tau2 <- maihda_a9_tau2(fit)

  # What the package returns: exp(eta), the mean for a stratum whose effect is 0.
  expect_equal(got, exp(eta), tolerance = 1e-12)

  # What "population average" would have meant: the lognormal-corrected marginal
  # mean, which quadrature confirms to be exactly exp(eta + tau2/2).
  marginal <- maihda_a9_marginal(exp, eta, tau2)
  expect_equal(marginal, exp(eta + tau2 / 2), tolerance = 1e-8)

  # The two are materially different, in the direction and by the factor the help
  # page now states -- so the old name was not a harmless synonym.
  expect_gt(tau2, 0.4)
  expect_equal(marginal / got, exp(tau2 / 2), tolerance = 1e-8)
  expect_gt(marginal / got, 1.3)
  expect_false(isTRUE(all.equal(got, marginal, tolerance = 1e-3)))
})


test_that("the package's own marginal count mean does integrate over the effects", {
  skip_on_cran()
  d <- maihda_a9_data()
  st <- make_strata(d, vars = c("gender", "race"))
  fit <- fit_maihda(cnt ~ x + (1 | stratum), data = st$data, family = "poisson")

  nd <- maihda_a9_unseen_row(st$data)
  eta <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "link", allow_new_levels = TRUE))
  tau2 <- maihda_a9_tau2(fit)

  # The latent-scale VPC path already carries the correction the prediction path
  # does not, and calls it a marginal mean; the two names must stay distinct.
  v <- maihda_latent_re_variance_rows(
    reformulas::findbars(stats::formula(fit$model)),
    lme4::VarCorr(fit$model),
    st$data[1, , drop = FALSE]
  )
  expect_equal(as.numeric(v), tau2, tolerance = 1e-10)
  expect_equal(exp(eta + as.numeric(v) / 2),
               maihda_a9_marginal(exp, eta, tau2), tolerance = 1e-8)
})


test_that("the Gaussian identity link makes conditional and marginal coincide", {
  skip_on_cran()
  d <- maihda_a9_data()
  st <- make_strata(d, vars = c("gender", "race"))
  fit <- fit_maihda(y ~ x + (1 | stratum), data = st$data, family = "gaussian")

  nd <- maihda_a9_unseen_row(st$data)
  eta <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "link", allow_new_levels = TRUE))
  got <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "response", allow_new_levels = TRUE))
  tau2 <- maihda_a9_tau2(fit)

  expect_equal(got, eta, tolerance = 0)
  expect_equal(got, maihda_a9_marginal(identity, eta, tau2), tolerance = 1e-8)
})


test_that("an unseen-stratum probability is attenuated relative to the marginal", {
  skip_on_cran()
  d <- maihda_a9_data()
  st <- make_strata(d, vars = c("gender", "race"))
  fit <- fit_maihda(b ~ x + (1 | stratum), data = st$data, family = "binomial")

  nd <- maihda_a9_unseen_row(st$data)
  eta <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "link", allow_new_levels = TRUE))
  got <- as.numeric(predict_maihda(fit, newdata = nd, type = "individual",
                                   scale = "response", allow_new_levels = TRUE))
  tau2 <- maihda_a9_tau2(fit)
  marginal <- maihda_a9_marginal(stats::plogis, eta, tau2)

  expect_equal(got, stats::plogis(eta), tolerance = 1e-12)
  # Returned value above 0.5, so the marginal probability sits strictly between
  # 0.5 and it: the mislabel is not confined to the log link.
  expect_gt(got, 0.5)
  expect_lt(marginal, got)
  expect_gt(marginal, 0.5)
  expect_false(isTRUE(all.equal(got, marginal, tolerance = 1e-3)))
})


test_that("the allow_new_levels hint names the zero random effect, not an average", {
  # No fit and no filesystem, so this one runs on CRAN too -- it is the most
  # direct guard on the finding, and the house pattern is to skip only the
  # model-fitting blocks.
  hint <- maihda_new_levels_hint()

  expect_true(grepl("random effect set to zero", hint, fixed = TRUE))
  expect_true(grepl("not a response-scale population average", hint, fixed = TRUE))
  expect_true(grepl("?predict_maihda", hint, fixed = TRUE))
  # Spelling-agnostic ban, the same rule the help-page block uses: strip the one
  # legitimate mention (the disclaimer) and nothing resembling the old name may
  # remain. Banning only the exact pre-fix string would miss a reworded relapse.
  kept <- gsub("not a response-scale population average", "", hint, fixed = TRUE)
  expect_false(grepl("population average", kept, fixed = TRUE))
  expect_false(grepl("population-average", kept, fixed = TRUE))
})


test_that("every unseen-stratum error path carries the corrected hint", {
  skip_on_cran()
  d <- maihda_a9_data()
  st <- make_strata(d, vars = c("gender", "race"))
  fit <- fit_maihda(cnt ~ x + (1 | stratum), data = st$data, family = "poisson")

  msg_of <- function(expr) {
    conditionMessage(tryCatch(expr, error = identity))
  }

  # (a) rebuilt label for an unseen intersection
  unseen_combo <- msg_of(predict_maihda(fit, newdata = maihda_a9_unseen_row(st$data)))
  # (b) a supplied 'stratum' the model never saw
  nd_named <- st$data[1, c("cnt", "x", "stratum"), drop = FALSE]
  nd_named$stratum <- "NOT_A_STRATUM"
  unseen_named <- msg_of(predict_maihda(fit, newdata = nd_named))
  # (c) a missing (NA) stratum
  nd_na <- nd_named
  nd_na$stratum <- NA_character_
  na_stratum <- msg_of(predict_maihda(fit, newdata = nd_na))
  # (d) a row whose stratum-defining dimension is missing
  nd_dim <- st$data[1, , drop = FALSE]
  nd_dim$stratum <- NULL
  nd_dim$race <- NA
  na_dim <- msg_of(predict_maihda(fit, newdata = nd_dim))

  for (m in list(unseen_combo, unseen_named, na_stratum, na_dim)) {
    expect_true(grepl("allow_new_levels = TRUE", m, fixed = TRUE))
    expect_true(grepl("random effect set to zero", m, fixed = TRUE))
    expect_true(grepl("not a response-scale population average", m, fixed = TRUE))
    kept <- gsub("not a response-scale population average", "", m, fixed = TRUE)
    expect_false(grepl("population average", kept, fixed = TRUE))
    expect_false(grepl("population-average", kept, fixed = TRUE))
  }

  # A stratum-level prediction has no random effect to fall back on, so it must
  # still offer no hint at all -- the fix must not have widened the escape.
  strata_msg <- msg_of(predict_maihda(fit, newdata = nd_named, type = "strata"))
  expect_true(grepl("not present in the fitted model", strata_msg, fixed = TRUE))
  expect_false(grepl("allow_new_levels", strata_msg, fixed = TRUE))
})


test_that("the shipped help pages no longer call the zero-effect value an average", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree (test_local / load_all);
  # R CMD check copies only tests/, so there is nothing to read there.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")

  # Whitespace-NORMALISED, because Rd re-wraps: the pre-fix phrase survived a
  # line-based grep of the sources by sitting across a line break ("...population"
  # / "average..."), which is exactly how one occurrence was missed on the first
  # sweep of this pass. The only legitimate mention is the disclaimer that the
  # value is NOT one, so that phrase is removed before the ban is applied.
  for (f in c("predict_maihda.Rd", "predict.maihda_model.Rd",
              "maihda_wemix_linpred.Rd", "maihda_clmm_linpred.Rd")) {
    rd <- gsub("[[:space:]]+", " ",
               paste(readLines(file.path(man, f), warn = FALSE), collapse = " "))
    kept <- gsub("a response-scale population average", "", rd, fixed = TRUE)
    expect_false(grepl("population average", kept, fixed = TRUE), info = f)
    expect_false(grepl("population-average", kept, fixed = TRUE), info = f)
  }

  pm <- paste(readLines(file.path(man, "predict_maihda.Rd"), warn = FALSE),
              collapse = "\n")
  expect_true(grepl("at a zero random effect", pm, fixed = TRUE))
  expect_true(grepl("not} a response-scale population average", pm, fixed = TRUE))
})
