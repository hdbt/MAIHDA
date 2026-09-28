# Audit 2026-09-28f: the second follow-up to F02 (test-audit-2026-09-28d.R).
# fit_maihda(engine = "brms") accepted every brms response addition term, but the
# package models only two of them: trials() and weights(). Measured on real Stan fits
# (a Poisson with exposure 1-10, and Gaussians with known truth):
#   * rate(expo) enters the likelihood outside the linear predictor brms reports
#     (poisson_log_lpmf(Y | mu + log_denom)), so the count VPC was evaluated at the rate
#     per unit of exposure (lambda 0.54 instead of 2.86): 0.180 where the same model
#     written offset(log(expo)) gives 0.435 (lme4 0.384).
#   * se(s) without sigma = TRUE fixes the residual SD at 0: VPC 1.000 [1.000, 1.000].
#     With sigma = TRUE it was a true-score VPC (0.357 against a true 0.36) the package
#     never defines.
#   * cens(): the latent VPC held (0.353), but maihda_describe() and the observed axis
#     of plot_obs_vs_shrunken() took censored values as exact -- a mean of 0.955 against
#     a true 1.117, plotted beside latent predictions of 1.115.
#   * trunc(), mi() and subset() were accepted as well.
# Per the user: refuse every addition term but trials() and weights(), pointing rate()
# to `+ offset(log(expo))`, which the package reads correctly -- the same model for a
# Poisson; for a negative binomial brms's rate() also multiplies the shape by the
# exposure (shape .* denom in its Stan code) and the offset does not, which the message
# says. Fits saved before the change are left alone. The cumulative path already
# refused every addition term with its own message, which stays.

a28f_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

a28f_data <- function() {
  set.seed(1)
  d <- data.frame(g1 = sample(c("a", "b", "c"), 120, TRUE),
                  g2 = sample(c("u", "v"), 120, TRUE), x = stats::rnorm(120))
  d$stratum <- as.integer(factor(paste(d$g1, d$g2)))
  d$g <- stats::rnorm(120)
  d$s <- stats::runif(120, 0.2, 1)
  d$cnt <- stats::rpois(120, 2)
  d$expo <- stats::runif(120, 1, 5)
  d$cc <- ifelse(d$g > 1, "right", "none")
  d$n <- 10L
  d$y <- stats::rbinom(120, 10, 0.4)
  d$yb <- stats::rbinom(120, 1, 0.5)
  d$w <- stats::runif(120, 0.5, 2)
  d$sw <- stats::runif(120, 0.5, 3)
  d$grp <- rep(c("A", "B"), 60)
  d$ord <- factor(sample(1:4, 120, TRUE), ordered = TRUE)
  d$ymis <- d$g
  d$ymis[3] <- NA
  d
}

# A brms fit_maihda() call that never compiles: brm(empty = TRUE).
a28f_fit <- function(f, family, data = a28f_data(), ...) {
  a28f_quiet(fit_maihda(f, data = data, engine = "brms", family = family,
                        empty = TRUE, ...))
}

a28f_refusal <- "supports the trials\\(\\) and weights\\(\\) addition terms only"

# The refusal an expression raises, or "(accepted)" when it returns, so an assertion
# on the message fails rather than errors when there is none.
a28f_message <- function(expr) {
  tryCatch({
    expr
    "(accepted)"
  }, error = function(e) conditionMessage(e))
}

test_that("engine = 'brms' refuses the addition terms the package does not model", {
  skip_if_not_installed("brms")
  expect_error(a28f_fit(g | se(s) ~ x + (1 | stratum), "gaussian"), a28f_refusal)
  expect_error(a28f_fit(g | se(s, sigma = TRUE) ~ x + (1 | stratum), "gaussian"),
               "carries se\\(\\)")
  expect_error(a28f_fit(cnt | rate(expo) ~ x + (1 | stratum), "poisson"), a28f_refusal)
  expect_error(a28f_fit(cnt | rate(expo) ~ x + (1 | stratum), "negbinomial"),
               "carries rate\\(\\)")
  expect_error(a28f_fit(g | cens(cc) ~ x + (1 | stratum), "gaussian"), "carries cens\\(\\)")
  expect_error(a28f_fit(g | trunc(lb = -3) ~ x + (1 | stratum), "gaussian"),
               "carries trunc\\(\\)")
  expect_error(a28f_fit(ymis | mi() ~ x + (1 | stratum), "gaussian"), "carries mi\\(\\)")
  expect_error(a28f_fit(g | subset(yb) ~ x + (1 | stratum), "gaussian"),
               "carries subset\\(\\)")
  # Beside an allowed term, and with sampling weights (which add weights(.maihda_sw)).
  expect_error(a28f_fit(cnt | weights(w) + rate(expo) ~ x + (1 | stratum), "poisson"),
               "carries rate\\(\\)")
  expect_error(a28f_fit(cnt | rate(expo) ~ x + (1 | stratum), "poisson",
                        sampling_weights = "sw"), "carries rate\\(\\)")
})

test_that("the rate() refusal gives the offset() spelling, and when it is the same model", {
  skip_if_not_installed("brms")
  msg <- a28f_message(a28f_fit(cnt | rate(expo) ~ x + (1 | stratum), "poisson"))
  expect_match(msg, "+ offset(log(expo))", fixed = TRUE)
  expect_match(msg, "remove rate(expo)", fixed = TRUE)
  # The offset is the same model only for a Poisson: brms's negative-binomial rate()
  # also multiplies the shape by the exposure, and the message says so.
  expect_match(msg, "For a Poisson model that is the same model", fixed = TRUE)
  expect_match(msg, "multiplies the shape by the exposure", fixed = TRUE)
  # A non-syntactic name keeps its backticks (brms itself then refuses such a name
  # in any spelling, offset() included).
  d <- a28f_data()
  d$`my expo` <- d$expo
  msg2 <- a28f_message(a28f_fit(cnt | rate(`my expo`) ~ x + (1 | stratum), "poisson",
                                data = d))
  expect_match(msg2, "offset(log(`my expo`))", fixed = TRUE)
  # The recipe is accepted.
  expect_s3_class(a28f_fit(cnt ~ x + (1 | stratum) + offset(log(expo)), "poisson"),
                  "maihda_model")
})

test_that("controls: trials(), weights() and the ordinal path are unchanged", {
  skip_if_not_installed("brms")
  expect_s3_class(a28f_fit(y | trials(n) ~ x + (1 | stratum), "binomial"), "maihda_model")
  expect_s3_class(a28f_fit(g | weights(w) ~ x + (1 | stratum), "gaussian"), "maihda_model")
  expect_s3_class(a28f_fit(y | trials(n) + weights(w) ~ x + (1 | stratum), "binomial"),
                  "maihda_model")
  expect_s3_class(a28f_fit(y | trials(n) ~ x + (1 | stratum), "binomial",
                           sampling_weights = "sw"), "maihda_model")
  expect_s3_class(a28f_fit(g ~ x + (1 | stratum), "gaussian"), "maihda_model")
  expect_error(a28f_fit(ord | thres(3) ~ x + (1 | stratum), "ordinal"),
               "needs a single outcome column")
})

test_that("maihda() and compare_maihda_groups() refuse once, up front", {
  skip_if_not_installed("brms")
  d <- a28f_data()
  d$stratum <- NULL
  expect_error(a28f_quiet(maihda(cnt | rate(expo) ~ x + (1 | g1:g2), data = d,
                                 engine = "brms", family = "poisson", empty = TRUE)),
               a28f_refusal)
  w <- character(0)
  err <- a28f_message(withCallingHandlers(
    compare_maihda_groups(cnt | rate(expo) ~ x + (1 | g1:g2), data = d, group = "grp",
                          engine = "brms", family = "poisson", empty = TRUE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")))
  expect_match(err, a28f_refusal)
  expect_length(w, 0L)
})

test_that("the addition terms are read off a formula or a brmsformula", {
  skip_if_not_installed("brms")
  expect_identical(maihda_brms_addition_terms(y | trials(n) + weights(w) ~ x),
                   c("trials", "weights"))
  expect_identical(maihda_brms_addition_terms(brms::bf(y | rate(e) ~ x)), "rate")
  expect_identical(maihda_brms_addition_terms(y | se(s, sigma = TRUE) + cens(c) ~ x),
                   c("se", "cens"))
  expect_identical(maihda_brms_addition_terms(y ~ x), character(0))
  expect_identical(maihda_brms_addition_terms(~ x), character(0))
  expect_invisible(maihda_brms_check_addition_terms(y | trials(n) ~ x))
  expect_invisible(maihda_brms_check_addition_terms(y ~ x))
  expect_error(maihda_brms_check_addition_terms(y | vreal(z) ~ x), "carries vreal\\(\\)")
})
