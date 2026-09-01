# Negative-binomial MAIHDA: the family-label normalization that makes lme4's
# theta-embedding "Negative Binomial(<theta>)" comparable across fits, the
# glmer.nb routing in fit_maihda(), the Nakagawa et al. (2017) latent-scale
# level-1 variance log1p(1/mu + 1/theta), and the brms 'shape' plumbing.

make_nb_data <- function(seed = 4242, n = 700, theta = 1.5) {
  set.seed(seed)
  d <- data.frame(
    gender = sample(c("F", "M"), n, replace = TRUE),
    race = sample(c("A", "B", "C"), n, replace = TRUE),
    edu = sample(c("low", "high"), n, replace = TRUE),
    age = rnorm(n),
    stringsAsFactors = FALSE
  )
  stratum <- interaction(d$gender, d$race, d$edu, drop = TRUE)
  u <- rnorm(nlevels(stratum), sd = 0.4)[stratum]
  d$y <- rnbinom(n, mu = exp(1.2 + 0.2 * (d$gender == "M") + 0.1 * d$age + u),
                 size = theta)
  d
}

fit_nb <- function(d, formula = y ~ age + (1 | gender:race:edu)) {
  suppressWarnings(fit_maihda(formula, data = d, family = "negbinomial"))
}

# ---- family-name normalization (no fit required) -----------------------------

test_that("maihda_normalize_family_name canonicalises engine-specific labels", {
  expect_identical(maihda_normalize_family_name("Negative Binomial(2.34)"),
                   "negbinomial")
  expect_identical(maihda_normalize_family_name("Negative Binomial(0.7213)"),
                   "negbinomial")
  expect_identical(maihda_normalize_family_name("ordinal"), "cumulative")

  # Identity for everything already canonical (or unknown).
  expect_identical(maihda_normalize_family_name("gaussian"), "gaussian")
  expect_identical(maihda_normalize_family_name("binomial"), "binomial")
  expect_identical(maihda_normalize_family_name("poisson"), "poisson")
  expect_identical(maihda_normalize_family_name("negbinomial"), "negbinomial")
  expect_identical(maihda_normalize_family_name("Gamma"), "Gamma")

  # Degenerate inputs pass through untouched.
  expect_null(maihda_normalize_family_name(NULL))
  expect_identical(maihda_normalize_family_name(NA_character_), NA_character_)
  expect_identical(maihda_normalize_family_name(1), 1)
})

test_that("maihda_model_family_key normalizes and falls back to the stored family", {
  # stats::family() is undefined for the placeholder fit object, so the key
  # must come from the wrapper-recorded family. An ESTIMATED theta (the canonical
  # "negbinomial" marker glmer.nb/brms record) normalizes to "negbinomial(log)".
  m <- structure(
    list(model = structure(list(), class = "no_such_fit"),
         family = list(family = "negbinomial", link = "log")),
    class = "maihda_model"
  )
  expect_identical(maihda_model_family_key(m), "negbinomial(log)")

  # A FIXED, user-specified theta (the MASS "Negative Binomial(<theta>)" label
  # the wrapper stores for a glmer(family = MASS::negative.binomial(theta)) fit)
  # stays in the key: two different fixed thetas are different specifications.
  m$family <- list(family = "Negative Binomial(2.5)", link = "log")
  expect_identical(maihda_model_family_key(m), "Negative Binomial(2.5)(log)")
  m$family <- list(family = "Negative Binomial(10)", link = "log")
  expect_identical(maihda_model_family_key(m), "Negative Binomial(10)(log)")

  m$family <- NULL
  expect_identical(maihda_model_family_key(m), "NA(NA)")
})

test_that("calculate_pcv rejects two fixed-theta NB fits with different thetas", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  d <- make_nb_data()
  m1 <- suppressWarnings(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               family = MASS::negative.binomial(1)))
  m10 <- suppressWarnings(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               family = MASS::negative.binomial(10)))

  # Different fixed dispersion assumptions are not comparable: the family/link
  # check must fire on the theta-bearing key rather than silently returning a PCV.
  expect_error(calculate_pcv(m1, m10), "same model family")

  # The SAME fixed theta is comparable (here a null/adjusted pair).
  m1_adj <- suppressWarnings(
    fit_maihda(y ~ age + gender + race + edu + (1 | gender:race:edu), data = d,
               family = MASS::negative.binomial(1)))
  pcv <- calculate_pcv(m1, m1_adj)
  expect_true(is.finite(pcv$pcv))
})

test_that("the family string switch accepts negbinomial and rejects non-log links", {
  d <- make_nb_data(n = 60)
  expect_error(
    fit_maihda(y ~ (1 | gender:race), data = d, family = "negbinomiall"),
    "Unsupported family"
  )
  skip_if_not_installed("MASS")
  expect_error(
    fit_maihda(y ~ (1 | gender:race), data = d,
               family = MASS::negative.binomial(2, link = "identity")),
    "log link"
  )
})

# ---- lme4 glmer.nb fit, variance math, summary -------------------------------

test_that("fit_maihda(family = 'negbinomial') fits via glmer.nb with the canonical marker", {
  skip_on_cran()

  d <- make_nb_data()
  m <- fit_nb(d)

  expect_s4_class(m$model, "glmerMod")
  expect_identical(m$engine, "lme4")
  expect_identical(m$family, list(family = "negbinomial", link = "log"))
  expect_identical(m$diagnostics$engine, "lme4")

  # The raw lme4 label embeds theta; the normalized view is canonical.
  expect_match(stats::family(m$model)$family, "^Negative Binomial\\(")
  expect_identical(maihda_family(m$model)$family, "negbinomial")

  theta <- maihda_negbin_theta_lme4(m$model)
  expect_true(is.finite(theta) && theta > 0)

  # Nakagawa, Johnson & Schielzeth (2017) level-1 variance at the MARGINAL
  # expected count -- the single global mean of lambda_i = exp(x_i'beta + v/2),
  # taken before the log1p transform (see test-audit-2026-07-30.R) -- vs hand
  # computation.
  rv <- maihda_residual_variance_lme4(m$model)
  eta <- as.numeric(stats::predict(m$model, re.form = NA, type = "link"))
  v <- as.numeric(lme4::VarCorr(m$model)$stratum[1, 1])
  lambda <- pmax(exp(eta + v / 2), .Machine$double.eps)
  expect_equal(rv, log1p(1 / mean(lambda) + 1 / theta))
  # The 1/theta overdispersion term makes it strictly larger than the Poisson
  # approximation on the same marginal mean.
  expect_gt(rv, log1p(1 / mean(lambda)))

  s <- suppressWarnings(summary(m))
  expect_true(s$vpc$estimate > 0 && s$vpc$estimate < 1)
  resid_row <- s$variance_components$variance[
    s$variance_components$component == "Within-stratum (residual)"]
  expect_equal(resid_row, rv)
  # Discriminatory accuracy stays binomial-only.
  expect_null(s$discriminatory_accuracy)
})

test_that("a fixed-theta negative.binomial family object fits via glmer and recovers theta", {
  skip_on_cran()
  skip_if_not_installed("MASS")

  d <- make_nb_data()
  m <- suppressWarnings(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               family = MASS::negative.binomial(2))
  )
  expect_s4_class(m$model, "glmerMod")
  # getME() reads theta off the negative.binomial family for any glmerMod, so
  # the user's fixed theta is recovered exactly.
  expect_equal(maihda_negbin_theta_lme4(m$model), 2)

  s <- suppressWarnings(summary(m))
  expect_true(s$vpc$estimate > 0 && s$vpc$estimate < 1)
})

test_that("maihda_negbin_theta_lme4 falls back to parsing the family label", {
  skip_if_not_installed("MASS")

  # A non-merMod fit (plain glm) has no getME(); the helper parses theta out of
  # the "Negative Binomial(<theta>)" label instead.
  d <- data.frame(y = rpois(50, 3), x = rnorm(50))
  g <- suppressWarnings(stats::glm(y ~ x, data = d,
                                   family = MASS::negative.binomial(2.5)))
  expect_equal(maihda_negbin_theta_lme4(g), 2.5)

  # Unrecoverable: no getME, no NB family label.
  g2 <- stats::glm(y ~ x, data = d, family = stats::poisson())
  expect_error(maihda_negbin_theta_lme4(g2), "Could not recover")
})

# ---- PCV / workflow ----------------------------------------------------------

test_that("calculate_pcv works across two glmer.nb fits despite differing thetas", {
  skip_on_cran()

  d <- make_nb_data()
  null_m <- fit_nb(d, y ~ age + (1 | gender:race:edu))
  adj_m <- fit_nb(d, y ~ age + gender + race + edu + (1 | gender:race:edu))

  # The raw labels embed each fit's own theta estimate; pre-normalization this
  # made the family/link equality check fail between a null/adjusted pair.
  pcv <- calculate_pcv(null_m, adj_m)
  expect_true(is.finite(pcv$pcv))
  expect_true(pcv$var_model1 > 0)
})

test_that("calculate_pcv still rejects a negbinomial vs poisson pair", {
  skip_on_cran()

  d <- make_nb_data()
  m_nb <- fit_nb(d, y ~ age + (1 | gender:race:edu))
  m_pois <- suppressWarnings(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d, family = "poisson"))
  expect_error(calculate_pcv(m_nb, m_pois), "same model family")
})

test_that("maihda() runs the two-model negbinomial decomposition end-to-end", {
  skip_on_cran()

  d <- make_nb_data()
  a <- suppressWarnings(
    maihda(y ~ age + gender + race + edu + (1 | gender:race:edu),
           data = d, family = "negbinomial")
  )
  expect_s3_class(a, "maihda_analysis")
  expect_identical(a$mode, "two-model")
  expect_true(is.finite(a$pcv$pcv))
  expect_true(a$summary$vpc$estimate > 0 && a$summary$vpc$estimate < 1)
  expect_output(print(a), "negbinomial")
})

test_that("bootstrap VPC intervals work for a glmer.nb fit (theta held fixed)", {
  skip_on_cran()

  d <- make_nb_data()
  m <- fit_nb(d)
  s <- suppressWarnings(summary(m, bootstrap = TRUE, n_boot = 10))
  expect_true(is.finite(s$vpc$ci_lower))
  expect_true(is.finite(s$vpc$ci_upper))
  expect_true(s$vpc$ci_lower <= s$vpc$ci_upper)
})

# ---- rejections ---------------------------------------------------------------

test_that("negbinomial is rejected where it is not defined", {
  skip_on_cran()

  d <- make_nb_data()
  d$w <- runif(nrow(d), 0.5, 2)

  expect_error(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               family = "negbinomial", engine = "wemix", sampling_weights = "w"),
    "gaussian\\(identity\\) and binomial\\(logit\\)"
  )

  m <- fit_nb(d)
  expect_error(maihda_discriminatory_accuracy(m), "negbinomial")
  expect_error(maihda_mor(m), "negbinomial")
  expect_error(maihda_vpc_response(m), "negbinomial")
})

# ---- brms plumbing (Stan-free) -------------------------------------------------

test_that("fit_maihda routes family = 'negbinomial' to brms::negbinomial()", {
  skip_if_not_installed("brms")

  d <- make_nb_data(n = 80)
  captured <- NULL
  local_mocked_bindings(
    brm = function(formula, data, family, ...) {
      captured <<- list(formula = formula, data = data, family = family)
      structure(list(), class = "brmsfit")
    },
    .package = "brms"
  )

  m <- fit_maihda(y ~ age + (1 | gender:race), data = d,
                  family = "negbinomial", engine = "brms")
  expect_identical(captured$family$family, "negbinomial")
  expect_identical(captured$family$link, "log")
  expect_identical(maihda_normalize_family_name(m$family$family), "negbinomial")
})

test_that("a fixed-theta family object is rejected for brms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("MASS")

  d <- make_nb_data(n = 80)
  expect_error(
    fit_maihda(y ~ age + (1 | gender:race), data = d,
               family = MASS::negative.binomial(2), engine = "brms"),
    "only\\s+supported by engine"
  )
})

test_that("brms negbinomial summary returns a draws-based VPC", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true). The shape /
  # residual-variance arithmetic is covered Stan-free elsewhere.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- make_nb_data(n = 600)
  # Same defect the 2026-08-31 audit found in the ordinal twin of this block: at
  # chains = 2, iter = 500 the chains had not mixed, and "vpc in (0, 1)" plus
  # "the interval is finite" pass on an unmixed posterior just as happily as on a
  # converged one. iter = 2000 (1000 warmup) with adapt_delta above the default
  # clears both halves of the package's bar; the seed keeps it reproducible.
  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ age + (1 | gender:race:edu), data = d,
               family = "negbinomial", engine = "brms",
               chains = 2, iter = 2000, warmup = 1000, refresh = 0, seed = 1,
               control = list(adapt_delta = 0.95))))
  # Fixture guard -- see the note in test-ordinal-engine.R.
  expect_true(isTRUE(m$diagnostics$converged))

  # summary() DELIBERATELY warns on this fixture and always has: mean(y) is 3.3,
  # but the effective count after the negative-binomial overdispersion is 1.15 --
  # below the 2 above which Nakagawa, Johnson & Schielzeth's (2017) lognormal,
  # delta and trigamma approximations agree. Here they give 0.626 / 0.870 / 1.347,
  # a factor of two apart, so the low-count guard added by the 2026-08-30 pass is
  # right to fire. It fires identically on the old chains = 2, iter = 500 fit, so
  # it is not an artefact of this retune -- but the block used to let it escape
  # unasserted, which is the same "a warning nobody owns" pattern this audit is
  # about. Pin it instead: suppressing it would hide the guard regressing.
  expect_warning(s <- summary(m), "approximations agree only above 2")
  expect_true(s$vpc$estimate > 0 && s$vpc$estimate < 1)
  expect_true(is.finite(s$vpc$ci_lower) && is.finite(s$vpc$ci_upper))
  # Internal consistency the block was missing. No recovery assertion on the VPC
  # itself: a count VPC's denominator is the level-1 approximation evaluated at
  # the package's own marginal lambda convention, so pinning a hand-computed
  # target here would test the package against its own formula rather than
  # against the data-generating truth.
  expect_true(s$vpc$estimate >= s$vpc$ci_lower && s$vpc$estimate <= s$vpc$ci_upper)
})

# ---- plots --------------------------------------------------------------------

test_that("the plot layer treats a glmer.nb fit as a count model", {
  skip_on_cran()

  d <- make_nb_data()
  m <- fit_nb(d)
  s <- suppressWarnings(summary(m))

  # The normalized family name routes the deviation panels to the count branch
  # (pre-normalization the raw label matched nothing).
  expect_identical(maihda_prediction_panel_auto_type(m$model), "poisson")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_s3_class(plot(m, type = "vpc", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "predicted", summary_obj = s), "ggplot")
  expect_s3_class(plot(m, type = "effect_decomp", summary_obj = s), "ggplot")
  expect_no_error(suppressWarnings(plot_prediction_deviation_panels(m)))
})
