# Audit 2026-09-25 -- na.exclude padding has to be gone BEFORE lme4 simulates.
#
# FINDING (external audit, F03). maihda_simulate_lme4() called stats::simulate() and
# then stripped the na.exclude padding from the RESULT. But lme4:::.simulateFun() also
# READS the fit before it draws anything:
#
#   weights <- weights(object)                    # padded under na.exclude
#   ...
#   val <- sfun(object, nsim = 1, ftd = ..., wts = weights)
#
# weights.merMod() pads through napredict(), so `wts` is the ORIGINAL-row vector --
# longer than the n fitted values it is paired with, and NA-carrying. Every GLMM
# simulator in lme4:::simfunList opens with a scalar test on it:
#
#   binomial                    if (any(wts %% 1 != 0)) stop(...)
#   poisson, negative.binomial  if (any(wts != 1))      warning(...)
#
# For an UNWEIGHTED fit those reduce to any(c(FALSE, ..., NA)) == NA, so the call dies
# inside lme4 with "missing value where TRUE/FALSE needed" BEFORE any draw exists --
# and a result that never came back cannot be post-processed.
#
# any() SHORT-CIRCUITS, which decides exactly who breaks (all measured on the pre-fix
# tree): a WEIGHTED Poisson or negative binomial survives, because a real weight != 1
# returns TRUE before the NA is reached, whereas the UNWEIGHTED ones die. The binomial
# has no such escape -- integer weights are all FALSE under `%% 1 != 0`, so the NA
# decides it whether the fit is weighted or not.
#
# So on an na.exclude fit, every bootstrap that simulates from the FITTED model --
# summary(bootstrap = TRUE), calculate_pcv(), the crossed-dimensions VPC, the
# longitudinal VPC(t) band -- was dead outright for ANY binomial and for an unweighted
# Poisson or negative binomial. A cbind() binomial would additionally have drawn
# rbinom(n, size = wts) from padded trial counts even without the NA.
#
# df_method = "bootstrap" and pcv_importance() were NOT affected: both simulate from a
# complete-case refit that records no na.action, so there was no padding to read. Both
# were measured, here or on HEAD, rather than assumed.
#
# Gaussian was immune: .simulateFun()'s LMM branch never touches `weights`. That is
# exactly why the 2026-09-20 pass, which exercised the bootstrap on a Gaussian fit,
# did not see this.
#
# CONFIRMED. The fix simulates from maihda_lme4_unpadded_fit(model) -- the same fit
# with the model frame's "na.action" re-classed from "exclude" to "omit", which is the
# one bit that makes napredict()/naresid() pad rather than act as the identity. The
# draws then come back identical(), row names and attributes included, to the SAME fit
# under na.action = na.omit. The unpadding of the result is kept as the standing
# guarantee it always was; it is simply a no-op now.

audit0925_data <- function(n = 480, seed = 7) {
  set.seed(seed)
  d <- data.frame(sex = factor(sample(c("F", "M"), n, TRUE)),
                  edu = factor(sample(c("low", "mid", "high"), n, TRUE)),
                  age = factor(sample(c("y", "o"), n, TRUE)),
                  x   = rnorm(n))
  st <- interaction(d$sex, d$edu, d$age, sep = "_", drop = TRUE)
  u  <- rnorm(nlevels(st), 0, 0.7)[as.integer(st)]
  d$y_g <- 1 + 0.5 * d$x + u + stats::rnorm(n)
  d$y_b <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * d$x + u))
  d$y_p <- stats::rpois(n, exp(0.6 + 0.3 * d$x + u))
  d$y_n <- stats::rnbinom(n, mu = exp(0.6 + 0.3 * d$x + u), size = 1.6)
  d$trials <- sample(5:12, n, TRUE)
  d$succ <- stats::rbinom(n, d$trials, stats::plogis(-0.2 + 0.5 * d$x + u))
  d$fail <- d$trials - d$succ
  d$pw <- sample(c(1, 2, 3), n, TRUE)
  d$x[c(5, 100, 300)] <- NA          # 3 dropped rows, 477 analytic
  d
}

# The families that reach lme4's GLMM branch, plus the Gaussian control. `w = TRUE`
# adds integer prior weights, and on the pre-fix tree that CHANGED the outcome for
# Poisson: its any(wts != 1) short-circuits to TRUE on a real weight before it reaches
# the NA, so `pois_w` survived where `poisson` died. The binomial's modulo test has no
# such escape (integer weights are all FALSE), so `binom_w` died like `binomial`.
audit0925_specs <- function() {
  list(
    gaussian = list(f = y_g ~ x + (1 | sex:edu:age), fam = "gaussian", w = FALSE),
    binomial = list(f = y_b ~ x + (1 | sex:edu:age), fam = stats::binomial(), w = FALSE),
    poisson  = list(f = y_p ~ x + (1 | sex:edu:age), fam = stats::poisson(), w = FALSE),
    binom_cb = list(f = cbind(succ, fail) ~ x + (1 | sex:edu:age),
                    fam = stats::binomial(), w = FALSE),
    gauss_w  = list(f = y_g ~ x + (1 | sex:edu:age), fam = "gaussian", w = TRUE),
    pois_w   = list(f = y_p ~ x + (1 | sex:edu:age), fam = stats::poisson(), w = TRUE),
    binom_w  = list(f = y_b ~ x + (1 | sex:edu:age), fam = stats::binomial(), w = TRUE)
  )
}

# Fitting seven families twice is not cheap, so build once and reuse across blocks.
audit0925_cache <- new.env(parent = emptyenv())

audit0925_fits <- function() {
  if (!is.null(audit0925_cache$fits)) {
    return(audit0925_cache$fits)
  }
  d <- audit0925_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  out <- lapply(audit0925_specs(), function(s) {
    mk <- function(na) {
      if (s$w) q(fit_maihda(s$f, data = d, family = s$fam, weights = pw, na.action = na))
      else     q(fit_maihda(s$f, data = d, family = s$fam, na.action = na))
    }
    list(exclude = mk(stats::na.exclude), omit = mk(stats::na.omit))
  })
  audit0925_cache$fits <- out
  out
}

audit0925_draws <- function(m, nsim = 3, seed = 99) {
  set.seed(seed)
  s <- suppressWarnings(suppressMessages(MAIHDA:::maihda_simulate_lme4(m, nsim = nsim)))
  # Columns may be plain vectors or 2-column matrices (a cbind() binomial); flatten to
  # one comparable numeric block either way.
  unname(as.matrix(as.data.frame(lapply(s, as.vector))))
}


# -- the premise -------------------------------------------------------------

test_that("an na.exclude fit hands lme4 a padded, NA-carrying prior-weight vector", {
  skip_on_cran()
  f <- audit0925_fits()$poisson
  wx <- stats::weights(f$exclude$model, type = "prior")
  wo <- stats::weights(f$omit$model, type = "prior")

  expect_identical(nrow(f$exclude$model@frame), 477L)
  expect_identical(nrow(f$omit$model@frame), 477L)
  expect_identical(length(wx), 480L)              # padded back to the input rows
  expect_identical(sum(is.na(wx)), 3L)
  expect_identical(length(wo), 477L)              # na.omit does not pad
  expect_identical(sum(is.na(wo)), 0L)
  # The package's own accessor already stripped it -- the padding only escaped through
  # accessors read INSIDE lme4.
  expect_identical(length(MAIHDA:::maihda_fit_rows_weights(f$exclude$model)), 477L)

  # And this is the arithmetic that aborts lme4's simulators on an UNWEIGHTED fit.
  expect_true(is.na(any(wx %% 1 != 0)))           # binomial's test
  expect_true(is.na(any(wx != 1)))                # poisson / negbin's test
  # ... and an NA condition is what aborts lme4. Caught rather than matched: the
  # message is translated, so its text cannot be asserted on.
  expect_s3_class(tryCatch(if (any(wx != 1)) TRUE, error = function(e) e), "error")

  # The asymmetry that decides the blast radius, pinned as arithmetic: any() SHORT-
  # CIRCUITS, so a real weight != 1 makes the Poisson/negbin test TRUE before it ever
  # reaches the padded NA -- which is why a WEIGHTED Poisson survived on the pre-fix
  # tree while the unweighted one died. The binomial's modulo test has no such escape:
  # integer weights are all FALSE, so the NA decides it either way.
  ww <- stats::weights(audit0925_fits()$pois_w$exclude$model, type = "prior")
  expect_identical(length(ww), 480L)
  expect_identical(sum(is.na(ww)), 3L)
  expect_false(all(stats::na.omit(ww) == 1))      # the fixture really is weighted
  expect_true(any(ww != 1))                       # TRUE, not NA -- the escape
  expect_true(is.na(any(ww %% 1 != 0)))           # no escape for the binomial test

  # The two fits really are the same fit, so nothing below can be an estimation
  # difference dressed up as an NA-handling one.
  expect_equal(lme4::fixef(f$exclude$model), lme4::fixef(f$omit$model))
  expect_identical(nrow(f$exclude$data), nrow(f$omit$data))
})


# -- the helper --------------------------------------------------------------

test_that("maihda_lme4_unpadded_fit() unpads, and is a no-op otherwise", {
  skip_on_cran()
  f <- audit0925_fits()$binomial
  u <- MAIHDA:::maihda_lme4_unpadded_fit(f$exclude$model)

  # Every accessor lme4 reads off the returned fit is now on the analytic rows.
  expect_identical(length(stats::weights(u, type = "prior")), 477L)
  expect_identical(sum(is.na(stats::weights(u, type = "prior"))), 0L)
  expect_identical(length(stats::fitted(u)), 477L)
  expect_identical(length(stats::predict(u)), 477L)
  expect_identical(nrow(stats::model.frame(u)), 477L)
  # The rows themselves are untouched -- only the na.action's class changed.
  expect_identical(class(attr(u@frame, "na.action")), "omit")
  expect_equal(unname(as.integer(attr(u@frame, "na.action"))),
               unname(as.integer(attr(f$exclude$model@frame, "na.action"))))
  expect_equal(lme4::fixef(u), lme4::fixef(f$exclude$model))
  expect_equal(unname(stats::sigma(u)), unname(stats::sigma(f$exclude$model)))
  expect_null(MAIHDA:::maihda_na_exclude_rows(u))

  # A COPY was modified: the caller's fit keeps its own na.exclude semantics.
  expect_identical(class(attr(f$exclude$model@frame, "na.action")), "exclude")
  expect_identical(length(stats::weights(f$exclude$model, type = "prior")), 480L)

  # No-op on a fit that does not pad: na.omit, and no missing data at all.
  expect_identical(MAIHDA:::maihda_lme4_unpadded_fit(f$omit$model), f$omit$model)
  d <- audit0925_data()
  d$x[is.na(d$x)] <- 0
  k <- suppressWarnings(suppressMessages(
    fit_maihda(y_g ~ x + (1 | sex:edu:age), data = d, na.action = stats::na.exclude)))
  expect_null(attr(k$model@frame, "na.action"))
  expect_identical(MAIHDA:::maihda_lme4_unpadded_fit(k$model), k$model)
})


# -- the defect --------------------------------------------------------------

test_that("maihda_simulate_lme4() draws from every family under na.exclude", {
  skip_on_cran()
  # THE BUG: `binomial`, `poisson`, `binom_cb` and `binom_w` stopped inside lme4 with
  # "missing value where TRUE/FALSE needed" -- no draw, no bootstrap. The message is
  # translated, so the assertion is structural: the draws have to exist and be on the
  # analytic rows.
  #
  # `gaussian`, `gauss_w` and `pois_w` are CONTROLS here -- measured OK on the pre-fix
  # tree. Gaussian never reaches the weight test at all, and a WEIGHTED Poisson escapes
  # it because any(wts != 1) short-circuits to TRUE on a real weight before reaching the
  # padded NA (pinned directly in the premise block). A weighted BINOMIAL still dies:
  # its any(wts %% 1 != 0) sees only FALSEs and the NA.
  fits <- audit0925_fits()
  for (nm in names(fits)) {
    sim <- suppressWarnings(suppressMessages(
      MAIHDA:::maihda_simulate_lme4(fits[[nm]]$exclude$model, nsim = 3)))
    expect_identical(NROW(sim), 477L, info = nm)
    expect_identical(length(sim), 3L, info = nm)
    expect_false(anyNA(unlist(sim, use.names = FALSE)), info = nm)
  }
})

test_that("the negative-binomial simulator is reached under na.exclude too", {
  skip_on_cran()
  # glmer.nb() is slow, so it sits in its own block rather than in the family loop.
  d <- audit0925_data()
  q <- function(e) suppressWarnings(suppressMessages(e))
  f <- lapply(list(exclude = stats::na.exclude, omit = stats::na.omit), function(na)
    q(fit_maihda(y_n ~ x + (1 | sex:edu:age), data = d, family = "negbinomial",
                 na.action = na)))
  expect_equal(lme4::fixef(f$exclude$model), lme4::fixef(f$omit$model))
  sim <- q(MAIHDA:::maihda_simulate_lme4(f$exclude$model, nsim = 3))
  expect_identical(NROW(sim), 477L)
  expect_equal(audit0925_draws(f$exclude$model), audit0925_draws(f$omit$model))
})

test_that("na.exclude draws are identical to the same fit's na.omit draws", {
  skip_on_cran()
  # Not "agree to within Monte-Carlo error": the same seed has to consume the same RNG
  # stream and produce the same numbers, because it is the same fit on the same rows.
  fits <- audit0925_fits()
  for (nm in names(fits)) {
    expect_identical(audit0925_draws(fits[[nm]]$exclude$model),
                     audit0925_draws(fits[[nm]]$omit$model), info = nm)
  }

  # And the whole simulate() object matches -- values, row names and attributes -- so
  # nothing downstream can tell the two NA actions apart.
  m <- fits$poisson
  strip <- function(s) { attr(s, "seed") <- NULL; s }
  set.seed(99)
  a <- strip(suppressWarnings(stats::simulate(
    MAIHDA:::maihda_lme4_unpadded_fit(m$exclude$model), nsim = 2)))
  set.seed(99)
  b <- strip(suppressWarnings(stats::simulate(m$omit$model, nsim = 2)))
  expect_identical(a, b)
})


# -- the public entry points -------------------------------------------------

test_that("summary(bootstrap = TRUE) agrees under na.exclude for every family", {
  skip_on_cran()
  fits <- audit0925_fits()
  for (nm in names(fits)) {
    set.seed(3)
    be <- suppressWarnings(summary(fits[[nm]]$exclude, bootstrap = TRUE, n_boot = 12))
    set.seed(3)
    bo <- suppressWarnings(summary(fits[[nm]]$omit, bootstrap = TRUE, n_boot = 12))
    expect_true(is.finite(be$vpc$ci_lower) && is.finite(be$vpc$ci_upper), info = nm)
    expect_equal(be$vpc, bo$vpc, info = nm)
    # Guarded, so the comparison cannot pass on two NULLs. ($decomposition is NULL on
    # a plain fit_maihda() summary -- it belongs to the two-model maihda() analysis --
    # so it is not compared here.)
    expect_s3_class(be$variance_components, "data.frame")
    expect_equal(be$variance_components, bo$variance_components, info = nm)
  }
})

test_that("df_method = 'bootstrap' agrees under na.exclude for every family", {
  skip_on_cran()
  # A CONTROL, not a tooth: this block passes on the pre-fix tree too. The draws come
  # from the restricted NULL fit, not the fitted model -- on y ~ x that is
  # maihda_refit_reduced()'s refit of the reduced formula on the fit's own analytic
  # frame (its first attempt, update() on the original call, re-admits the three NA
  # rows once x is dropped and is rejected for having 480 rows). That refit records no
  # na.action, so lme4 never saw a padded weight vector here and the path already
  # worked. Kept to show the fix does not perturb a bootstrap that was already correct.
  fits <- audit0925_fits()

  # Pin the reason, not just the outcome.
  m <- fits$binomial$exclude$model
  red <- MAIHDA:::maihda_refit_reduced(m, stats::update(stats::formula(m), . ~ . - x))
  expect_false(is.null(red))
  expect_identical(stats::nobs(red), stats::nobs(m))     # 477: the analytic rows
  expect_null(attr(red@frame, "na.action"))              # nothing to pad from
  for (nm in c("binomial", "poisson", "binom_cb", "pois_w")) {
    set.seed(5)
    pe <- suppressWarnings(summary(fits[[nm]]$exclude, df_method = "bootstrap",
                                   n_boot = 10))
    set.seed(5)
    po <- suppressWarnings(summary(fits[[nm]]$omit, df_method = "bootstrap",
                                   n_boot = 10))
    # The intercept has no null model to simulate from and is documented NA; every
    # other coefficient has to carry a real bootstrap p-value.
    slope <- pe$fixed_effects$term != "(Intercept)"
    expect_true(any(slope), info = nm)
    expect_true(all(is.finite(pe$fixed_effects$p_value[slope])), info = nm)
    expect_equal(pe$fixed_effects, po$fixed_effects, info = nm)
  }
})

test_that("calculate_pcv(bootstrap = TRUE) agrees under na.exclude", {
  skip_on_cran()
  # The NA goes in the RESPONSE so the null and the adjusted model drop the same rows
  # (calculate_pcv() refuses two different analytic samples, and rightly so).
  d <- audit0925_data()
  d$x[is.na(d$x)] <- 0
  d$y_p[c(5, 100, 300)] <- NA
  q <- function(e) suppressWarnings(suppressMessages(e))
  pcv <- lapply(list(exclude = stats::na.exclude, omit = stats::na.omit), function(na) {
    m0 <- q(fit_maihda(y_p ~ (1 | sex:edu:age), data = d, family = stats::poisson(),
                       na.action = na))
    m1 <- q(fit_maihda(y_p ~ x + (1 | sex:edu:age), data = d, family = stats::poisson(),
                       na.action = na))
    set.seed(4)
    q(calculate_pcv(m0, m1, bootstrap = TRUE, n_boot = 12))
  })
  expect_true(is.finite(pcv$exclude$ci_lower) && is.finite(pcv$exclude$ci_upper))
  expect_equal(pcv$exclude$pcv, pcv$omit$pcv)
  expect_equal(pcv$exclude$ci_lower, pcv$omit$ci_lower)
  expect_equal(pcv$exclude$ci_upper, pcv$omit$ci_upper)
})

test_that("the longitudinal count VPC(t) band survives na.exclude", {
  skip_on_cran()
  set.seed(31)
  nid <- 220; tt <- 4
  ld <- data.frame(id = rep(seq_len(nid), each = tt),
                   time = rep(0:(tt - 1), times = nid),
                   sex = factor(rep(sample(c("F", "M"), nid, TRUE), each = tt)),
                   edu = factor(rep(sample(c("low", "mid", "high"), nid, TRUE),
                                    each = tt)))
  lst <- interaction(ld$sex, ld$edu, sep = "_", drop = TRUE)
  ld$cnt <- stats::rpois(nrow(ld),
                         exp(1.1 + 0.15 * ld$time +
                               rnorm(nlevels(lst), 0, 0.4)[as.integer(lst)] +
                               rnorm(nid, 0, 0.3)[ld$id]))
  ld$cov <- stats::rnorm(nrow(ld))
  ld$cov[c(7, 300, 700)] <- NA

  q <- function(e) suppressWarnings(suppressMessages(e))
  s <- lapply(list(exclude = stats::na.exclude, omit = stats::na.omit), function(na) {
    m <- q(fit_maihda(cnt ~ cov + (1 | sex:edu), data = ld, family = stats::poisson(),
                      id = "id", time = "time", na.action = na))
    set.seed(8)
    q(summary(m, bootstrap = TRUE, n_boot = 10))   # the documented minimum
  })
  # The band lives at $longitudinal$vpc_t (time / estimate / lower / upper). Asserting
  # it is populated first, so the comparison below cannot pass on two NULLs.
  vt <- lapply(s, function(x) x$longitudinal$vpc_t)
  expect_s3_class(vt$exclude, "data.frame")
  expect_true(nrow(vt$exclude) > 0L)
  expect_true(any(is.finite(vt$exclude$lower)))
  expect_true(any(is.finite(vt$exclude$upper)))
  expect_equal(vt$exclude, vt$omit)
  expect_length(s$exclude$longitudinal$vpc_intercept_ci, 2L)
  expect_equal(s$exclude$longitudinal$vpc_intercept_ci,
               s$omit$longitudinal$vpc_intercept_ci)
  expect_true(is.finite(s$exclude$vpc$estimate))
  expect_equal(s$exclude$vpc, s$omit$vpc)
})
