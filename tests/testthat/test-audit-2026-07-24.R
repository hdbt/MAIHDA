# Regression tests for the crossed-dimensions double-fit audit finding.
#
# maihda(decomposition = "crossed-dimensions") builds and fits a SINGLE
# cross-classified model from the resolved strata/family metadata. It used to fit
# the supplied formula in full first and then discard that fit, keeping only its
# metadata -- roughly doubling the fitting cost (worst for brms: an extra compile
# + MCMC). maihda() now takes fit_maihda()'s metadata-only fast path for that
# preliminary pass, so the engine is fit exactly once. The two-model and
# longitudinal decompositions genuinely reuse the preliminary fit as the null or
# adjusted model and are unaffected (still two fits).

# ---- unit: metadata-only fast path matches a full fit's metadata ---------------

test_that("fit_maihda(.metadata_only=TRUE) returns resolved metadata without a fit", {
  skip_on_cran()
  set.seed(101)
  N   <- 240
  d1  <- factor(sample(c("m", "f"), N, replace = TRUE))
  d2  <- factor(sample(c("lo", "mid", "hi"), N, replace = TRUE))
  age <- round(stats::runif(N, 20, 70))
  dat <- data.frame(y = stats::rnorm(N), bin = stats::rbinom(N, 1, 0.4),
                    d1 = d1, d2 = d2, age = age)

  check <- function(formula) {
    full <- suppressMessages(suppressWarnings(fit_maihda(formula, dat)))
    meta <- suppressMessages(suppressWarnings(
      fit_maihda(formula, dat, .metadata_only = TRUE)))
    # The stub is deliberately NOT a fitted model, so a stray summary()/predict()
    # cannot silently treat it as one.
    expect_s3_class(meta, "maihda_model_metadata")
    expect_false(inherits(meta, "maihda_model"))
    expect_null(meta$model)
    expect_true(isTRUE(meta$metadata_only))
    # Every field maihda()'s crossed-dimensions path reads must match the full fit,
    # so the cross-classified model it builds from them is identical (no drift).
    expect_identical(deparse(meta$formula), deparse(full$formula))
    expect_identical(meta$strata_vars, full$strata_vars)
    expect_identical(meta$strata_autobin_info, full$strata_autobin_info)
    expect_identical(meta$original_data, full$original_data)
    expect_identical(meta$family$family, full$family$family)
    expect_identical(meta$family$link, full$family$link)
    expect_identical(meta$strata_info, full$strata_info)
  }

  check(y ~ d1 + d2 + (1 | d1:d2))     # gaussian shorthand
  check(y ~ 1 + (1 | d1:d2))            # null
  check(bin ~ d1 + d2 + (1 | d1:d2))   # binary auto-detected -> binomial (0/1 recode)
  check(y ~ d1 + age + (1 | d1:age))   # numeric dimension auto-binned to tertiles
})

# ---- integration: crossed-dimensions fits the engine exactly once --------------

# Count lme4 engine fits triggered while evaluating `expr` (evaluated lazily, after
# the trace is armed). The tracer is an EXPRESSION, not a function: trace() splices
# it into lmer()'s body, where a local closure name would not resolve -- so it keys
# off an option (base getOption()/options(), resolvable from any frame). on.exit
# restores the option and removes the trace even on error.
maihda_count_lmer_fits <- function(expr) {
  old <- getOption(".maihda_test_lmer_n")
  options(.maihda_test_lmer_n = 0L)
  on.exit(options(.maihda_test_lmer_n = old), add = TRUE)
  suppressMessages(trace(
    what = "lmer", where = asNamespace("lme4"), print = FALSE,
    tracer = quote(options(.maihda_test_lmer_n = getOption(".maihda_test_lmer_n") + 1L))))
  on.exit(suppressMessages(untrace(what = "lmer", where = asNamespace("lme4"))),
          add = TRUE)
  force(expr)
  getOption(".maihda_test_lmer_n")
}

test_that("crossed-dimensions fits the engine once; two-model still fits it twice", {
  skip_on_cran()
  set.seed(202)
  N   <- 300
  d1  <- factor(sample(c("m", "f"), N, replace = TRUE))
  d2  <- factor(sample(c("lo", "mid", "hi"), N, replace = TRUE))
  sid <- droplevels(interaction(d1, d2))
  u   <- stats::rnorm(nlevels(sid), 0, 0.5)[as.integer(sid)]
  y   <- 1 + 0.5 * (d1 == "m") + u + stats::rnorm(N)
  dat <- data.frame(y = y, d1 = d1, d2 = d2)

  # Was 2 (preliminary fit of the supplied formula + the cross-classified fit);
  # the preliminary fit is now metadata-only, so the engine runs once.
  cc_fits <- maihda_count_lmer_fits(
    suppressMessages(suppressWarnings(
      maihda(y ~ d1 + d2 + (1 | d1:d2), data = dat,
             decomposition = "crossed-dimensions", bootstrap = FALSE))))
  expect_identical(cc_fits, 1L)

  # The two-model path reuses the preliminary fit as the adjusted model and fits a
  # fresh null: still two engine fits, unchanged by the fast path.
  tm_fits <- maihda_count_lmer_fits(
    suppressMessages(suppressWarnings(
      maihda(y ~ d1 + d2 + (1 | d1:d2), data = dat, bootstrap = FALSE))))
  expect_identical(tm_fits, 2L)
})

# ---- integration: crossed-dimensions result is still correct -------------------

test_that("crossed-dimensions still returns a valid additive/interaction split", {
  skip_on_cran()
  set.seed(303)
  N   <- 300
  d1  <- factor(sample(c("m", "f"), N, replace = TRUE))
  d2  <- factor(sample(c("lo", "mid", "hi"), N, replace = TRUE))
  sid <- droplevels(interaction(d1, d2))
  u   <- stats::rnorm(nlevels(sid), 0, 0.5)[as.integer(sid)]
  y   <- 1 + 0.5 * (d1 == "m") + u + stats::rnorm(N)
  dat <- data.frame(y = y, d1 = d1, d2 = d2)

  cc <- suppressMessages(suppressWarnings(
    maihda(y ~ d1 + d2 + (1 | d1:d2), data = dat,
           decomposition = "crossed-dimensions", bootstrap = FALSE)))
  expect_identical(cc$mode, "crossed-dimensions")
  expect_false(is.null(cc$model$cc_info))
  dec <- cc$summary$decomposition
  expect_true(all(c("additive_var", "interaction_var") %in% names(dec)))
  expect_true(is.finite(dec$additive_var) && dec$additive_var >= 0)
  expect_true(is.finite(dec$interaction_var) && dec$interaction_var >= 0)
})
