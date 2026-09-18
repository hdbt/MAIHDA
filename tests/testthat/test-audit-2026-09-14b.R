# Audit 2026-09-14b: a partial spelling of a forwarded engine argument bypassed every
# check that reads it by its exact name -- CONFIRMED, and wider than the ordinal refusal
# it was reported against.
#
# R binds a supplied argument by exact name, then by unique PARTIAL name (to a formal
# before any `...`), then by position. fit_maihda(), maihda() and
# compare_maihda_groups() read subset / weights / offset by exact name -- for family
# detection, strata auto-binning, precision-weight normalisation and the per-engine
# refusals -- and maihda_fit_clmm() reads contrasts for its coding record, while
# ordinal::clmm() (all its formals precede `...`) and lme4's lmer() / glmer() (no `...`)
# bind the partial spellings. On this file's fixture, against the unmodified code:
#   * ordinal: subs = keep reached clmm() as subset past the refusal of 'subset' --
#     clmm fitted 301 of 600 rows while the package stored and predicted all 600, and
#     so did both of maihda()'s fits; weig = pw fitted a weighted model (coefficient
#     1.847282, equal to a direct weighted clmm() fit) that the package treated as
#     unweighted; a partial contrasts name on a design with an aliased column escaped
#     the fit-time coding record, and predictions were 0.155 off in probability.
#   * lme4: subse = keep turned family detection from binomial to gaussian (fit_maihda()
#     and maihda()) and auto-binned 4 strata instead of 6, and offs = off did the same
#     where missing offsets drop rows; compare_maihda_groups() fitted its groups as
#     Gaussian (residual variance 0.245 against the binomial 3.29); and a partially
#     named precision weight skipped the zero-weight normalisation, so an outcome whose
#     third level sat only on zero-weight rows switched to the ordinal engine and failed
#     inside clmm() ("non-conformable arguments") where weights = gave a binomial fit.
#     The same held on the glmer.nb route (strata binned on all rows) and for
#     longitudinal fits, whose time centring used all rows (intercept 0.078 against
#     0.338, on the data in their block).
#   * WeMix (resolved since the 2026-09-14 pass) and brms (unknown names fall into
#     brm()'s `...`, which rstan refuses: "passing unknown arguments") failed loudly.
#
# FIX: maihda_resolve_dot_names(), generalised from that WeMix-only resolver, applies
# R's own binding through match.call() to a call carrying the same argument names --
# and, when one name breaks the whole match (count_approximation travelling through
# maihda() is no lme4 formal), to each name alone. maihda_resolve_engine_dots() renames
# the spellings the engine would bind to subset / weights / offset where each of the
# three entry points evaluates its forwarded arguments, and maihda_fit_clmm() resolves
# against clmm() in full before its coding record. Every other name is left for R to
# bind at the engine call. Renaming two spellings of one argument alike would have let
# fit_maihda() pass only the first to the engine where R used to refuse the call, so
# one of the three supplied more than once is now an error -- which also closes the
# silent use of the first value when a name was simply given twice.

b14_data <- function() {
  set.seed(20260914)
  n <- 600
  d <- data.frame(gender = factor(sample(c("f", "m"), n, TRUE)),
                  race = factor(sample(c("a", "b", "c"), n, TRUE)),
                  x = stats::rnorm(n, 3),
                  grp = factor(sample(c("G1", "G2"), n, TRUE)),
                  f = factor(sample(c("lo", "hi"), n, TRUE), levels = c("lo", "hi")))
  cell <- as.integer(interaction(d$gender, d$race))
  u <- stats::rnorm(6)[cell]
  d$yo <- factor(cut(d$x + 0.8 * (d$f == "hi") + u + stats::rnorm(n), 4, labels = FALSE),
                 ordered = TRUE)
  d$y <- d$x + u + stats::rnorm(n)
  d$keep <- d$x > 3
  # Binary only on the kept rows, so family detection depends on the rows analysed.
  d$yb <- ifelse(d$keep, stats::rbinom(n, 1, 0.4), 2)
  d$pw <- stats::runif(n, 0.5, 2)
  d$x2 <- 2 * d$x
  d
}

b14_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

b14_ordinal <- function(formula, data, ...) {
  b14_quiet(fit_maihda(formula, data = data, engine = "ordinal", family = "ordinal", ...))
}

# clmm's own fitted probability of each observed category against the one implied by
# the package's rebuilt location.
b14_gap <- function(fit) {
  eta <- predict_maihda(fit, scale = "link")
  P <- maihda_ordinal_category_probs(eta, maihda_clmm_cutpoints(fit$model), fit$family$link)
  yi <- as.integer(fit$data$yo)
  max(abs(P[cbind(seq_along(yi), yi)] - fit$model$fitted.values))
}

test_that("the ordinal engine refuses partial spellings of subset and weights", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- b14_data()
  refused <- "Argument(s) not supported by engine = \"ordinal\": "
  for (nm in c("su", "subs", "subse")) {
    expect_error(do.call(b14_ordinal, c(list(yo ~ x + (1 | gender:race), d),
                                        stats::setNames(list(d$keep), nm))),
                 paste0(refused, "subset"), fixed = TRUE, info = nm)
  }
  for (nm in c("w", "weig")) {
    expect_error(do.call(b14_ordinal, c(list(yo ~ x + (1 | gender:race), d),
                                        stats::setNames(list(d$pw), nm))),
                 paste0(refused, "weights"), fixed = TRUE, info = nm)
  }
  expect_error(b14_quiet(maihda(yo ~ x + (1 | gender:race), data = d, engine = "ordinal",
                                family = "ordinal", subs = keep)),
               paste0(refused, "subset"), fixed = TRUE)
})

test_that("a partial contrasts name reaches the ordinal coding record", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- b14_data()
  cm <- matrix(c(-1, 1), ncol = 1, dimnames = list(c("lo", "hi"), "hi"))
  fit <- b14_ordinal(yo ~ f + x + x2 + (1 | gender:race), d, contra = list(f = cm))
  # The aliased x2 leaves clmm without its own contrasts record.
  expect_null(fit$model$contrasts)
  expect_lt(b14_gap(fit), 1e-12)
  fit_sum <- b14_ordinal(yo ~ f + x + x2 + (1 | gender:race), d,
                         contras = list(f = "contr.sum"))
  expect_lt(b14_gap(fit_sum), 1e-12)
})

test_that("lme4 treats a partial subset exactly as the full name", {
  skip_on_cran()
  d <- b14_data()
  exact <- b14_quiet(fit_maihda(yb ~ x + (1 | gender:race), data = d, subset = keep))
  partial <- b14_quiet(fit_maihda(yb ~ x + (1 | gender:race), data = d, subse = keep))
  expect_identical(exact$family$family, "binomial")
  expect_identical(partial$family$family, "binomial")
  expect_identical(lme4::fixef(partial$model), lme4::fixef(exact$model))

  strata_of <- function(fit) sort(unique(as.character(fit$data$stratum)))
  binned_exact <- b14_quiet(fit_maihda(y ~ 1 + (1 | gender:x), data = d, subset = keep))
  binned_partial <- b14_quiet(fit_maihda(y ~ 1 + (1 | gender:x), data = d, subse = keep))
  expect_identical(strata_of(binned_partial), strata_of(binned_exact))

  # lme4 drops the rows whose offset is missing, so family detection must too.
  d$off <- ifelse(d$keep, 0.1, NA_real_)
  offset_exact <- b14_quiet(fit_maihda(yb ~ x + (1 | gender:race), data = d, offset = off))
  offset_partial <- b14_quiet(fit_maihda(yb ~ x + (1 | gender:race), data = d, offs = off))
  expect_identical(offset_exact$family$family, "binomial")
  expect_identical(offset_partial$family$family, "binomial")

  # maihda() evaluates its own forwarded arguments; count_approximation, which it passes
  # on to fit_maihda(), is no lme4 formal and must not hide the partial subset.
  analysis <- b14_quiet(maihda(yb ~ x + (1 | gender:race), data = d, subse = keep,
                               count_approximation = "delta"))
  expect_identical(analysis$model$family$family, "binomial")
})

test_that("glmer.nb and longitudinal fits read a partial subset like the full name", {
  skip_on_cran()
  set.seed(4242)
  m <- 700
  nb <- data.frame(gender = sample(c("F", "M"), m, TRUE),
                   race = sample(c("A", "B", "C"), m, TRUE),
                   age = stats::rnorm(m))
  u <- stats::rnorm(6, sd = 0.4)[as.integer(interaction(nb$gender, nb$race))]
  nb$y <- stats::rnbinom(m, mu = exp(1.2 + 0.1 * nb$age + u), size = 1.5)
  nb$keep <- nb$age > -0.5
  nb$z <- stats::runif(m, 18, 80)
  # z is auto-binned into strata on the analytic rows.
  nb_exact <- b14_quiet(fit_maihda(y ~ age + (1 | gender:z), data = nb,
                                   family = "negbinomial", subset = keep))
  nb_partial <- b14_quiet(fit_maihda(y ~ age + (1 | gender:z), data = nb,
                                     family = "negbinomial", subse = keep))
  expect_identical(lme4::fixef(nb_partial$model), lme4::fixef(nb_exact$model))

  set.seed(5)
  people <- data.frame(pid = 1:150, gender = sample(c("F", "M"), 150, TRUE),
                       edu = sample(c("lo", "hi"), 150, TRUE))
  long <- people[rep(seq_len(nrow(people)), each = 4), ]
  long$wave <- rep(0:3, nrow(people))
  long$y <- stats::rnorm(nrow(long)) + 0.3 * long$wave
  long$keep <- long$wave > 0
  # Time is centred on the analytic rows, so the intercept moves with them.
  long_exact <- b14_quiet(fit_maihda(y ~ wave + (1 | gender:edu), data = long,
                                     id = "pid", time = "wave", subset = keep))
  long_partial <- b14_quiet(fit_maihda(y ~ wave + (1 | gender:edu), data = long,
                                       id = "pid", time = "wave", subse = keep))
  expect_identical(lme4::fixef(long_partial$model), lme4::fixef(long_exact$model))
  expect_identical(long_partial$longitudinal_info, long_exact$longitudinal_info)
})

test_that("subset, weights or offset supplied more than once is refused", {
  skip_on_cran()
  d <- b14_data()
  # Two spellings of one argument met R's own matching error at the engine call;
  # renamed alike they would reach it once, so both are refused up front.
  expect_error(b14_quiet(fit_maihda(yb ~ x + (1 | gender:race), data = d,
                                    subset = keep, subse = !keep)),
               "'subset' is supplied more than once (as 'subset' and 'subse')", fixed = TRUE)
  expect_error(b14_quiet(fit_maihda(yb ~ x + (1 | gender:race), data = d,
                                    subs = keep, subse = !keep)),
               "'subset' is supplied more than once", fixed = TRUE)
  expect_error(b14_quiet(fit_maihda(y ~ x + (1 | gender:race), data = d,
                                    weights = pw, weig = pw)),
               "'weights' is supplied more than once", fixed = TRUE)
  expect_error(b14_quiet(maihda(yb ~ x + (1 | gender:race), data = d,
                                subset = keep, subse = keep)),
               "'subset' is supplied more than once", fixed = TRUE)
  # A name given twice used to reach the engine with its first value only, silently.
  expect_error(b14_quiet(fit_maihda(y ~ x + (1 | gender:race), data = d,
                                    subset = keep, subset = !keep)),
               "'subset' is supplied more than once", fixed = TRUE)
})

test_that("a partially named precision weight is normalised like the full name", {
  skip_on_cran()
  d <- b14_data()
  # A third ordered level only on rows whose weight is 0: normalising those weights to
  # NA leaves a binary outcome, which an unnormalised weight does not.
  third <- which(!d$keep)[1:20]
  d$y3 <- factor(ifelse(d$keep, "lo", "hi"), levels = c("lo", "hi", "top"), ordered = TRUE)
  d$y3[sample(which(d$keep), 150)] <- "hi"
  d$y3[third] <- "top"
  d$pw[third] <- 0
  exact <- b14_quiet(fit_maihda(y3 ~ x + (1 | gender:race), data = d, weights = pw))
  partial <- b14_quiet(fit_maihda(y3 ~ x + (1 | gender:race), data = d, weig = pw))
  expect_identical(exact$family$family, "binomial")
  expect_identical(lme4::fixef(partial$model), lme4::fixef(exact$model))
  analysis <- b14_quiet(maihda(y3 ~ x + (1 | gender:race), data = d, weig = pw))
  expect_identical(analysis$model$family$family, "binomial")
})

test_that("compare_maihda_groups() applies a partial subset exactly as the full name", {
  skip_on_cran()
  d <- b14_data()
  exact <- b14_quiet(compare_maihda_groups(yb ~ x + (1 | gender:race), data = d,
                                           group = "grp", subset = keep))
  partial <- b14_quiet(compare_maihda_groups(yb ~ x + (1 | gender:race), data = d,
                                             group = "grp", subse = keep))
  num <- vapply(exact, is.numeric, NA)
  expect_true(any(num))
  expect_identical(as.data.frame(partial)[, num], as.data.frame(exact)[, num])
})

test_that("forwarded names resolve to what the engine binds, and nothing else", {
  skip_if_not_installed("ordinal")
  clmm_supplied <- c("formula", "data", "link", "Hess")
  full <- maihda_resolve_dot_names(
    list(subs = 1, weig = 2, contra = 3, cont = 4, thr = 5, offs = 6),
    ordinal::clmm, clmm_supplied)
  # cont binds contrasts and control alike; offs binds no formal and falls into `...`.
  expect_identical(names(full), c("subset", "weights", "contrasts", "cont", "threshold", "offs"))
  watched <- maihda_resolve_dot_names(
    list(subs = 1, weig = 2, contra = 3, thr = 5), ordinal::clmm, clmm_supplied,
    watch = c("subset", "weights", "offset"))
  expect_identical(names(watched), c("subset", "weights", "contra", "thr"))
  # lmer() has no `...`: a name it cannot bind breaks the whole match, but not the rest.
  mixed <- maihda_resolve_dot_names(list(count_approximation = "delta", subse = TRUE),
                                    lme4::lmer, c("formula", "data"))
  expect_identical(names(mixed), c("count_approximation", "subset"))

  expect_identical(names(maihda_resolve_engine_dots(list(subs = 1, thr = 2), "ordinal",
                                                    "ordinal")), c("subset", "thr"))
  expect_identical(names(maihda_resolve_engine_dots(list(o = 1), "lme4", "gaussian")),
                   "offset")
  expect_identical(names(maihda_resolve_engine_dots(list(subs = 1), "no-such-engine",
                                                    "gaussian")), "subs")
  if (requireNamespace("brms", quietly = TRUE)) {
    # brm() has no subset formal, so the spelling used to fall into its `...`, where
    # rstan refused it only after compiling; since audit 2026-09-17c a name the engine
    # cannot bind is renamed to the lme4-only argument it abbreviates and refused here.
    expect_identical(names(maihda_resolve_engine_dots(list(subs = 1), "brms", "gaussian")),
                     "subset")
    # A name brm() DOES bind keeps its own binding.
    expect_identical(names(maihda_resolve_engine_dots(list(w = 500), "brms", "gaussian")),
                     "w")
  }
})

test_that("a partial name the checks do not read binds as its full name did", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- b14_data()
  full <- b14_ordinal(yo ~ x + (1 | gender:race), d, threshold = "equidistant")
  partial <- b14_ordinal(yo ~ x + (1 | gender:race), d, thr = "equidistant")
  expect_identical(partial$model$beta, full$model$beta)
  expect_identical(partial$model$alpha, full$model$alpha)
})
