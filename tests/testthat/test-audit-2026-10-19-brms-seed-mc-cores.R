# Audit 2026-10-19 (brms-seed-mc-cores): set.seed() before an engine = "brms" fit
# reproduces it whether the chains run in this R session or in parallel (cores > 1).
#
# brms passes an unseeded fit's seed = NA to rstan, which draws the seed from R's
# random-number stream in the process that runs the chains. Serially that is this
# session, so set.seed() pins it; with parallel chains each chain runs in a process
# of its own whose stream starts afresh, so set.seed() alone pins nothing: two maihda(engine = "brms") calls after
# set.seed(11) under options(mc.cores = 2) gave different posteriors. fit_maihda()
# draws the seed itself when the chains run in parallel, from this session's stream,
# the number a serial fit's seed is; serial fits draw nothing extra.
#
# The first block pins the seed fit_maihda() hands brms::brm() (no Stan); the second
# fits twice with parallel chains when MAIHDA_TEST_BRMS=true.

# Evaluate `expr` after set.seed(seed) and put the caller's RNG state back, so this
# file leaves the random stream the later test files read as it found it (see
# setup-seed.R).
bsm_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

bsm_data <- function() {
  d <- data.frame(g1 = rep(c("a", "b"), each = 30),
                  g2 = rep(c("x", "y", "z"), times = 20),
                  x = seq_len(60) / 60)
  d$y <- 0.5 * d$x + as.numeric(factor(paste(d$g1, d$g2))) / 6 +
    rep(c(-0.3, 0.1, 0.2), 20)
  d
}

test_that("fit_maihda() hands brm() a seed from R's stream when the chains run in parallel", {
  skip_if_not_installed("brms")

  d <- bsm_data()
  captured <- NULL
  # seed, cores and future come before `...`, as in brms::brm(), so R binds a partial
  # spelling to them as it does at the real call; `extra` is whatever stays unbound.
  local_mocked_bindings(
    brm = function(formula, data, family, seed = NULL, cores = NULL, future = NULL, ...) {
      captured <<- list(seed = seed, extra = names(list(...)))
      stop("brm() reached")
    },
    .package = "brms"
  )
  # The seed brm() receives, the arguments it could not bind, and R's random-number
  # state right after the call.
  seed_sent <- function(...) {
    captured <<- NULL
    expect_error(fit_maihda(y ~ x + (1 | g1:g2), data = d, engine = "brms", ...),
                 "brm() reached", fixed = TRUE)
    list(seed = captured$seed, extra = captured$extra,
         state = get(".Random.seed", envir = globalenv()))
  }
  with_options <- function(opts, expr) {
    old <- options(opts)
    tryCatch(expr, finally = options(old))
  }
  next_draw <- bsm_with_seed(42, sample.int(.Machine$integer.max, 1L))
  state_after_draw <- bsm_with_seed(42, {
    sample.int(.Machine$integer.max, 1L)
    get(".Random.seed", envir = globalenv())
  })
  state_untouched <- bsm_with_seed(42, get(".Random.seed", envir = globalenv()))

  # Parallel chains and no seed: the seed is the stream's next number, taken here.
  s <- bsm_with_seed(42, seed_sent(cores = 2))
  expect_identical(s$seed, next_draw)
  expect_length(s$extra, 0)
  # identical() rather than expect_identical(): testthat's diff of two .Random.seed
  # vectors overflows integer arithmetic and errors instead of failing.
  expect_true(identical(s$state, state_after_draw))

  # options(mc.cores), brms's default for cores, counts the same, and an explicit
  # cores wins over it either way; seed = NA is no seed.
  expect_identical(with_options(list(mc.cores = 2), bsm_with_seed(42, seed_sent())$seed),
                   next_draw)
  expect_identical(with_options(list(mc.cores = 1),
                                bsm_with_seed(42, seed_sent(cores = 2))$seed), next_draw)
  expect_null(with_options(list(mc.cores = 2), bsm_with_seed(42, seed_sent(cores = 1))$seed))
  expect_identical(bsm_with_seed(42, seed_sent(cores = 2, seed = NA))$seed, next_draw)

  # Partial spellings count as brm() binds them: cor = 2 is cores; se = NA is no seed,
  # and the drawn one takes its place; se = 7 is the user's seed, passed on with
  # nothing drawn.
  expect_identical(bsm_with_seed(42, seed_sent(cor = 2))$seed, next_draw)
  s <- bsm_with_seed(42, seed_sent(cores = 2, se = NA))
  expect_identical(s$seed, next_draw)
  expect_length(s$extra, 0)
  s <- bsm_with_seed(42, seed_sent(cores = 2, se = 7))
  expect_identical(s$seed, 7)
  expect_length(s$extra, 0)
  expect_true(identical(s$state, state_untouched))

  # A seed given is passed on as it is, and nothing is drawn.
  s <- bsm_with_seed(42, seed_sent(cores = 2, seed = 7))
  expect_identical(s$seed, 7)
  expect_true(identical(s$state, state_untouched))

  # Chains in this session, or run as futures (which brms seeds from this session):
  # no seed is added and nothing is drawn here.
  for (args in list(list(), list(cores = 1), list(cores = 2, future = TRUE))) {
    s <- bsm_with_seed(42, do.call(seed_sent, args))
    expect_null(s$seed)
    expect_true(identical(s$state, state_untouched))
  }
})

test_that("set.seed() reproduces a brms fit whose chains run in parallel (opt-in Stan)", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # Short chains: the check is that two runs give the same draws, which holds
  # whether or not the chains have converged.
  d <- bsm_data()
  fit <- function() bsm_with_seed(11, suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | g1:g2), data = d, engine = "brms",
               chains = 2, cores = 2, iter = 300, refresh = 0))))
  a <- fit()
  b <- fit()
  expect_identical(as.matrix(a$model), as.matrix(b$model))
  expect_identical(as.numeric(a$model$fit@stan_args[[1]]$seed),
                   as.numeric(bsm_with_seed(11, sample.int(.Machine$integer.max, 1L))))
})
