# Audit 2026-10-12 (crossed-brms-few-levels): a crossed-dimensions fit estimates each
# stratum dimension's variance from that dimension's levels alone. With few levels an
# engine = "brms" fit leaves that variance to the prior on its standard deviation, and
# the VPC and the additive share go with it.
#
# Measured against a known truth under brms's default prior, three dimensions and 100
# Bernoulli rows a stratum: 2 x 2 x 2 strata give a VPC of 0.57 for a true 0.07 (lme4
# on the same data 0.07) and an additive share of 0.98 for a true 0.83 to 0.90 -- and
# of 0.91 where the dimensions have no additive effect at all (0.67 at 3 x 3 x 3, 0.31
# at 4 x 4 x 4); 5 x 5 x 5 strata still give a VPC of 0.18 for 0.08. With effects of
# variance about 0.1, the posterior median of a dimension's variance is about 5, 3, 2,
# 1.7 and 1.25 times the variance of its effects with 2, 3, 4, 5 and 10 levels. Under
# normal(0, 0.5), a fifth of the default prior's scale, a two-level dimension's median
# is 0.2 of what it is under the default and a ten-level one's 0.93. adapt_delta = 0.99
# takes the divergent transitions from about a hundred a fit to one or fewer on
# average and leaves the VPC where it was.
#
# The printed decomposition of such a fit therefore names its dimensions with fewer
# than ten levels. The blocks below pin the count, the rule and the words without
# Stan, and repeat them on real fits when MAIHDA_TEST_BRMS=true.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
cfl_with_seed <- function(seed, expr) {
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

# Printed lines without colour codes.
cfl_print <- function(expr) cli::ansi_strip(utils::capture.output(expr))

cfl_old_note <- c(
  "  Note: the additive share is the crossed-dimensions analogue of the PCV but",
  "  a different estimator; interpret the interaction share cautiously.")

# A crossed-dimensions maihda_model shell around a brmsfit stub that carries its own
# posterior draws and the fitted data.
cfl_stub <- function(data, dim_groups) {
  sd_seq <- function(from, to) seq(from, to, length.out = 101)
  draws <- data.frame(b_Intercept = 0, sd_stratum__Intercept = sd_seq(0.05, 0.45))
  for (g in unname(dim_groups)) draws[[paste0("sd_", g, "__Intercept")]] <- sd_seq(0.2, 1.4)
  fam <- list(family = "bernoulli", link = "logit")
  f <- y ~ 1
  structure(
    list(model = structure(list(formula = f, family = fam, draws = draws),
                           class = "brmsfit"),
         engine = "brms", formula = f, family = fam, data = data,
         cc_info = list(dim_groups = dim_groups, interaction_group = "stratum")),
    class = "maihda_model"
  )
}

cfl_mock_brms <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) model$draws,
    .env = env
  )
}

# 2 x 3 strata of 20 Gaussian rows, built without random numbers.
cfl_gaussian_data <- function() {
  g <- expand.grid(d1 = c("a", "b"), d2 = c("p", "q", "r"), stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 20), ]
  rownames(d) <- NULL
  cell <- rep(c(0, 0.5, 0.3, 0.9, 0.6, 1.4), each = 20)
  d$y <- cell + sin(seq_len(nrow(d)) * 12.9898)
  d
}

test_that("a dimension's levels are counted in the fitted data", {
  d <- data.frame(
    sex = c("f", "m", "f", NA, "m", "f"),
    edu = factor(c("lo", "mid", "hi", "lo", "mid", "hi"), levels = c("lo", "mid", "hi", "phd")),
    .maihda_dim_age = c(1L, 2L, 3L, 1L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  groups <- c(sex = "sex", edu = "edu", age = ".maihda_dim_age", gone = "not_a_column")
  lv <- maihda_cc_dim_levels(d, groups)
  # Named by the dimension, not by its grouping column; a missing value is no level;
  # a declared level without rows is no level; an absent column is not counted as 0.
  expect_identical(lv, c(sex = 2L, edu = 3L, age = 3L, gone = NA_integer_))
  expect_identical(maihda_cc_dim_levels(NULL, c(sex = "sex")), c(sex = NA_integer_))
})

test_that("few-level dimensions are named for a brms decomposition only, below ten levels", {
  expect_identical(.maihda_cc_few_levels, 10L)
  lv <- c(a = 2L, b = 9L, c = 10L, d = NA_integer_, e = 25L)
  expect_identical(maihda_cc_few_level_dims(list(method = "posterior", dim_levels = lv)),
                   c(a = 2L, b = 9L))
  # Every dimension with enough levels: nothing to name.
  expect_length(
    maihda_cc_few_level_dims(list(method = "posterior", dim_levels = c(a = 10L, b = 40L))),
    0L)
  # An lme4 decomposition has no prior to depend on, whatever it carries.
  expect_length(maihda_cc_few_level_dims(list(bootstrap = FALSE, dim_levels = lv)), 0L)
  expect_length(maihda_cc_few_level_dims(list(method = "bootstrap", dim_levels = lv)), 0L)
  # A brms decomposition without the counts names none.
  expect_length(maihda_cc_few_level_dims(list(method = "posterior")), 0L)
})

test_that("the note names the dimensions and says what follows", {
  expect_identical(maihda_cc_few_levels_note(integer(0)), character(0))

  one <- maihda_cc_few_levels_note(c(gender = 2L))
  expect_identical(paste(trimws(one), collapse = " "), paste(
    "Note: gender (2 levels) has few levels. With engine = \"brms\" the variance of",
    "such a dimension depends on the prior on its standard deviation (with two levels",
    "roughly in proportion to that prior's scale), and the VPC and the shares follow",
    "it. See ?maihda."))

  three <- maihda_cc_few_levels_note(c(d1 = 2L, d2 = 3L, d3 = 4L))
  expect_identical(paste(trimws(three), collapse = " "), paste(
    "Note: d1 (2 levels), d2 (3 levels) and d3 (4 levels) have few levels. With",
    "engine = \"brms\" the variance of such a dimension depends on the prior on its",
    "standard deviation (with two levels roughly in proportion to that prior's scale),",
    "and the VPC and the shares follow it. See ?maihda."))

  # The clause on two levels is a claim about two-level dimensions only.
  no_two <- maihda_cc_few_levels_note(c(race = 5L, education = 4L))
  expect_identical(paste(trimws(no_two), collapse = " "), paste(
    "Note: race (5 levels) and education (4 levels) have few levels. With",
    "engine = \"brms\" the variance of such a dimension depends on the prior on its",
    "standard deviation, and the VPC and the shares follow it. See ?maihda."))

  # One level reads as one level.
  expect_true(grepl("site (1 level) has few levels.",
                    paste(trimws(maihda_cc_few_levels_note(c(site = 1L))), collapse = " "),
                    fixed = TRUE))

  # Laid out like the note above it: indented by two, within the console width.
  for (note in list(one, three, no_two)) {
    expect_true(all(startsWith(note, "  ")))
    expect_false(any(startsWith(note, "   ")))
    expect_true(all(nchar(note) <= 77L))
  }
})

test_that("a brms crossed-dimensions summary carries the level counts and prints the note", {
  skip_if_not_installed("brms")
  cfl_mock_brms()
  dat <- data.frame(ga = rep(c("f", "m"), times = 12), gb = rep(letters[1:12], each = 2),
                    stringsAsFactors = FALSE)
  dat$stratum <- paste(dat$ga, dat$gb)
  m <- cfl_stub(dat, c(gender = "ga", region = "gb"))

  dc <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)$decomposition
  expect_identical(dc$dim_levels, c(gender = 2L, region = 12L))
  expect_identical(maihda_cc_few_level_dims(dc), c(gender = 2L))

  out <- cfl_print(maihda_print_cc_decomposition(dc))
  i <- which(out == cfl_old_note[2])
  expect_length(i, 1L)
  expect_identical(out[i - 1L], cfl_old_note[1])
  note <- maihda_cc_few_levels_note(c(gender = 2L))
  expect_identical(out[i + seq_along(note)], note)
  expect_identical(out[i + length(note) + 1L], "")
  expect_length(out, i + length(note) + 1L)
  expect_false(any(grepl("region (12 levels)", out, fixed = TRUE)))

  # Twelve and fifteen levels: the decomposition prints as it always did.
  dat2 <- data.frame(ga = rep(letters[1:12], times = 15),
                     gb = rep(LETTERS[1:15], each = 12), stringsAsFactors = FALSE)
  dat2$stratum <- paste(dat2$ga, dat2$gb)
  dc2 <- maihda_cc_summary_brms(cfl_stub(dat2, c(a = "ga", b = "gb")),
                                list(dim_groups = c(a = "ga", b = "gb"),
                                     interaction_group = "stratum"),
                                conf_level = 0.95)$decomposition
  expect_identical(dc2$dim_levels, c(a = 12L, b = 15L))
  out2 <- cfl_print(maihda_print_cc_decomposition(dc2))
  expect_false(any(grepl("few levels", out2, fixed = TRUE)))
  expect_identical(utils::tail(out2, 3L), c(cfl_old_note, ""))
})

test_that("print() of the analysis and of its summary show the note; an lme4 print has none", {
  d <- cfl_gaussian_data()
  a <- cfl_with_seed(1, suppressWarnings(suppressMessages(
    maihda(y ~ 1 + (1 | d1:d2), data = d, decomposition = "crossed-dimensions"))))
  # An lme4 decomposition is what it was: no level counts, no second note.
  expect_identical(names(a$decomposition),
                   c("additive_var", "interaction_var", "between_var", "within_var",
                     "additive_share", "interaction_share", "per_dim", "bootstrap"))
  out_a <- cfl_print(print(a))
  out_s <- cfl_print(print(a$summary))
  expect_false(any(grepl("few levels", c(out_a, out_s), fixed = TRUE)))

  # The counts are read from a real fit's frame under the grouping columns the crossed
  # model uses; an auto-binned numeric dimension is counted in its bins and named by
  # the variable.
  expect_identical(maihda_cc_dim_levels(a$model$data, a$model$cc_info$dim_groups),
                   c(d1 = 2L, d2 = 3L))
  d$age <- 20 + (seq_len(nrow(d)) * 7) %% 41
  ab <- cfl_with_seed(1, suppressWarnings(suppressMessages(
    maihda(y ~ 1 + (1 | d1:age), data = d, decomposition = "crossed-dimensions"))))
  expect_identical(unname(ab$model$cc_info$dim_groups), c("d1", ".maihda_dim_age"))
  expect_identical(maihda_cc_dim_levels(ab$model$data, ab$model$cc_info$dim_groups),
                   c(d1 = 2L, age = 3L))

  # The same analysis with the decomposition of a brms fit whose dimensions have two
  # and three levels: both print methods add the note directly below the old one, and
  # nothing else in either print moves.
  b <- a
  b$decomposition$method <- "posterior"
  b$decomposition$dim_levels <- c(d1 = 2L, d2 = 3L)
  b$summary$decomposition <- b$decomposition
  note <- maihda_cc_few_levels_note(c(d1 = 2L, d2 = 3L))
  for (pair in list(list(cfl_print(print(b)), out_a),
                    list(cfl_print(print(b$summary)), out_s))) {
    out_b <- pair[[1]]
    i <- which(out_b == cfl_old_note[2])
    expect_length(i, 1L)
    expect_identical(out_b[i + seq_along(note)], note)
    expect_identical(out_b[-(i + seq_along(note))], pair[[2]])
  }
})

test_that("the help page gives the number of levels the note uses", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  rd <- gsub("[[:space:]]+", " ",
             paste(readLines(file.path(man, "maihda.Rd"), warn = FALSE), collapse = " "))
  expect_true(grepl(sprintf("fewer than %d levels", .maihda_cc_few_levels), rd, fixed = TRUE))
})

test_that("real brms fits: the level counts, the printed note and a forwarded prior", {
  # Compiles two Stan models, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # 2 x 3 x 4 strata of 40 binomial trials each, additive effects plus a small stratum
  # effect. brms draws its sampler seed from R's stream, so the fits sit inside the
  # seeded block too.
  fits <- cfl_with_seed(20261013, {
    g <- expand.grid(d1 = c("a", "b"), d2 = c("c", "d", "e"), d3 = c("f", "g", "h", "i"),
                     stringsAsFactors = FALSE)
    eta <- -0.9 + c(0, 0.5)[match(g$d1, c("a", "b"))] +
      c(0, 0.3, 0.6)[match(g$d2, c("c", "d", "e"))] +
      c(0, 0.2, 0.4, 0.7)[match(g$d3, c("f", "g", "h", "i"))] + stats::rnorm(24, 0, 0.2)
    g$n <- 40L
    g$s <- stats::rbinom(24, 40, stats::plogis(eta))
    fit <- function(...) {
      suppressWarnings(suppressMessages(
        maihda(s | trials(n) ~ 1 + (1 | d1:d2:d3), data = g, engine = "brms",
               family = "binomial", decomposition = "crossed-dimensions",
               chains = 2, iter = 1000, refresh = 0, ...)))
    }
    list(default = fit(),
         tight = fit(prior = brms::set_prior("normal(0, 0.5)", class = "sd")))
  })
  a <- fits$default

  # The counts are the engine's own: brms's number of levels of each grouping factor.
  ng <- brms::ngrps(a$model$model)
  expect_identical(a$decomposition$dim_levels, c(d1 = 2L, d2 = 3L, d3 = 4L))
  expect_identical(unname(a$decomposition$dim_levels),
                   as.integer(unlist(ng[c("d1", "d2", "d3")])))
  expect_identical(a$summary$decomposition$dim_levels, a$decomposition$dim_levels)

  note <- maihda_cc_few_levels_note(c(d1 = 2L, d2 = 3L, d3 = 4L))
  for (out in list(cfl_print(print(a)), cfl_print(print(a$summary)))) {
    i <- which(out == cfl_old_note[2])
    expect_length(i, 1L)
    expect_identical(out[i + seq_along(note)], note)
  }

  # `prior` reaches brms::brm(): the fit holds it, and its Stan program uses it.
  ps <- as.data.frame(brms::prior_summary(fits$tight$model$model))
  expect_identical(ps$prior[ps$class == "sd" & ps$source == "user"], "normal(0, 0.5)")
  expect_true(grepl("normal_lpdf(sd_1 | 0, 0.5)", brms::stancode(fits$tight$model$model),
                    fixed = TRUE))
  expect_false(grepl("normal_lpdf(sd_1 | 0, 0.5)", brms::stancode(a$model$model),
                     fixed = TRUE))
  # The tighter prior pulls every dimension's variance, and the VPC, down. (Which
  # dimension it pulls most is not asserted: a default-prior fit of 2 x 1000 draws that
  # does not mix can reverse that order.)
  ratio <- fits$tight$decomposition$per_dim / a$decomposition$per_dim
  expect_true(all(ratio < 1))
  expect_lt(fits$tight$summary$vpc$estimate, a$summary$vpc$estimate)
})
