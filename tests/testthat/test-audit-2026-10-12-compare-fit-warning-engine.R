# Audit 2026-10-12 (packet 2026-10-12-compare-fit-warning-engine): the fit-quality
# warning of compare_maihda_groups() -- CONFIRMED. Wording only: every number, row and
# column is as it was, and the same groups are named.
#
# compare_maihda_groups() names, in one warning, the groups whose fit is singular and
# those whose fit did not converge. Three things were wrong with its text.
#   1. It opened "lme4 reported fit problems (" on every engine. The groups are
#      collected from each fit's diagnostics whatever fitted it: a wemix comparison
#      with one group of identical strata (between-stratum variance 6.4e-23), an
#      ordinal one and a brms one were all told that lme4 reported the problem.
#   2. It always ended "A singular fit means the between-stratum variance is at the
#      boundary (often VPC = 0); a non-converged fit may be unreliable." An lme4 fit
#      is singular when ANY variance component is zero: with exactly additive strata
#      a crossed-dimensions comparison named both groups at VPCs of 0.51 and 0.27 (the
#      interaction term is the zero), and a context without variance named both at
#      0.32 and 0.37. And a warning that listed no singular fit explained one.
#   3. On brms the second entry stands for a largest R-hat above 1.01 or a divergent
#      transition (maihda_fit_diagnostics()), which "did not converge" does not say
#      for a fit whose chains mixed.
#
# Per the user: no engine in the opening; each entry is explained only when it is
# listed; one sentence about a singular fit that holds on every engine and
# decomposition; and a brms entry that says what was flagged.
#
# No block draws a random number: the fixtures are built from fixed sequences, and
# the opt-in Stan block puts the caller's seed back.

fq12_observe <- function(expr) {
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res = res, warns = warns)
}

# The fit-quality warnings among `warns`.
fq12_quality <- function(warns) warns[grepl("fit problems", warns, fixed = TRUE)]

fq12_keep_seed <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  force(expr)
}

fq12_open <- "compare_maihda_groups(): fit problems ("
fq12_singular <- paste0("A singular fit has a variance component estimated at zero; ",
                        "where that is the between-stratum variance, var_between and ",
                        "the VPC are 0.")
fq12_nonconv <- "A fit that did not converge may be unreliable."
fq12_brms <- "The posterior summaries of these brms fits may be unreliable."
fq12_close <- " Interpret these groups' VPC/ICC values with caution."

fq12_spread <- seq(-1.45, 1.45, length.out = 30)
fq12_flat <- rep(0, 6)
fq12_s2 <- c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)
fq12_s3 <- c(0.7, -0.6, 0.1, 0.9, -0.8, -0.3)

# One group of gender x ses strata, 30 rows each: a stratum's outcomes are its shift
# plus the same fixed spread, so shifts of zero give identical strata. `w` holds unit
# sampling weights.
fq12_group <- function(country, shift) {
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 30), ]
  d$y <- rep(shift, each = 30) + rep(fq12_spread, times = nrow(g))
  d$country <- country
  d$w <- 1
  rownames(d) <- NULL
  d
}

# The same rows with the six strata as one column and the dimensions dropped: a
# single stratum dimension has no decomposition, so each group is one fit.
fq12_one_column <- function(d) {
  d$stratum <- paste(d$gender, d$ses, sep = "_")
  d[, c("country", "stratum", "y")]
}

test_that("the warning text names no engine and explains only what it lists", {
  one <- paste0(fq12_open, "singular fit: g1, g4). ", fq12_singular, fq12_close)
  two <- paste0(fq12_open, "did not converge: g2). ", fq12_nonconv, fq12_close)
  both <- paste0(fq12_open, "singular fit: g1; did not converge: g1, g2). ",
                 fq12_singular, " ", fq12_nonconv, fq12_close)
  for (engine in c("lme4", "wemix", "ordinal")) {
    expect_identical(
      maihda_group_fit_quality_message(engine, c("g1", "g4"), character(0)), one)
    expect_identical(
      maihda_group_fit_quality_message(engine, character(0), "g2"), two)
    expect_identical(
      maihda_group_fit_quality_message(engine, "g1", c("g1", "g2")), both)
  }

  # A brms comparison lists the groups under what the package checks there.
  flagged <- paste0(fq12_open,
                    "largest R-hat above 1.01 or divergent transitions: g1, g3). ",
                    fq12_brms, fq12_close)
  expect_identical(
    maihda_group_fit_quality_message("brms", character(0), c("g1", "g3")), flagged)

  # Nothing listed, nothing said.
  for (engine in c("lme4", "wemix", "ordinal", "brms")) {
    expect_null(maihda_group_fit_quality_message(engine, character(0), character(0)))
  }

  # No text holds what the old one claimed.
  texts <- c(one, two, both, flagged)
  expect_false(any(grepl("lme4", texts, fixed = TRUE)))
  expect_false(any(grepl("reported", texts, fixed = TRUE)))
  expect_false(any(grepl("is at the boundary", texts, fixed = TRUE)))
  expect_false(grepl("did not converge", flagged, fixed = TRUE))
  # An entry that is not listed is not explained.
  expect_false(grepl("singular", two, fixed = TRUE))
  expect_false(grepl("converge", one, fixed = TRUE))
})

test_that("an lme4 comparison with a flat group gets the engine-free sentence", {
  d <- rbind(fq12_group("g1", fq12_flat), fq12_group("g2", fq12_s2))
  obs <- fq12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country"))
  expect_lt(obs$res$var_between[obs$res$group == "g1"], 1e-8)
  fq <- fq12_quality(obs$warns)
  expect_identical(fq, maihda_group_fit_quality_message("lme4", "g1", character(0)))
  expect_false(grepl("lme4", fq, fixed = TRUE))
  expect_false(grepl("g2", fq, fixed = TRUE))
  # The fit completed, so the status cells are as they were.
  expect_identical(obs$res$status, c("ok", "ok"))
})

test_that("a wemix comparison is not told that lme4 reported its singular fit", {
  skip_on_cran()
  skip_if_not_installed("WeMix")
  d <- rbind(fq12_group("g1", fq12_flat), fq12_group("g2", fq12_s2))
  obs <- fq12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country", engine = "wemix",
                                            sampling_weights = "w"))
  # The premise: WeMix puts the flat group's between-stratum variance at zero.
  expect_lt(obs$res$var_between[obs$res$group == "g1"], 1e-8)
  expect_gt(obs$res$var_between[obs$res$group == "g2"], 0.1)
  fq <- fq12_quality(obs$warns)
  expect_identical(fq, maihda_group_fit_quality_message("wemix", "g1", character(0)))
  expect_false(grepl("lme4", fq, fixed = TRUE))
  expect_identical(obs$res$status, c("ok", "ok"))
})

test_that("an optimiser cut short is listed without the singular-fit sentence", {
  # Neither group is flat.
  d <- fq12_one_column(rbind(fq12_group("g2", fq12_s2), fq12_group("g3", fq12_s3)))
  obs <- fq12_observe(compare_maihda_groups(
    y ~ (1 | stratum), data = d, group = "country",
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 3))))
  fq <- fq12_quality(obs$warns)
  expect_identical(
    fq, maihda_group_fit_quality_message("lme4", character(0), c("g2", "g3")))
  expect_false(grepl("singular", fq, fixed = TRUE))
  # With the optimiser left alone the same comparison raises no such warning.
  clean <- fq12_observe(compare_maihda_groups(y ~ (1 | stratum), data = d,
                                              group = "country"))
  expect_length(fq12_quality(clean$warns), 0L)
})

test_that("the engine that was used chooses the words, group by group", {
  d <- fq12_one_column(rbind(fq12_group("g1", fq12_s2), fq12_group("g2", fq12_s3),
                             fq12_group("g3", rev(fq12_s2))))
  # Each group's fit is an lme4 fit whose diagnostics record is replaced, so the
  # groups listed are chosen here and no sampler runs.
  real_fit <- fit_maihda
  with_record <- function(singular, converged) {
    function(formula, data, ...) {
      m <- real_fit(formula, data, engine = "lme4")
      g <- data$country[1]
      m$diagnostics$singular <- singular[[g]]
      m$diagnostics$converged <- converged[[g]]
      m
    }
  }
  compare <- function(engine) {
    fq12_observe(compare_maihda_groups(y ~ (1 | stratum), data = d,
                                       group = "country", engine = engine))
  }

  # brms: never singular, two groups flagged by the sampler diagnostics.
  local_mocked_bindings(fit_maihda = with_record(
    singular = list(g1 = NA, g2 = NA, g3 = NA),
    converged = list(g1 = FALSE, g2 = TRUE, g3 = FALSE)))
  obs <- compare("brms")
  fq <- fq12_quality(obs$warns)
  expect_identical(
    fq, maihda_group_fit_quality_message("brms", character(0), c("g1", "g3")))
  expect_true(grepl("largest R-hat above 1.01 or divergent transitions: g1, g3",
                    fq, fixed = TRUE))
  expect_false(grepl("lme4", fq, fixed = TRUE))
  expect_false(grepl("did not converge", fq, fixed = TRUE))
  expect_identical(obs$res$status, rep("ok", 3))
  # The same records under lme4 read "did not converge".
  expect_identical(
    fq12_quality(compare("lme4")$warns),
    maihda_group_fit_quality_message("lme4", character(0), c("g1", "g3")))

  # A group can be in both lists.
  local_mocked_bindings(fit_maihda = with_record(
    singular = list(g1 = TRUE, g2 = FALSE, g3 = FALSE),
    converged = list(g1 = FALSE, g2 = FALSE, g3 = TRUE)))
  expect_identical(
    fq12_quality(compare("lme4")$warns),
    maihda_group_fit_quality_message("lme4", "g1", c("g1", "g2")))

  # A record that is unknown (NA) or clean lists nothing.
  local_mocked_bindings(fit_maihda = with_record(
    singular = list(g1 = NA, g2 = FALSE, g3 = NA),
    converged = list(g1 = NA, g2 = TRUE, g3 = NA)))
  expect_length(fq12_quality(compare("brms")$warns), 0L)
  expect_length(fq12_quality(compare("lme4")$warns), 0L)
})

test_that("what the brms entry says is what the diagnostics check", {
  skip_if_not_installed("brms")
  # The record maihda_fit_diagnostics() makes for a brmsfit with these R-hat values
  # and divergence indicators, without fitting anything.
  record <- function(rhat, divergent = 0) {
    testthat::local_mocked_bindings(
      rhat = function(x, ...) rhat,
      nuts_params = function(x, ...) {
        data.frame(Parameter = "divergent__", Value = divergent)
      },
      .package = "brms")
    maihda_fit_diagnostics(structure(list(), class = "brmsfit"))
  }
  # "largest R-hat above 1.01": the largest of all, strictly above.
  expect_true(record(c(a = 1.002, b = 1.01))$converged)
  expect_false(record(c(a = 1.002, b = 1.0101))$converged)
  # "or divergent transitions": one is enough, with the chains mixed.
  expect_false(record(c(a = 1.002), divergent = c(0, 0, 1, 0))$converged)
  # A brms fit is never recorded as singular.
  expect_true(is.na(record(c(a = 1.2), divergent = 1)$singular))
  expect_true(is.na(record(c(a = 1.002))$singular))
})

test_that("real brms fits: the groups the diagnostics flag are the groups listed", {
  # Compiles two Stan models, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- fq12_one_column(rbind(fq12_group("g1", fq12_s2), fq12_group("g3", fq12_s3)))
  # Sixty iterations a chain cannot mix: the run is short on purpose. Each fit's own
  # record is kept as it is made, one fit a group, in the order of the groups.
  real_record <- maihda_fit_diagnostics
  records <- list()
  local_mocked_bindings(maihda_fit_diagnostics = function(model, ...) {
    r <- real_record(model, ...)
    records[[length(records) + 1L]] <<- r
    r
  })
  obs <- fq12_keep_seed(fq12_observe(compare_maihda_groups(
    y ~ (1 | stratum), data = d, group = "country", engine = "brms",
    chains = 2, iter = 60, warmup = 30, seed = 1, refresh = 0)))

  expect_length(records, 2L)
  expect_identical(vapply(records, function(r) r$engine, character(1)),
                   c("brms", "brms"))
  flagged <- c("g1", "g3")[vapply(records, function(r) isFALSE(r$converged),
                                  logical(1))]
  # The run is too short to pass the diagnostics; if this fails, shorten it further.
  expect_gt(length(flagged), 0L)
  fq <- fq12_quality(obs$warns)
  expect_identical(fq,
                   maihda_group_fit_quality_message("brms", character(0), flagged))
  expect_false(grepl("lme4", fq, fixed = TRUE))
  expect_false(any(vapply(records, function(r) isTRUE(r$singular), logical(1))))
})
