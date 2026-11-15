# The weekly Stan job on a new runner image (batch 2026-11-15-brms-weekly).
#
# F03. A seeded Stan fixture fit is the same fit on one platform only, so the blocks
#   that guard such a fit with the package's convergence verdict go through
#   stan_fixture() (helper-stan-tests.R): up to three seeds, the first fit that meets
#   the bar, one expectation that fails with every seed's reason when none does, and
#   R's random-number state left as it was found. Pinned here on stand-in fits.
# F04. The workflow checks that every Stan-gated block ran with
#   stan_gated_not_run(): a gated block that skipped for any reason, or that is not in
#   the results at all. A skip's reason names only the skip that fired, so a check on
#   reasons misses a gated block that skipped because a second package is missing.
#   Pinned here on a small test directory run through testthat.
# F01 (three blocks that used an optional package without a skip) has no block here:
#   which packages a block reaches cannot be read from its source, so the check is a
#   run of the suite against a library without the package, which the weekly job is.
#
# No block of this file may leave R's random-number state changed (setup-seed.R).

wk_keep_stream <- function(expr) {
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

wk_with_env <- function(name, value, expr) {
  old <- Sys.getenv(name, unset = NA)
  on.exit({
    if (is.na(old)) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(old), name))
  }, add = TRUE)
  if (is.na(value)) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(value), name))
  force(expr)
}

# A stand-in for a fitting function: the k-th call returns a fit whose verdict is
# verdicts[[k]], after drawing from R's stream as a sampler does. `calls` records
# the seed each call was given and the first number it drew.
wk_fits <- function(verdicts, messages = NULL) {
  calls <- new.env()
  calls$seed <- numeric(0)
  calls$draw <- numeric(0)
  fit <- function(seed) {
    k <- length(calls$seed) + 1L
    draw <- stats::runif(3)[1]
    calls$seed <- c(calls$seed, seed)
    calls$draw <- c(calls$draw, draw)
    list(seed = seed, draw = draw,
         diagnostics = list(converged = verdicts[[k]],
                            messages = if (is.null(messages)) character(0) else messages[[k]]))
  }
  list(fit = fit, calls = calls)
}

wk_rhat <- "Max Rhat = 1.023 (> 1.01): the chains may not have converged."
wk_div <- "4 divergent transition(s): the posterior may be biased; consider increasing adapt_delta."

test_that("stan_fixture_try() keeps the first fit that meets the bar, seeds in order", {
  wk_keep_stream({
    # Third seed.
    f <- wk_fits(list(FALSE, FALSE, TRUE), list(wk_rhat, c(wk_rhat, wk_div), character(0)))
    res <- stan_fixture_try(f$fit, c(11, 12, 13))
    expect_identical(f$calls$seed, c(11, 12, 13))
    expect_identical(res$fit$seed, 13)
    expect_identical(res$tried$seed, c(11, 12, 13))
    expect_identical(res$tried$converged, c(FALSE, FALSE, TRUE))
    expect_identical(res$tried$reason, c(wk_rhat, paste(wk_rhat, wk_div), ""))
    expect_true(all(is.finite(res$tried$seconds)) && all(res$tried$seconds >= 0))

    # First seed: no further fit.
    f <- wk_fits(list(TRUE, TRUE, TRUE))
    res <- stan_fixture_try(f$fit, c(11, 12, 13))
    expect_identical(f$calls$seed, 11)
    expect_identical(res$fit$seed, 11)
    expect_identical(res$tried$converged, TRUE)

    # None: the last fit comes back, every seed with its reason.
    f <- wk_fits(list(FALSE, FALSE, FALSE), list(wk_div, wk_rhat, wk_div))
    res <- stan_fixture_try(f$fit, c(11, 12, 13))
    expect_identical(f$calls$seed, c(11, 12, 13))
    expect_identical(res$fit$seed, 13)
    expect_identical(res$tried$converged, c(FALSE, FALSE, FALSE))
    expect_identical(res$tried$reason, c(wk_div, wk_rhat, wk_div))
  })
})

test_that("only a verdict of TRUE meets the bar", {
  wk_keep_stream({
    # NA: Rhat could not be read. Not a pass.
    f <- wk_fits(list(NA, TRUE))
    res <- stan_fixture_try(f$fit, c(1, 2))
    expect_identical(res$fit$seed, 2)
    expect_identical(res$tried$converged, c(FALSE, TRUE))
    expect_match(res$tried$reason[1], "Rhat could not be read", fixed = TRUE)

    # No diagnostics at all.
    res <- stan_fixture_try(function(seed) list(seed = seed), 1)
    expect_identical(res$tried$converged, FALSE)
    expect_identical(res$tried$reason, "the fit carries no diagnostics")

    # FALSE without a message.
    f <- wk_fits(list(FALSE))
    expect_identical(stan_fixture_try(f$fit, 1)$tried$reason,
                     "converged is FALSE, with no message")

    # The verdict is read where `diagnostics` says it is.
    nested <- function(seed) {
      list(seed = seed, model = list(diagnostics = list(converged = seed == 8,
                                                        messages = wk_div)))
    }
    res <- stan_fixture_try(nested, c(7, 8), diagnostics = function(x) x$model$diagnostics)
    expect_identical(res$fit$seed, 8)
    expect_identical(res$tried$converged, c(FALSE, TRUE))
    expect_identical(stan_fixture_try(nested, c(7, 8))$tried$converged, c(FALSE, FALSE))
  })
})

test_that("one to three distinct seeds, no more", {
  wk_keep_stream({
    f <- wk_fits(list(TRUE, TRUE, TRUE, TRUE))
    expect_error(stan_fixture_try(f$fit, c(1, 2, 3, 4)))
    expect_error(stan_fixture_try(f$fit, numeric(0)))
    expect_error(stan_fixture_try(f$fit, c(1, 1)))
    expect_error(stan_fixture_try(f$fit, c(1, NA)))
    expect_error(stan_fixture_try(f$fit, "1"))
    expect_identical(f$calls$seed, numeric(0))
  })
})

test_that("each fit starts from its seed, and the stream is put back", {
  wk_keep_stream({
    first_draw <- function(seed) {
      set.seed(seed)
      stats::runif(3)[1]
    }
    want <- vapply(c(5, 6, 7), first_draw, numeric(1))

    for (verdicts in list(list(TRUE), list(FALSE, FALSE, TRUE), list(FALSE, FALSE, FALSE))) {
      set.seed(99)
      stats::runif(10)
      before <- get(".Random.seed", envir = globalenv())
      f <- wk_fits(verdicts)
      res <- stan_fixture_try(f$fit, c(5, 6, 7))
      # The state after the call does not depend on how many seeds were tried.
      # (identical() inside expect_true(): the comparison expect_identical() prints
      # for two unequal states overflows.)
      expect_true(identical(get(".Random.seed", envir = globalenv()), before))
      # A fit's draws depend on its seed only.
      expect_identical(f$calls$draw, want[seq_along(f$calls$seed)])
    }

    # Also when the fitting function stops.
    set.seed(99)
    before <- get(".Random.seed", envir = globalenv())
    expect_error(stan_fixture_try(function(seed) {
      stats::runif(5)
      stop("the sampler stopped")
    }, c(5, 6)), "the sampler stopped")
    expect_true(identical(get(".Random.seed", envir = globalenv()), before))

    # And when there was no state yet.
    rm(".Random.seed", envir = globalenv())
    f <- wk_fits(list(FALSE, TRUE))
    stan_fixture_try(f$fit, c(5, 6))
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  })
})

test_that("stan_fixture() records one expectation and returns the fit it kept", {
  # The stand-in fits of this block are kept out of a log the caller of the suite
  # may have asked for.
  wk_keep_stream(wk_with_env("MAIHDA_STAN_FIXTURE_LOG", NA, {
    f <- wk_fits(list(FALSE, TRUE), list(wk_div, character(0)))
    m <- NULL
    expect_success(m <- stan_fixture("demo", f$fit, c(1, 2, 3)))
    expect_identical(m$seed, 2)
    expect_identical(f$calls$seed, c(1, 2))

    # No seed meets the bar: one failure, which gives every seed's reason, and the
    # last fit comes back so that the block goes on.
    f <- wk_fits(list(FALSE, FALSE, FALSE), list(wk_div, wk_rhat, c(wk_rhat, wk_div)))
    m <- NULL
    expect_failure(m <- stan_fixture("demo", f$fit, c(1, 2, 20260831)),
                   "Stan fixture 'demo': none of the 3 seed(s) tried", fixed = TRUE)
    expect_identical(m$seed, 20260831)
    f <- wk_fits(list(FALSE, FALSE, FALSE), list(wk_div, wk_rhat, c(wk_rhat, wk_div)))
    msg <- tryCatch(stan_fixture("demo", f$fit, c(1, 2, 20260831)),
                    expectation_failure = function(e) conditionMessage(e))
    expect_identical(msg, paste0(
      "Stan fixture 'demo': none of the 3 seed(s) tried gave a fit that meets the ",
      "convergence bar (largest Rhat at most 1.01 and no divergent transition).\n",
      "  seed 1: ", wk_div, "\n",
      "  seed 2: ", wk_rhat, "\n",
      "  seed 20260831: ", wk_rhat, " ", wk_div))

    # A block that applies a bar of its own has it named in the failure.
    f <- wk_fits(list(FALSE), list("Max Rhat = 1.014 (null model) and 1.004 (adjusted model)."))
    msg <- tryCatch(
      stan_fixture("pair", f$fit, 7, bar = "both models' largest Rhat below 1.01"),
      expectation_failure = function(e) conditionMessage(e))
    expect_identical(msg, paste0(
      "Stan fixture 'pair': none of the 1 seed(s) tried gave a fit that meets the ",
      "convergence bar (both models' largest Rhat below 1.01).\n",
      "  seed 7: Max Rhat = 1.014 (null model) and 1.004 (adjusted model)."))

    # A name and a bar are required.
    f <- wk_fits(list(TRUE))
    expect_error(stan_fixture(NA_character_, f$fit, 1))
    expect_error(stan_fixture("", f$fit, 1))
    expect_error(stan_fixture("demo", f$fit, 1, bar = ""))
    expect_identical(f$calls$seed, numeric(0))
  }))
})

test_that("stan_fixture() appends one line per fit to MAIHDA_STAN_FIXTURE_LOG", {
  wk_keep_stream({
    log <- tempfile("stan-fixture-log-")
    on.exit(unlink(log), add = TRUE)
    wk_with_env("MAIHDA_STAN_FIXTURE_LOG", log, {
      f <- wk_fits(list(FALSE, TRUE), list(wk_div, character(0)))
      expect_success(stan_fixture("first", f$fit, c(1, 2, 3)))
      f <- wk_fits(list(TRUE))
      expect_success(stan_fixture("second", f$fit, 100000))
    })
    got <- utils::read.delim(log, header = FALSE, quote = "", stringsAsFactors = FALSE,
                             col.names = c("name", "seed", "converged", "seconds", "reason"))
    expect_identical(got$name, c("first", "first", "second"))
    expect_identical(as.numeric(got$seed), c(1, 2, 100000))
    expect_match(readLines(log)[3], "^second\t100000\tTRUE\t[0-9]+\t$")
    expect_identical(got$converged, c(FALSE, TRUE, TRUE))
    expect_identical(got$reason[1], wk_div)

    # Unset: nothing is written.
    unlink(log)
    wk_with_env("MAIHDA_STAN_FIXTURE_LOG", NA, {
      f <- wk_fits(list(TRUE))
      expect_success(stan_fixture("third", f$fit, 1))
    })
    expect_false(file.exists(log))
  })
})

# ---- the Stan-gated blocks -------------------------------------------------------

# The gate's name, put together here so that no block of this file names it: the
# blocks below write it into test files of their own.
wk_gate <- paste0("MAIHDA_TEST", "_BRMS")

# A test directory with gated and ungated blocks. "notThere" stands for an optional
# package that is not installed.
wk_gated_dir <- function() {
  dir <- tempfile("stan-gate-")
  dir.create(dir)
  gate <- sprintf("skip_if(Sys.getenv(\"%s\") != \"true\", \"opt-in; set %s=true\")",
                  wk_gate, wk_gate)
  missing <- "skip_if_not_installed(\"aPackageThatIsNotThere\")"
  writeLines(c(
    "test_that(\"gated, runs\", {", gate, "expect_true(TRUE)", "})",
    "test_that(\"gated, a second package is missing\", {", gate, missing,
    "expect_true(TRUE)", "})",
    "test_that(\"not gated, a package is missing\", {", missing, "expect_true(TRUE)", "})",
    "testthat::test_that(desc = \"gated, qualified call and named arguments\", code = {",
    gate, "expect_true(TRUE)", "})",
    "test_that(\"not gated, runs\", {", "expect_true(TRUE)", "})"),
    file.path(dir, "test-a.R"))
  writeLines(c(
    sprintf("opt_in <- function() %s", gate),
    "test_that(\"gated through a function of the file\", {", "opt_in()",
    "expect_true(TRUE)", "})",
    "test_that(paste(\"gated,\", \"computed description\"), {", gate,
    "expect_true(TRUE)", "})"),
    file.path(dir, "test-b.R"))
  # Names the gate, and is not a test file.
  writeLines(sprintf("gate_name <- \"%s\"", wk_gate), file.path(dir, "helper-gate.R"))
  dir
}

# Run a test directory with the gate open ("true") or closed (NA: unset).
wk_run_dir <- function(dir, gate_value) {
  wk_with_env(wk_gate, gate_value, {
    testthat::test_dir(dir, reporter = "silent", stop_on_failure = FALSE,
                       stop_on_warning = FALSE,
                       env = new.env(parent = asNamespace("testthat")))
  })
}

test_that("stan_gated_blocks() lists the blocks whose code names the gate", {
  dir <- wk_gated_dir()
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  blocks <- stan_gated_blocks(dir)
  expect_identical(names(blocks), c("file", "test"))
  expect_identical(blocks$file, c("test-a.R", "test-a.R", "test-a.R", "test-b.R"))
  expect_identical(blocks$test, c("gated, runs", "gated, a second package is missing",
                                  "gated, qualified call and named arguments",
                                  NA_character_))

  # A directory without one.
  empty <- tempfile("stan-gate-none-")
  dir.create(empty)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  writeLines(c("test_that(\"plain\", {", "expect_true(TRUE)", "})"),
             file.path(empty, "test-c.R"))
  none <- stan_gated_blocks(empty)
  expect_identical(nrow(none), 0L)
  expect_identical(names(none), c("file", "test"))
})

test_that("the suite's own gated blocks are found, each with a description to match", {
  blocks <- stan_gated_blocks(test_path())
  expect_gt(nrow(blocks), 0L)
  expect_false(anyNA(blocks$test))
  expect_identical(anyDuplicated(blocks), 0L)
  # The two blocks the weekly workflow names must be among them.
  expect_true(all(c("fit_maihda fits a binary outcome with brms via bernoulli()",
                    "brms summary returns a posterior credible interval for the VPC/ICC")
                  %in% blocks$test))
  # The six fixture blocks are gated, and this file adds none.
  expect_true(all(c("test-contextual.R", "test-cross-classified.R", "test-tidiers.R",
                    "test-longitudinal.R", "test-negbinomial.R", "test-ordinal-engine.R")
                  %in% blocks$file))
  expect_false("test-audit-2026-11-15-brms-weekly-F01-F03-F04.R" %in% blocks$file)
})

test_that("a gated block that skipped for a second package is listed", {
  skip_on_cran()
  wk_keep_stream({
    dir <- wk_gated_dir()
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)

    # Gate open, as in the weekly job.
    res <- wk_run_dir(dir, "true")
    df <- as.data.frame(res)
    expect_identical(sum(df$failed), 0L)
    expect_false(any(df$error))
    reasons <- unlist(lapply(res, function(el) {
      vapply(Filter(function(r) inherits(r, "expectation_skip"), el$results),
             conditionMessage, character(1))
    }))
    # Two blocks skipped, and neither reason names the gate: a check on reasons
    # sees nothing.
    expect_identical(length(reasons), 2L)
    expect_false(any(grepl(wk_gate, reasons, fixed = TRUE)))

    not_run <- stan_gated_not_run(res, dir)
    expect_identical(names(not_run), c("file", "test", "why"))
    expect_identical(not_run$file, c("test-a.R", "test-b.R"))
    expect_identical(not_run$test, c("gated, a second package is missing", NA_character_))
    expect_match(not_run$why[1], "^skipped: .*aPackageThatIsNotThere")
    expect_identical(not_run$why[2], "its description is not a string constant")

    # A gated block that is not in the results at all.
    ran <- vapply(res, function(el) identical(el$test, "gated, runs"), logical(1))
    expect_identical(sum(ran), 1L)
    absent <- stan_gated_not_run(res[!ran], dir)
    expect_true(any(absent$test %in% "gated, runs" & absent$why == "did not run"))

    # Gate closed: every gated block is listed, the one gated through a function of
    # its file by its reason.
    res <- wk_run_dir(dir, NA)
    closed <- stan_gated_not_run(res, dir)
    expect_setequal(closed$test[!is.na(closed$test)],
                    c("gated, runs", "gated, a second package is missing",
                      "gated, qualified call and named arguments",
                      "gated through a function of the file",
                      "gated, computed description"))
    expect_false(any(c("not gated, a package is missing", "not gated, runs") %in% closed$test))
  })
})

test_that("with every gated block run, nothing is listed", {
  skip_on_cran()
  wk_keep_stream({
    dir <- tempfile("stan-gate-ok-")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    writeLines(c(
      "test_that(\"gated, runs\", {",
      sprintf("skip_if(Sys.getenv(\"%s\") != \"true\", \"opt-in\")", wk_gate),
      "expect_true(TRUE)", "})",
      "test_that(\"not gated, a package is missing\", {",
      "skip_if_not_installed(\"aPackageThatIsNotThere\")", "expect_true(TRUE)", "})"),
      file.path(dir, "test-a.R"))
    not_run <- stan_gated_not_run(wk_run_dir(dir, "true"), dir)
    expect_identical(nrow(not_run), 0L)
    expect_identical(names(not_run), c("file", "test", "why"))
  })
})
