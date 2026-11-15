# Stan-backed test blocks.
#
# The blocks that compile and sample a Stan model are opt-in: each starts with
# skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true", ...), and the weekly brms-tests
# workflow runs the suite with that gate open. Two things live here: the fit of a
# Stan fixture that tries further seeds, and the list of gated blocks that the
# workflow uses to check that every one of them ran.

# ---- Stan fixtures -------------------------------------------------------------
#
# Some gated blocks read a posterior they fit themselves, and guard the fit with the
# package's own verdict (maihda_fit_diagnostics(), brms branch: largest Rhat at most
# 1.01 and no divergent transition), so that a fit that has not converged fails the
# block instead of passing its assertions on noise. A seeded Stan fit is the same fit
# on one platform only: another compiler, Stan release or runner image draws other
# numbers from the same seed, and a fit that met the bar can then miss it by one
# divergent transition.
#
# stan_fixture() fits with up to three seeds, in the order given, and returns the
# first fit that meets the bar. It records one expectation: a success for such a fit,
# or, when no seed gives one, a failure that names each seed's reason, and the last
# fit is returned. The bar is the package's unless the block says otherwise; trying
# another seed never relaxes it.
#
# The settings and seeds in the blocks were measured over seeds 1 to 20 on GitHub's
# ubuntu-26.04 runner (2026-11-15). Two things carried over from one fit to the
# next: the largest Rhat is governed by the effective sample size per chain, so a
# fit that misses the bar on Rhat needs longer chains, and more chains of the same
# length do not help; a fit that misses it on divergent transitions needs a higher
# adapt_delta, and every further draw is one more chance of one.
#
#   name         a short label for the failure message and the log (see below)
#   fit          function(seed) returning the fitted object
#   seeds        one to three distinct seeds
#   diagnostics  function(object) returning the list that holds `converged` and
#                `messages`; the default reads the package's verdict
#   bar          the bar `diagnostics` applies, in words, for the failure message
#
# A block that deliberately asks for less than the package's bar (the growth-model
# pair of test-longitudinal.R asks for the Rhat half only, and says why) passes its
# own `diagnostics` and `bar`.
#
# R's random-number state is set from the seed before each fit, so a fit depends on
# its seed and not on what ran before it, also where the seed reaches the sampler
# only through R's stream. On exit the state is put back to what it was before the
# call: the stream the rest of the suite reads (see setup-seed.R) depends neither on
# how many seeds were tried nor on the sampler settings.
#
# When the environment variable MAIHDA_STAN_FIXTURE_LOG names a file, one
# tab-separated line per fit is appended to it (name, seed, converged, seconds,
# reason). The weekly workflow prints that file, so a fixture that has come to need
# its second or third seed shows in a run that is still green.

stan_fixture_reason <- function(diag) {
  if (is.null(diag)) return("the fit carries no diagnostics")
  if (length(diag$messages)) return(paste(diag$messages, collapse = " "))
  if (identical(diag$converged, FALSE)) return("converged is FALSE, with no message")
  "no verdict (converged is not TRUE: Rhat could not be read)"
}

# The fits and their verdicts, with no expectation recorded: list(fit, tried), where
# `tried` has one row per seed used (seed, converged, reason, seconds).
stan_fixture_try <- function(fit, seeds, diagnostics = function(x) x$diagnostics) {
  stopifnot(is.function(fit), is.function(diagnostics), is.numeric(seeds),
            length(seeds) >= 1L, length(seeds) <= 3L, !anyNA(seeds),
            anyDuplicated(seeds) == 0L)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)

  tried <- vector("list", length(seeds))
  for (i in seq_along(seeds)) {
    set.seed(seeds[i])
    start <- proc.time()[["elapsed"]]
    kept <- fit(seeds[i])
    diag <- diagnostics(kept)
    ok <- isTRUE(diag$converged)
    tried[[i]] <- data.frame(
      seed = seeds[i], converged = ok,
      reason = if (ok) "" else stan_fixture_reason(diag),
      seconds = proc.time()[["elapsed"]] - start,
      stringsAsFactors = FALSE)
    if (ok) break
  }
  list(fit = kept, tried = do.call(rbind, tried[seq_len(i)]))
}

stan_fixture <- function(name, fit, seeds, diagnostics = function(x) x$diagnostics,
                         bar = "largest Rhat at most 1.01 and no divergent transition") {
  stopifnot(is.character(name), length(name) == 1L, !is.na(name), nzchar(name),
            is.character(bar), length(bar) == 1L, !is.na(bar), nzchar(bar))
  res <- stan_fixture_try(fit, seeds, diagnostics)
  tried <- res$tried
  seed_label <- format(tried$seed, scientific = FALSE, trim = TRUE)

  log <- Sys.getenv("MAIHDA_STAN_FIXTURE_LOG")
  if (nzchar(log)) {
    cat(sprintf("%s\t%s\t%s\t%.0f\t%s", name, seed_label, tried$converged,
                tried$seconds, tried$reason),
        file = log, sep = "\n", append = TRUE)
  }

  testthat::expect(
    tried$converged[nrow(tried)],
    sprintf(paste0("Stan fixture '%s': none of the %d seed(s) tried gave a fit that ",
                   "meets the convergence bar (%s).\n%s"),
            name, nrow(tried), bar,
            paste0("  seed ", seed_label, ": ", tried$reason, collapse = "\n")))
  res$fit
}

# ---- Stan-gated blocks ---------------------------------------------------------
#
# stan_gated_blocks(dir): the test_that() blocks of the test files in `dir` whose
# code names MAIHDA_TEST_BRMS, as a data frame (file, test). A block whose
# description is not a string constant has test = NA.
stan_gated_blocks <- function(dir) {
  files <- sort(list.files(dir, pattern = "^test.*\\.[rR]$"))
  found <- list()
  visit <- function(e, file) {
    if (!is.call(e) && !is.expression(e)) return(invisible())
    if (is.call(e)) {
      fn <- e[[1L]]
      plain <- is.symbol(fn) && identical(as.character(fn), "test_that")
      qualified <- is.call(fn) && length(fn) == 3L &&
        as.character(fn[[1L]])[1L] %in% c("::", ":::") &&
        identical(as.character(fn[[3L]]), "test_that")
      if (plain || qualified) {
        mc <- tryCatch(match.call(function(desc, code) NULL, e), error = function(err) NULL)
        if (!is.null(mc$code) &&
            any(grepl("MAIHDA_TEST_BRMS", deparse(mc$code), fixed = TRUE))) {
          desc <- mc$desc
          found[[length(found) + 1L]] <<- data.frame(
            file = file,
            test = if (is.character(desc) && length(desc) == 1L) desc else NA_character_,
            stringsAsFactors = FALSE)
        }
      }
    }
    # Go on only where a test_that() call lies further in.
    if (!"test_that" %in% all.names(e)) return(invisible())
    for (i in seq_along(e)) {
      if (!identical(e[[i]], quote(expr = ))) visit(e[[i]], file)
    }
  }
  for (file in files) {
    path <- file.path(dir, file)
    # A file that never names the gate holds no gated block: not parsed.
    if (!any(grepl("MAIHDA_TEST_BRMS", readLines(path, warn = FALSE), fixed = TRUE))) next
    visit(parse(path, keep.source = FALSE), file)
  }
  if (!length(found)) {
    return(data.frame(file = character(0), test = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, found)
}

# stan_gated_not_run(results, dir): the gated blocks that did not run in `results`
# (what testthat::test_local() returns), as a data frame (file, test, why). A gated
# block counts as not run when it skipped, whatever the reason: the reason names only
# the skip that fired, so a block that skipped because a second package is missing
# does not mention the gate. It also counts when it is absent from the results (a
# skip at the top of its file, a description that cannot be matched). A skip whose
# reason names MAIHDA_TEST_BRMS is listed too, wherever it was raised. With the gate
# open the result must have no rows.
stan_gated_not_run <- function(results, dir) {
  blocks <- stan_gated_blocks(dir)
  key <- function(file, test) paste(file, test, sep = "\r")
  gated <- key(blocks$file, blocks$test)
  seen <- character(0)
  out <- list()
  add <- function(file, test, why) {
    out[[length(out) + 1L]] <<- data.frame(file = file, test = test, why = why,
                                           stringsAsFactors = FALSE)
  }
  for (el in results) {
    k <- key(el$file, el$test)
    seen <- c(seen, k)
    for (r in el$results) {
      if (!inherits(r, "expectation_skip")) next
      reason <- sub("^Reason: ", "", conditionMessage(r))
      if (k %in% gated || grepl("MAIHDA_TEST_BRMS", reason, fixed = TRUE)) {
        add(el$file, el$test, paste0("skipped: ", reason))
      }
    }
  }
  absent <- blocks[!gated %in% seen, , drop = FALSE]
  for (i in seq_len(nrow(absent))) {
    add(absent$file[i], absent$test[i],
        if (is.na(absent$test[i])) "its description is not a string constant" else "did not run")
  }
  if (!length(out)) {
    return(data.frame(file = character(0), test = character(0), why = character(0),
                      stringsAsFactors = FALSE))
  }
  unique(do.call(rbind, out))
}
