# Slow test blocks.
#
# Bootstrap and simulation regressions and the heaviest model fits take seconds
# each and make up much of the suite's run time. A block that starts with
# skip_slow_test() is skipped when the environment variable MAIHDA_SKIP_SLOW_TESTS
# is "true". R-CMD-check.yaml sets it on four of the five R CMD check legs and in
# the job that runs the test files reading the package tree, so each such block
# runs on the ubuntu-latest (release) leg only. Unset, as in local runs, covr and
# on CRAN, every block runs.
#
# setup-seed.R fixes one random stream for the whole suite, and some blocks draw
# from it without a seed of their own. A block that ran on one leg and was skipped
# on another would hand those later blocks different numbers on the two legs. So
# skip_slow_test() also puts R's random-number state back when the block exits:
# a gated block moves the stream neither when it runs nor when it is skipped. Call
# it first in the block, before anything draws a random number. A later on.exit()
# in the block needs add = TRUE, or it replaces that restore.
skip_slow_test <- function() {
  if (identical(Sys.getenv("MAIHDA_SKIP_SLOW_TESTS"), "true")) {
    testthat::skip("slow test (MAIHDA_SKIP_SLOW_TESTS is \"true\")")
  }
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  restore <- function() {
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }
  # Registered on the calling block, as on.exit(restore(), add = TRUE) written
  # there would be.
  do.call(base::on.exit, list(as.call(list(restore)), TRUE, FALSE),
          envir = parent.frame())
  invisible()
}
