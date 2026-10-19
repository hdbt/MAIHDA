# Audit 2026-10-18 (rank-deficient-strata): the longitudinal vignette said that a
# rank-deficient stratum growth block is rare when the stratum slope variance is not near
# zero and usually evidence that the strata do not diverge; print() on a longitudinal PCV
# offered stratum_slope = FALSE as the remedy, and ?fit_maihda said to use it when the
# block collapses. Measured with lme4 through fit_maihda() (time | id) + (time | stratum)
# fits (stratum intercept SD 0.5, person intercept SD 0.7, person slope SD 0.25, residual
# SD 0.6; 100 data sets per design):
# - without stratum slope variance the block is rank-deficient in 56-77% of fits in every
#   design (8 to 32 strata, 5 or 20 people per stratum, 3 or 5 waves);
# - with a stratum slope variance of 0.04 it is rank-deficient in 46 / 28 / 7% of fits
#   with 8 / 16 / 32 strata of 5 people seen 3 times, 31 / 6 / 0% with 5 waves, and 11%
#   with 8 strata of 20 people; such fits understate the slope variance (0.028 for a
#   realised 0.037 with 8 strata, 0.018 for 0.033 with 32);
# - the AIC comparison keeps that slope variance in 22% of fits with 8 strata of 5 people
#   seen 3 times (74% with 32 strata), and a slope variance of zero in 2-8%;
# - stratum_slope = FALSE on such strata moves their slope variance to the person level
#   (+0.03 to +0.04); the VPC at the last wave is 0.18 for a true 0.27 (3 waves) and 0.10
#   for 0.33 (5 waves), falling in 83-100% of data sets where the true VPC rises; with
#   slopes 80% and intercepts 30% additive, its PCV_intercept is 0.13 to 0.17 too high
#   (median) against -0.04 to +0.02 with the stratum slope.
# The note, the vignette and the help now say so; no number, default or computation
# changed. No random numbers are drawn, so the suite's random stream is left as it was.

# A rank 1 linear block with stratum deviations s[1] + s[2] * t.
rank1_stratum_block <- function(s) {
  matrix(c(s[1]^2, s[1] * s[2], s[1] * s[2], s[2]^2), 2, 2,
         dimnames = list(c("(Intercept)", "wave"), c("(Intercept)", "wave")))
}

# A maihda_long_pcv as maihda_longitudinal_pcv() builds it, from two stratum blocks.
rank_note_stub <- function(Sn, Sa, null_rd, adjusted_rd, times = 0:4) {
  vn <- maihda_var_at_time(Sn, times)
  va <- maihda_var_at_time(Sa, times)
  structure(
    list(pcv_intercept = (vn[1] - va[1]) / vn[1],
         pcv_slope = if (null_rd) NA_real_ else 1 - Sa[2, 2] / Sn[2, 2],
         var_baseline_null = vn[1], var_baseline_adjusted = va[1],
         var_slope_null = Sn[2, 2], var_slope_adjusted = Sa[2, 2],
         ref_time = times[1], time_center = 0,
         pcv_t = data.frame(time = times, var_null = vn, var_adjusted = va,
                            pcv = (vn - va) / vn),
         Sigma_stratum_null = Sn, Sigma_stratum_adjusted = Sa,
         time = "wave", time_degree = 1L, ml_refit = FALSE,
         estimation = "fitted", estimation_used = "fitted",
         null_at_boundary = FALSE, adjusted_at_boundary = FALSE,
         null_rank_deficient = null_rd, null_slope_at_boundary = null_rd,
         adjusted_rank_deficient = adjusted_rd, stratum_slope = TRUE,
         variance_summary = "estimate", Sigma_summary = "estimate"),
    class = "maihda_long_pcv")
}

rank_note_lines <- function(x) cli::ansi_strip(utils::capture.output(print(x)))
rank_note_text <- function(x) gsub("\\s+", " ", paste(rank_note_lines(x), collapse = " "))

test_that("the rank-deficiency note does not read the block as strata that do not diverge", {
  new_sentences <- paste(
    "A rank-deficient block is common with few strata even when the strata diverge over",
    "time, so it does not show that they do not; it is rarer with more strata, more people",
    "per stratum or more measurement occasions. stratum_slope = FALSE models a",
    "time-constant between-stratum variance; where the strata do diverge, it moves their",
    "slope differences to the person level.")
  full <- matrix(c(0.2, 0.01, 0.01, 0.02), 2, 2,
                 dimnames = list(c("(Intercept)", "wave"), c("(Intercept)", "wave")))
  Sn <- rank1_stratum_block(c(0.4, 0.1))
  Sa <- rank1_stratum_block(c(0.2, 0.05))
  for (flags in list(c(TRUE, TRUE), c(TRUE, FALSE), c(FALSE, TRUE))) {
    x <- rank_note_stub(if (flags[1]) Sn else full, if (flags[2]) Sa else full / 4,
                        null_rd = flags[1], adjusted_rd = flags[2])
    out <- rank_note_text(x)
    info <- paste(c("null", "adjusted")[flags], collapse = " + ")
    expect_match(out, "RANK-DEFICIENT", fixed = TRUE, info = info)
    expect_match(out, new_sentences, fixed = TRUE, info = info)
    expect_false(grepl("needs either more measurement occasions", out, fixed = TRUE),
                 info = info)
    expect_false(grepl("says so explicitly", out, fixed = TRUE), info = info)
    # The paragraph keeps the width of the lines around it.
    lines <- rank_note_lines(x)
    from <- grep("are still reported", lines, fixed = TRUE)
    expect_length(from, 1L)
    expect_true(all(nchar(lines[from:length(lines)]) <= 80L), info = info)
  }
  # Without a rank-deficient block there is no such paragraph.
  healthy <- rank_note_text(rank_note_stub(full, full / 4, FALSE, FALSE))
  expect_false(grepl("common with few strata", healthy, fixed = TRUE))
  expect_false(grepl("stratum_slope = FALSE models", healthy, fixed = TRUE))
})

test_that("the vignette and ?fit_maihda no longer call the block evidence of no divergence", {
  skip_on_cran()
  root <- testthat::test_path("..", "..")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(file.path(root, "man")), "the package tree is not available")
  flat <- function(...) {
    p <- file.path(root, ...)
    if (!file.exists(p)) return(NA_character_)
    gsub("[[:space:]]+", " ", paste(readLines(p, warn = FALSE), collapse = " "))
  }
  rd <- flat("man", "fit_maihda.Rd")
  expect_false(grepl("routinely collapses", rd, fixed = TRUE))
  expect_false(grepl("Use it when the stratum slope variance is at the singularity boundary",
                     rd, fixed = TRUE))
  expect_true(grepl("It assumes that the strata do not diverge over time.", rd, fixed = TRUE))
  expect_true(grepl(paste("produce such fits when the strata diverge too, and where they do,",
                          "\\code{FALSE} moves their slope differences into the",
                          "person-level slope variance."), rd, fixed = TRUE))

  vignette <- flat("vignettes", "longitudinal.Rmd")
  if (!is.na(vignette)) {
    expect_false(grepl("rare when it is not", vignette, fixed = TRUE))
    expect_false(grepl("is usually evidence that the strata do not", vignette, fixed = TRUE))
    expect_true(grepl("Where the strata do not diverge it is common in any design.",
                      vignette, fixed = TRUE))
    expect_true(grepl(paste("With few strata, a rank-deficient block is therefore weak",
                            "evidence that the strata do not diverge over time."),
                      vignette, fixed = TRUE))
    expect_true(grepl("With few strata this comparison has little power", vignette,
                      fixed = TRUE))
    expect_true(grepl(paste("the constraint moves their slope differences into the",
                            "individual-level slope variance"), vignette, fixed = TRUE))
  }
})
