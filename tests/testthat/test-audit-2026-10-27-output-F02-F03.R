# Audit 2026-10-27 (output batch, packet P2): print() of a maihda() analysis.
#
# F02  A note in the PCV block of the two-model print opens and closes with one blank
#      line. Two stacked up after the maximum-likelihood note when the
#      discriminatory-accuracy block followed, and after the boundary note when the
#      maximum-likelihood note or that block followed.
# F03  The two-model print and print.maihda_long_pcv() name the variances the PCV
#      compares, in calculate_pcv()'s words, wherever no other line of the print names
#      them (the user's decision, 2026-10-27).

# Leaves R's random-number state as it found it (see setup-seed.R).
p2_keep_rng <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  force(expr)
}

p2_lines <- function(x) cli::ansi_strip(utils::capture.output(print(x)))

# A binary outcome on 3 x 3 strata of 40, built without random numbers.
p2_binomial <- function(...) {
  g <- expand.grid(a = factor(1:3), b = factor(1:3), rep = 1:40)
  ia <- as.integer(g$a)
  ib <- as.integer(g$b)
  eta <- -0.4 + 0.5 * (ia - 2) + 0.3 * (ib - 2) + 0.6 * sin(3.1 * ia * ib)
  g$y <- as.integer(((g$rep * 37) %% 40) / 40 < stats::plogis(eta))
  p2_keep_rng(suppressWarnings(suppressMessages(
    maihda(y ~ a + b + (1 | a:b), data = g, family = "binomial", ...))))
}

# A Gaussian outcome on 3 x 3 strata of 20, built without random numbers.
p2_gaussian <- function(...) {
  g <- expand.grid(a = factor(1:3), b = factor(1:3), rep = 1:20)
  ia <- as.integer(g$a)
  ib <- as.integer(g$b)
  g$y <- 1 + 0.6 * ia - 0.4 * ib + 0.5 * sin(2.3 * ia * ib) + sin(1.7 * seq_len(nrow(g)))
  p2_keep_rng(suppressWarnings(suppressMessages(
    maihda(y ~ a + b + (1 | a:b), data = g, ...))))
}

# A longitudinal PCV built by hand, as maihda_longitudinal_pcv() returns it for an
# lmer growth pair under the default; `...` replaces fields, and a NULL drops one.
p2_long_pcv <- function(...) {
  S <- matrix(c(0.4, 0.01, 0.01, 0.04), 2, 2,
              dimnames = list(c("(Intercept)", "wave"), c("(Intercept)", "wave")))
  x <- list(pcv_intercept = 0.9, pcv_slope = 0.85,
            var_baseline_null = 0.4, var_baseline_adjusted = 0.04,
            var_slope_null = 0.04, var_slope_adjusted = 0.006,
            ref_time = 0, time_center = 0,
            pcv_t = data.frame(time = 0:2, var_null = 0.4, var_adjusted = 0.04,
                               pcv = 0.9),
            Sigma_stratum_null = S, Sigma_stratum_adjusted = S / 10,
            time = "wave", time_degree = 1L, ml_refit = FALSE, estimation = "fitted",
            estimation_used = "fitted", ml_variances = FALSE, null_at_boundary = FALSE,
            adjusted_at_boundary = FALSE, null_rank_deficient = FALSE,
            null_slope_at_boundary = FALSE, adjusted_rank_deficient = FALSE,
            stratum_slope = TRUE, variance_summary = "estimate",
            Sigma_summary = "estimate")
  args <- list(...)
  for (nm in names(args)) x[nm] <- list(args[[nm]])
  x <- x[!vapply(x, is.null, logical(1))]
  structure(x, class = "maihda_long_pcv")
}

reml_label <- "as fitted (REML for Gaussian lmer, matching summary())"

test_that("each PCV note of the two-model print has one blank line either side", {
  a <- p2_binomial()
  expect_true(isTRUE(a$pcv$ml_variances))
  expect_false(is.null(a$summary$discriminatory_accuracy))

  grid <- expand.grid(boundary = c(FALSE, TRUE), ml = c(FALSE, TRUE),
                      da = c(FALSE, TRUE), group = c(FALSE, TRUE))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    info <- paste(names(g), unlist(g), sep = " = ", collapse = ", ")
    x <- a
    x$pcv$adjusted_at_boundary <- g$boundary
    x$pcv$ml_variances <- g$ml
    if (!g$da) x$summary$discriminatory_accuracy <- NULL
    if (g$group) {
      x$groups <- data.frame(group = c("g1", "g2"), n = c(10L, 12L))
      x$group_var <- "c"
    }
    out <- p2_lines(x)
    blank <- out == ""
    expect_false(any(blank[-1] & blank[-length(blank)]), info = info)

    starts <- grep("^Note: ", out)
    expect_length(starts, g$boundary + g$ml)
    last <- NA_integer_
    for (s in starts) {
      expect_identical(out[s - 1], "", info = info)
      last <- s
      while (last < length(out) && nzchar(out[last + 1])) last <- last + 1
      expect_identical(out[last + 1], "", info = info)
    }
    if (length(starts) > 0) {
      # After the last note: the discriminatory-accuracy block, else the Strata line.
      expect_match(out[last + 2],
                   if (g$da) "^Discriminatory accuracy" else "^Strata: ", info = info)
    }
  }
})

test_that("the two-model print names the basis in calculate_pcv()'s words", {
  a <- p2_gaussian()
  expect_identical(a$pcv$estimation_used, "fitted")
  out <- p2_lines(a)
  i <- grep("^Between-stratum variance: ", out)
  expect_length(i, 1)
  expect_identical(out[i + 1], paste0("  Variance basis: ", reml_label))
  expect_identical(sum(grepl("Variance basis", out, fixed = TRUE)), 1L)

  a_ml <- p2_gaussian(estimation = "ML")
  expect_identical(a_ml$pcv$estimation_used, "ML")
  out_ml <- p2_lines(a_ml)
  expect_identical(out_ml[grep("^Between-stratum variance: ", out_ml) + 1],
                   "  Variance basis: ML-refit (maximum-likelihood variances; see ?calculate_pcv)")

  # Every basis reads as print(calculate_pcv()) words it, wrapped within 84 characters.
  for (basis in c("fitted", "ML", "mixed", "posterior")) {
    x <- a
    x$pcv$estimation_used <- basis
    out_b <- p2_lines(x)
    first <- grep("^  Variance basis: ", out_b)
    expect_length(first, 1)
    end <- first
    while (grepl("^    \\S", out_b[end + 1])) end <- end + 1
    expect_lte(max(nchar(out_b[first:end])), 84)
    expect_identical(paste(trimws(out_b[first:end]), collapse = " "),
                     grep("^Variance basis: ", p2_lines(x$pcv), value = TRUE),
                     info = basis)
  }
  expect_gt(length(strwrap(paste("Variance basis:",
                                 MAIHDA:::maihda_pcv_basis_label("mixed")),
                           width = 84)), 1)

  # An object stored without `estimation` does not record its basis: no line.
  old <- a
  old$pcv$estimation <- NULL
  old$pcv$estimation_used <- NULL
  expect_false(any(grepl("Variance basis", p2_lines(old), fixed = TRUE)))
})

test_that("a glmer pair's maximum-likelihood note or REML-type line stands alone", {
  a <- p2_binomial()
  out <- p2_lines(a)
  expect_true(any(grepl("Note: these are maximum-likelihood variances", out,
                        fixed = TRUE)))
  expect_false(any(grepl("Variance basis", out, fixed = TRUE)))

  r <- p2_binomial(estimation = "REML")
  expect_identical(r$pcv$estimation_used, "REML")
  out_r <- p2_lines(r)
  expect_true(any(grepl("REML-type variances (estimation = \"REML\")", out_r,
                        fixed = TRUE)))
  expect_false(any(grepl("Variance basis", out_r, fixed = TRUE)))
  expect_false(any(grepl("Note: these are maximum-likelihood variances", out_r,
                         fixed = TRUE)))
})

test_that("print.maihda_long_pcv() names the basis where no note does", {
  # Default lmer growth pair: after the variance lines, before the closing paragraph.
  out <- p2_lines(p2_long_pcv())
  i <- grep("multiplicative/interaction part", out, fixed = TRUE)
  expect_identical(out[i + 1], paste0("Variance basis: ", reml_label))
  expect_identical(out[i + 2], "")
  expect_match(out[i + 3], "^The PCV is the share")

  # A glmer growth pair, whose variances are maximum-likelihood estimates.
  out_glmer <- p2_lines(p2_long_pcv(ml_variances = TRUE))
  expect_true(any(out_glmer == "Variance basis: as fitted (maximum likelihood)"))
  expect_false(any(grepl("REML for Gaussian lmer", out_glmer, fixed = TRUE)))

  # estimation = "ML" with nothing to refit (a glmer growth pair).
  out_ml <- p2_lines(p2_long_pcv(estimation = "ML", estimation_used = "ML",
                                 ml_variances = TRUE))
  expect_true(any(out_ml == paste0("Variance basis: ML-refit (maximum-likelihood ",
                                   "variances; see ?calculate_pcv)")))

  # An intercept-only stratum block: right after PCV_intercept, before its note.
  S1 <- matrix(0.4, 1, 1, dimnames = list("(Intercept)", "(Intercept)"))
  out_int <- p2_lines(p2_long_pcv(Sigma_stratum_null = S1,
                                  Sigma_stratum_adjusted = S1 / 10,
                                  stratum_slope = FALSE, pcv_slope = NA_real_,
                                  var_slope_null = NA_real_,
                                  var_slope_adjusted = NA_real_))
  j <- grep("PCV_intercept:", out_int, fixed = TRUE)
  expect_identical(out_int[j + 1], paste0("Variance basis: ", reml_label))
  expect_identical(out_int[j + 2], "")
  expect_match(out_int[j + 3], "random intercept only", fixed = TRUE)
})

test_that("print.maihda_long_pcv() leaves the basis to a note that names it", {
  # ML refit, mixed and brms: each note states the basis, so no line is added.
  cases <- list(
    ml_refit = list(p2_long_pcv(estimation = "ML", estimation_used = "ML",
                                ml_refit = TRUE, ml_variances = TRUE),
                    "refitted with maximum likelihood"),
    mixed = list(p2_long_pcv(estimation = "ML", estimation_used = "mixed",
                             ml_variances = TRUE),
                 "kept its REML fit"),
    brms = list(p2_long_pcv(estimation_used = "posterior",
                            variance_summary = "posterior median"),
                "posterior median of its per-draw value"))
  for (nm in names(cases)) {
    out <- p2_lines(cases[[nm]][[1]])
    expect_true(any(grepl(cases[[nm]][[2]], out, fixed = TRUE)), info = nm)
    expect_false(any(grepl("Variance basis", out, fixed = TRUE)), info = nm)
  }

  # A brms object stored without variance_summary holds posterior means: no
  # posterior-median claim on any line.
  out_brms_old <- p2_lines(p2_long_pcv(estimation_used = "posterior",
                                       variance_summary = NULL))
  expect_false(any(grepl("posterior", out_brms_old, fixed = TRUE)))

  # An object stored without estimation_used does not record its basis.
  out_old <- p2_lines(p2_long_pcv(estimation = NULL, estimation_used = NULL,
                                  ml_variances = NULL))
  expect_false(any(grepl("Variance basis", out_old, fixed = TRUE)))
  out_old_ml <- p2_lines(p2_long_pcv(estimation = NULL, estimation_used = NULL,
                                     ml_variances = NULL, ml_refit = TRUE))
  expect_false(any(grepl("Variance basis", out_old_ml, fixed = TRUE)))
  expect_true(any(grepl("refitted with maximum likelihood", out_old_ml, fixed = TRUE)))
})

test_that("maihda_longitudinal_pcv() records whether the variances are maximum likelihood", {
  d <- subset(maihda_long_data, id %in% unique(maihda_long_data$id)[1:150])
  a <- p2_keep_rng(suppressWarnings(suppressMessages(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education), data = d, id = "id",
           time = "wave", decomposition = "longitudinal"))))
  expect_false(a$pcv$ml_variances)
  out <- p2_lines(a)
  expect_identical(sum(grepl("Variance basis", out, fixed = TRUE)), 1L)
  expect_true(any(out == paste0("Variance basis: ", reml_label)))

  ml <- MAIHDA:::maihda_longitudinal_pcv(a$model, a$model_adjusted, estimation = "ML")
  expect_true(ml$ml_refit)
  expect_true(ml$ml_variances)
  expect_false(any(grepl("Variance basis", p2_lines(ml), fixed = TRUE)))
})
