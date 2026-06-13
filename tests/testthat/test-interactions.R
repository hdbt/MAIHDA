# Tests for maihda_interactions() -- the "which strata show significant
# interaction" diagnostic -- and the highlight_interactions plot option.

# Balanced 3x3 strata with a KNOWN, orthogonal (zero row/column margin) 2x2
# interaction block planted in the top-left: (1,1)=+d, (2,2)=+d, (1,2)=-d,
# (2,1)=-d, every other cell 0. Zero margins keep the additive main-effects fit
# unbiased, so the interaction BLUPs recover the planted pattern: the 4 block cells
# are clearly non-zero (flagged) and the other 5 are ~0 (not flagged) -- a clean,
# deterministic recovery. d = 0 plants no interaction (a near-singular adjusted
# fit). Additive A/B effects are large, so they sit in the adjusted model's fixed
# part and the stratum random effect isolates the interaction.
maihda_interaction_data <- function(seed = 123, n_per = 60, d = 2.5, sd = 1) {
  set.seed(seed)
  combos <- expand.grid(A = 1:3, B = 1:3)
  inter <- function(ai, bi) {
    if (ai == 1 && bi == 1) d
    else if (ai == 2 && bi == 2) d
    else if (ai == 1 && bi == 2) -d
    else if (ai == 2 && bi == 1) -d
    else 0
  }
  do.call(rbind, lapply(seq_len(nrow(combos)), function(k) {
    ai <- combos$A[k]
    bi <- combos$B[k]
    mu <- 2 * ai + 1.5 * bi + inter(ai, bi)
    data.frame(A = factor(ai), B = factor(bi),
               y = mu + stats::rnorm(n_per, 0, sd))
  }))
}

# The adjusted analysis (A, B as fixed main effects; A:B random intercept).
maihda_interaction_analysis <- function(...) {
  d <- maihda_interaction_data(...)
  suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = d)
  ))
}

test_that("maihda_interactions returns a classed table with the documented columns", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)

  expect_s3_class(mi, "maihda_interactions")
  expect_s3_class(mi, "data.frame")
  for (col in c("stratum", "label", "n", "interaction", "se", "lower", "upper",
                "p_value", "flagged", "direction")) {
    expect_true(col %in% names(mi))
  }
  # adjust = "none" => no p_adjusted column
  expect_false("p_adjusted" %in% names(mi))

  # Attributes are consistent with the flagged column.
  expect_equal(attr(mi, "n_strata"), nrow(mi))
  expect_equal(attr(mi, "n_flagged"), sum(mi$flagged))
  expect_identical(attr(mi, "scale"), "link")
  expect_identical(attr(mi, "engine"), "lme4")
  # 9 balanced strata of 60.
  expect_equal(nrow(mi), 9L)
  expect_true(all(mi$n == 60))
})

test_that("maihda_interactions recovers the planted orthogonal interaction block", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)

  # Exactly the 4 block cells are flagged; the 5 zero-interaction cells are not.
  expect_equal(attr(mi, "n_flagged"), 4L)

  # Rows are sorted flagged-first by |interaction|, so the top 4 are the block.
  expect_true(all(mi$flagged[1:4]))
  expect_false(any(mi$flagged[5:9]))

  # The block is sign-balanced (two +d, two -d cells).
  expect_equal(sum(mi$flagged & mi$direction == "above"), 2L)
  expect_equal(sum(mi$flagged & mi$direction == "below"), 2L)

  # Clean separation: the flagged interactions dwarf the unflagged ones.
  expect_gt(min(abs(mi$interaction[1:4])), max(abs(mi$interaction[5:9])))
})

test_that("multiplicity correction is monotone: bonferroni subset of BH subset of none", {
  a <- maihda_interaction_analysis()
  none <- maihda_interactions(a, adjust = "none")
  bh   <- maihda_interactions(a, adjust = "BH")
  bon  <- maihda_interactions(a, adjust = "bonferroni")

  flagged_set <- function(x) sort(x$stratum[x$flagged %in% TRUE])
  expect_true(all(flagged_set(bon) %in% flagged_set(bh)))
  expect_true(all(flagged_set(bh) %in% flagged_set(none)))

  # An adjustment adds the p_adjusted column.
  expect_true("p_adjusted" %in% names(bh))
})

test_that("a higher conf_level flags a subset of a lower one", {
  a <- maihda_interaction_analysis()
  wide   <- maihda_interactions(a, conf_level = 0.99)
  narrow <- maihda_interactions(a, conf_level = 0.90)

  flagged_set <- function(x) sort(x$stratum[x$flagged %in% TRUE])
  expect_true(all(flagged_set(wide) %in% flagged_set(narrow)))
})

test_that("an unknown adjust method errors", {
  a <- maihda_interaction_analysis()
  expect_error(maihda_interactions(a, adjust = "not-a-method"))
})

test_that("a bare null model warns about total deviation; an analysis does not", {
  d <- maihda_interaction_data()
  null_mod <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ (1 | A:B), data = d)
  ))
  expect_warning(maihda_interactions(null_mod), "null model")

  a <- maihda_interaction_analysis()
  # No guardrail warning when reading from a maihda() analysis (adjusted model).
  expect_warning(maihda_interactions(a), regexp = NA)
})

test_that("crossed-dimensions analyses use the interaction RE with no guardrail warning", {
  d <- maihda_interaction_data()
  cc <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = d, decomposition = "crossed-dimensions")
  ))
  mi <- suppressWarnings(maihda_interactions(cc))
  expect_s3_class(mi, "maihda_interactions")
  expect_identical(attr(mi, "model_type"), "crossed-dimensions")
  expect_warning(maihda_interactions(cc), regexp = NA)
})

test_that("print reports the flagged count and is exploratory", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)
  expect_output(print(mi), "strata flagged")
  # The default (no correction) steers toward BH.
  expect_output(print(mi), "BH")
})

test_that("a singular/boundary fit flags nothing and prints without error", {
  # Purely additive outcome (no interaction) -> between-stratum interaction
  # variance near zero; the BLUP SEs collapse, so nothing should flag.
  dat <- maihda_interaction_data(d = 0)
  a <- suppressMessages(suppressWarnings(maihda(y ~ A + B + (1 | A:B), data = dat)))
  mi <- maihda_interactions(a)
  expect_s3_class(mi, "maihda_interactions")
  expect_output(print(mi))
  if (isTRUE(attr(mi, "singular"))) {
    expect_equal(attr(mi, "n_flagged"), 0L)
    expect_output(print(mi), "singular")
  }
})

# ---- plot highlight option --------------------------------------------------

test_that("highlight_interactions = FALSE is unchanged; TRUE / object add a ring layer", {
  a <- maihda_interaction_analysis()

  p_off   <- plot(a, type = "effect_decomp")
  p_false <- plot(a, type = "effect_decomp", highlight_interactions = FALSE)
  p_on    <- plot(a, type = "effect_decomp", highlight_interactions = TRUE)

  expect_s3_class(p_off, "ggplot")
  expect_s3_class(p_on, "ggplot")
  # FALSE is identical to the default (no extra layer); TRUE adds the ring layer.
  expect_equal(length(p_false$layers), length(p_off$layers))
  expect_gt(length(p_on$layers), length(p_off$layers))
})

test_that("a precomputed maihda_interactions object can drive the highlight", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, adjust = "BH")
  p <- plot(a, type = "effect_decomp", highlight_interactions = mi)
  expect_s3_class(p, "ggplot")

  # Works on the predicted view too (routed to the null model, flags reused).
  p2 <- plot(a, type = "predicted", highlight_interactions = mi)
  expect_s3_class(p2, "ggplot")
})

test_that("an invalid highlight_interactions argument errors", {
  a <- maihda_interaction_analysis()
  expect_error(plot(a, type = "effect_decomp", highlight_interactions = "yes"),
               "maihda_interactions")
})

# ---- brms: exact posterior tail ---------------------------------------------

test_that("brms uses the exact posterior tail and ignores adjust", {
  # Compiles a Stan model, so OPT-IN (set MAIHDA_TEST_BRMS=true), matching the
  # other brms Stan tests -- the main R-CMD-check runners have brms installed but
  # no Boost/Stan toolchain, so skip_on_cran() is not enough.
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- maihda_interaction_data(n_per = 40)
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = d, engine = "brms",
           chains = 1, iter = 500, refresh = 0, seed = 1)
  ))

  mi <- suppressWarnings(maihda_interactions(a))
  expect_s3_class(mi, "maihda_interactions")
  # Bayesian columns: probability of direction, exact interval; no frequentist p.
  expect_true("pd" %in% names(mi))
  expect_false("p_value" %in% names(mi))
  expect_true(all(mi$pd >= 0 & mi$pd <= 1, na.rm = TRUE))
  expect_true(all(is.finite(mi$lower) & is.finite(mi$upper)))

  # adjust is inert for brms and says so.
  expect_message(maihda_interactions(a, adjust = "BH"), "ignored for brms")
  mi_bh <- suppressMessages(maihda_interactions(a, adjust = "BH"))
  expect_false("p_adjusted" %in% names(mi_bh))
})
