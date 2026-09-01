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

maihda_effect_decomp_label_data <- function(plot) {
  label_layers <- plot$layers[vapply(
    plot$layers,
    function(layer) inherits(layer$geom, "GeomLabelRepel"),
    logical(1)
  )]
  expect_length(label_layers, 1L)
  label_layers[[1]]$data
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
  # default adjust = "BH" => p_adjusted column present
  expect_true("p_adjusted" %in% names(mi))

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
  mi <- maihda_interactions(a, adjust = "none")  # pure recovery, independent of FDR

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

test_that("a non-syntactic dimension name does not trigger a false null-model warning", {
  # Regression: the adjusted-term detector compared RAW dimension names to terms()
  # labels, which backtick-quote a non-syntactic name like `gender var`, so a
  # fully-specified adjusted model was misread as a null model and warned.
  set.seed(7)
  n <- 300
  d <- data.frame(
    y = stats::rnorm(n),
    `gender var` = sample(c("F", "M"), n, TRUE),
    race = sample(c("A", "B", "C"), n, TRUE),
    check.names = FALSE
  )
  # Both additive main effects ARE in the fixed part -> a genuine adjusted model.
  adj <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ `gender var` + race + (1 | `gender var`:race), data = d)))
  expect_warning(maihda_interactions(adj), regexp = NA)

  # The true null model (no main effects) still warns.
  null_mod <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ (1 | `gender var`:race), data = d)))
  expect_warning(maihda_interactions(null_mod), "null model")
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
  # The default correction is now BH (FDR), named in the evidence line.
  expect_output(print(mi), "BH")
  # The frequentist p-values are labelled "conservative": a Wald tail on a shrunken
  # BLUP's conditional SE, stochastically large under a true null, so the BH flag
  # under-flags (not a calibrated p-value BH could repair).
  expect_output(print(mi), "conservative p-values")
})

# Subsetting (head / `[` / dplyr verbs) keeps the maihda_interactions class but, in
# many R/dplyr versions, drops the attributes the print method reads (n_flagged,
# n_strata, conf_level, ...). Whether they survive is version-dependent, so the
# tests reproduce the attribute-less state deterministically rather than relying on
# head() to drop it.
maihda_strip_interaction_attrs <- function(x) {
  for (nm in c("n_flagged", "n_strata", "conf_level", "adjust", "engine",
               "model_type", "rope", "singular")) {
    attr(x, nm) <- NULL
  }
  x
}

test_that("printing a maihda_interactions object stripped of its attributes does not error", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)

  bare <- maihda_strip_interaction_attrs(mi)
  expect_s3_class(bare, "maihda_interactions")     # still classed...
  expect_null(attr(bare, "n_flagged"))             # ...but the attribute is gone (the crash precondition)

  expect_no_error(print(bare))                     # previously errored on `if (NULL > 0)`
  # Counts are recomputed from the rows themselves.
  expect_output(print(bare),
                paste0(sum(bare$flagged), " of ", nrow(bare), " strata flagged"))

  # A stripped 3-row subset reports its own recomputed counts.
  bare3 <- maihda_strip_interaction_attrs(mi[1:3, ])
  expect_no_error(print(bare3))
  expect_output(print(bare3), paste0(sum(bare3$flagged), " of 3 strata flagged"))

  # Robust even when the flagged column itself was selected away.
  expect_no_error(print(maihda_strip_interaction_attrs(mi[, c("stratum", "label", "interaction")])))

  # The one-line summary helper shares the same guard.
  expect_no_error(MAIHDA:::maihda_print_interactions_line(bare))
})

test_that("subsetting a maihda_interactions object then printing never errors", {
  # Belt-and-suspenders over the real subsetting operations, whatever the running
  # R/dplyr version does with the attributes.
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)
  expect_no_error(print(head(mi, 3)))
  expect_no_error(print(mi[1:3, ]))
  expect_no_error(print(mi[mi$flagged, ]))
  expect_no_error(print(mi[, c("stratum", "label", "interaction")]))
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

test_that("highlight_interactions = FALSE is unchanged; TRUE applies the fade highlight", {
  a <- maihda_interaction_analysis()

  p_off   <- plot(a, type = "effect_decomp")
  p_false <- plot(a, type = "effect_decomp", highlight_interactions = FALSE)
  p_on    <- plot(a, type = "effect_decomp", highlight_interactions = TRUE)

  expect_s3_class(p_off, "ggplot")
  expect_s3_class(p_on, "ggplot")
  # Highlighting no longer adds a ring layer; it focuses by contrast (fade the
  # non-flagged), which keeps the layer count constant but adds the discrete alpha
  # scale and gives the segments more than one distinct opacity.
  expect_equal(length(p_false$layers), length(p_off$layers))
  expect_equal(length(p_on$layers), length(p_off$layers))
  expect_equal(length(p_false$scales$scales), length(p_off$scales$scales))
  expect_gt(length(p_on$scales$scales), length(p_off$scales$scales))

  seg_alpha_off <- unique(ggplot2::ggplot_build(p_off)$data[[2]]$alpha)
  seg_alpha_on  <- unique(ggplot2::ggplot_build(p_on)$data[[2]]$alpha)
  expect_length(seg_alpha_off, 1L)         # uniform opacity when not highlighting
  expect_gt(length(seg_alpha_on), 1L)       # flagged (solid) vs non-flagged (dimmed)
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

test_that("plot labels can follow multiplicity-adjusted interaction flags", {
  dat <- maihda_interaction_data(n_per = 30, d = 0.55)
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = dat)))

  none <- maihda_interactions(a, adjust = "none")
  bh <- maihda_interactions(a, adjust = "BH")
  expect_gt(attr(none, "n_flagged"), attr(bh, "n_flagged"))
  expect_gt(attr(bh, "n_flagged"), 0L)

  hl <- maihda_resolve_analysis_highlight(a, "BH")
  expect_s3_class(hl, "maihda_interactions")
  expect_identical(attr(hl, "adjust"), "BH")
  expect_equal(sort(hl$stratum[hl$flagged]), sort(bh$stratum[bh$flagged]))

  p <- plot(a, type = "effect_decomp", highlight_interactions = "BH")
  label_data <- maihda_effect_decomp_label_data(p)
  expect_equal(sort(as.character(label_data$stratum)),
               sort(as.character(bh$stratum[bh$flagged])))
  expect_true(all(label_data$.maihda_flag))

  dat_null <- maihda_interaction_data(n_per = 30, d = 0)
  a_null <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = dat_null)))
  bh_null <- maihda_interactions(a_null, adjust = "BH")
  expect_equal(attr(bh_null, "n_flagged"), 0L)

  p_null <- plot(a_null, type = "effect_decomp", highlight_interactions = "BH")
  label_data_null <- maihda_effect_decomp_label_data(p_null)
  expect_equal(nrow(label_data_null), 0L)
})

test_that("an invalid highlight_interactions argument errors", {
  a <- maihda_interaction_analysis()
  expect_error(plot(a, type = "effect_decomp", highlight_interactions = "yes"),
               "maihda_interactions")
})

# ---- only_flagged: filter to the flagged strata -----------------------------

test_that("only_flagged restricts the predicted view to the flagged strata", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)
  flagged <- as.character(mi$stratum[mi$flagged])
  expect_gt(length(flagged), 0L)

  p_all  <- plot(a, type = "predicted")
  p_only <- plot(a, type = "predicted", only_flagged = TRUE)

  expect_s3_class(p_only, "ggplot")
  # exactly the flagged strata, nothing else
  expect_setequal(as.character(p_only$data$stratum), flagged)
  expect_equal(nrow(p_only$data), length(flagged))
  expect_lt(nrow(p_only$data), nrow(p_all$data))
  # the cap can no longer hide a flagged stratum: every flagged point is drawn
  expect_true(all(p_only$data$.maihda_flag))
  # the caption names the screen honestly
  expect_match(p_only$labels$caption, "flagged strata")
  expect_match(p_only$labels$caption, "BH-adjusted")
})

test_that("only_flagged restricts the obs_vs_shrunken view too", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)
  flagged <- as.character(mi$stratum[mi$flagged])

  p <- plot(a, type = "obs_vs_shrunken", only_flagged = TRUE)
  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(p$data$stratum), flagged)
  expect_match(p$labels$caption, "flagged strata")
})

test_that("only_flagged returns a captioned empty panel when nothing is flagged", {
  d0 <- maihda_interaction_data(d = 0)
  a0 <- suppressMessages(suppressWarnings(maihda(y ~ A + B + (1 | A:B), data = d0)))
  expect_equal(attr(maihda_interactions(a0), "n_flagged"), 0L)

  p <- plot(a0, type = "predicted", only_flagged = TRUE)
  expect_s3_class(p, "ggplot")
  # the title is preserved and the panel explains why it is empty
  expect_match(p$labels$title, "Predicted Subgroup Values")
  bdata <- ggplot2::ggplot_build(p)$data[[1]]
  expect_true(any(grepl("No strata flagged", bdata$label)))
})

# ---- flag-aware truncation: the cap never drops a flagged stratum -----------

test_that("the n_strata cap is flag-aware: flagged strata survive truncation", {
  a <- maihda_interaction_analysis()          # 9 strata, 4 flagged
  mi <- maihda_interactions(a)
  flagged <- as.character(mi$stratum[mi$flagged])

  # cap BELOW the flagged count, highlighting on but not only_flagged
  p <- plot(a, type = "predicted", n_strata = 2, highlight_interactions = TRUE)
  shown <- as.character(p$data$stratum)
  # every flagged stratum is kept even though n_strata = 2
  expect_true(all(flagged %in% shown))
  expect_gte(nrow(p$data), length(flagged))
  expect_match(p$labels$caption, "exceeded to keep every flagged stratum")
})

test_that("flag-aware truncation fills the remaining slots with non-flagged strata", {
  a <- maihda_interaction_analysis()          # 9 strata, 4 flagged
  mi <- maihda_interactions(a)
  flagged <- as.character(mi$stratum[mi$flagged])

  p <- plot(a, type = "predicted", n_strata = 6, highlight_interactions = TRUE)
  shown <- as.character(p$data$stratum)
  expect_true(all(flagged %in% shown))        # all flagged kept
  expect_equal(length(shown), 6L)             # filled up to the cap
  expect_match(p$labels$caption, "plus the first")
})

test_that("the unhighlighted cap is unchanged (plain head, captioned)", {
  set.seed(202)
  dat <- data.frame(stratum = rep(1:20, each = 10),
                    age = rnorm(200), outcome = rnorm(200))
  m <- fit_maihda(outcome ~ age + (1 | stratum), data = dat, engine = "lme4")
  p <- plot(m, type = "predicted", n_strata = 10)
  expect_equal(nrow(p$data), 10L)
  expect_match(p$labels$caption, "first 10 of 20")
})

# ---- effect_decomp stays highlighted in context -----------------------------

test_that("only_flagged does not filter effect_decomp but says so", {
  a <- maihda_interaction_analysis()
  expect_message(
    p <- plot(a, type = "effect_decomp", only_flagged = TRUE),
    "only_flagged"
  )
  expect_s3_class(p, "ggplot")
})

test_that("an invalid only_flagged argument errors", {
  a <- maihda_interaction_analysis()
  expect_error(plot(a, type = "predicted", only_flagged = "yes"),
               "single TRUE or FALSE")
})

# ---- select composes with the flag-aware cap --------------------------------

test_that("select governs the flag-aware fill while every flagged stratum is kept", {
  a <- maihda_interaction_analysis()          # 9 strata, 4 flagged
  mi <- maihda_interactions(a)
  flagged <- as.character(mi$stratum[mi$flagged])

  p <- plot(a, type = "predicted", n_strata = 6,
            highlight_interactions = TRUE, select = "deviation")
  shown <- as.character(p$data$stratum)
  expect_true(all(flagged %in% shown))        # flagged always survive
  expect_equal(length(shown), 6L)             # filled to the cap
  expect_match(p$labels$caption, "most extreme")
})

test_that("select picks which flagged strata survive an only_flagged cap", {
  a <- maihda_interaction_analysis()          # 4 flagged
  p <- plot(a, type = "predicted", only_flagged = TRUE,
            n_strata = 2, select = "deviation")
  expect_equal(nrow(p$data), 2L)
  expect_match(p$labels$caption, "most extreme of 4 flagged")
})

# ---- probability of direction (Stan-free) -----------------------------------

test_that("maihda_pd is the conventional probability of direction in [0.5, 1]", {
  # A strong NEGATIVE effect has pd ~ 1 (its direction is near-certain), NOT ~ 0.
  # The old code reported mean(d > 0), which would read ~ 0 here -- the bug.
  d_neg <- rnorm(4000, mean = -3, sd = 1)
  expect_gt(MAIHDA:::maihda_pd(d_neg), 0.99)
  expect_equal(MAIHDA:::maihda_pd(d_neg),
               max(mean(d_neg > 0), mean(d_neg < 0)))

  # A strong positive effect is symmetric (also ~ 1); a null effect ~ 0.5.
  expect_gt(MAIHDA:::maihda_pd(rnorm(4000, mean = 3)), 0.99)
  expect_equal(MAIHDA:::maihda_pd(rnorm(20000, mean = 0)), 0.5, tolerance = 0.05)

  # Always in [0.5, 1]; NA for an all-non-finite vector.
  expect_true(MAIHDA:::maihda_pd(rnorm(500, 0.3)) >= 0.5)
  expect_true(is.na(MAIHDA:::maihda_pd(c(NA_real_, NaN, Inf))))
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
  # seed = does not reach brms::brm() from maihda() (it is the response-scale VPC
  # simulation seed, and maihda() now warns when combined with engine = "brms");
  # set.seed() pins the sampler instead.
  set.seed(20260831)
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = d, engine = "brms",
           chains = 1, iter = 500, refresh = 0)
  ))

  mi <- suppressWarnings(maihda_interactions(a))
  expect_s3_class(mi, "maihda_interactions")
  # Bayesian columns: probability of direction, exact interval; no frequentist p.
  expect_true("pd" %in% names(mi))
  expect_false("p_value" %in% names(mi))
  # Conventional probability of direction: always in [0.5, 1], never the mislabelled
  # P(BLUP > 0) that would dip below 0.5 for strata with negative interactions.
  expect_true(all(mi$pd >= 0.5 & mi$pd <= 1, na.rm = TRUE))
  expect_true(all(is.finite(mi$lower) & is.finite(mi$upper)))

  # adjust is inert for brms and says so.
  expect_message(maihda_interactions(a, adjust = "BH"), "ignored for brms")
  mi_bh <- suppressMessages(maihda_interactions(a, adjust = "BH"))
  expect_false("p_adjusted" %in% names(mi_bh))
})

# --- the interaction diagnostic built into maihda() / fit_maihda() ------------

test_that("maihda() attaches the interaction diagnostic by default", {
  a <- maihda_interaction_analysis()
  expect_s3_class(a$interactions, "maihda_interactions")
  expect_identical(attr(a$interactions, "adjust"), "BH")  # FDR default for an all-strata scan
  # identical to calling the diagnostic directly on the analysis (no recompute drift)
  direct <- maihda_interactions(a)
  expect_equal(attr(a$interactions, "n_flagged"), attr(direct, "n_flagged"))
  expect_equal(a$interactions$interaction, direct$interaction)
})

test_that("interactions = FALSE skips the diagnostic", {
  d <- maihda_interaction_data()
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = d, interactions = FALSE)))
  expect_null(a$interactions)
})

test_that("interactions = 'none' overrides the FDR default", {
  d <- maihda_interaction_data()
  a <- suppressMessages(suppressWarnings(
    maihda(y ~ A + B + (1 | A:B), data = d, interactions = "none")))
  expect_s3_class(a$interactions, "maihda_interactions")
  expect_identical(attr(a$interactions, "adjust"), "none")
  expect_false("p_adjusted" %in% names(a$interactions))
})

test_that("an invalid interactions argument errors", {
  d <- maihda_interaction_data()
  expect_error(
    suppressMessages(suppressWarnings(
      maihda(y ~ A + B + (1 | A:B), data = d, interactions = 1L))),
    "TRUE, FALSE, or a multiple-comparison method")
})

test_that("plot(highlight_interactions = TRUE) reuses the stored diagnostic", {
  a <- maihda_interaction_analysis()
  hl <- maihda_resolve_analysis_highlight(a, TRUE)
  expect_identical(hl, a$interactions)
})

test_that("fit_maihda(interactions = ) is opt-in and parallels maihda()", {
  d <- maihda_interaction_data()
  # default off
  m_off <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ A + B + (1 | A:B), data = d)))
  expect_null(m_off$interactions)
  # opt-in, on the adjusted model
  m_on <- suppressMessages(suppressWarnings(
    fit_maihda(y ~ A + B + (1 | A:B), data = d, interactions = TRUE)))
  expect_s3_class(m_on$interactions, "maihda_interactions")
  # on a null model it warns (the stratum RE is the total deviation, not the
  # pure interaction the diagnostic claims)
  expect_warning(
    fit_maihda(y ~ 1 + (1 | A:B), data = d, interactions = TRUE),
    "looks like a null model")
})

# ---- equivalence / ROPE reading ---------------------------------------------

test_that("rope adds an equivalence decision column", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, rope = 0.5)

  expect_true("decision" %in% names(mi))
  expect_true(all(mi$decision %in% c("relevant", "negligible", "inconclusive") |
                    is.na(mi$decision)))
  expect_equal(attr(mi, "rope"), c(-0.5, 0.5))
  # rows are sorted by |interaction|; the strongest planted block cell sits well
  # outside a small ROPE, so it is 'relevant'.
  expect_identical(mi$decision[1], "relevant")
  # a near-zero stratum's interval sits inside the ROPE -> 'negligible'
  expect_true(any(mi$decision == "negligible", na.rm = TRUE))
})

test_that("a single-number rope is the symmetric region; c(lo, hi) passes through", {
  a <- maihda_interaction_analysis()
  expect_equal(attr(maihda_interactions(a, rope = 0.4), "rope"), c(-0.4, 0.4))
  expect_equal(attr(maihda_interactions(a, rope = c(-0.2, 0.6)), "rope"), c(-0.2, 0.6))
})

test_that("invalid rope arguments error", {
  a <- maihda_interaction_analysis()
  expect_error(maihda_interactions(a, rope = -1), "positive")
  expect_error(maihda_interactions(a, rope = c(0.5, -0.5)), "lower < upper")
  expect_error(maihda_interactions(a, rope = c(1, 2, 3)), "length")
  expect_error(maihda_interactions(a, rope = "wide"), "NULL")
})

test_that("rope = NULL (default) adds no decision column and prints normally", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a)
  expect_false("decision" %in% names(mi))
  expect_null(attr(mi, "rope"))
})

test_that("rope print reports the equivalence breakdown", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, rope = 0.5)
  expect_output(print(mi), "Equivalence vs ROPE")
})

# ---- highlight_by = "rope": highlight the ROPE-relevant strata ---------------

test_that("highlight_by = 'rope' highlights exactly the decision == 'relevant' strata", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, rope = 0.5)
  relevant <- sort(as.character(mi$stratum[mi$decision == "relevant"]))
  expect_gt(length(relevant), 0L)

  # predicted view, analysis path
  p <- plot(a, type = "predicted", highlight_interactions = TRUE,
            highlight_by = "rope", rope = 0.5)
  expect_setequal(sort(as.character(p$data$stratum[p$data$.maihda_flag])), relevant)

  # obs_vs_shrunken applies the same set (routed to the null model, flags reused)
  po <- plot(a, type = "obs_vs_shrunken", highlight_interactions = TRUE,
             highlight_by = "rope", rope = 0.5)
  expect_setequal(sort(as.character(po$data$stratum[po$data$.maihda_flag])), relevant)

  # model path (the adjusted model directly) gives the same set
  pm <- plot(a$model_adjusted, type = "predicted", highlight_interactions = TRUE,
             highlight_by = "rope", rope = 0.5)
  expect_setequal(sort(as.character(pm$data$stratum[pm$data$.maihda_flag])), relevant)
})

test_that("highlight_by = 'rope' marks the ROPE-relevant strata on effect_decomp", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, rope = 0.5)
  relevant <- sort(as.character(mi$stratum[mi$decision == "relevant"]))

  p <- plot(a, type = "effect_decomp", highlight_interactions = TRUE,
            highlight_by = "rope", rope = 0.5)
  ld <- maihda_effect_decomp_label_data(p)
  expect_setequal(sort(as.character(ld$stratum)), relevant)
  expect_true(all(ld$.maihda_flag))
})

test_that("highlight_by switches the set: a wide ROPE flags nothing where the zero-test flags four", {
  a <- maihda_interaction_analysis()
  # The zero-centred flag marks the 4 planted block cells...
  expect_equal(sum(maihda_interactions(a)$flagged), 4L)
  # ...but a wide equivalence region (|effect| > 2.4) classifies none "relevant".
  mi_rope <- maihda_interactions(a, rope = 2.4)
  expect_equal(sum(mi_rope$decision == "relevant", na.rm = TRUE), 0L)

  # Flag mode still highlights four; rope mode highlights none.
  p_flag <- plot(a, type = "predicted", highlight_interactions = TRUE)
  expect_equal(sum(p_flag$data$.maihda_flag), 4L)

  # only_flagged + a wide ROPE -> a captioned empty panel with ROPE-specific wording
  # (the flag-mode panel would have kept four strata).
  p_rope <- plot(a, type = "predicted", highlight_by = "rope", rope = 2.4,
                 only_flagged = TRUE)
  expect_s3_class(p_rope, "ggplot")
  expect_match(p_rope$labels$title, "Predicted Subgroup Values")
  bdata <- ggplot2::ggplot_build(p_rope)$data[[1]]
  expect_true(any(grepl("ROPE-relevant", bdata$label)))
})

test_that("only_flagged + highlight_by = 'rope' restricts to the ROPE-relevant strata", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, rope = 0.5)
  relevant <- sort(as.character(mi$stratum[mi$decision == "relevant"]))

  p <- plot(a, type = "predicted", highlight_by = "rope", rope = 0.5, only_flagged = TRUE)
  expect_setequal(as.character(p$data$stratum), relevant)
  expect_true(all(p$data$.maihda_flag))
  # the caption names the ROPE screen honestly, not "flagged"/"BH-adjusted"
  expect_match(p$labels$caption, "ROPE-relevant strata")
  expect_match(p$labels$caption, "latent \\(link\\) scale")
})

test_that("a precomputed maihda_interactions(rope=) drives the ROPE highlight without recomputing", {
  a <- maihda_interaction_analysis()
  mi <- maihda_interactions(a, rope = 0.5)
  relevant <- sort(as.character(mi$stratum[mi$decision == "relevant"]))

  # No `rope` passed here: the decision column already rides on the object.
  p <- plot(a, type = "predicted", highlight_interactions = mi, highlight_by = "rope")
  expect_s3_class(p, "ggplot")
  expect_setequal(sort(as.character(p$data$stratum[p$data$.maihda_flag])), relevant)

  p2 <- plot(a, type = "effect_decomp", highlight_interactions = mi, highlight_by = "rope")
  expect_s3_class(p2, "ggplot")
})

test_that("highlight_by = 'rope' without a usable rope errors with an actionable message", {
  a <- maihda_interaction_analysis()
  # analysis path (single type -> the error surfaces rather than being caught)
  expect_error(
    plot(a, type = "predicted", highlight_interactions = TRUE, highlight_by = "rope"),
    "no 'decision' column")
  # model path
  expect_error(
    plot(a$model_adjusted, type = "predicted", highlight_interactions = TRUE,
         highlight_by = "rope"),
    "rope")
  # a precomputed screen built WITHOUT a rope is rejected the same way
  mi_norope <- maihda_interactions(a)
  expect_error(
    plot(a, type = "predicted", highlight_interactions = mi_norope, highlight_by = "rope"),
    "no 'decision' column")
})

test_that("highlight_by defaults to 'flag' and is validated", {
  a <- maihda_interaction_analysis()
  d_default <- plot(a, type = "predicted", highlight_interactions = TRUE)
  d_flag <- plot(a, type = "predicted", highlight_interactions = TRUE, highlight_by = "flag")
  expect_equal(sort(as.character(d_default$data$stratum[d_default$data$.maihda_flag])),
               sort(as.character(d_flag$data$stratum[d_flag$data$.maihda_flag])))
  expect_error(
    plot(a, type = "predicted", highlight_interactions = TRUE, highlight_by = "nope"))
  expect_error(
    plot(a$model_adjusted, type = "predicted", highlight_by = "nope"))
})

# --- link-scale note (2026-08-29) -------------------------------------------
# The stratum BLUP is a departure in LINK units. On a logit fit that makes it a
# log-odds (multiplicative) interaction -- a different quantity from the
# probability-scale one Evans et al. (2024, sec. 2.5.1) define as pi_j - pi^A_j --
# so the result carries the link and print() names it. On an identity link the BLUP
# is already in the outcome's units and nothing extra is said.

test_that("maihda_interactions() records the model's link", {
  a <- maihda_interaction_analysis()
  expect_identical(attr(maihda_interactions(a), "link"), "identity")

  d <- maihda_interaction_data()
  d$yb <- as.integer(d$y > stats::median(d$y))
  m <- suppressMessages(suppressWarnings(
    fit_maihda(yb ~ A + B + (1 | A:B), data = d, family = "binomial")
  ))
  expect_identical(attr(maihda_interactions(m), "link"), "logit")
})

test_that("print() names the link-scale departure only for a non-identity link", {
  a <- maihda_interaction_analysis()
  gaussian_out <- paste(utils::capture.output(print(maihda_interactions(a))),
                        collapse = "\n")
  expect_false(grepl("MULTIPLICATIVE", gaussian_out, fixed = TRUE))

  d <- maihda_interaction_data()
  d$yb <- as.integer(d$y > stats::median(d$y))
  m <- suppressMessages(suppressWarnings(
    fit_maihda(yb ~ A + B + (1 | A:B), data = d, family = "binomial")
  ))
  binomial_out <- paste(utils::capture.output(print(maihda_interactions(m))),
                        collapse = "\n")
  expect_true(grepl("log-odds departure (logit link)", binomial_out, fixed = TRUE))
  expect_true(grepl("additivity on the link", binomial_out, fixed = TRUE))
  expect_true(grepl("not in risks", binomial_out, fixed = TRUE))
  # The surviving departure is called "risks", never "the probability scale": that
  # phrase would collide with the probability-POINTS restatement two lines below and
  # reintroduce the very conflation the note exists to prevent.
  expect_false(grepl("probability scale", binomial_out, fixed = TRUE))
  # Additivity is claimed for the LINK scale, never for the odds -- when log-odds
  # add, odds multiply, so "odds combine additively" would be an error on screen.
  expect_false(grepl("odds combine", binomial_out, ignore.case = TRUE))
  # Nor may the note claim risks NEVER add. The risk-difference interaction is
  # positive below the curve's midpoint and negative above it, so it crosses zero:
  # with effects a and b it vanishes exactly at an intercept of -(a + b)/2, where the
  # two comparisons straddle the midpoint symmetrically. Pinned below.
  expect_false(grepl("never add", binomial_out, fixed = TRUE))
})

test_that("risks CAN add exactly under a logit link (so the note must not say never)", {
  rd_int <- function(b0, a, b) {
    (stats::plogis(b0 + a + b) - stats::plogis(b0 + a)) -
      (stats::plogis(b0 + b) - stats::plogis(b0))
  }
  a <- 0.7; b <- 0.7

  # Vanishes exactly at the symmetric straddle ...
  expect_equal(rd_int(-(a + b) / 2, a, b), 0)
  # ... and is non-zero, with opposite signs, on either side of it.
  expect_gt(rd_int(-2, a, b), 0)
  expect_lt(rd_int(0, a, b), 0)
})

test_that("the link note is silent when the link is unknown", {
  expect_null(maihda_interaction_link_note(NA_character_))
  expect_null(maihda_interaction_link_note(character(0)))
  expect_null(maihda_interaction_link_note("identity"))
  expect_match(maihda_interaction_link_note("probit"), "not in risks")
  expect_match(maihda_interaction_link_note("log"), "not in counts")
})

# --- response-scale interaction, pi^B_j (2026-08-29) ------------------------
# Evans et al. (2024, sec. 2.5.1): pi^B_j = pi_j - pi^A_j, the stratum's total
# predicted outcome minus the outcome implied by its additive main effects alone.
# The map from the BLUP is strictly increasing through zero, so the evidence
# (flags, direction, p-values) must be invariant while the units change.

maihda_binomial_interaction_model <- function() {
  d <- maihda_interaction_data()
  d$yb <- as.integer(d$y > stats::median(d$y))
  suppressMessages(suppressWarnings(
    fit_maihda(yb ~ A + B + (1 | A:B), data = d, family = "binomial")
  ))
}

test_that("scale = 'response' reproduces the Evans et al. (2024) definition", {
  m <- maihda_binomial_interaction_model()
  link <- maihda_interactions(m)
  resp <- maihda_interactions(m, scale = "response")

  # pi^B_j computed straight from the definition, averaged over the stratum's own
  # rows (which reduces to Evans' constant-x_j form for this dimension-only model).
  eta <- stats::predict(m$model, re.form = NA, type = "link")
  key <- as.character(m$data$stratum)
  hand <- vapply(as.character(resp$stratum), function(s) {
    e <- eta[key == s]
    u <- link$interaction[match(s, as.character(link$stratum))]
    mean(stats::plogis(e + u) - stats::plogis(e))
  }, numeric(1))
  expect_equal(resp$interaction, unname(hand))

  # ... and the interval endpoints are the same map applied to the same endpoints.
  hand_lo <- vapply(as.character(resp$stratum), function(s) {
    e <- eta[key == s]
    lo <- link$lower[match(s, as.character(link$stratum))]
    mean(stats::plogis(e + lo) - stats::plogis(e))
  }, numeric(1))
  expect_equal(resp$lower, unname(hand_lo))
})

test_that("the response scale changes the units, not the evidence", {
  m <- maihda_binomial_interaction_model()
  link <- maihda_interactions(m)
  resp <- maihda_interactions(m, scale = "response")

  i <- match(as.character(resp$stratum), as.character(link$stratum))
  expect_identical(resp$flagged, link$flagged[i])
  expect_identical(resp$direction, link$direction[i])
  expect_equal(resp$p_value, link$p_value[i])
  expect_equal(resp$p_adjusted, link$p_adjusted[i])
  expect_identical(sign(resp$interaction), sign(link$interaction[i]))
  # A probability difference is bounded by 1; the log-odds departure here is not.
  expect_true(all(abs(resp$interaction) < 1))
  expect_true(max(abs(link$interaction)) > 1)

  # se is link-scale only: the response interval is asymmetric about the estimate.
  expect_true("se" %in% names(link))
  expect_false("se" %in% names(resp))
  expect_identical(attr(resp, "scale"), "response")
  expect_false(isTRUE(all.equal(resp$interaction - resp$lower,
                                resp$upper - resp$interaction)))
})

test_that("an identity link returns the same numbers on either scale", {
  a <- maihda_interaction_analysis()
  link <- maihda_interactions(a)
  resp <- maihda_interactions(a, scale = "response")
  i <- match(as.character(resp$stratum), as.character(link$stratum))
  expect_equal(resp$interaction, link$interaction[i])
  expect_equal(resp$lower, link$lower[i])
  expect_equal(resp$upper, link$upper[i])
})

test_that("the ROPE is read on the requested scale", {
  m <- maihda_binomial_interaction_model()
  # 0.02 probability points is a tight region the log-odds BLUPs would clear
  # trivially, so the two scales must classify differently.
  resp <- maihda_interactions(m, scale = "response", rope = 0.02)
  link <- maihda_interactions(m, rope = 0.02)
  expect_identical(attr(resp, "rope"), c(-0.02, 0.02))
  expect_false(identical(sort(resp$decision), sort(link$decision)))
  # ... and the plot caption names the scale it was screened on.
  expect_match(maihda_highlight_screen_label(
    structure(1, highlight_by = "rope", rope = c(-0.02, 0.02), scale = "response")),
    "response (outcome) scale", fixed = TRUE)
})

test_that("print() names the response scale and the quantity", {
  m <- maihda_binomial_interaction_model()
  out <- paste(utils::capture.output(print(maihda_interactions(m, scale = "response"))),
               collapse = "\n")
  expect_match(out, "interaction on the response (outcome) scale", fixed = TRUE)
  expect_match(out, "pi_j - pi^A_j", fixed = TRUE)
  expect_false(grepl("MULTIPLICATIVE", out, fixed = TRUE))
})

test_that("scale is validated", {
  a <- maihda_interaction_analysis()
  expect_error(maihda_interactions(a, scale = "probability"))
})

# --- printed outcome-scale column (2026-08-29) ------------------------------
# print() shows pi^B_j beside the link-scale BLUP on a non-identity link, so the
# reader never sees the log-odds figure alone. It is DISPLAY ONLY: the returned
# columns must be exactly what they were before, with scale = "response" the way to
# get the quantity (and its interval) as data.

test_that("print() shows the outcome-scale column without changing the object", {
  m <- maihda_binomial_interaction_model()
  x <- maihda_interactions(m)

  expect_identical(
    names(x),
    c("stratum", "label", "n", "interaction", "se", "lower", "upper",
      "p_value", "p_adjusted", "flagged", "direction"))
  expect_false("prob_diff" %in% names(x))

  out <- paste(utils::capture.output(print(x)), collapse = "\n")
  expect_match(out, "prob_diff", fixed = TRUE)
  expect_match(out, "probability points", fixed = TRUE)
})

test_that("the printed column is the same number scale = 'response' returns", {
  m <- maihda_binomial_interaction_model()
  x <- maihda_interactions(m)
  resp <- maihda_interactions(m, scale = "response")

  map <- attr(x, "response_interaction")
  expect_false(is.null(map))
  expect_equal(unname(map[as.character(resp$stratum)]), resp$interaction)
})

test_that("no outcome-scale column on an identity link or already-response result", {
  a <- maihda_interaction_analysis()
  expect_null(attr(maihda_interactions(a), "response_interaction"))
  expect_false(grepl("resp_diff|prob_diff",
                     paste(utils::capture.output(print(maihda_interactions(a))),
                           collapse = "\n")))

  m <- maihda_binomial_interaction_model()
  # Already reported on the response scale: the estimate IS pi^B, so no second copy.
  expect_null(attr(maihda_interactions(m, scale = "response"), "response_interaction"))
  resp_out <- paste(
    utils::capture.output(print(maihda_interactions(m, scale = "response"))),
    collapse = "\n")
  expect_false(grepl("prob_diff", resp_out, fixed = TRUE))
})

test_that("the display helper is keyed by stratum and degrades quietly", {
  df <- data.frame(stratum = c("b", "a"), interaction = c(2, 1), n = c(5L, 6L),
                   stringsAsFactors = FALSE)
  map <- c(a = 0.1, b = 0.2)

  out <- maihda_interactions_add_response_col(df, map, "probability")
  # Matched by name, not by position, so a reordered frame still lines up ...
  expect_identical(out$prob_diff, c(0.2, 0.1))
  # ... and it sits next to the quantity it restates.
  expect_identical(names(out), c("stratum", "interaction", "prob_diff", "n"))

  expect_identical(maihda_interactions_add_response_col(df, NULL, "probability"), df)
  expect_identical(
    maihda_interactions_add_response_col(df, c(z = 1), "probability"), df)
  expect_identical(
    maihda_interactions_add_response_col(df[0, ], map, "probability"), df[0, ])
  expect_identical(
    names(maihda_interactions_add_response_col(df, map, "count")),
    c("stratum", "interaction", "count_diff", "n"))
  expect_identical(
    names(maihda_interactions_add_response_col(df, map, "score")),
    c("stratum", "interaction", "score_diff", "n"))
})


# --- response-scale labelling follows the FAMILY, not the link (2026-08-30) ---
# A cumulative (ordinal) fit is LOGIT-linked but its response scale is the expected
# category score, not a probability. Labelling off the link alone called that
# "prob_diff ... probability points", which was simply wrong.

test_that("the response kind is read from the family, not the link alone", {
  expect_identical(maihda_response_kind_from_link("logit"), "probability")
  expect_identical(maihda_response_kind_from_link("log"), "count")
  expect_identical(maihda_response_kind_from_link("wibble"), "response")

  expect_identical(maihda_interaction_response_header("probability"), "prob_diff")
  expect_identical(maihda_interaction_response_header("count"), "count_diff")
  expect_identical(maihda_interaction_response_header("score"), "score_diff")
  expect_identical(maihda_interaction_response_header("anything else"), "resp_diff")
})

test_that("an ordinal fit is never labelled as probabilities", {
  skip_if_not_installed("ordinal")
  d <- maihda_interaction_data()
  d$ord <- factor(cut(d$y, 3), ordered = TRUE)
  m <- suppressMessages(suppressWarnings(
    fit_maihda(ord ~ A + B + (1 | A:B), data = d)
  ))
  skip_if_not(identical(m$engine, "ordinal"))

  x <- suppressWarnings(maihda_interactions(m))
  expect_identical(attr(x, "response_kind"), "score")
  # The DEPARTURE is still log-odds -- that part the link does settle.
  expect_identical(attr(x, "link"), "logit")

  out <- paste(utils::capture.output(print(x)), collapse = "
")
  expect_match(out, "log-odds departure (logit link)", fixed = TRUE)
  expect_match(out, "not in expected scores", fixed = TRUE)
  expect_false(grepl("prob_diff", out, fixed = TRUE))
  expect_false(grepl("probability points", out, fixed = TRUE))
  expect_false(grepl("not in risks", out, fixed = TRUE))
  # pi is Evans et al.'s probability notation and has no business on a score scale.
  expect_false(grepl("pi_j - pi^A_j", out, fixed = TRUE))
})

test_that("the pi notation appears only where the response really is a probability", {
  probability <- maihda_interaction_link_note("logit", has_response = TRUE,
                                              kind = "probability")
  expect_match(probability, "pi_j - pi^A_j", fixed = TRUE)
  for (k in c("count", "score", "response")) {
    note <- maihda_interaction_link_note("logit", has_response = TRUE, kind = k)
    expect_false(grepl("pi_j - pi^A_j", note, fixed = TRUE))
  }
  expect_match(maihda_interaction_link_note("logit", "response", kind = "score"),
               "expected-score difference", fixed = TRUE)
  expect_match(maihda_interaction_link_note("logit", "response", kind = "count"),
               "expected-count difference", fixed = TRUE)
})
