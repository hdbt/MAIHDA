# Audit 2026-10-13 (packet 2026-10-13-interactions-singular-note): when print() of a
# maihda_interactions() table says that the stratum standard errors collapse --
# CONFIRMED, and wider than reported. No number, flag or column changes.
#
# F1. The note "singular/boundary fit (between-stratum variance ~ 0); the BLUP SEs
#     collapse toward zero, so 'no flag' is not evidence of no interaction." printed
#     whenever the fit was singular in ANY random-effect term. Those standard errors
#     collapse only when the term the stratum effects come from is at zero. With a
#     zero in another term the note was false: beside standard errors of 5.4 to 5.6
#     (maihda_country_data, Japan, crossed-dimensions, gender at zero), and beside 16
#     of 20 strata FLAGGED (a two-model analysis whose only zero is a context). For a
#     crossed-dimensions fit whose interaction is at zero the bracket was wrong: the
#     between-stratum variance there is the dimension terms plus the interaction
#     (1034.5 for Finland).
#
# Per the user: the note prints only when the term the interaction is read from is at
# zero by the fit's record ('stratum', or the interaction term of a crossed-dimensions
# fit), and a record without a list of terms keeps it; a crossed-dimensions fit says
# "interaction variance"; attr(, "singular") keeps meaning a zero in any term, and the
# new attribute interaction_at_boundary records the note's condition.
#
# No block draws a random number: the fixtures are built from fixed sequences.

is13_quiet <- function(expr) suppressMessages(suppressWarnings(expr))

is13_head_crossed <- paste0("singular/boundary fit (interaction variance ~ 0); ",
                            "the BLUP SEs")
is13_head_other <- paste0("singular/boundary fit (between-stratum variance ~ 0); ",
                          "the BLUP SEs")
is13_tail <- "collapse toward zero, so 'no flag' is not evidence of no interaction."

is13_print <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")

# A two-row interactions table with the attributes print() reads.
is13_table <- function(model_type, singular, at_boundary) {
  out <- data.frame(stratum = c("1", "2"), label = c("a1 x b1", "a1 x b2"),
                    n = c(30L, 30L), interaction = c(0, 0), se = c(0, 0),
                    lower = c(0, 0), upper = c(0, 0),
                    p_value = c(NA_real_, NA_real_), p_adjusted = c(NA_real_, NA_real_),
                    flagged = c(FALSE, FALSE), direction = c("above", "above"),
                    stringsAsFactors = FALSE)
  attr(out, "conf_level") <- 0.95
  attr(out, "adjust") <- "BH"
  attr(out, "engine") <- "lme4"
  attr(out, "model_type") <- model_type
  attr(out, "n_strata") <- 2L
  attr(out, "n_flagged") <- 0L
  attr(out, "scale") <- "link"
  attr(out, "link") <- "identity"
  attr(out, "singular") <- singular
  if (!is.null(at_boundary)) {
    attr(out, "interaction_at_boundary") <- at_boundary
  }
  class(out) <- c("maihda_interactions", "data.frame")
  out
}

is13_spread <- seq(-1.45, 1.45, length.out = 30)

# Strata on the grid `g`: a stratum's outcomes are its mean plus the same fixed values.
is13_cells <- function(g, mu, within = is13_spread) {
  d <- g[rep(seq_len(nrow(g)), each = length(within)), , drop = FALSE]
  d$y <- rep(mu, each = length(within)) + rep(within, times = nrow(g))
  rownames(d) <- NULL
  d
}
is13_g23 <- function() {
  expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
              stringsAsFactors = FALSE)
}
is13_g45 <- function() {
  expand.grid(d1 = paste0("a", 1:4), d2 = paste0("b", 1:5), stringsAsFactors = FALSE)
}
# Means of the 4 x 5 strata that are exactly additive in the two dimensions.
is13_mu_additive <- function(g = is13_g45()) {
  c(-0.9, -0.2, 0.3, 0.8)[match(g$d1, paste0("a", 1:4))] +
    c(-0.7, -0.3, 0, 0.4, 0.6)[match(g$d2, paste0("b", 1:5))]
}
# A departure from additivity that sums to zero over each dimension.
is13_mu_departure <- function(g = is13_g45()) {
  0.5 * c(1, -1, 1, -1)[match(g$d1, paste0("a", 1:4))] *
    c(1, -1, 0, 1, -1)[match(g$d2, paste0("b", 1:5))]
}
# Six schools of five rows in every stratum, the same five values in each school: a
# school's mean is its stratum's mean, so the school variance is exactly zero.
is13_with_school <- function(g, mu) {
  d <- is13_cells(g, mu, within = rep(c(-1.4, -0.7, 0, 0.7, 1.4), times = 6))
  d$school <- rep(rep(paste0("s", 1:6), each = 5), times = nrow(g))
  d
}

test_that("the note is printed for the term the interaction is read from", {
  # At zero: the bracket names the interaction variance in a crossed-dimensions fit
  # and the between-stratum variance otherwise.
  crossed <- is13_print(is13_table("crossed-dimensions", TRUE, TRUE))
  expect_match(crossed, is13_head_crossed, fixed = TRUE)
  expect_match(crossed, is13_tail, fixed = TRUE)
  expect_false(grepl("between-stratum variance", crossed, fixed = TRUE))
  for (type in c("adjusted (two-model)", "adjusted")) {
    other <- is13_print(is13_table(type, TRUE, TRUE))
    expect_match(other, is13_head_other, fixed = TRUE, info = type)
    expect_match(other, is13_tail, fixed = TRUE, info = type)
  }

  # A zero elsewhere in a singular fit: no note.
  for (type in c("crossed-dimensions", "adjusted (two-model)", "adjusted")) {
    expect_false(grepl("singular/boundary", is13_print(is13_table(type, TRUE, FALSE)),
                       fixed = TRUE), info = type)
  }
  expect_false(grepl("singular/boundary",
                     is13_print(is13_table("adjusted", FALSE, FALSE)), fixed = TRUE))

  # A table made without interaction_at_boundary falls back to `singular`.
  expect_match(is13_print(is13_table("crossed-dimensions", TRUE, NULL)),
               is13_head_crossed, fixed = TRUE)
  expect_match(is13_print(is13_table("adjusted (two-model)", TRUE, NULL)),
               is13_head_other, fixed = TRUE)
  expect_false(grepl("singular/boundary",
                     is13_print(is13_table("adjusted", FALSE, NULL)), fixed = TRUE))
  # A table whose attributes were dropped (head(), `[`): no note.
  bare <- is13_table("adjusted", TRUE, TRUE)
  attr(bare, "singular") <- NULL
  attr(bare, "interaction_at_boundary") <- NULL
  expect_false(grepl("singular/boundary", is13_print(bare), fixed = TRUE))
})

test_that("a crossed-dimensions table notes a zero interaction and nothing else", {
  crossed <- function(d, ...) {
    is13_quiet(maihda(y ~ (1 | gender:ses), data = d,
                      decomposition = "crossed-dimensions", ...))
  }

  # Both dimensions at zero, the interaction not: four of six strata are flagged.
  shifts <- crossed(is13_cells(is13_g23(), c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)))
  mi <- shifts$interactions
  expect_setequal(shifts$model$diagnostics$singular_terms, c("gender", "ses"))
  expect_true(attr(mi, "singular"))
  expect_false(attr(mi, "interaction_at_boundary"))
  expect_equal(attr(mi, "n_flagged"), 4L)
  expect_true(all(mi$se > 0.15))
  expect_false(grepl("singular/boundary", is13_print(mi), fixed = TRUE))

  # The interaction and gender at zero: the standard errors are 0.
  ses_only <- crossed(is13_cells(is13_g23(), c(-0.5, -0.5, 0, 0, 0.5, 0.5)))
  mi <- ses_only$interactions
  expect_setequal(ses_only$model$diagnostics$singular_terms, c("stratum", "gender"))
  expect_true(attr(mi, "interaction_at_boundary"))
  expect_true(all(mi$se == 0))
  expect_gt(ses_only$decomposition$between_var, 0.2)
  out <- is13_print(mi)
  expect_match(out, is13_head_crossed, fixed = TRUE)
  expect_match(out, is13_tail, fixed = TRUE)
  expect_false(grepl("between-stratum variance ~ 0", out, fixed = TRUE))

  # Every crossed term at zero: the note, about the interaction.
  flat <- crossed(is13_cells(is13_g23(), rep(0, 6)))
  expect_setequal(flat$model$diagnostics$singular_terms, c("stratum", "gender", "ses"))
  expect_true(attr(flat$interactions, "interaction_at_boundary"))
  expect_match(is13_print(flat$interactions), is13_head_crossed, fixed = TRUE)

  # Exactly additive 4 x 5 strata: the interaction at zero.
  additive <- is13_quiet(maihda(y ~ (1 | d1:d2), data = is13_cells(is13_g45(),
                                                                    is13_mu_additive()),
                                decomposition = "crossed-dimensions"))
  expect_identical(additive$model$diagnostics$singular_terms, "stratum")
  expect_true(attr(additive$interactions, "interaction_at_boundary"))
  expect_match(is13_print(additive$interactions), is13_head_crossed, fixed = TRUE)

  # A context at zero beside positive crossed terms: no note.
  d <- is13_with_school(is13_g45(), is13_mu_additive() + is13_mu_departure())
  ctx <- is13_quiet(maihda(y ~ (1 | d1:d2), data = d, context = "school",
                           decomposition = "crossed-dimensions"))
  expect_identical(ctx$model$diagnostics$singular_terms, "school")
  expect_true(attr(ctx$interactions, "singular"))
  expect_false(attr(ctx$interactions, "interaction_at_boundary"))
  expect_true(all(ctx$interactions$se > 0.3))
  expect_false(grepl("singular/boundary", is13_print(ctx$interactions), fixed = TRUE))
})

test_that("a two-model table notes the adjusted model's zero stratum variance only", {
  # A context at zero in both models: 16 of 20 strata are flagged, and no note.
  d <- is13_with_school(is13_g45(), is13_mu_additive() + is13_mu_departure())
  a <- is13_quiet(maihda(y ~ (1 | d1:d2), data = d, context = "school"))
  mi <- a$interactions
  expect_identical(attr(mi, "model_type"), "adjusted (two-model)")
  expect_identical(a$model_adjusted$diagnostics$singular_terms, "school")
  expect_true(attr(mi, "singular"))
  expect_false(attr(mi, "interaction_at_boundary"))
  expect_equal(attr(mi, "n_flagged"), 16L)
  expect_equal(unique(round(mi$se, 4)), 0.1743)
  expect_false(grepl("singular/boundary", is13_print(mi), fixed = TRUE))

  # The same adjusted model given to maihda_interactions() directly.
  fit <- is13_quiet(fit_maihda(y ~ d1 + d2 + (1 | d1:d2), data = d, context = "school"))
  mf <- is13_quiet(maihda_interactions(fit))
  expect_identical(attr(mf, "model_type"), "adjusted")
  expect_false(attr(mf, "interaction_at_boundary"))
  expect_equal(attr(mf, "n_flagged"), 16L)
  expect_false(grepl("singular/boundary", is13_print(mf), fixed = TRUE))

  # A record without a list of terms counts as the stratum term at zero, and a
  # record that is not singular (a brms fit's) as not.
  unread <- fit
  unread$diagnostics$singular_terms <- character(0)
  mu <- is13_quiet(maihda_interactions(unread))
  expect_true(attr(mu, "interaction_at_boundary"))
  expect_match(is13_print(mu), is13_head_other, fixed = TRUE)
  never <- fit
  never$diagnostics$singular <- NA
  never$diagnostics$singular_terms <- character(0)
  mn <- is13_quiet(maihda_interactions(never))
  expect_false(attr(mn, "singular"))
  expect_false(attr(mn, "interaction_at_boundary"))

  # Exactly additive strata: the adjusted model's stratum variance is the zero, the
  # standard errors are 0, and the note keeps its words.
  add <- is13_quiet(maihda(y ~ (1 | d1:d2), data = is13_cells(is13_g45(),
                                                               is13_mu_additive())))
  expect_identical(add$model_adjusted$diagnostics$singular_terms, "stratum")
  expect_true(attr(add$interactions, "interaction_at_boundary"))
  expect_true(all(add$interactions$se == 0))
  out <- is13_print(add$interactions)
  expect_match(out, is13_head_other, fixed = TRUE)
  expect_match(out, is13_tail, fixed = TRUE)
})

test_that("wemix and ordinal tables keep the note at a zero stratum variance", {
  skip_on_cran()
  d <- is13_cells(is13_g45(), is13_mu_additive())

  skip_if_not_installed("WeMix")
  dw <- d
  dw$w <- 1
  w <- is13_quiet(maihda(y ~ (1 | d1:d2), data = dw, engine = "wemix",
                         sampling_weights = "w"))
  expect_true(isTRUE(w$model_adjusted$diagnostics$singular))
  expect_length(w$model_adjusted$diagnostics$singular_terms, 0L)
  expect_true(attr(w$interactions, "interaction_at_boundary"))
  expect_match(is13_print(w$interactions), is13_head_other, fixed = TRUE)

  skip_if_not_installed("ordinal")
  do <- d
  do$y <- cut(do$y, breaks = c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
  o <- is13_quiet(maihda(y ~ (1 | d1:d2), data = do, engine = "ordinal"))
  expect_true(isTRUE(o$model_adjusted$diagnostics$singular))
  expect_true(attr(o$interactions, "interaction_at_boundary"))
  expect_match(is13_print(o$interactions), is13_head_other, fixed = TRUE)
})

test_that("the package's crossed-dimensions example notes Finland and not Japan", {
  # Which term of a fit to real data ends at zero is the optimiser's to say.
  skip_on_cran()
  pisa <- function(country) {
    is13_quiet(maihda(math ~ 1 + (1 | gender:ses),
                      data = maihda_country_data[maihda_country_data$country == country, ],
                      decomposition = "crossed-dimensions"))
  }
  jap <- pisa("Japan")
  expect_identical(jap$model$diagnostics$singular_terms, "gender")
  expect_true(all(jap$interactions$se > 5))
  expect_false(attr(jap$interactions, "interaction_at_boundary"))
  expect_false(grepl("singular/boundary", is13_print(jap$interactions), fixed = TRUE))

  fin <- pisa("Finland")
  expect_identical(fin$model$diagnostics$singular_terms, "stratum")
  expect_true(all(fin$interactions$se == 0))
  expect_gt(fin$decomposition$between_var, 1000)
  out <- is13_print(fin$interactions)
  expect_match(out, is13_head_crossed, fixed = TRUE)
  expect_match(out, is13_tail, fixed = TRUE)
})
