# Audit 2026-10-13 (packet 2026-10-13-crossed-singular-note): what the "Singular fit"
# block of print() says for a decomposition = "crossed-dimensions" fit -- CONFIRMED.
# Words of a printed note; no number and no returned object changes.
#
# F1. The block chose its closing sentence from whether the 'stratum' random-effect
#     term was at zero. In the crossed-dimensions model that term is the INTERACTION,
#     and the between-stratum variance is the dimension terms plus the interaction.
#     (a) Exactly additive strata, the interaction at zero: "The between-stratum
#         variance and any VPC/PCV derived from it may be unreliable." above a VPC of
#         0.51 (maihda_country_data: four of six countries, at VPCs of 0.14 to 0.18).
#     (b) A dimension at zero: "The 'stratum' block is NOT at the boundary, so the
#         between-stratum variance and the VPC/PCV are unaffected by this.", although
#         that dimension's variance is part of the additive and of the between-stratum
#         variance (the fit of vignette("cross_classified"), Japan).
#
# Per the user: the block stays where it is, on print() of the analysis, the fit and
# its summary and on the Shiny card. For a crossed-dimensions fit with a zero in some
# of its dimension and interaction terms it says what each zero is -- strata estimated
# as additive; an additive variance without that dimension, shares to be read with
# caution -- and where the between-stratum variance and the VPC come from. A zero in
# all of those terms is a between-stratum variance of zero and keeps the usual
# sentence, as does a zero in a context alone and every fit that is not a
# crossed-dimensions one.
#
# No block draws a random number: the fixtures are built from fixed sequences.

cs13_keep_seed <- function(expr) {
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

cs13_quiet <- function(expr) suppressMessages(suppressWarnings(expr))

# The sentences of the block, as maihda_format_fit_diagnostics() returns them.
cs13_cause <- paste0("  (a variance estimated at ~0, or a perfect correlation between ",
                     "an intercept and a slope).")
cs13_unreliable <- paste0("  The between-stratum variance and any VPC/PCV derived ",
                          "from it may be unreliable.")
cs13_unaffected <- c(
  "  The 'stratum' block is NOT at the boundary, so the between-stratum variance and",
  "  the VPC/PCV are unaffected by this.")
cs13_interaction <- c(
  "  'stratum' is the interaction term of this crossed-dimensions model: the strata are",
  "  estimated as additive (interaction share 0%).")
cs13_closing <- paste0("  The between-stratum variance and the VPC come from the terms ",
                       "not at the boundary.")
cs13_one_block <- function(term) {
  sprintf("Singular fit: the random-effects block for '%s' is at (or near) the boundary",
          term)
}
cs13_blocks <- function(listed) {
  sprintf("Singular fit: the random-effects blocks for %s are at (or near) the boundary",
          listed)
}

# A singular lme4 record with these terms at zero.
cs13_record <- function(terms) {
  list(engine = "lme4", singular = TRUE, singular_terms = terms,
       stratum_singular = if (length(terms) == 0L) NA else "stratum" %in% terms,
       converged = TRUE, messages = character(0), adequacy = NULL)
}
cs13_cc <- list(dim_groups = c(gender = "gender", ses = "ses"),
                interaction_group = "stratum", dim_labels = c("gender", "ses"))

# Lines as one text, for sentences whose line breaks depend on the names in them.
cs13_text <- function(lines) paste(trimws(lines), collapse = " ")

cs13_spread <- seq(-1.45, 1.45, length.out = 30)

# Strata on the grid `g`: a stratum's outcomes are its mean plus the same fixed values.
cs13_cells <- function(g, mu, within = cs13_spread) {
  d <- g[rep(seq_len(nrow(g)), each = length(within)), , drop = FALSE]
  d$y <- rep(mu, each = length(within)) + rep(within, times = nrow(g))
  rownames(d) <- NULL
  d
}
cs13_g23 <- function() {
  expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
              stringsAsFactors = FALSE)
}
cs13_g45 <- function() {
  expand.grid(d1 = paste0("a", 1:4), d2 = paste0("b", 1:5), stringsAsFactors = FALSE)
}
# Means of the 4 x 5 strata that are exactly additive in the two dimensions.
cs13_mu_additive <- function(g = cs13_g45()) {
  c(-0.9, -0.2, 0.3, 0.8)[match(g$d1, paste0("a", 1:4))] +
    c(-0.7, -0.3, 0, 0.4, 0.6)[match(g$d2, paste0("b", 1:5))]
}
# A departure from additivity that sums to zero over each dimension.
cs13_mu_departure <- function(g = cs13_g45()) {
  0.5 * c(1, -1, 1, -1)[match(g$d1, paste0("a", 1:4))] *
    c(1, -1, 0, 1, -1)[match(g$d2, paste0("b", 1:5))]
}
cs13_shifts <- c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)

# The same strata with six schools of five rows in every stratum, the same five values
# in each school: a school's mean is its stratum's mean, so the school variance is zero.
cs13_with_school <- function(g, mu) {
  d <- cs13_cells(g, mu, within = rep(c(-1.4, -0.7, 0, 0.7, 1.4), times = 6))
  d$school <- rep(rep(paste0("s", 1:6), each = 5), times = nrow(g))
  d
}

cs13_crossed <- function(formula, data, ...) {
  cs13_quiet(maihda(formula, data = data, decomposition = "crossed-dimensions", ...))
}

# The "Fit diagnostics" block of a print, up to the blank line that ends it.
cs13_printed <- function(x) {
  out <- utils::capture.output(print(x))
  i <- grep("^Fit diagnostics", out)
  if (length(i) == 0L) {
    return(character(0))
  }
  j <- i[1L]
  while (j < length(out) && nzchar(trimws(out[j + 1L]))) j <- j + 1L
  out[i[1L]:j]
}
# The block a print makes of these lines.
cs13_as_printed <- function(lines) c("Fit diagnostics:", paste0("  ", lines))

# One variance per random-effect term, read off the fitted lme4 model.
cs13_variances <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  stats::setNames(vc$vcov, ifelse(is.na(vc$grp), "Residual", vc$grp))
}

# The Shiny app's objects, sourced as test-run_app.R does it.
cs13_source_app <- function() {
  app_file <- system.file("shiny", "app.R", package = "MAIHDA")
  if (app_file == "") {
    app_file <- file.path("inst", "shiny", "app.R")
  }
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  app_env <- new.env(parent = globalenv())
  suppressPackageStartupMessages(suppressWarnings(sys.source(app_file, envir = app_env)))
  app_env
}

test_that("a crossed-dimensions block says what each zero is", {
  fmt <- function(terms, cc = cs13_cc) {
    maihda_format_fit_diagnostics(cs13_record(terms), cc_info = cc)
  }

  # The interaction alone: strata estimated as additive.
  expect_identical(fmt("stratum"),
                   c(cs13_one_block("stratum"), cs13_cause, cs13_interaction,
                     cs13_closing))

  # One dimension: the additive variance has nothing from it.
  expect_identical(
    fmt("gender"),
    c(cs13_one_block("gender"), cs13_cause,
      "  'gender' is a stratum dimension of this crossed-dimensions model: the additive variance",
      "  has nothing from it, so read the additive and interaction shares with caution.",
      cs13_closing))

  # Every dimension, the interaction not: named in the order of the first line.
  expect_identical(
    fmt(c("ses", "gender")),
    c(cs13_blocks("'ses' and 'gender'"), cs13_cause,
      "  'ses' and 'gender' are stratum dimensions of this crossed-dimensions model: the",
      "  additive variance has nothing from them, so read the additive and interaction shares",
      "  with caution.",
      cs13_closing))

  # The interaction and a dimension: one sentence each, the model named once.
  both <- fmt(c("stratum", "gender"))
  expect_identical(
    both,
    c(cs13_blocks("'stratum' and 'gender'"), cs13_cause, cs13_interaction,
      "  'gender' is a stratum dimension: the additive variance has nothing from it, so read the",
      "  additive and interaction shares with caution.",
      cs13_closing))
  expect_length(gregexpr("crossed-dimensions model", cs13_text(both), fixed = TRUE)[[1L]],
                1L)
  # The interaction sentence comes first whatever the order of the terms.
  expect_identical(fmt(c("gender", "stratum"))[-1L], both[-1L])

  # None of these says that the between-stratum variance is zero or untouched.
  for (terms in list("stratum", "gender", c("ses", "gender"), c("stratum", "gender"))) {
    txt <- cs13_text(fmt(terms))
    expect_false(grepl("may be unreliable", txt, fixed = TRUE), info = toString(terms))
    expect_false(grepl("unaffected", txt, fixed = TRUE), info = toString(terms))
    expect_false(grepl("NOT at the boundary", txt, fixed = TRUE), info = toString(terms))
  }
})

test_that("the usual sentences serve a zero in every crossed term or in none", {
  fmt <- function(terms, cc = cs13_cc) {
    maihda_format_fit_diagnostics(cs13_record(terms), cc_info = cc)
  }
  # Every dimension and the interaction: the between-stratum variance is the zero.
  expect_identical(fmt(c("stratum", "ses", "gender")),
                   c(cs13_blocks("'stratum', 'ses' and 'gender'"), cs13_cause,
                     cs13_unreliable))
  expect_identical(fmt(c("school", "stratum", "ses", "gender")),
                   c(cs13_blocks("'school', 'stratum', 'ses' and 'gender'"), cs13_cause,
                     cs13_unreliable))
  # A context alone: no crossed term is at zero.
  expect_identical(fmt("school"),
                   c(cs13_one_block("school"), cs13_cause, cs13_unaffected))

  # A context next to crossed terms: the crossed terms get their sentences, and the
  # context is named on the first line only.
  with_ctx <- fmt(c("stratum", "school"))
  expect_identical(with_ctx,
                   c(cs13_blocks("'stratum' and 'school'"), cs13_cause, cs13_interaction,
                     cs13_closing))
  expect_length(grep("school", with_ctx, fixed = TRUE), 1L)
  expect_match(cs13_text(fmt(c("school", "ses"))),
               "'ses' is a stratum dimension of this crossed-dimensions model:",
               fixed = TRUE)

  # Every set of terms at zero: the block speaks of a between-stratum variance at zero
  # exactly when the group comparison would name the fit for one
  # (maihda_between_variance_singular()), and of an untouched one exactly when no
  # crossed term is at zero.
  all_terms <- c("gender", "ses", "stratum", "school")
  crossed <- c("gender", "ses", "stratum")
  for (k in seq_len(15L)) {
    terms <- all_terms[bitwAnd(k, c(1L, 2L, 4L, 8L)) > 0L]
    lines <- fmt(terms)
    between_zero <- maihda_between_variance_singular(cs13_record(terms), crossed)
    expect_identical(between_zero, all(crossed %in% terms), info = toString(terms))
    expect_identical(cs13_unreliable %in% lines, between_zero, info = toString(terms))
    expect_identical(all(cs13_unaffected %in% lines), !any(crossed %in% terms),
                     info = toString(terms))
    expect_identical(cs13_closing %in% lines,
                     any(crossed %in% terms) && !between_zero, info = toString(terms))
    expect_identical(all(cs13_interaction %in% lines),
                     "stratum" %in% terms && !between_zero, info = toString(terms))
    expect_identical(grepl("the additive variance has nothing from", cs13_text(lines),
                           fixed = TRUE),
                     any(c("gender", "ses") %in% terms) && !between_zero,
                     info = toString(terms))
  }

  # The terms could not be read: the conservative sentence, with or without the tag.
  unread <- c("Singular fit: at least one variance component is estimated at (or near) zero.",
              cs13_unreliable)
  expect_identical(fmt(character(0)), unread)
  expect_identical(maihda_format_fit_diagnostics(cs13_record(character(0))), unread)

  # Not singular, or not known to be (a brms record): no block.
  quiet <- cs13_record("stratum")
  quiet$singular <- FALSE
  expect_length(maihda_format_fit_diagnostics(quiet, cc_info = cs13_cc), 0L)
  quiet$singular <- NA
  quiet$engine <- "brms"
  expect_length(maihda_format_fit_diagnostics(quiet, cc_info = cs13_cc), 0L)
  expect_null(maihda_crossed_singular_lines(cs13_record("stratum"), NULL))
})

test_that("a fit that is not a crossed-dimensions one keeps its block", {
  plain <- function(terms, ...) maihda_format_fit_diagnostics(cs13_record(terms), ...)

  # The stratum term IS the between-stratum variance there.
  expect_identical(plain("stratum"),
                   c(cs13_one_block("stratum"), cs13_cause, cs13_unreliable))
  expect_identical(plain(c("id", "stratum")),
                   c(cs13_blocks("'id' and 'stratum'"), cs13_cause, cs13_unreliable))
  expect_identical(plain("id"), c(cs13_one_block("id"), cs13_cause, cs13_unaffected))
  # The names of a crossed-dimensions model's terms mean nothing without its tag.
  expect_identical(plain("gender"),
                   c(cs13_one_block("gender"), cs13_cause, cs13_unaffected))
  expect_identical(plain("stratum", cc_info = NULL), plain("stratum"))
  # The adjusted model of a two-model print: a zero stratum variance stays unprinted.
  expect_length(plain("stratum", suppress_stratum_singular = TRUE), 0L)

  # Convergence messages follow the block as before, with the tag as well.
  nc <- cs13_record("stratum")
  nc$converged <- FALSE
  nc$messages <- "Model failed to converge with max|grad| = 0.0071"
  tail_lines <- c("Convergence warnings reported by lme4:",
                  "  - Model failed to converge with max|grad| = 0.0071")
  expect_identical(maihda_format_fit_diagnostics(nc),
                   c(cs13_one_block("stratum"), cs13_cause, cs13_unreliable, tail_lines))
  expect_identical(maihda_format_fit_diagnostics(nc, cc_info = cs13_cc),
                   c(cs13_one_block("stratum"), cs13_cause, cs13_interaction, cs13_closing,
                     tail_lines))
})

test_that("the sentences name terms as the first line does and stay within the block's width", {
  # An auto-binned dimension's term is its grouping factor.
  cc3 <- list(dim_groups = c(Gender = "Gender", Age = ".maihda_dim_Age", Race = "Race"),
              interaction_group = "stratum")
  one <- maihda_format_fit_diagnostics(cs13_record(".maihda_dim_Age"), cc_info = cc3)
  expect_identical(one[1L], cs13_one_block(".maihda_dim_Age"))
  expect_identical(
    cs13_text(one[-(1:2)]),
    paste("'.maihda_dim_Age' is a stratum dimension of this crossed-dimensions model:",
          "the additive variance has nothing from it, so read the additive and",
          "interaction shares with caution. The between-stratum variance and the VPC",
          "come from the terms not at the boundary."))

  three <- maihda_format_fit_diagnostics(
    cs13_record(c("Race", ".maihda_dim_Age", "Gender")), cc_info = cc3)
  expect_identical(
    cs13_text(three[-(1:2)]),
    paste("'Race', '.maihda_dim_Age' and 'Gender' are stratum dimensions of this",
          "crossed-dimensions model: the additive variance has nothing from them, so",
          "read the additive and interaction shares with caution. The between-stratum",
          "variance and the VPC come from the terms not at the boundary."))

  # No sentence line is longer than the line every block already has, and each
  # sentence starts a line.
  for (lines in list(one, three,
                     maihda_format_fit_diagnostics(cs13_record(c("stratum", "Race")),
                                                   cc_info = cc3))) {
    expect_lte(max(nchar(lines[-(1:2)])), nchar(cs13_cause))
    expect_true(cs13_closing %in% lines)
    expect_true(all(startsWith(lines[-1L], "  ")))
  }
  starts <- maihda_crossed_singular_lines(cs13_record(c("stratum", "Race")), cc3)
  expect_true(startsWith(starts[1L], "  'stratum' is the interaction term"))
  expect_length(grep("^  'Race' is a stratum dimension:", starts), 1L)
})

test_that("print() of a crossed-dimensions fit with additive strata says so", {
  skip_on_cran()
  d <- cs13_cells(cs13_g45(), cs13_mu_additive())
  a <- cs13_crossed(y ~ (1 | d1:d2), d)

  # The premise: the interaction term alone is at zero, and the record is as it was.
  rec <- a$model$diagnostics
  expect_identical(rec$singular_terms, "stratum")
  expect_true(rec$stratum_singular)
  expect_length(setdiff(names(rec), c("engine", "singular", "singular_terms",
                                      "stratum_singular", "converged", "messages",
                                      "adequacy")), 0L)

  # The numbers: the between-stratum variance is the dimensions plus the interaction.
  v <- cs13_variances(a$model$model)
  expect_equal(v[["stratum"]], 0)
  expect_equal(a$decomposition$additive_var, v[["d1"]] + v[["d2"]], tolerance = 1e-12)
  expect_equal(a$decomposition$additive_var, 0.7903, tolerance = 1e-4)
  expect_equal(a$decomposition$between_var, a$decomposition$additive_var)
  expect_equal(a$summary$vpc$estimate, (v[["d1"]] + v[["d2"]]) / sum(v), tolerance = 1e-12)
  expect_equal(a$summary$vpc$estimate, 0.5100, tolerance = 1e-4)

  want <- cs13_as_printed(c(cs13_one_block("stratum"), cs13_cause, cs13_interaction,
                            cs13_closing))
  expect_identical(cs13_printed(a), want)
  expect_identical(cs13_printed(a$model), want)
  expect_identical(cs13_printed(a$summary), want)
  expect_identical(cs13_printed(cs13_quiet(summary(a$model))), want)
  out <- paste(utils::capture.output(print(a)), collapse = "\n")
  expect_false(grepl("may be unreliable", out, fixed = TRUE))
  expect_match(out, "VPC/ICC: 0.5100", fixed = TRUE)
  expect_match(out, "Interaction share of between-strata variance: 0.0%", fixed = TRUE)

  # The same model without the crossed-dimensions tag: the dimension terms are other
  # random effects there, the VPC is the 'stratum' term's and the block is the usual one.
  d$stratum <- paste(d$d1, d$d2, sep = "_")
  untagged <- cs13_quiet(fit_maihda(y ~ (1 | d1) + (1 | d2) + (1 | stratum), data = d))
  expect_null(untagged$cc_info)
  expect_equal(cs13_quiet(summary(untagged))$vpc$estimate, 0)
  expect_identical(cs13_printed(untagged),
                   cs13_as_printed(c(cs13_one_block("stratum"), cs13_cause,
                                     cs13_unreliable)))

  # The two-model print of the same rows: no block for the adjusted model, and the note
  # at the PCV.
  two <- paste(utils::capture.output(print(cs13_quiet(
    maihda(y ~ (1 | d1:d2), data = d[c("d1", "d2", "y")], interactions = FALSE)))),
    collapse = "\n")
  expect_false(grepl("Fit diagnostics", two, fixed = TRUE))
  expect_match(two, "pinned near 100", fixed = TRUE)
})

test_that("print() of a crossed-dimensions fit with a dimension at zero points to the shares", {
  skip_on_cran()
  # Strata that are not additive at all: both dimensions at zero.
  d <- cs13_cells(cs13_g23(), cs13_shifts)
  a <- cs13_crossed(y ~ (1 | gender:ses), d)
  expect_setequal(a$model$diagnostics$singular_terms, c("gender", "ses"))
  expect_false(a$model$diagnostics$stratum_singular)
  v <- cs13_variances(a$model$model)
  expect_equal(a$decomposition$additive_var, 0)
  expect_equal(a$decomposition$interaction_var, 0.3622, tolerance = 1e-4)
  expect_equal(a$summary$vpc$estimate, v[["stratum"]] / sum(v), tolerance = 1e-12)
  expect_equal(a$summary$vpc$estimate, 0.3185, tolerance = 1e-4)
  listed <- paste0("'", a$model$diagnostics$singular_terms, "'", collapse = " and ")
  want <- cs13_as_printed(c(
    cs13_blocks(listed), cs13_cause,
    paste0("  ", strwrap(paste0(
      listed, " are stratum dimensions of this crossed-dimensions model: the additive ",
      "variance has nothing from them, so read the additive and interaction shares ",
      "with caution."), width = 88)),
    cs13_closing))
  expect_identical(cs13_printed(a), want)
  expect_identical(cs13_printed(a$model), want)
  expect_identical(cs13_printed(a$summary), want)
  expect_false(any(grepl("unaffected", cs13_printed(a), fixed = TRUE)))

  # One dimension with an effect and nothing else: the interaction and the other
  # dimension are at zero.
  d1 <- cs13_cells(cs13_g23(), c(-0.5, -0.5, 0, 0, 0.5, 0.5))
  b <- cs13_crossed(y ~ (1 | gender:ses), d1)
  expect_setequal(b$model$diagnostics$singular_terms, c("stratum", "gender"))
  expect_equal(b$summary$vpc$estimate, 0.2375, tolerance = 1e-4)
  block <- cs13_printed(b)
  expect_identical(block[4:5], paste0("  ", cs13_interaction))
  expect_match(cs13_text(block), paste("'gender' is a stratum dimension: the additive",
                                       "variance has nothing from it, so read the additive",
                                       "and interaction shares with caution."),
               fixed = TRUE)
  expect_identical(block[length(block)], paste0("  ", cs13_closing))
  expect_identical(cs13_printed(b$summary), block)

  # Identical strata: every term at zero, a VPC of 0 and the usual sentence.
  flat <- cs13_crossed(y ~ (1 | gender:ses), cs13_cells(cs13_g23(), rep(0, 6)))
  expect_setequal(flat$model$diagnostics$singular_terms, c("stratum", "gender", "ses"))
  expect_equal(flat$summary$vpc$estimate, 0)
  block <- cs13_printed(flat)
  expect_identical(block[3:4], paste0("  ", c(cs13_cause, cs13_unreliable)))
  expect_length(block, 4L)
})

test_that("a context at zero beside the crossed terms is told from them", {
  skip_on_cran()
  # The context alone: every crossed term is positive.
  d <- cs13_with_school(cs13_g45(), cs13_mu_additive() + cs13_mu_departure())
  a <- cs13_crossed(y ~ (1 | d1:d2), d, context = "school")
  expect_identical(a$model$diagnostics$singular_terms, "school")
  expect_gt(a$decomposition$interaction_var, 0.2)
  want <- cs13_as_printed(c(cs13_one_block("school"), cs13_cause, cs13_unaffected))
  expect_identical(cs13_printed(a), want)
  expect_identical(cs13_printed(a$summary), want)

  # The context and the interaction.
  d <- cs13_with_school(cs13_g45(), cs13_mu_additive())
  b <- cs13_crossed(y ~ (1 | d1:d2), d, context = "school")
  expect_setequal(b$model$diagnostics$singular_terms, c("stratum", "school"))
  expect_gt(b$summary$vpc$estimate, 0.4)
  block <- cs13_printed(b)
  expect_identical(block[-(1:2)],
                   paste0("  ", c(cs13_cause, cs13_interaction, cs13_closing)))
  expect_identical(cs13_printed(b$model), block)

  # The same context on a fit that is not a crossed-dimensions one: unchanged.
  plain <- cs13_quiet(fit_maihda(y ~ (1 | d1:d2), data = d, context = "school"))
  expect_identical(plain$diagnostics$singular_terms, "school")
  expect_identical(cs13_printed(plain), want)
})

test_that("the package's crossed-dimensions example prints what each country's zero is", {
  # Which term of a fit to real data ends at zero is the optimiser's to say.
  skip_on_cran()
  pisa <- function(country) {
    cs13_crossed(math ~ 1 + (1 | gender:ses),
                 maihda_country_data[maihda_country_data$country == country, ])
  }
  fin <- pisa("Finland")
  expect_identical(fin$model$diagnostics$singular_terms, "stratum")
  expect_equal(fin$summary$vpc$estimate, 0.1403, tolerance = 1e-3)
  expect_identical(cs13_printed(fin),
                   cs13_as_printed(c(cs13_one_block("stratum"), cs13_cause,
                                     cs13_interaction, cs13_closing)))

  jap <- pisa("Japan")
  expect_identical(jap$model$diagnostics$singular_terms, "gender")
  expect_equal(jap$summary$vpc$estimate, 0.1622, tolerance = 1e-3)
  block <- cs13_printed(jap)
  expect_match(cs13_text(block), "'gender' is a stratum dimension of this crossed-dimensions model:",
               fixed = TRUE)
  expect_false(any(grepl("unaffected", block, fixed = TRUE)))
  expect_identical(block[length(block)], paste0("  ", cs13_closing))
})

test_that("the Shiny card shows the same sentences", {
  skip_on_cran()
  for (pkg in maihda_app_required_packages()) skip_if_not_installed(pkg)
  d <- cs13_cells(cs13_g45(), cs13_mu_additive())
  res <- cs13_quiet(maihda_app_fit_models(d, "y", c("d1", "d2"),
                                          decomposition = "crossed-dimensions"))
  expect_identical(res$model$diagnostics$singular_terms, "stratum")
  expect_false(is.null(res$model$cc_info))
  expect_identical(
    maihda_format_fit_diagnostics(res$model$diagnostics, cc_info = res$model$cc_info),
    c(cs13_one_block("stratum"), cs13_cause, cs13_interaction, cs13_closing))

  cs13_keep_seed({
    app_env <- cs13_source_app()
    shiny::testServer(app_env$server, {
      session$setInputs(dataset = "pisa")
      model_results(res$model)
      summary_results(res$summary_obj)
      null_summary_results(res$summary_obj)
      decomposition_results(res$decomposition)
      decomposition_mode("crossed-dimensions")
      fitted_family("gaussian")
      pcv_results(NULL)
      stepwise_results(NULL)
      session$flushReact()

      card <- paste(unlist(output$model_summary_ui), collapse = " ")
      expect_match(card, "is the interaction term of this crossed-dimensions model: the strata are",
                   fixed = TRUE)
      expect_match(card, "come from the terms not at the boundary.", fixed = TRUE)
      expect_false(grepl("may be unreliable", card, fixed = TRUE))
    })
  })
})
