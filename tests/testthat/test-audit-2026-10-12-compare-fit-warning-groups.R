# Audit 2026-10-12 (packet 2026-10-12-compare-fit-warning-groups): WHICH groups the
# fit-quality warning of compare_maihda_groups() names -- both findings CONFIRMED. No
# number, row, column or status cell changes.
#
# F1. A group was named under "singular fit", and its VPC put under "with caution",
#     whenever lme4::isSingular() was TRUE for its fit, whichever variance component
#     was the zero. With exactly additive strata a crossed-dimensions comparison named
#     both groups at VPCs of 0.51 and 0.27 (the zero is the interaction, which the
#     two-model path records as pcv_status = "singular" without a warning); a context
#     without variance named both at 0.32 and 0.37; and the package's own example,
#     maihda_country_data with decomposition = "crossed-dimensions", named all six
#     countries at VPCs of 0.08 to 0.18.
# F2. Only each group's own model reached the warning. A group whose ADJUSTED model
#     did not converge kept pcv_status = "ok", was named nowhere, and its PCV was
#     reported like any other: lme4's own warnings name no group, and WeMix raises
#     none. On maihda_health_data under default control (Obese by Race, three
#     dimensions) two of five groups were in that state.
#
# Per the user: a group is named as a singular fit only when the zero is its
# between-stratum variance (in the crossed-dimensions model every dimension and the
# interaction), so that var_between and the VPC are 0; a group whose adjusted model did
# not converge is named in the same warning, which says that it is the adjusted model
# and points to the PCV; the same on every engine, a brms adjusted model being named
# for a largest R-hat above 1.01 or divergent transitions. pcv_status keeps its three
# values.
#
# No block draws a random number: the fixtures are built from fixed sequences, and
# the opt-in Stan block puts the caller's seed back.

fg12_observe <- function(expr) {
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
fg12_quality <- function(warns) warns[grepl("fit problems", warns, fixed = TRUE)]

fg12_keep_seed <- function(expr) {
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

fg12_open <- "compare_maihda_groups(): fit problems ("
fg12_singular <- paste0("A singular fit listed here has its between-stratum variance ",
                        "estimated at zero, so var_between and the VPC are 0.")
fg12_nonconv <- "A fit that did not converge may be unreliable."
fg12_brms <- "The posterior summaries of these brms fits may be unreliable."
fg12_flag <- "largest R-hat above 1.01 or divergent transitions: "
fg12_vpc <- " Interpret these groups' VPC/ICC values with caution."
fg12_pcv <- " Interpret these groups' PCV values with caution."
fg12_both <- paste0(" Interpret these groups' VPC/ICC values with caution; for an ",
                    "adjusted model, the PCV instead.")

fg12_spread <- seq(-1.45, 1.45, length.out = 30)
fg12_flat <- rep(0, 6)
fg12_s2 <- c(-0.9, 0.2, 0.5, -0.4, 0.8, -0.2)
fg12_s3 <- c(0.7, -0.6, 0.1, 0.9, -0.8, -0.3)

# One group of gender x ses strata, 30 rows each: a stratum's outcomes are its shift
# plus the same fixed spread, so shifts of zero give identical strata. Twelve schools
# cut across every stratum alike and carry no effect; `w` holds unit weights.
fg12_group <- function(country, shift) {
  g <- expand.grid(gender = c("F", "M"), ses = c("low", "mid", "high"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 30), ]
  d$y <- rep(shift, each = 30) + rep(fg12_spread, times = nrow(g))
  d$country <- country
  d$school <- paste0("s", rep(rep(1:12, length.out = 30), times = nrow(g)))
  d$w <- 1
  rownames(d) <- NULL
  d
}

# One group of 4 x 5 strata whose means are exactly additive in the two dimensions.
fg12_additive_group <- function(country, scale) {
  a <- c(-0.9, -0.2, 0.3, 0.8) * scale
  b <- c(-0.7, -0.3, 0, 0.4, 0.6) * scale
  g <- expand.grid(d1 = paste0("a", 1:4), d2 = paste0("b", 1:5),
                   stringsAsFactors = FALSE)
  mu <- a[match(g$d1, paste0("a", 1:4))] + b[match(g$d2, paste0("b", 1:5))]
  d <- g[rep(seq_len(nrow(g)), each = 30), ]
  d$y <- rep(mu, each = 30) + rep(fg12_spread, times = nrow(g))
  d$country <- country
  rownames(d) <- NULL
  d
}

# One group of 3 x 4 strata with a binary outcome: `k` of a stratum's 30 rows are 1.
fg12_binary_group <- function(country, k) {
  g <- expand.grid(gender = paste0("g", 1:3), race = paste0("r", 1:4),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 30), ]
  d$yb <- unlist(lapply(k, function(s) rep(c(1, 0), c(s, 30 - s))))
  d$country <- country
  rownames(d) <- NULL
  d
}

# TRUE for the adjusted model of a two-model decomposition: the fixtures have no
# covariate, so only that model has fixed-effect terms.
fg12_is_adjusted <- function(formula) {
  length(attr(stats::terms(maihda_nobars(formula)), "term.labels")) > 0L
}

# A stand-in for fit_maihda() that fits as usual and keeps every fit's diagnostics
# record, in the order the comparison made them.
fg12_recorder <- function(real_fit) {
  log <- new.env()
  log$fits <- list()
  list(
    log = log,
    fit = function(formula, data, ...) {
      m <- real_fit(formula, data, ...)
      log$fits[[length(log$fits) + 1L]] <- list(
        group = as.character(data$country[1]),
        adjusted = fg12_is_adjusted(formula),
        diagnostics = m$diagnostics)
      m
    })
}

# The groups among `fits` (own models, or adjusted ones) whose record satisfies `test`.
fg12_groups <- function(fits, adjusted, test) {
  keep <- vapply(fits, function(f) identical(f$adjusted, adjusted) &&
                   isTRUE(test(f$diagnostics)), logical(1))
  vapply(fits[keep], function(f) f$group, character(1))
}

test_that("a singular fit counts only when the between-stratum variance is the zero", {
  rec <- function(singular, terms = character(0)) {
    list(singular = singular, singular_terms = terms)
  }
  # Not singular, or not known to be: nothing to name.
  expect_false(maihda_between_variance_singular(rec(FALSE)))
  expect_false(maihda_between_variance_singular(rec(NA)))
  expect_false(maihda_between_variance_singular(NULL))
  expect_false(maihda_between_variance_singular(rec(FALSE, "stratum")))
  expect_false(maihda_between_variance_singular(rec(NA, "stratum")))

  # One stratum term: its variance is the between-stratum variance.
  expect_true(maihda_between_variance_singular(rec(TRUE, "stratum")))
  expect_true(maihda_between_variance_singular(rec(TRUE, c("school", "stratum"))))
  expect_false(maihda_between_variance_singular(rec(TRUE, "school")))
  expect_false(maihda_between_variance_singular(rec(TRUE, c("school", "clinic"))))

  # Crossed-dimensions: every dimension and the interaction together.
  cc <- c("gender", "ses", "stratum")
  expect_true(maihda_between_variance_singular(rec(TRUE, c("stratum", "ses", "gender")),
                                               cc))
  expect_true(maihda_between_variance_singular(
    rec(TRUE, c("school", "stratum", "ses", "gender")), cc))
  # Additive strata: the interaction alone.
  expect_false(maihda_between_variance_singular(rec(TRUE, "stratum"), cc))
  # No additive part: the dimensions alone, or one of them.
  expect_false(maihda_between_variance_singular(rec(TRUE, c("ses", "gender")), cc))
  expect_false(maihda_between_variance_singular(rec(TRUE, "gender"), cc))
  expect_false(maihda_between_variance_singular(rec(TRUE, c("stratum", "gender")), cc))
  expect_false(maihda_between_variance_singular(rec(TRUE, c("school", "ses", "gender")),
                                                cc))
  expect_false(maihda_between_variance_singular(rec(TRUE, "school"), cc))

  # A record without a list of terms (wemix and ordinal, which flag the stratum
  # variance itself; an lme4 fit whose terms could not be read) is taken as one.
  expect_true(maihda_between_variance_singular(rec(TRUE)))
  expect_true(maihda_between_variance_singular(rec(TRUE), cc))
})

test_that("the warning text names the adjusted model and points to the right number", {
  for (engine in c("lme4", "wemix", "ordinal")) {
    # Only adjusted models: the PCV, and no word about the VPC.
    adjusted_only <- paste0(fg12_open, "adjusted model did not converge: A, B). ",
                            fg12_nonconv, fg12_pcv)
    expect_identical(
      maihda_group_fit_quality_message(engine, character(0), character(0),
                                       c("A", "B")),
      adjusted_only)
    expect_false(grepl("VPC", adjusted_only, fixed = TRUE))
    expect_false(grepl("singular", adjusted_only, fixed = TRUE))

    # A group's own model and an adjusted one: each entry says which model it is,
    # the explanation is given once, and the closing sentence tells the two apart.
    expect_identical(
      maihda_group_fit_quality_message(engine, character(0), "A", c("A", "B")),
      paste0(fg12_open, "did not converge: A; adjusted model did not converge: A, B). ",
             fg12_nonconv, fg12_both))
    expect_identical(
      maihda_group_fit_quality_message(engine, "g1", character(0), "g2"),
      paste0(fg12_open, "singular fit: g1; adjusted model did not converge: g2). ",
             fg12_singular, " ", fg12_nonconv, fg12_both))
    expect_identical(
      maihda_group_fit_quality_message(engine, c("g1", "g4"), "g2", "g3"),
      paste0(fg12_open, "singular fit: g1, g4; did not converge: g2; ",
             "adjusted model did not converge: g3). ",
             fg12_singular, " ", fg12_nonconv, fg12_both))

    # Without an adjusted model the closing sentence is about the VPC alone, and a
    # call that passes no fourth argument is the same call.
    own_only <- paste0(fg12_open, "singular fit: g1; did not converge: g2). ",
                       fg12_singular, " ", fg12_nonconv, fg12_vpc)
    expect_identical(
      maihda_group_fit_quality_message(engine, "g1", "g2", character(0)), own_only)
    expect_identical(maihda_group_fit_quality_message(engine, "g1", "g2"), own_only)
    expect_false(grepl("PCV", own_only, fixed = TRUE))
    expect_false(grepl("adjusted", own_only, fixed = TRUE))
  }

  # brms: both entries say what the package checks there.
  expect_identical(
    maihda_group_fit_quality_message("brms", character(0), character(0), "g1"),
    paste0(fg12_open, "adjusted model, ", fg12_flag, "g1). ", fg12_brms, fg12_pcv))
  flagged <- maihda_group_fit_quality_message("brms", character(0), "g1",
                                              c("g1", "g2", "g3"))
  expect_identical(
    flagged,
    paste0(fg12_open, fg12_flag, "g1; adjusted model, ", fg12_flag, "g1, g2, g3). ",
           fg12_brms, fg12_both))
  expect_false(grepl("did not converge", flagged, fixed = TRUE))

  # Nothing listed, nothing said.
  for (engine in c("lme4", "wemix", "ordinal", "brms")) {
    expect_null(maihda_group_fit_quality_message(engine, character(0), character(0),
                                                 character(0)))
  }

  # A singular fit is explained as what it is named for, on every engine.
  for (engine in c("lme4", "wemix", "ordinal", "brms")) {
    one <- maihda_group_fit_quality_message(engine, "g1", character(0))
    expect_identical(one, paste0(fg12_open, "singular fit: g1). ", fg12_singular,
                                 fg12_vpc))
    expect_false(grepl("where that is", one, fixed = TRUE))
    expect_false(grepl("a variance component", one, fixed = TRUE))
  }
})

test_that("a zero interaction, dimension or context variance names no group", {
  real_fit <- fit_maihda

  # Exactly additive strata, crossed-dimensions: the interaction is the zero.
  additive <- rbind(fg12_additive_group("g1", 1), fg12_additive_group("g2", 0.6))
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | d1:d2), data = additive,
                                            group = "country",
                                            decomposition = "crossed-dimensions"))
  # The premise: both fits ARE singular, in the interaction term alone.
  expect_length(rec$log$fits, 2L)
  for (f in rec$log$fits) {
    expect_true(isTRUE(f$diagnostics$singular), info = f$group)
    expect_identical(f$diagnostics$singular_terms, "stratum", info = f$group)
  }
  expect_true(all(obs$res$var_interaction < 1e-8))
  expect_equal(obs$res$vpc, c(0.5100, 0.2674), tolerance = 1e-3)
  expect_equal(obs$res$additive_share, c(1, 1))
  expect_length(fg12_quality(obs$warns), 0L)
  expect_identical(obs$res$status, c("ok", "ok"))

  # The two-model path on the same rows: silent as well, with the state in pcv_status.
  two <- fg12_observe(compare_maihda_groups(y ~ (1 | d1:d2), data = additive,
                                            group = "country"))
  expect_length(fg12_quality(two$warns), 0L)
  expect_identical(two$res$pcv_status, c("singular", "singular"))

  # Strata that are not additive at all, crossed-dimensions: the dimensions are the zero.
  shifted <- rbind(fg12_group("g2", fg12_s2), fg12_group("g3", fg12_s3))
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = shifted,
                                            group = "country",
                                            decomposition = "crossed-dimensions"))
  expect_length(rec$log$fits, 2L)
  expect_true(all(vapply(rec$log$fits, function(f) isTRUE(f$diagnostics$singular),
                         logical(1))))
  expect_setequal(rec$log$fits[[1L]]$diagnostics$singular_terms, c("gender", "ses"))
  expect_identical(rec$log$fits[[2L]]$diagnostics$singular_terms, "gender")
  expect_gt(min(obs$res$var_interaction), 0.3)
  expect_gt(min(obs$res$vpc), 0.3)
  expect_length(fg12_quality(obs$warns), 0L)

  # A context without variance: the stratum variance and the VPC do not depend on it.
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = shifted,
                                            group = "country", context = "school"))
  own <- Filter(function(f) !f$adjusted, rec$log$fits)
  expect_length(own, 2L)
  for (f in own) {
    expect_true(isTRUE(f$diagnostics$singular), info = f$group)
    expect_identical(f$diagnostics$singular_terms, "school", info = f$group)
    expect_false(f$diagnostics$stratum_singular, info = f$group)
  }
  expect_true(all(obs$res$var_context < 1e-8))
  expect_equal(obs$res$vpc, c(0.3185, 0.3695), tolerance = 1e-3)
  expect_length(fg12_quality(obs$warns), 0L)

  # The same in the crossed-dimensions model with the context.
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = shifted,
                                            group = "country", context = "school",
                                            decomposition = "crossed-dimensions"))
  expect_gt(min(obs$res$vpc), 0.3)
  expect_length(fg12_quality(obs$warns), 0L)

  # A further random effect without variance, supplied in the formula.
  with_strata <- make_strata(shifted, vars = c("gender", "ses"))$data
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | stratum) + (1 | school),
                                            data = with_strata, group = "country"))
  expect_equal(obs$res$var_other, c(0, 0))
  expect_gt(min(obs$res$vpc), 0.3)
  expect_length(fg12_quality(obs$warns), 0L)
})

test_that("the package's own crossed-dimensions example raises no fit-quality warning", {
  real_fit <- fit_maihda
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(math ~ 1 + (1 | gender:ses),
                                            data = maihda_country_data,
                                            group = "country",
                                            decomposition = "crossed-dimensions"))
  # Every country's fit is singular, none in all three terms, and no VPC is zero.
  expect_length(rec$log$fits, 6L)
  for (f in rec$log$fits) {
    expect_true(isTRUE(f$diagnostics$singular), info = f$group)
    expect_false(all(c("gender", "ses", "stratum") %in% f$diagnostics$singular_terms),
                 info = f$group)
  }
  expect_gt(min(obs$res$vpc), 0.05)
  expect_gt(min(obs$res$var_between), 100)
  expect_length(fg12_quality(obs$warns), 0L)
  expect_identical(obs$res$status, rep("ok", 6))
})

test_that("a between-stratum variance of zero is still named, in every mode", {
  real_fit <- fit_maihda
  d <- rbind(fg12_group("g1", fg12_flat), fg12_group("g2", fg12_s2))
  named <- paste0(fg12_open, "singular fit: g1). ", fg12_singular, fg12_vpc)

  # Two-model: g1's stratum variance is the zero.
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country"))
  expect_lt(obs$res$var_between[obs$res$group == "g1"], 1e-8)
  expect_identical(fg12_quality(obs$warns), named)

  # Crossed-dimensions: all three terms of g1 are zero. g2's two dimensions are zero as
  # well, its interaction is not, and g2 is not named.
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country",
                                            decomposition = "crossed-dimensions"))
  expect_setequal(rec$log$fits[[1L]]$diagnostics$singular_terms,
                  c("gender", "ses", "stratum"))
  expect_setequal(rec$log$fits[[2L]]$diagnostics$singular_terms, c("gender", "ses"))
  expect_lt(obs$res$var_between[obs$res$group == "g1"], 1e-8)
  expect_gt(obs$res$var_between[obs$res$group == "g2"], 0.3)
  expect_identical(fg12_quality(obs$warns), named)

  # With a context that has no variance in either group: g1 for its stratum variance,
  # g2 not at all.
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                            group = "country", context = "school"))
  own <- Filter(function(f) !f$adjusted, rec$log$fits)
  expect_setequal(own[[1L]]$diagnostics$singular_terms, c("school", "stratum"))
  expect_identical(own[[2L]]$diagnostics$singular_terms, "school")
  expect_identical(fg12_quality(obs$warns), named)
  expect_identical(obs$res$status, c("ok", "ok"))
})

test_that("a group whose adjusted model did not converge is named for that model", {
  d <- rbind(fg12_group("g2", fg12_s2), fg12_group("g3", fg12_s3))
  clean <- fg12_observe(compare_maihda_groups(y ~ (1 | gender:ses), data = d,
                                              group = "country"))
  expect_length(clean$warns, 0L)
  expect_identical(clean$res$pcv_status, c("ok", "ok"))

  # Each fit is an lme4 fit whose `converged` is replaced: `own` for a group's own
  # model, `adjusted` for the adjusted model of its decomposition. No sampler and no
  # WeMix fit runs; the engine passed to the comparison chooses the words.
  real_fit <- fit_maihda
  with_records <- function(own, adjusted) {
    function(formula, data, ...) {
      m <- real_fit(formula, data, engine = "lme4")
      g <- data$country[1]
      m$diagnostics$converged <-
        if (fg12_is_adjusted(formula)) adjusted[[g]] else own[[g]]
      m
    }
  }
  compare <- function(engine, data = d, ...) {
    args <- list(y ~ (1 | gender:ses), data = data, group = "country",
                 engine = engine, ...)
    if (identical(engine, "wemix")) args$sampling_weights <- "w"
    fg12_observe(do.call(compare_maihda_groups, args))
  }

  for (engine in c("lme4", "wemix", "brms")) {
    # Only g3's adjusted model.
    local_mocked_bindings(fit_maihda = with_records(
      own = list(g2 = TRUE, g3 = TRUE), adjusted = list(g2 = TRUE, g3 = FALSE)))
    obs <- compare(engine)
    expect_identical(
      fg12_quality(obs$warns),
      maihda_group_fit_quality_message(engine, character(0), character(0), "g3"),
      info = engine)
    expect_true(grepl("adjusted model", fg12_quality(obs$warns), fixed = TRUE),
                info = engine)
    expect_false(grepl("g2", fg12_quality(obs$warns), fixed = TRUE), info = engine)
    # The table is what it is without the record: the PCV is reported, as "ok".
    expect_identical(obs$res$pcv_status, c("ok", "ok"), info = engine)
    expect_identical(obs$res$status, c("ok", "ok"), info = engine)
    expect_identical(obs$res$pcv, clean$res$pcv, info = engine)
    expect_identical(obs$res$vpc, clean$res$vpc, info = engine)

    # A group's own model and adjusted models, in the order of the groups.
    local_mocked_bindings(fit_maihda = with_records(
      own = list(g2 = FALSE, g3 = TRUE), adjusted = list(g2 = FALSE, g3 = FALSE)))
    expect_identical(
      fg12_quality(compare(engine)$warns),
      maihda_group_fit_quality_message(engine, character(0), "g2", c("g2", "g3")),
      info = engine)

    # An adjusted record that is unknown (NA) or clean lists nothing.
    local_mocked_bindings(fit_maihda = with_records(
      own = list(g2 = TRUE, g3 = NA), adjusted = list(g2 = NA, g3 = TRUE)))
    expect_length(fg12_quality(compare(engine)$warns), 0L)
  }

  # The words differ by engine, the groups do not.
  local_mocked_bindings(fit_maihda = with_records(
    own = list(g2 = TRUE, g3 = TRUE), adjusted = list(g2 = FALSE, g3 = FALSE)))
  expect_true(grepl("adjusted model did not converge: g2, g3",
                    fg12_quality(compare("lme4")$warns), fixed = TRUE))
  expect_true(grepl(paste0("adjusted model, ", fg12_flag, "g2, g3"),
                    fg12_quality(compare("brms")$warns), fixed = TRUE))

  # A group without a PCV is not named for its adjusted model: g1's strata are
  # identical, its decomposition fails, and the warning for that names it.
  flat <- rbind(fg12_group("g1", fg12_flat), fg12_group("g2", fg12_s2))
  local_mocked_bindings(fit_maihda = with_records(
    own = list(g1 = TRUE, g2 = TRUE), adjusted = list(g1 = FALSE, g2 = FALSE)))
  obs <- compare("lme4", data = flat)
  expect_identical(obs$res$pcv_status, c("failed", "ok"))
  expect_identical(
    fg12_quality(obs$warns),
    maihda_group_fit_quality_message("lme4", "g1", character(0), "g2"))

  # One model a group, so no adjusted model to name: a single stratum dimension, and
  # the crossed-dimensions model.
  one_column <- d
  one_column$stratum <- paste(one_column$gender, one_column$ses, sep = "_")
  one_column <- one_column[, c("country", "stratum", "y")]
  local_mocked_bindings(fit_maihda = with_records(
    own = list(g2 = FALSE, g3 = TRUE), adjusted = list(g2 = FALSE, g3 = FALSE)))
  obs <- fg12_observe(compare_maihda_groups(y ~ (1 | stratum), data = one_column,
                                            group = "country"))
  expect_false("pcv_status" %in% names(obs$res))
  expect_identical(
    fg12_quality(obs$warns),
    maihda_group_fit_quality_message("lme4", character(0), "g2"))
  obs <- compare("lme4", decomposition = "crossed-dimensions")
  expect_identical(
    fg12_quality(obs$warns),
    maihda_group_fit_quality_message("lme4", character(0), "g2"))
})

test_that("an optimiser cap that stops only the adjusted model names the group", {
  skip_on_cran()
  # 3 x 4 strata a group, a binary outcome. With bobyqa in both stages glmer needs
  # under 30 evaluations for a null model here and about 180 or more for an adjusted
  # one, so a cap of 60 stops the adjusted models alone.
  d <- rbind(
    fg12_binary_group("A", c(4, 9, 14, 20, 25, 7, 12, 18, 22, 27, 6, 16)),
    fg12_binary_group("B", c(21, 5, 11, 26, 8, 17, 3, 24, 13, 19, 28, 10)))
  real_fit <- fit_maihda
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_observe(compare_maihda_groups(
    yb ~ (1 | gender:race), data = d, group = "country", family = "binomial",
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 60))))

  expect_length(rec$log$fits, 4L)
  not_converged <- function(r) isFALSE(r$converged)
  own <- fg12_groups(rec$log$fits, adjusted = FALSE, not_converged)
  adjusted <- fg12_groups(rec$log$fits, adjusted = TRUE, not_converged)
  # The premise; if a platform's optimiser differs, move the cap.
  expect_identical(own, character(0))
  expect_identical(adjusted, c("A", "B"))

  # lme4 warns as well, about the cap and the gradient, and names no group.
  fq <- fg12_quality(obs$warns)
  expect_identical(
    fq, maihda_group_fit_quality_message("lme4", character(0), own, adjusted))
  expect_identical(
    fq, paste0(fg12_open, "adjusted model did not converge: A, B). ", fg12_nonconv,
               fg12_pcv))
  others <- setdiff(obs$warns, fq)
  expect_gt(length(others), 0L)
  expect_false(any(grepl("\\b(A|B)\\b", others)))
  # The PCVs are reported, as before.
  expect_true(all(is.finite(obs$res$pcv)))
  expect_identical(obs$res$pcv_status, c("ok", "ok"))
  expect_identical(obs$res$status, c("ok", "ok"))
})

test_that("real brms fits: an adjusted model the diagnostics flag names its group", {
  # Compiles four Stan models, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  d <- rbind(fg12_group("g2", fg12_s2), fg12_group("g3", fg12_s3))
  # Sixty iterations a chain cannot mix: the run is short on purpose.
  real_fit <- fit_maihda
  rec <- fg12_recorder(real_fit)
  local_mocked_bindings(fit_maihda = rec$fit)
  obs <- fg12_keep_seed(fg12_observe(compare_maihda_groups(
    y ~ (1 | gender:ses), data = d, group = "country", engine = "brms",
    chains = 2, iter = 60, warmup = 30, seed = 1, refresh = 0)))

  expect_length(rec$log$fits, 4L)
  expect_true(all(vapply(rec$log$fits, function(f) identical(f$diagnostics$engine, "brms"),
                         logical(1))))
  not_converged <- function(r) isFALSE(r$converged)
  own <- fg12_groups(rec$log$fits, adjusted = FALSE, not_converged)
  adjusted <- fg12_groups(rec$log$fits, adjusted = TRUE, not_converged)
  # The run is too short to pass the diagnostics; if this fails, shorten it further.
  expect_gt(length(adjusted), 0L)
  # A PCV was computed for every group, so every flagged adjusted model is listed.
  expect_true(all(obs$res$pcv_status %in% c("ok", "singular")))
  fq <- fg12_quality(obs$warns)
  expect_identical(
    fq, maihda_group_fit_quality_message("brms", character(0), own, adjusted))
  expect_true(grepl(paste0("adjusted model, ", fg12_flag), fq, fixed = TRUE))
  expect_false(grepl("did not converge", fq, fixed = TRUE))
  expect_false(grepl("singular", fq, fixed = TRUE))
})
