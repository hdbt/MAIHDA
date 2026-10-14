# Audit 2026-10-14 (context-few-levels-brms): a `context` enters every model as a
# random intercept, so its variance is estimated from its levels alone. With few levels
# an engine = "brms" fit leaves that variance to the prior on its standard deviation,
# and the context share and the VPC depend on it; the package advised engine = "brms"
# for exactly that case.
#
# Measured against a known truth (8 strata crossed with C contexts, 100 binomial rows
# a cell, context effects of variance 0.10, brms's default prior): the posterior median
# of the context variance is 0.58 / 0.32 / 0.18 / 0.12 with 2 / 3 / 5 / 10 contexts
# (lme4 on the same data 0.07 / 0.08 / 0.10 / 0.09), the context share 0.13 / 0.08 /
# 0.045 / 0.03 for a true 0.027, and with no context effect at all two contexts still
# get a share of about 0.02 to 0.03. A Gaussian outcome behaves the same way. Under
# normal(0, 0.5) the two-context median is about a quarter of the default's, and
# adapt_delta = 0.99 moves it toward the exact posterior, not toward the truth.
#
# The help pages, maihda_describe(), the group-comparison warning, the README and the
# vignette therefore say what each engine does with few levels instead of recommending
# brms, and the printed partition of a brms fit names its contexts with fewer than ten
# levels. The blocks below pin the count, the rule and the words without Stan, and
# repeat them on real fits when MAIHDA_TEST_BRMS=true.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
cfc_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# Printed lines without colour codes.
cfc_print <- function(expr) cli::ansi_strip(utils::capture.output(expr))

cfc_old_note <- c(
  "  Note: the headline VPC/ICC is the between-stratum share conditional on",
  "  the context random effect(s). The context share is the between-context",
  "  component of the unexplained variance.")

# 2 x 2 strata crossed with three sites, 15 Gaussian rows a cell, built without random
# numbers.
cfc_gaussian_data <- function() {
  g <- expand.grid(d1 = c("a", "b"), d2 = c("p", "q"), site = c("s1", "s2", "s3"),
                   stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 15), ]
  rownames(d) <- NULL
  d$y <- 0.6 * (d$d1 == "b") + 0.4 * (d$d2 == "q") +
    c(s1 = -0.3, s2 = 0, s3 = 0.3)[d$site] + sin(seq_len(nrow(d)) * 12.9898)
  d
}

# A brmsfit stub carrying its own posterior draws, and the maihda_model around it.
cfc_brms_object <- function(data, draws, family = list(family = "gaussian",
                                                        link = "identity")) {
  list(model = structure(list(formula = y ~ 1 + (1 | stratum) + (1 | site),
                              family = family, draws = draws), class = "brmsfit"),
       engine = "brms", data = data)
}

# Three regions seeing 3, 9 and 10 sites, 300 rows each, every site of a region
# present; drawn under a fixed seed.
cfc_compare_data <- function() {
  cfc_with_seed(11, {
    sites <- c(A = 3L, B = 9L, C = 10L)
    d <- do.call(rbind, lapply(names(sites), function(r) {
      data.frame(region = r, site = paste0(r, rep_len(seq_len(sites[[r]]), 300)),
                 stringsAsFactors = FALSE)
    }))
    n <- nrow(d)
    d$g1 <- sample(c("m", "f"), n, replace = TRUE)
    d$g2 <- sample(c("lo", "hi"), n, replace = TRUE)
    d$x <- stats::rnorm(n)
    d$y <- 0.3 * d$x + stats::rnorm(n)
    d
  })
}
cfc_weak_cells <- c("A x site (3 levels)", "B x site (9 levels)")

cfc_draws <- function() {
  data.frame(sd_stratum__Intercept = seq(0.5, 0.7, length.out = 101),
             sd_site__Intercept = seq(0.2, 0.9, length.out = 101),
             sigma = seq(0.9, 1.1, length.out = 101))
}

test_that("a context's levels are counted in the fitted data", {
  d <- data.frame(
    site = c("s1", "s2", NA, "s1", "s3"),
    region = factor(c("n", "s", "n", "s", "n"), levels = c("n", "s", "e", "w")),
    stringsAsFactors = FALSE)
  # A missing value is no level, a declared level without rows is no level, and an
  # absent column is not counted as 0.
  expect_identical(maihda_context_levels(d, c("site", "region", "gone")),
                   c(site = 3L, region = 2L, gone = NA_integer_))
  expect_identical(maihda_context_levels(NULL, "site"), c(site = NA_integer_))
})

test_that("contexts are named for a brms partition only, below ten levels", {
  expect_identical(.maihda_context_few_levels, 10L)
  lv <- c(a = 2L, b = 9L, c = 10L, d = NA_integer_, e = 30L)
  expect_identical(maihda_context_few_levels(list(method = "posterior", context_levels = lv)),
                   c(a = 2L, b = 9L))
  expect_length(maihda_context_few_levels(
    list(method = "posterior", context_levels = c(a = 10L, b = 12L))), 0L)
  # An lme4 partition has no prior to depend on, whatever it carries; a bootstrapped
  # one neither.
  expect_length(maihda_context_few_levels(list(bootstrap = FALSE, context_levels = lv)), 0L)
  expect_length(maihda_context_few_levels(list(method = "bootstrap", context_levels = lv)), 0L)
  # A brms partition without the counts names none.
  expect_length(maihda_context_few_levels(list(method = "posterior")), 0L)
})

test_that("the note names the contexts and says what follows", {
  expect_identical(maihda_context_few_levels_note(integer(0)), character(0))

  one <- maihda_context_few_levels_note(c(ctx = 3L))
  expect_identical(paste(trimws(one), collapse = " "), paste(
    "Note: context 'ctx' has 3 levels. With engine = \"brms\" the variance of such a",
    "context depends on the prior on its standard deviation, and the context share",
    "and the VPC depend on it. See ?maihda."))

  # The clause on two levels is a claim about two-level contexts only.
  two <- maihda_context_few_levels_note(c(country = 2L))
  expect_identical(paste(trimws(two), collapse = " "), paste(
    "Note: context 'country' has 2 levels. With engine = \"brms\" the variance of such",
    "a context depends on the prior on its standard deviation (with two levels roughly",
    "in proportion to that prior's scale), and the context share and the VPC depend on",
    "it. See ?maihda."))

  several <- maihda_context_few_levels_note(c(school = 4L, region = 2L, area = 7L))
  expect_true(startsWith(paste(trimws(several), collapse = " "), paste(
    "Note: contexts 'school', 'region' and 'area' have 4, 2 and 7 levels. With",
    "engine = \"brms\"")))
  expect_true(grepl("with two levels roughly",
                    paste(trimws(several), collapse = " "), fixed = TRUE))
  same <- maihda_context_few_levels_note(c(school = 4L, region = 4L))
  expect_true(startsWith(paste(trimws(same), collapse = " "),
                         "Note: contexts 'school' and 'region' have 4 levels each."))
  expect_false(grepl("with two levels", paste(trimws(same), collapse = " "), fixed = TRUE))
  expect_true(grepl("context 'site' has 1 level.",
                    paste(trimws(maihda_context_few_levels_note(c(site = 1L))), collapse = " "),
                    fixed = TRUE))

  # Laid out like the note above it: indented by two, within the console width.
  for (note in list(one, two, several, same)) {
    expect_true(all(startsWith(note, "  ")))
    expect_false(any(startsWith(note, "   ")))
    expect_true(all(nchar(note) <= 77L))
  }
})

test_that("a brms contextual summary carries the level counts and its partition names the context", {
  skip_if_not_installed("brms")
  draws <- cfc_draws()
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) model$draws)

  d3 <- data.frame(site = c("s1", "s2", "s3", NA, "s1"), stringsAsFactors = FALSE)
  res <- maihda_context_summary_brms(cfc_brms_object(d3, draws), list(context_vars = "site"),
                                     conf_level = 0.95)
  ctx <- res$context
  expect_identical(ctx$context_levels, c(site = 3L))
  out <- cfc_print(maihda_print_context_partition(ctx))
  i <- which(out == cfc_old_note[3])
  expect_length(i, 1L)
  expect_identical(out[i - 2:1], cfc_old_note[1:2])
  note <- maihda_context_few_levels_note(c(site = 3L))
  expect_identical(out[i + seq_along(note)], note)
  expect_identical(out[i + length(note) + 1L], "")
  expect_length(out, i + length(note) + 1L)

  # Twelve sites: the partition prints as it always did, ending in the old note.
  d12 <- data.frame(site = paste0("s", 1:12), stringsAsFactors = FALSE)
  ctx12 <- maihda_context_summary_brms(cfc_brms_object(d12, draws), list(context_vars = "site"),
                                       conf_level = 0.95)$context
  expect_identical(ctx12$context_levels, c(site = 12L))
  out12 <- cfc_print(maihda_print_context_partition(ctx12))
  expect_false(any(grepl("levels. With engine", out12, fixed = TRUE)))
  expect_identical(utils::tail(out12, 4L), c(cfc_old_note, ""))
  # Apart from the note, the two partitions print the same lines.
  expect_identical(out[-(i + seq_along(note))], out12)
})

test_that("a brms crossed-dimensions summary with a context carries the counts too", {
  skip_if_not_installed("brms")
  local_mocked_bindings(maihda_posterior_draws_brms = function(model) model$draws)
  dat <- expand.grid(ga = c("f", "m"), gb = c("x", "y", "z"), site = c("s1", "s2"),
                     stringsAsFactors = FALSE)
  dat <- dat[rep(seq_len(nrow(dat)), each = 3), ]
  dat$stratum <- paste(dat$ga, dat$gb)
  sd_seq <- function(from, to) seq(from, to, length.out = 101)
  draws <- data.frame(b_Intercept = 0, sd_stratum__Intercept = sd_seq(0.05, 0.45),
                      sd_ga__Intercept = sd_seq(0.2, 1.4), sd_gb__Intercept = sd_seq(0.2, 1.2),
                      sd_site__Intercept = sd_seq(0.1, 1.5))
  fam <- list(family = "bernoulli", link = "logit")
  m <- structure(
    list(model = structure(list(formula = y ~ 1, family = fam, draws = draws),
                           class = "brmsfit"),
         engine = "brms", formula = y ~ 1, family = fam, data = dat,
         cc_info = list(dim_groups = c(gender = "ga", region = "gb"),
                        interaction_group = "stratum"),
         context_info = list(context_vars = "site")),
    class = "maihda_model")
  res <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)
  expect_identical(res$context$context_levels, c(site = 2L))
  expect_identical(res$context$method, "posterior")
  out <- cfc_print(maihda_print_context_partition(res$context))
  note <- maihda_context_few_levels_note(c(site = 2L))
  i <- which(out == cfc_old_note[3])
  expect_identical(out[i + seq_along(note)], note)
})

test_that("print() of an analysis and of its summary show the note; an lme4 print has none", {
  d <- cfc_gaussian_data()
  a <- cfc_with_seed(1, suppressWarnings(suppressMessages(
    maihda(y ~ 1 + (1 | d1:d2), data = d, context = "site"))))
  # An lme4 partition is what it was: no level counts, no note.
  expect_null(a$summary$context$context_levels)
  out_a <- cfc_print(print(a))
  out_s <- cfc_print(print(a$summary))
  expect_false(any(grepl("levels. With engine", c(out_a, out_s), fixed = TRUE)))
  # The counts are read from a real fit's frame.
  expect_identical(maihda_context_levels(a$model$data, "site"), c(site = 3L))

  # The same analysis with the partition of a brms fit: the analysis prints the note
  # directly below its context share line, the summary below the partition's note, and
  # nothing else in either print moves.
  b <- a
  b$summary$context$method <- "posterior"
  b$summary$context$context_levels <- c(site = 3L)
  note <- maihda_context_few_levels_note(c(site = 3L))
  out_b <- cfc_print(print(b))
  i <- grep("^Context share \\(null\\):", out_b)
  expect_length(i, 1L)
  expect_identical(out_b[i + seq_along(note)], note)
  expect_identical(out_b[-(i + seq_along(note))], out_a)
  out_bs <- cfc_print(print(b$summary))
  j <- which(out_bs == cfc_old_note[3])
  expect_length(j, 1L)
  expect_identical(out_bs[j + seq_along(note)], note)
  expect_identical(out_bs[-(j + seq_along(note))], out_s)
})

test_that("the group-comparison warning says what the fitting engine does", {
  cells <- c("A x site (3 levels)", "B x site (2 levels)")
  w_lme4 <- maihda_weak_context_message("lme4", cells)
  w_brms <- maihda_weak_context_message("brms", cells)
  expect_identical(w_lme4, paste0(
    "compare_maihda_groups(): few context levels within group(s) (A x site (3 levels); ",
    "B x site (2 levels)). A context variance estimated from fewer than 10 levels is ",
    "imprecise, and at zero (a singular fit) when the levels differ little; interpret ",
    "those groups' context partition with caution."))
  expect_identical(w_brms, paste0(
    "compare_maihda_groups(): few context levels within group(s) (A x site (3 levels); ",
    "B x site (2 levels)). With fewer than 10 levels a context's variance depends on the ",
    "prior on its standard deviation, and the context share with it; interpret those ",
    "groups' context partition with caution."))
  # Neither engine is sent to the other.
  expect_false(any(grepl("engine =", c(w_lme4, w_brms), fixed = TRUE)))
  expect_null(maihda_weak_context_message("lme4", character(0)))

  # End to end on lme4: regions with 3 and 9 sites are named, the one with 10 is not.
  w <- capture_warnings(compare_maihda_groups(y ~ x + (1 | g1:g2), data = cfc_compare_data(),
                                              group = "region", context = "site"))
  expect_true(maihda_weak_context_message("lme4", cfc_weak_cells) %in% w)
  expect_false(maihda_weak_context_message("brms", cfc_weak_cells) %in% w)
})

test_that("a brms comparison's warning is worded for brms", {
  # Each group's fit is an lme4 fit made inside a brms comparison, so the comparison's
  # own engine reaches the warning without a sampler.
  real_fit <- fit_maihda
  local_mocked_bindings(fit_maihda = function(formula, data, ...) {
    dots <- list(...)
    dots$engine <- "lme4"
    do.call(real_fit, c(list(formula, data), dots))
  })
  w <- capture_warnings(compare_maihda_groups(y ~ x + (1 | g1:g2), data = cfc_compare_data(),
                                              group = "region", context = "site",
                                              engine = "brms"))
  expect_true(maihda_weak_context_message("brms", cfc_weak_cells) %in% w)
  expect_false(maihda_weak_context_message("lme4", cfc_weak_cells) %in% w)
})

test_that("maihda_describe() says what each engine does with a few-level context", {
  data(maihda_country_data, envir = environment())
  desc <- maihda_describe(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
                          context = "country")
  msg <- desc$warnings$message[desc$warnings$check == "context_few_levels"]
  expect_identical(msg, paste0(
    "Context 'country' has 6 levels; with fewer than ~10 levels the context variance is ",
    "imprecise: lme4 may estimate it at zero (a singular fit), and a brms fit's posterior ",
    "depends on the prior on its standard deviation."))
  expect_false(grepl("Consider engine", msg, fixed = TRUE))
})

test_that("the help pages, README and vignettes no longer send few-level contexts to brms", {
  skip_on_cran()
  root <- testthat::test_path("..", "..")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(file.path(root, "man")), "the package tree is not available")
  flat <- function(...) {
    gsub("[[:space:]]+", " ", paste(readLines(file.path(root, ...), warn = FALSE),
                                    collapse = " "))
  }
  fit_rd <- flat("man", "fit_maihda.Rd")
  maihda_rd <- flat("man", "maihda.Rd")
  expect_false(grepl("the \\code{brms} engine handles this better", fit_rd, fixed = TRUE))
  expect_false(grepl("weakly identifies its variance (consider", maihda_rd, fixed = TRUE))
  for (rd in list(fit_rd, maihda_rd)) {
    expect_true(grepl("identifies its variance weakly with either engine", rd, fixed = TRUE))
    expect_true(grepl("depends on the prior on its standard deviation", rd, fixed = TRUE))
  }
  # ?maihda states the number of levels the printed note uses.
  expect_true(grepl(sprintf("a \\code{brms} fit that have fewer than %d levels, and the contextual",
                            .maihda_context_few_levels), maihda_rd, fixed = TRUE))
  readme <- file.path(root, "README.md")
  if (file.exists(readme)) {
    expect_false(grepl("prefer many-level contexts or `engine = \"brms\"`", flat("README.md"),
                       fixed = TRUE))
  }
  cc <- file.path(root, "vignettes", "cross_classified.Rmd")
  if (file.exists(cc)) {
    v <- flat("vignettes", "cross_classified.Rmd")
    expect_false(grepl("whose priors regularise the variance", v, fixed = TRUE))
    expect_true(grepl("names the contexts with fewer than 10 levels", v, fixed = TRUE))
  }
  sp <- file.path(root, "vignettes", "bayesian_sparse_maihda.Rmd")
  if (file.exists(sp)) {
    expect_true(grepl("Refitted under brms's default prior instead",
                      flat("vignettes", "bayesian_sparse_maihda.Rmd"), fixed = TRUE))
  }
})

test_that("real brms fits: the level counts, the printed note and a forwarded prior", {
  # Compiles two Stan models, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # 2 x 2 strata crossed with two contexts, 50 binomial trials a cell. brms draws its
  # sampler seed from R's stream, so the fits sit inside the seeded block too.
  fits <- cfc_with_seed(20261014, {
    g <- expand.grid(d1 = c("a", "b"), d2 = c("p", "q"), ctx = c("k1", "k2"),
                     stringsAsFactors = FALSE)
    eta <- -0.6 + 0.5 * (g$d1 == "b") + 0.4 * (g$d2 == "q") +
      ifelse(g$ctx == "k2", 0.45, 0) + stats::rnorm(nrow(g), 0, 0.1)
    g$n <- 50L
    g$s <- stats::rbinom(nrow(g), 50L, stats::plogis(eta))
    fit <- function(...) {
      suppressWarnings(suppressMessages(
        fit_maihda(s | trials(n) ~ 1 + (1 | d1:d2), data = g, context = "ctx",
                   engine = "brms", family = "binomial", chains = 2, iter = 1000,
                   refresh = 0, ...)))
    }
    list(default = fit(),
         tight = fit(prior = brms::set_prior("normal(0, 0.5)", class = "sd")))
  })
  m <- fits$default
  # summary() passes on brms's own warnings about the sampler (divergent transitions).
  s <- suppressWarnings(summary(m))

  # The counts are the engine's own: brms's number of levels of the grouping factor.
  expect_identical(s$context$context_levels, c(ctx = 2L))
  expect_identical(unname(s$context$context_levels),
                   as.integer(brms::ngrps(m$model)[["ctx"]]))
  note <- maihda_context_few_levels_note(c(ctx = 2L))
  out <- cfc_print(suppressWarnings(print(s)))
  i <- which(out == cfc_old_note[3])
  expect_length(i, 1L)
  expect_identical(out[i + seq_along(note)], note)

  # `prior` reaches brms::brm(), and with two contexts it moves the context variance
  # roughly in proportion to its scale: a fifth of the default scale, well under half
  # the default posterior median.
  ps <- as.data.frame(brms::prior_summary(fits$tight$model))
  expect_identical(ps$prior[ps$class == "sd" & ps$source == "user"], "normal(0, 0.5)")
  st <- suppressWarnings(summary(fits$tight))
  expect_lt(st$context$per_context[["ctx"]], 0.5 * s$context$per_context[["ctx"]])
  expect_lt(st$context$vpc_context_total, s$context$vpc_context_total)
})
