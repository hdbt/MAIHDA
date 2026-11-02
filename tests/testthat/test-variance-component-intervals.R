# Intervals for the variance components of summary() (GitHub issue #74).
#
# Each row of the variance-components table that is an estimated variance or covariance
# carries an interval at conf_level in `lower` / `upper`:
#   lme4, bootstrap = TRUE: the alpha / 2 and 1 - alpha / 2 quantiles (type 7) of the
#     row's value over the parametric-bootstrap refits that form the VPC interval
#     (stats::simulate() from the fit, lme4::refit() per draw), so no refit is added;
#   brms: the same quantiles of the posterior draws whose median (a growth model: mean)
#     the row reports.
# No interval for the Total row, nor for a level-1 variance the model does not estimate:
# the latent pi^2 / 3 of a logit link and a count model's approximation at the marginal
# expected count. tidy(component = "variance") always has conf.low / conf.high (NA where
# there is none). maihda_table()'s between-stratum variance and SD rows carry the
# interval of that variance (the SD row its square roots): the table row for a
# single-stratum or contextual fit, decomposition$between_var_ci (dimensions plus
# interaction, per replicate) for a crossed-dimensions fit, and
# longitudinal$var_stratum_ref_ci (a(t0)' Sigma a(t0) per replicate) for a growth fit.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
vci_with_seed <- function(seed, expr) {
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

# The value of `expr` and the messages of the warnings it raised.
vci_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

vci_q <- function(v, level = 0.95) {
  a <- 1 - level
  stats::quantile(v, c(a / 2, 1 - a / 2), names = FALSE)
}

# Twelve unbalanced strata (gender x race x edu), a covariate, precision weights, a
# six-level context, and Gaussian, binary and count outcomes.
vci_data <- function() {
  vci_with_seed(20261101, {
    d <- expand.grid(gender = c("f", "m"), race = c("a", "b", "c"), edu = c("lo", "hi"),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    d <- d[rep(seq_len(nrow(d)), c(40, 12, 35, 18, 45, 10, 30, 22, 38, 15, 28, 25)), ]
    cell <- interaction(d$gender, d$race, d$edu, drop = TRUE)
    d$x <- stats::rnorm(nrow(d))
    d$w <- stats::runif(nrow(d), 0.5, 2)
    d$site <- sample(paste0("s", 1:6), nrow(d), replace = TRUE)
    u <- stats::rnorm(nlevels(cell), 0, 0.7)[cell]
    v <- stats::rnorm(6, 0, 0.5)[match(d$site, paste0("s", 1:6))]
    d$y <- 1 + 0.3 * d$x + u + v + stats::rnorm(nrow(d))
    d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(-0.3 + u))
    d$cnt <- stats::rpois(nrow(d), exp(0.3 + 0.5 * u))
    rownames(d) <- NULL
    d
  })
}

# Every variance of a refit, named by its grouping factor ("Residual" for sigma^2).
vci_refit_variances <- function(r) {
  vc <- as.data.frame(lme4::VarCorr(r))
  vc <- vc[is.na(vc$var2), ]
  stats::setNames(vc$vcov, vc$grp)
}

# The package's draw sequence written out: stats::simulate() from the fit after
# set.seed(seed), lme4::refit() per draw. One row of variances per refit.
vci_ref_refits <- function(fit, n_boot, seed) {
  nm <- names(vci_refit_variances(fit))
  vci_with_seed(seed, {
    sims <- stats::simulate(fit, nsim = n_boot)
    v <- vapply(seq_len(n_boot), function(i) {
      unname(vci_refit_variances(suppressWarnings(suppressMessages(
        lme4::refit(fit, sims[[i]])))))
    }, numeric(length(nm)))
    matrix(v, nrow = n_boot, byrow = TRUE, dimnames = list(NULL, nm))
  })
}

# summary(bootstrap = TRUE) under a seed; returns the summary and its warnings.
vci_boot_summary <- function(m, n_boot, seed, ...) {
  vci_with_seed(seed, vci_warnings(summary(m, bootstrap = TRUE, n_boot = n_boot, ...)))
}

vci_row <- function(tab, label) tab[tab$component == label, , drop = FALSE]

# ---- the pure helpers --------------------------------------------------------------

test_that("maihda_variance_draw_intervals() returns type-7 quantiles per component", {
  a <- c(0.5, 0.1, 0.9, 0.3, 0.7, 0.2, 0.8, 0.4, 0.6, 1.0, 0.05)
  b <- rev(seq(2, 3, length.out = 11))
  iv <- maihda_variance_draw_intervals(list(first = a, second = b), 0.95)
  expect_identical(iv$component, c("first", "second"))
  expect_equal(iv$lower, c(stats::quantile(a, 0.025, names = FALSE),
                           stats::quantile(b, 0.025, names = FALSE)), tolerance = 1e-15)
  expect_equal(iv$upper, c(stats::quantile(a, 0.975, names = FALSE),
                           stats::quantile(b, 0.975, names = FALSE)), tolerance = 1e-15)
  # The type-7 quantile interpolates: 0.025 of 11 values sits a quarter of the way
  # from the smallest to the second smallest.
  expect_equal(iv$lower[1], 0.05 + 0.25 * (0.1 - 0.05), tolerance = 1e-15)

  iv90 <- maihda_variance_draw_intervals(list(first = a), 0.90)
  expect_equal(c(iv90$lower, iv90$upper), vci_q(a, 0.90), tolerance = 1e-15)

  # A length-1 element is recycled; a replicate that is not finite in one element is
  # dropped from every element.
  iv_s <- maihda_variance_draw_intervals(list(first = a, fixed = 3), 0.95)
  expect_identical(c(iv_s$lower[2], iv_s$upper[2]), c(3, 3))
  a_na <- a
  a_na[3] <- NA
  iv_na <- maihda_variance_draw_intervals(list(first = a_na, second = b), 0.95)
  expect_equal(c(iv_na$lower[2], iv_na$upper[2]), vci_q(b[-3]), tolerance = 1e-15)
  expect_equal(c(iv_na$lower[1], iv_na$upper[1]), vci_q(a[-3]), tolerance = 1e-15)
  iv_none <- maihda_variance_draw_intervals(list(first = c(NA, Inf)), 0.95)
  expect_true(is.na(iv_none$lower) && is.na(iv_none$upper))
})

test_that("maihda_add_variance_intervals() matches rows by label and keeps the table", {
  tab <- maihda_context_components_table(0.4, c(site = 0.2), 0, 1.5)
  draws <- list(`Context: site` = c(0.1, 0.2, 0.3), `Between-stratum (random)` = c(1, 2, 3))
  out <- maihda_add_variance_intervals(tab, draws, 0.5)
  expect_identical(attr(out, "kind"), "contextual")
  expect_identical(attr(out, "n_contexts"), 1L)
  expect_identical(out[c("component", "variance", "sd", "proportion")],
                   tab[c("component", "variance", "sd", "proportion")])
  expect_equal(unlist(vci_row(out, "Between-stratum (random)")[c("lower", "upper")]),
               c(lower = 1.5, upper = 2.5))
  expect_equal(unlist(vci_row(out, "Context: site")[c("lower", "upper")]),
               c(lower = 0.15, upper = 0.25))
  # Rows without replicates -- the residual left out, the Total -- are NA.
  expect_true(all(is.na(unlist(vci_row(out, "Within-stratum (residual)")[c("lower", "upper")]))))
  expect_true(all(is.na(unlist(vci_row(out, "Total")[c("lower", "upper")]))))
})

test_that("only a Gaussian level-1 variance gets replicates", {
  fam <- function(f, l) structure(list(family = structure(list(family = f, link = l),
                                                          class = "family")),
                                  class = "brmsfit")
  expect_true(maihda_level1_is_estimated(fam("gaussian", "identity")))
  for (f in list(c("bernoulli", "logit"), c("binomial", "probit"), c("cumulative", "logit"),
                 c("poisson", "log"), c("negbinomial", "log"))) {
    expect_false(maihda_level1_is_estimated(fam(f[1], f[2])), info = f[1])
  }

  m <- matrix(1:6, 2, 3, dimnames = list(NULL, c("Between-stratum (random)",
                                                  "Within-stratum (residual)",
                                                  "Context: site")))
  kept <- maihda_component_draw_list(m, TRUE)
  expect_identical(names(kept), colnames(m))
  expect_identical(kept[["Context: site"]], c(5L, 6L))
  expect_identical(names(maihda_component_draw_list(m, FALSE)),
                   c("Between-stratum (random)", "Context: site"))
  lng <- list(`Between-stratum: intercept (time = 0)` = 1, `Within (residual)` = 2)
  expect_identical(names(maihda_component_draw_list(lng, FALSE)),
                   "Between-stratum: intercept (time = 0)")
  std <- maihda_standard_variance_draws(1, 2, 3, FALSE)
  expect_identical(names(std), c("Between-stratum (random)", "Other random effects"))
  expect_identical(maihda_standard_variance_draws(1, 2, 3, TRUE)[["Within-stratum (residual)"]], 3)
})

test_that("a growth block's entries come in the components table's row order", {
  s1 <- matrix(0.7, 1, 1)
  s2 <- matrix(c(0.7, -0.1, -0.1, 0.2), 2, 2)
  s3 <- matrix(c(0.9, 0.1, -0.05, 0.1, 0.3, 0.02, -0.05, 0.02, 0.04), 3, 3)
  for (S in list(s1, s2, s3)) {
    tab <- maihda_longitudinal_components_table(S, s3, 0.5, "wave", "id")
    expect_equal(unname(c(maihda_lng_block_entries(S), maihda_lng_block_entries(s3), 0.5)),
                 tab$variance, tolerance = 0)
  }
})

test_that("the print note names the interval's kind and level on one line", {
  tab <- maihda_variance_components_table(0.4, 0, 1)
  s <- list(variance_components = tab, conf_level = 0.9,
            vpc = list(method = "bootstrap"))
  expect_null(maihda_variance_interval_label(s))
  s$variance_components <- maihda_add_variance_intervals(
    tab, list(`Between-stratum (random)` = c(0.1, 0.5)), 0.9)
  expect_identical(maihda_variance_interval_label(s),
                   "(Bootstrap 90% percentile intervals of the model's variance parameters)")
  s$vpc$method <- "posterior"
  s$conf_level <- 0.999
  expect_identical(maihda_variance_interval_label(s),
                   "(Posterior 99.9% credible intervals of the model's variance parameters)")
})

# ---- lme4: the single-stratum bootstrap -------------------------------------------

test_that("an lme4 bootstrap gives each estimated variance the percentile of its refits", {
  d <- vci_data()
  m <- fit_maihda(y ~ x + (1 | gender:race:edu), data = d)
  s0 <- summary(m)
  expect_false(any(c("lower", "upper") %in% names(s0$variance_components)))

  b <- vci_boot_summary(m, 40, 7)
  expect_identical(length(b$warnings), 1L)          # n_boot = 40 is low (and only that)
  expect_match(b$warnings, "n_boot = 40 is low", fixed = TRUE)
  s <- b$value
  vc <- s$variance_components
  ref <- vci_ref_refits(m$model, 40, 7)
  expect_equal(unlist(vci_row(vc, "Between-stratum (random)")[c("lower", "upper")]),
               c(lower = vci_q(ref[, "stratum"])[1], upper = vci_q(ref[, "stratum"])[2]),
               tolerance = 1e-10)
  expect_equal(unlist(vci_row(vc, "Within-stratum (residual)")[c("lower", "upper")]),
               c(lower = vci_q(ref[, "Residual"])[1], upper = vci_q(ref[, "Residual"])[2]),
               tolerance = 1e-10)
  expect_true(all(is.na(unlist(vci_row(vc, "Total")[c("lower", "upper")]))))
  # The same refits as the VPC interval, and nothing else in the summary moved.
  share <- ref[, "stratum"] / (ref[, "stratum"] + ref[, "Residual"])
  expect_equal(c(s$vpc$ci_lower, s$vpc$ci_upper), vci_q(share), tolerance = 1e-10)
  expect_identical(vc[c("component", "variance", "sd", "proportion")],
                   s0$variance_components)
  expect_identical(s$vpc$estimate, s0$vpc$estimate)

  tv <- tidy(s, component = "variance")
  expect_identical(names(tv), c("component", "variance", "sd", "proportion",
                                "conf.low", "conf.high"))
  expect_identical(tv$conf.low, vc$lower)
  expect_identical(tv$conf.high, vc$upper)
  out <- utils::capture.output(print(s))
  expect_identical(sum(grepl("(Bootstrap 95% percentile intervals of the model's variance parameters)",
                             out, fixed = TRUE)), 1L)

  # conf_level sets the level of these intervals as of the VPC's.
  b90 <- vci_boot_summary(m, 40, 7, conf_level = 0.9)
  vc90 <- b90$value$variance_components
  expect_equal(vc90$lower[1:2], c(vci_q(ref[, "stratum"], 0.9)[1],
                                  vci_q(ref[, "Residual"], 0.9)[1]), tolerance = 1e-10)
  expect_equal(vc90$upper[1:2], c(vci_q(ref[, "stratum"], 0.9)[2],
                                  vci_q(ref[, "Residual"], 0.9)[2]), tolerance = 1e-10)
})

test_that("a weighted lme4 fit's residual interval is that of mean(sigma^2 / w)", {
  d <- vci_data()
  m <- fit_maihda(y ~ x + (1 | gender:race:edu), data = d, weights = w)
  b <- vci_boot_summary(m, 40, 8)
  expect_identical(length(b$warnings), 1L)
  vc <- b$value$variance_components
  # The package draws each residual at sigma / sqrt(w_i) around a residual-free draw.
  fit <- m$model
  w <- stats::weights(fit)
  ref <- vci_with_seed(8, {
    sims <- stats::simulate(fit, nsim = 40, cond.sim = FALSE)
    sd_i <- stats::sigma(fit) / sqrt(w)
    for (i in seq_len(ncol(sims))) sims[[i]] <- sims[[i]] + stats::rnorm(length(w), 0, sd_i)
    t(vapply(seq_len(40), function(i) {
      r <- suppressWarnings(suppressMessages(lme4::refit(fit, sims[[i]])))
      c(vci_refit_variances(r)[["stratum"]], mean(stats::sigma(r)^2 / w))
    }, numeric(2)))
  })
  expect_equal(c(vc$lower[2], vc$upper[2]), vci_q(ref[, 2]), tolerance = 1e-10)
  expect_equal(c(vc$lower[1], vc$upper[1]), vci_q(ref[, 1]), tolerance = 1e-10)
  # The point estimate is the same mean, at the fit.
  expect_equal(vc$variance[2], mean(stats::sigma(fit)^2 / w), tolerance = 1e-10)
})

test_that("binomial and Poisson lme4 bootstraps leave the level-1 row without an interval", {
  skip_slow_test()
  d <- vci_data()
  mb <- fit_maihda(yb ~ 1 + (1 | gender:race:edu), data = d, family = "binomial")
  bb <- vci_boot_summary(mb, 30, 9)
  expect_identical(sum(grepl("n_boot = 30 is low", bb$warnings, fixed = TRUE)), 1L)
  vc <- bb$value$variance_components
  ref <- vci_ref_refits(mb$model, 30, 9)
  expect_equal(c(vc$lower[1], vc$upper[1]), vci_q(ref[, "stratum"]), tolerance = 1e-10)
  expect_identical(vc$component[2], "Within-stratum (residual)")
  expect_true(is.na(vc$lower[2]) && is.na(vc$upper[2]))
  expect_true(is.na(vc$lower[3]) && is.na(vc$upper[3]))
  tv <- tidy(bb$value, component = "variance")
  expect_identical(is.na(tv$conf.low), c(FALSE, TRUE, TRUE))

  mp <- fit_maihda(cnt ~ 1 + (1 | gender:race:edu), data = d, family = "poisson")
  bp <- vci_boot_summary(mp, 30, 10)
  expect_identical(sum(grepl("n_boot = 30 is low", bp$warnings, fixed = TRUE)), 1L)
  vcp <- bp$value$variance_components
  refp <- vci_ref_refits(mp$model, 30, 10)
  expect_equal(c(vcp$lower[1], vcp$upper[1]), vci_q(refp[, "stratum"]), tolerance = 1e-10)
  expect_true(is.na(vcp$lower[2]) && is.na(vcp$upper[2]))
})

test_that("without an interval the table is unchanged and tidy() reports NA", {
  d <- vci_data()
  m <- fit_maihda(y ~ x + (1 | gender:race:edu), data = d)
  s <- summary(m)
  expect_identical(names(s$variance_components),
                   c("component", "variance", "sd", "proportion"))
  tv <- tidy(m, component = "variance")
  expect_identical(names(tv), c("component", "variance", "sd", "proportion",
                                "conf.low", "conf.high"))
  expect_true(all(is.na(tv$conf.low)) && all(is.na(tv$conf.high)))
  expect_identical(tv$variance, s$variance_components$variance)
  out <- utils::capture.output(print(s))
  expect_false(any(grepl("intervals of the model's variance parameters", out, fixed = TRUE)))
  tab <- maihda_table(m)$models
  expect_true(all(is.na(unlist(tab[tab$statistic %in% c("Between-stratum variance",
                                                         "Between-stratum SD"),
                                   c("estimate_lower", "estimate_upper")]))))
})

test_that("the WeMix and ordinal engines report the same point estimates with NA intervals", {
  skip_on_cran()
  d <- vci_data()
  d$yord <- factor(cut(d$y, stats::quantile(d$y, c(0, 1 / 3, 2 / 3, 1)),
                       include.lowest = TRUE, labels = 1:3), ordered = TRUE)
  fits <- list()
  if (requireNamespace("WeMix", quietly = TRUE)) {
    fits$wemix <- suppressWarnings(suppressMessages(
      fit_maihda(y ~ x + (1 | gender:race:edu), data = d, engine = "wemix",
                 sampling_weights = "w")))
  }
  if (requireNamespace("ordinal", quietly = TRUE)) {
    fits$ordinal <- suppressWarnings(suppressMessages(
      fit_maihda(yord ~ x + (1 | gender:race:edu), data = d, family = "ordinal")))
  }
  skip_if(length(fits) == 0, "neither WeMix nor ordinal is installed")
  for (nm in names(fits)) {
    s <- summary(fits[[nm]])
    expect_identical(names(s$variance_components),
                     c("component", "variance", "sd", "proportion"), info = nm)
    tv <- tidy(s, component = "variance")
    expect_identical(names(tv), c("component", "variance", "sd", "proportion",
                                  "conf.low", "conf.high"), info = nm)
    expect_identical(tv$variance, s$variance_components$variance, info = nm)
    expect_true(all(is.na(tv$conf.low)) && all(is.na(tv$conf.high)), info = nm)
  }
})

# ---- lme4: contextual, crossed-dimensions and growth designs -----------------------

test_that("a contextual lme4 bootstrap gives the stratum and context rows their refits' interval", {
  d <- vci_data()
  m <- fit_maihda(y ~ x + (1 | gender:race:edu), data = d, context = "site")
  b <- vci_boot_summary(m, 40, 11)
  expect_identical(length(b$warnings), 1L)
  vc <- b$value$variance_components
  expect_identical(attr(vc, "kind"), "contextual")
  ref <- vci_ref_refits(m$model, 40, 11)
  for (lab in c("Between-stratum (random)", "Context: site", "Within-stratum (residual)")) {
    col <- c(`Between-stratum (random)` = "stratum", `Context: site` = "site",
             `Within-stratum (residual)` = "Residual")[[lab]]
    expect_equal(unlist(vci_row(vc, lab)[c("lower", "upper")]),
                 c(lower = vci_q(ref[, col])[1], upper = vci_q(ref[, col])[2]),
                 tolerance = 1e-10, info = lab)
  }
  # maihda_table()'s between-stratum rows read the stratum row; the SD is its root.
  st <- maihda_collect_model_stats(b$value)
  row <- vci_row(vc, "Between-stratum (random)")
  expect_equal(st[["Between-stratum variance"]], c(row$variance, row$lower, row$upper),
               tolerance = 1e-12)
  expect_equal(st[["Between-stratum SD"]], sqrt(c(row$variance, row$lower, row$upper)),
               tolerance = 1e-12)
  # tidy() keeps the table's attributes beside the new columns.
  tv <- tidy(b$value, component = "variance")
  expect_identical(attr(tv, "kind"), "contextual")
  expect_identical(attr(tv, "n_contexts"), attr(vc, "n_contexts"))
  expect_identical(tv$conf.low, vc$lower)
})

test_that("a crossed-dimensions lme4 bootstrap fills each row and the between-strata interval", {
  skip_slow_test()
  d <- vci_data()
  for (ctx in list(NULL, "site")) {
    info <- if (is.null(ctx)) "no context" else "context"
    cc <- suppressMessages(maihda(y ~ x + (1 | gender:race:edu), data = d, context = ctx,
                                  decomposition = "crossed-dimensions"))
    b <- vci_boot_summary(cc$model, 40, 12)
    expect_identical(sum(grepl("n_boot = 40 is low", b$warnings, fixed = TRUE)), 1L,
                     info = info)
    s <- b$value
    vc <- s$variance_components
    dims <- cc$model$cc_info$dim_groups
    ref <- vci_ref_refits(cc$model$model, 40, 12)
    for (k in names(dims)) {
      expect_equal(unlist(vci_row(vc, paste0("Additive: ", k))[c("lower", "upper")]),
                   c(lower = vci_q(ref[, dims[[k]]])[1], upper = vci_q(ref[, dims[[k]]])[2]),
                   tolerance = 1e-10, info = paste(info, k))
    }
    ig <- cc$model$cc_info$interaction_group
    expect_equal(unlist(vci_row(vc, "Intersectional interaction")[c("lower", "upper")]),
                 c(lower = vci_q(ref[, ig])[1], upper = vci_q(ref[, ig])[2]),
                 tolerance = 1e-10, info = info)
    if (!is.null(ctx)) {
      expect_equal(unlist(vci_row(vc, "Context: site")[c("lower", "upper")]),
                   c(lower = vci_q(ref[, "site"])[1], upper = vci_q(ref[, "site"])[2]),
                   tolerance = 1e-10, info = info)
    }
    expect_equal(unlist(vci_row(vc, "Within-stratum (residual)")[c("lower", "upper")]),
                 c(lower = vci_q(ref[, "Residual"])[1], upper = vci_q(ref[, "Residual"])[2]),
                 tolerance = 1e-10, info = info)
    between <- rowSums(ref[, c(unname(dims), ig), drop = FALSE])
    expect_equal(s$decomposition$between_var_ci, vci_q(between), tolerance = 1e-10,
                 info = info)

    st <- maihda_collect_model_stats(s)
    expect_equal(st[["Between-stratum variance"]],
                 c(s$decomposition$between_var, vci_q(between)), tolerance = 1e-10, info = info)
    expect_equal(st[["Between-stratum SD"]],
                 sqrt(c(s$decomposition$between_var, vci_q(between))), tolerance = 1e-10,
                 info = info)
  }
  # Without the bootstrap the decomposition has no such interval and the table none.
  s0 <- summary(cc$model)
  expect_null(s0$decomposition$between_var_ci)
  expect_true(all(is.na(maihda_collect_model_stats(s0)[["Between-stratum variance"]][2:3])))
})

test_that("a growth lme4 bootstrap fills every block entry and the reference-time variance", {
  skip_on_cran()
  skip_slow_test()
  dl <- vci_with_seed(20261102, {
    g <- expand.grid(a = c("p", "q"), b = c("r", "s", "t"), KEEP.OUT.ATTRS = FALSE,
                     stringsAsFactors = FALSE)
    people <- g[rep(seq_len(nrow(g)), each = 8), ]
    people$id <- sprintf("i%02d", seq_len(nrow(people)))
    st <- interaction(people$a, people$b, drop = TRUE)
    u0 <- stats::rnorm(nlevels(st), 0, 0.6)[st]
    u1 <- stats::rnorm(nlevels(st), 0, 0.15)[st]
    e0 <- stats::rnorm(nrow(people), 0, 0.8)
    e1 <- stats::rnorm(nrow(people), 0, 0.2)
    long <- people[rep(seq_len(nrow(people)), each = 4), ]
    long$wave <- rep(0:3, nrow(people))
    k <- rep(seq_len(nrow(people)), each = 4)
    long$y <- 2 + (u0 + e0)[k] + (0.3 + u1 + e1)[k] * long$wave + stats::rnorm(nrow(long), 0, 0.5)
    long
  })
  m <- suppressMessages(fit_maihda(y ~ wave + (1 | a:b), data = dl, id = "id", time = "wave"))
  n_boot <- 20
  b <- vci_boot_summary(m, n_boot, 13)
  expect_identical(sum(grepl("n_boot = 20 is low", b$warnings, fixed = TRUE)), 1L)
  s <- b$value
  vc <- s$variance_components
  lng <- m$longitudinal_info
  expect_true(all(is.finite(vc$lower)) && all(is.finite(vc$upper)))

  # Replay the package's refits: its draws, refitted as fit_maihda() fits the model.
  spec <- maihda_lme4_growth_spec(m$model, lng)
  tt <- maihda_lng_time_term(lng)
  sc <- maihda_lng_time_scale(lng)
  ref_c <- lng$ref_time - maihda_lng_time_center(lng)
  ref <- vci_with_seed(13, {
    sims <- maihda_simulate_lme4(m$model, nsim = n_boot)
    t(vapply(seq_len(n_boot), function(i) {
      bm <- suppressWarnings(suppressMessages(maihda_refit_growth_draw(m$model, sims[[i]], spec)))
      Ss <- maihda_re_block_lme4(bm, "stratum", tt, maihda_lng_stratum_degree(lng), scale = sc)
      Si <- maihda_re_block_lme4(bm, lng$id, tt, lng$time_degree, scale = sc)
      c(maihda_lng_block_entries(Ss), maihda_lng_block_entries(Si), stats::sigma(bm)^2,
        maihda_var_at_time(Ss, ref_c))
    }, numeric(nrow(vc) + 1L)))
  })
  for (j in seq_len(nrow(vc))) {
    expect_equal(c(vc$lower[j], vc$upper[j]), vci_q(ref[, j]), tolerance = 1e-8,
                 info = vc$component[j])
  }
  expect_true(any(grepl("covariance", vc$component)))
  expect_equal(s$longitudinal$var_stratum_ref_ci, vci_q(ref[, nrow(vc) + 1L]),
               tolerance = 1e-8)
  st <- maihda_collect_model_stats(s)
  expect_equal(st[["Between-stratum variance"]][2:3], s$longitudinal$var_stratum_ref_ci,
               tolerance = 1e-12)
  tv <- tidy(s, component = "variance")
  expect_identical(names(tv), c("component", "variance", "sd", "conf.low", "conf.high"))
  # Without the bootstrap the reference-time interval is NA.
  expect_identical(summary(m)$longitudinal$var_stratum_ref_ci, c(NA_real_, NA_real_))
})

test_that("a bootstrapped maihda() analysis carries the intervals to tidy() and maihda_table()", {
  skip_slow_test()
  d <- vci_data()
  res <- vci_with_seed(14, vci_warnings(suppressMessages(
    maihda(y ~ x + gender + race + edu + (1 | gender:race:edu), data = d,
           bootstrap = TRUE, n_boot = 40))))
  expect_true(any(grepl("n_boot = 40 is low", res$warnings, fixed = TRUE)))
  a <- res$value
  for (which in c("null", "adjusted")) {
    s <- if (which == "null") a$summary else a$summary_adjusted
    tv <- tidy(a, component = "variance", which = which)
    expect_identical(tv$conf.low, s$variance_components$lower, info = which)
    expect_true(is.finite(tv$conf.low[1]) && is.finite(tv$conf.high[1]), info = which)
  }
  tab <- maihda_table(a)$models
  bv <- tab[tab$statistic == "Between-stratum variance", ]
  bsd <- tab[tab$statistic == "Between-stratum SD", ]
  row <- vci_row(a$summary$variance_components, "Between-stratum (random)")
  expect_equal(c(bv$null_lower, bv$null_upper), c(row$lower, row$upper), tolerance = 1e-12)
  expect_equal(c(bsd$null_lower, bsd$null_upper), sqrt(c(row$lower, row$upper)),
               tolerance = 1e-12)
  row_a <- vci_row(a$summary_adjusted$variance_components, "Between-stratum (random)")
  expect_equal(c(bv$adjusted_lower, bv$adjusted_upper), c(row_a$lower, row_a$upper),
               tolerance = 1e-12)
  out <- utils::capture.output(print(maihda_table(a)))
  expect_true(any(grepl(sprintf("[%.3f, %.3f]", row$lower, row$upper), out, fixed = TRUE)))
})

test_that("the Shiny app attaches the worker's variance intervals to its summary", {
  d <- vci_data()
  m <- fit_maihda(y ~ x + (1 | gender:race:edu), data = d)
  ci <- vci_with_seed(15, vci_warnings(bootstrap_vpc(m$model, m$data, m$formula, 40, 0.95)))
  expect_identical(length(ci$warnings), 0L)
  s_app <- maihda_app_attach_vpc_ci(summary(m), ci$value)
  s_boot <- vci_boot_summary(m, 40, 15)$value
  expect_equal(s_app$variance_components, s_boot$variance_components, tolerance = 1e-12)
  expect_equal(c(s_app$vpc$ci_lower, s_app$vpc$ci_upper),
               c(s_boot$vpc$ci_lower, s_boot$vpc$ci_upper), tolerance = 1e-12)
})

# ---- brms, on posterior draws (no Stan) ---------------------------------------------

# Skewed SD draws built without random numbers; the shuffles keep the components from
# rising together.
vci_n <- 401L
vci_sd_s <- c(seq(0.30, 0.80, length.out = 400), 3)[order(sin(seq_len(401) * 12.9898))]
vci_sd_o <- c(seq(0.05, 0.40, length.out = 400), 2)[order(cos(seq_len(401) * 78.233))]
vci_sd_c <- c(seq(0.10, 0.60, length.out = 400), 2.5)[order(sin(seq_len(401) * 3.7))]
vci_sd_d <- c(seq(0.20, 0.90, length.out = 400), 4)[order(cos(seq_len(401) * 5.1))]
vci_sigma <- seq(0.8, 1.2, length.out = 401)[order(sin(seq_len(401) * 9.3))]
vci_cor <- seq(-0.6, 0.4, length.out = 401)[order(cos(seq_len(401) * 2.2))]

vci_stub <- function(draws, formula = y ~ 1 + (1 | stratum),
                     family = list(family = "gaussian", link = "identity"), ...) {
  structure(
    c(list(model = structure(list(formula = formula, family = family, draws = draws),
                             class = "brmsfit"),
           engine = "brms", formula = formula, family = family),
      list(...)),
    class = "maihda_model")
}

vci_mock_brms <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) model$draws,
    maihda_stratum_ranef_brms = function(model, group = "stratum", ...) {
      data.frame(stratum = c("1", "2"), stratum_id = 1:2, random_effect = c(-0.1, 0.1),
                 se = 0.1, lower_95 = c(-0.3, -0.1), upper_95 = c(0.1, 0.3))
    },
    maihda_discriminatory_accuracy = function(object, ...) NULL,
    .env = env)
  local_mocked_bindings(
    VarCorr = function(x, ...) {
      sd_cols <- grep("^sd_", names(x$draws), value = TRUE)
      out <- lapply(sd_cols, function(cn) {
        list(sd = matrix(c(stats::median(x$draws[[cn]]), 0, 0, 0), nrow = 1,
                         dimnames = list("Intercept", c("Estimate", "Est.Error",
                                                         "Q2.5", "Q97.5"))))
      })
      stats::setNames(out, sub("__Intercept$", "", sub("^sd_", "", sd_cols)))
    },
    fixef = function(object, ...) {
      matrix(c(0, 0.1, -0.2, 0.2), nrow = 1,
             dimnames = list("Intercept", c("Estimate", "Est.Error", "Q2.5", "Q97.5")))
    },
    .package = "brms", .env = env)
  local_mocked_s3_method("summary", "brmsfit", function(object, ...) NULL, frame = env)
}

test_that("a brms summary gives each estimated variance its credible interval", {
  skip_if_not_installed("brms")
  skip_if_not(exists("local_mocked_s3_method", envir = asNamespace("testthat")))
  vci_mock_brms()
  m <- vci_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = vci_sd_s,
                           sd_site__Intercept = vci_sd_o, sigma = vci_sigma),
                formula = y ~ 1 + (1 | stratum) + (1 | site))
  s <- summary(m)
  vc <- s$variance_components
  expect_equal(unlist(vci_row(vc, "Between-stratum (random)")[c("lower", "upper")]),
               c(lower = vci_q(vci_sd_s^2)[1], upper = vci_q(vci_sd_s^2)[2]), tolerance = 1e-12)
  expect_equal(unlist(vci_row(vc, "Other random effects")[c("lower", "upper")]),
               c(lower = vci_q(vci_sd_o^2)[1], upper = vci_q(vci_sd_o^2)[2]), tolerance = 1e-12)
  expect_equal(unlist(vci_row(vc, "Within-stratum (residual)")[c("lower", "upper")]),
               c(lower = vci_q(vci_sigma^2)[1], upper = vci_q(vci_sigma^2)[2]), tolerance = 1e-12)
  expect_true(all(is.na(unlist(vci_row(vc, "Total")[c("lower", "upper")]))))
  # The rows stay the posterior medians.
  expect_identical(vci_row(vc, "Between-stratum (random)")$variance, stats::median(vci_sd_s^2))
  out <- utils::capture.output(print(s))
  expect_identical(sum(grepl("(Posterior 95% credible intervals of the model's variance parameters)",
                             out, fixed = TRUE)), 1L)
  s90 <- summary(m, conf_level = 0.9)
  expect_equal(s90$variance_components$lower[1], vci_q(vci_sd_s^2, 0.9)[1], tolerance = 1e-12)

  # A logit model: the level-1 row has no interval.
  mb <- vci_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = vci_sd_s),
                 family = list(family = "bernoulli", link = "logit"))
  vcb <- summary(mb)$variance_components
  expect_equal(c(vcb$lower[1], vcb$upper[1]), vci_q(vci_sd_s^2), tolerance = 1e-12)
  expect_identical(vcb$variance[2], pi^2 / 3)
  expect_true(is.na(vcb$lower[2]) && is.na(vcb$upper[2]))
  expect_identical(tidy(summary(mb), component = "variance")$conf.high, vcb$upper)
})

test_that("contextual and crossed-dimensions brms summaries interval every component", {
  skip_if_not_installed("brms")
  vci_mock_brms()
  m <- vci_stub(data.frame(b_Intercept = 0, sd_stratum__Intercept = vci_sd_s,
                           sd_site__Intercept = vci_sd_c, sigma = vci_sigma),
                formula = y ~ 1 + (1 | stratum) + (1 | site))
  res <- maihda_context_summary_brms(m, list(context_vars = "site"), 0.95)
  vc <- res$variance_components
  expect_equal(unlist(vci_row(vc, "Context: site")[c("lower", "upper")]),
               c(lower = vci_q(vci_sd_c^2)[1], upper = vci_q(vci_sd_c^2)[2]), tolerance = 1e-12)
  expect_equal(unlist(vci_row(vc, "Between-stratum (random)")[c("lower", "upper")]),
               c(lower = vci_q(vci_sd_s^2)[1], upper = vci_q(vci_sd_s^2)[2]), tolerance = 1e-12)
  expect_equal(unlist(vci_row(vc, "Within-stratum (residual)")[c("lower", "upper")]),
               c(lower = vci_q(vci_sigma^2)[1], upper = vci_q(vci_sigma^2)[2]), tolerance = 1e-12)
  expect_true(is.na(vci_row(vc, "Total")$lower))

  mc <- vci_stub(
    data.frame(b_Intercept = 0, sd_ga__Intercept = vci_sd_d, sd_gb__Intercept = vci_sd_c,
               sd_stratum__Intercept = vci_sd_o, sd_site__Intercept = vci_sd_s),
    formula = y ~ 1 + (1 | ga) + (1 | gb) + (1 | stratum) + (1 | site),
    family = list(family = "bernoulli", link = "logit"),
    context_info = list(context_vars = "site"))
  cc <- list(dim_groups = c(a = "ga", b = "gb"), interaction_group = "stratum")
  rc <- maihda_cc_summary_brms(mc, cc, conf_level = 0.95)
  vcc <- rc$variance_components
  want <- list(`Additive: a` = vci_sd_d^2, `Additive: b` = vci_sd_c^2,
               `Intersectional interaction` = vci_sd_o^2, `Context: site` = vci_sd_s^2)
  for (lab in names(want)) {
    expect_equal(unlist(vci_row(vcc, lab)[c("lower", "upper")]),
                 c(lower = vci_q(want[[lab]])[1], upper = vci_q(want[[lab]])[2]),
                 tolerance = 1e-12, info = lab)
  }
  expect_true(is.na(vci_row(vcc, "Within-stratum (residual)")$lower))
  between <- vci_sd_d^2 + vci_sd_c^2 + vci_sd_o^2
  expect_equal(rc$decomposition$between_var_ci, vci_q(between), tolerance = 1e-12)
  # Not the sum of the rows' bounds: the components do not rise together.
  expect_gt(abs(rc$decomposition$between_var_ci[2] -
                  sum(vcc$upper[vcc$component %in% c("Additive: a", "Additive: b",
                                                     "Intersectional interaction")])), 0.01)
})

test_that("a growth brms summary intervals each block entry and the reference-time variance", {
  skip_if_not_installed("brms")
  vci_mock_brms()
  sd_s0 <- vci_sd_s; sd_s1 <- vci_sd_o / 4; sd_i0 <- vci_sd_c; sd_i1 <- vci_sd_d / 5
  draws <- data.frame(b_Intercept = 0,
                      sd_stratum__Intercept = sd_s0, sd_stratum__wave = sd_s1,
                      cor_stratum__Intercept__wave = vci_cor,
                      sd_id__Intercept = sd_i0, sd_id__wave = sd_i1,
                      cor_id__Intercept__wave = rev(vci_cor), sigma = vci_sigma)
  m <- vci_stub(draws, formula = y ~ wave + (1 + wave | stratum) + (1 + wave | id),
                longitudinal_info = list(time = "wave", id = "id", time_degree = 1L,
                                         ref_time = 0, stratum_slope = TRUE),
                data = data.frame(wave = 0:3))
  res <- maihda_longitudinal_summary_brms(m, 0.95)
  vc <- res$variance_components
  want <- list(sd_s0^2, sd_s1^2, vci_cor * sd_s0 * sd_s1, sd_i0^2, sd_i1^2,
               rev(vci_cor) * sd_i0 * sd_i1, vci_sigma^2)
  expect_identical(nrow(vc), length(want))
  for (j in seq_along(want)) {
    expect_equal(c(vc$lower[j], vc$upper[j]), vci_q(want[[j]]), tolerance = 1e-12,
                 info = vc$component[j])
  }
  expect_match(vc$component[3], "covariance")
  # At the reference time 0 the between-stratum variance is the intercept variance.
  expect_equal(res$longitudinal$var_stratum_ref_ci, vci_q(sd_s0^2), tolerance = 1e-12)
  # The rows stay the posterior means.
  expect_equal(vc$variance[1], mean(sd_s0^2), tolerance = 1e-12)

  # No stratum slope: a one-entry stratum block.
  m0 <- vci_stub(draws[, setdiff(names(draws), c("sd_stratum__wave",
                                                 "cor_stratum__Intercept__wave"))],
                 formula = y ~ wave + (1 | stratum) + (1 + wave | id),
                 longitudinal_info = list(time = "wave", id = "id", time_degree = 1L,
                                          ref_time = 0, stratum_slope = FALSE),
                 data = data.frame(wave = 0:3))
  vc0 <- maihda_longitudinal_summary_brms(m0, 0.95)$variance_components
  expect_identical(nrow(vc0), 5L)
  expect_equal(c(vc0$lower[1], vc0$upper[1]), vci_q(sd_s0^2), tolerance = 1e-12)
  expect_equal(c(vc0$lower[5], vc0$upper[5]), vci_q(vci_sigma^2), tolerance = 1e-12)
})

# ---- brms, a real fit (opt-in) -------------------------------------------------------

test_that("a brms fit's intervals are the quantiles of its own variance draws", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")
  d <- vci_data()
  # The sampler also draws from R's random stream, which this file leaves as it found it.
  m <- vci_with_seed(74, suppressWarnings(suppressMessages(
    fit_maihda(y ~ x + (1 | gender:race:edu), data = d, engine = "brms",
               chains = 2, iter = 1500, warmup = 750, refresh = 0, seed = 74,
               control = list(adapt_delta = 0.95)))))
  s <- summary(m)
  vc <- s$variance_components
  dr <- as.data.frame(brms::as_draws_df(m$model))
  v_s <- dr[["sd_stratum__Intercept"]]^2
  v_e <- dr[["sigma"]]^2
  expect_equal(c(vc$lower[1], vc$upper[1]), vci_q(v_s), tolerance = 1e-10)
  expect_equal(c(vc$lower[2], vc$upper[2]), vci_q(v_e), tolerance = 1e-10)
  expect_true(is.na(vc$lower[3]))
  tab <- maihda_table(m)$models
  bv <- tab[tab$statistic == "Between-stratum variance", ]
  expect_equal(c(bv$estimate_lower, bv$estimate_upper), vci_q(v_s), tolerance = 1e-10)
})
