# Audit 2026-09-10 -- the longitudinal trajectory plot drew every stratum on ONE
# fixed trajectory, suppressing the additive between-stratum differences.
#
# plot_stratum_trajectories() computed a single fixed-part curve at the mean/modal
# covariate profile and added each stratum's random intercept & slope to it. The
# profile loop that built that curve replaced EVERY non-time fixed variable --
# including the stratum-defining dimensions. In an adjusted growth model the dimension
# main effects and their dim:time interactions are part of the plotted stratum's own
# prediction, so freezing them at the modal profile drew all strata on the modal
# stratum's curve.
#
# On maihda_long_data (12 strata, adjusted longitudinal Gaussian) the plotted baselines
# spanned 4.510-4.939 (SD 0.137) where the correct predictions span 4.059-5.971
# (SD 0.637), a max discrepancy of 2.32 over the time grid. The mechanism was exact:
# plotted_j(t) - correct_j(t) = (X_modal(t) - X_j(t))' beta to 8.9e-16, and the one
# line with zero error was the stratum that IS the modal profile.
#
# maihda_longitudinal_fixed_trajectory() now takes `strata` and returns one column per
# stratum, seeding each block of grid rows from that stratum's own fitted rows so the
# dimension columns follow the stratum while every OTHER covariate stays at the shared
# reference profile. `strata = NULL` keeps the old single population trajectory, so a
# null growth model -- whose fixed part carries no dimension terms -- is untouched.

skip_if_not_installed("lme4")

fit_q <- function(...) suppressMessages(suppressWarnings(fit_maihda(...)))

# Fixed part at each stratum's own dimensions + that stratum's random effects, WITHOUT
# the person-level effects: the estimand the view reports, from lme4's own predict().
native_traj <- function(model, strata, grid, re_form, tweak = identity) {
  d <- model$data
  tv <- model$longitudinal_info$time
  unlist(lapply(as.character(strata), function(j) {
    rj <- which(as.character(d$stratum) == j)[1]
    nd <- d[rep(rj, length(grid)), , drop = FALSE]
    nd[[tv]] <- grid
    as.numeric(stats::predict(model$model, newdata = tweak(nd), re.form = re_form,
                              allow.new.levels = TRUE))
  }))
}

# What the plot produced BEFORE the fix: one population trajectory + each stratum's
# random deviation. Used as the confinement reference for the paths that must not move.
prefix_values <- function(m, s) {
  grid <- s$longitudinal$time_grid
  gcen <- grid - MAIHDA:::maihda_lng_time_center(m$longitudinal_info)
  e0 <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid)
  re <- MAIHDA:::maihda_longitudinal_stratum_re(m)
  unlist(lapply(seq_len(nrow(re)), function(i) {
    co <- re$coef[[i]]
    e0 + vapply(gcen, function(t) sum(co * t^(0:(length(co) - 1))), numeric(1))
  }))
}

# ---- the headline defect ----------------------------------------------------

test_that("adjusted longitudinal trajectories carry each stratum's own dimensions", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal")))
  m <- a$model_adjusted
  # The adjusted growth formula really does put the dimensions AND their time
  # interactions in the fixed part -- without that there is nothing to suppress.
  lbl <- attr(stats::terms(MAIHDA:::maihda_nobars(m$formula)), "term.labels")
  expect_true(all(c("gender", "ethnicity", "education") %in% lbl))
  expect_true(any(grepl("^wave:", lbl)))

  # The invariant that makes seeding a grid row per stratum safe: the dimension
  # columns really are constant within a stratum, so WHICH fitted row is picked
  # cannot matter. Stated as a comment in the helper; asserted here.
  fv <- all.vars(MAIHDA:::maihda_nobars(m$formula)[[3]])
  dcols <- MAIHDA:::maihda_longitudinal_dimension_columns(m, fv, names(m$data))
  expect_setequal(dcols, c("gender", "ethnicity", "education"))
  for (v in dcols) {
    expect_equal(max(tapply(as.character(m$data[[v]]),
                            as.character(m$data$stratum),
                            function(x) length(unique(x)))), 1)
  }

  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  strata <- unique(as.character(pd$stratum))
  expect_equal(pd$value, native_traj(m, strata, grid, ~ (wave | stratum)),
               tolerance = 1e-10)

  # The symptom itself: the plotted baselines must span the real between-stratum
  # range, not the collapsed one the modal-profile curve produced (SD 0.137).
  b <- pd$value[pd$time == grid[1]]
  expect_gt(stats::sd(b), 0.5)
  expect_gt(diff(range(b)), 1.5)
})

test_that("a null growth model's trajectories are unchanged by the per-stratum split", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  m <- fit_q(wellbeing ~ wave + (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave")
  s <- suppressMessages(summary(m))
  # No dimension reaches the fixed part, so the per-stratum grid must reproduce the
  # single population trajectory + random deviation EXACTLY.
  expect_equal(MAIHDA:::maihda_longitudinal_dimension_columns(
    m, all.vars(MAIHDA:::maihda_nobars(m$formula)[[3]]), names(m$data)),
    character(0))
  expect_equal(plot_stratum_trajectories(m, s)$data$value, prefix_values(m, s),
               tolerance = 0)
})

test_that("a covariate-only adjustment keeps its reference profile", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  # `age` is a covariate, not a stratum dimension: it stays at the shared mean, so
  # this fit is also bit-identical to the pre-fix construction.
  m <- fit_q(wellbeing ~ wave + age + (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave")
  s <- suppressMessages(summary(m))
  expect_equal(plot_stratum_trajectories(m, s)$data$value, prefix_values(m, s),
               tolerance = 0)
})

test_that("strata = NULL still returns the population reference trajectory", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal")))
  m <- a$model_adjusted
  grid <- suppressMessages(summary(m))$longitudinal$time_grid
  eN <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid)
  expect_true(is.null(dim(eN)))
  expect_length(eN, length(grid))
  # ... at the MODAL dimension profile, exactly as before: on this branch the
  # dimensions are treated as ordinary covariates.
  mf <- m$data
  nd <- mf[rep(1L, length(grid)), , drop = FALSE]
  for (v in c("gender", "ethnicity", "education")) {
    nd[[v]] <- names(sort(table(mf[[v]]), decreasing = TRUE))[1]
  }
  nd$wave <- grid
  expect_equal(eN, as.numeric(stats::predict(m$model, newdata = nd, re.form = NA,
                                             allow.new.levels = TRUE)),
               tolerance = 1e-10)
})

# ---- pairing, shapes and guards ---------------------------------------------

test_that("each returned column belongs to the stratum that asked for it", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal")))
  m <- a$model_adjusted
  grid <- suppressMessages(summary(m))$longitudinal$time_grid
  st <- as.character(MAIHDA:::maihda_longitudinal_stratum_re(m)$stratum)
  E <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid, strata = st)
  expect_equal(dim(E), c(length(grid), length(st)))
  expect_identical(colnames(E), st)
  # Columns genuinely differ -- otherwise the permutation check below has no teeth.
  expect_gt(max(apply(E, 1, stats::sd)), 0.1)
  Er <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid, strata = rev(st))
  expect_equal(unname(Er), unname(E[, rev(seq_along(st)), drop = FALSE]))
  # A single stratum returns that stratum's own column, not the population curve.
  E1 <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid, strata = st[3])
  expect_equal(dim(E1), c(length(grid), 1L))
  expect_equal(as.numeric(E1), unname(E[, 3]))
})

test_that("the n_strata cap keeps the stratum/column pairing", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
           data = maihda_long_data, id = "id", time = "wave",
           decomposition = "longitudinal")))
  m <- a$model_adjusted
  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  for (sel in c("order", "deviation")) {
    pd <- plot(m, type = "trajectories", n_strata = 4L, select = sel)$data
    kept <- unique(as.character(pd$stratum))
    expect_length(kept, 4L)
    expect_equal(pd$value, native_traj(m, kept, grid, ~ (wave | stratum)),
                 tolerance = 1e-10)
  }
})

test_that("a stratum with no fitted row errors instead of falling back silently", {
  skip_on_cran()
  d <- data.frame(stratum = rep(c("a", "b"), each = 3), x = 1:6)
  expect_equal(MAIHDA:::maihda_stratum_seed_rows(d, c("b", "a")), c(4L, 1L))
  expect_error(MAIHDA:::maihda_stratum_seed_rows(d, c("a", "zz")), "No fitted rows")
  expect_error(MAIHDA:::maihda_stratum_seed_rows(d["x"], "a"), "No 'stratum' column")
})

test_that("a mis-sized fixed-part return is rejected, not recycled into the grid", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  m <- fit_q(wellbeing ~ wave + gender + ethnicity + education +
               (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave")
  grid <- suppressMessages(summary(m))$longitudinal$time_grid
  st <- as.character(MAIHDA:::maihda_longitudinal_stratum_re(m)$stratum)[1:3]
  # matrix() recycles silently, which would mis-pair every stratum with a
  # plausible-looking curve; the length check must turn that into an error.
  testthat::local_mocked_bindings(
    maihda_lme4_fixed_link = function(model, newdata, offset = NULL, ...) {
      rep(0, nrow(newdata) - 1L)
    })
  expect_error(
    MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid, strata = st),
    "cannot align them with the strata")
})

# ---- the dimension-column rule ----------------------------------------------

test_that("maihda_adjusted_term_names matches the terms maihda_adjusted_terms builds", {
  d <- data.frame(g = c("a", "b", "a", "b"), z = c(1, 2, 3, 4))
  ab <- list(z = list(breaks = c(0, 2, 4), labels = c("lo", "hi")))
  expect_identical(MAIHDA:::maihda_adjusted_term_names(c("g", "z"), ab),
                   MAIHDA:::maihda_adjusted_terms(c("g", "z"), ab, d)$terms)
  expect_identical(MAIHDA:::maihda_adjusted_term_names(c("g", "z"), NULL),
                   MAIHDA:::maihda_adjusted_terms(c("g", "z"), NULL, d)$terms)
  expect_identical(MAIHDA:::maihda_adjusted_term_names(character(0), NULL),
                   character(0))
})

test_that("an auto-binned numeric dimension follows the stratum", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  set.seed(4)
  d <- maihda_long_data
  d$income <- rep(stats::rnorm(length(unique(d$id)), 50, 12), each = 5)
  a <- suppressMessages(suppressWarnings(
    maihda(wellbeing ~ wave + (1 | gender:income), data = d, id = "id",
           time = "wave", decomposition = "longitudinal")))
  m <- a$model_adjusted
  # The adjusted model enters the binned tertile factor, not raw income.
  expect_true(".maihda_dim_income" %in%
                all.vars(MAIHDA:::maihda_nobars(m$formula)[[3]]))
  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  kept <- unique(as.character(pd$stratum))
  expect_equal(pd$value, native_traj(m, kept, grid, ~ (wave | stratum)),
               tolerance = 1e-10)
  expect_gt(stats::sd(pd$value[pd$time == grid[1]]), 0.01)
})

test_that("a raw numeric dimension is held at ITS stratum's mean", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  set.seed(4)
  d <- maihda_long_data
  d$income <- rep(stats::rnorm(length(unique(d$id)), 50, 12), each = 5)
  # A hand-written adjusted formula may enter an auto-binned dimension as a linear
  # term. It is not constant within a stratum, so a single seed row would be
  # arbitrary -- the stratum's own mean is the representative value.
  m <- fit_q(wellbeing ~ wave + gender + income + wave:income + (1 | gender:income),
             data = d, id = "id", time = "wave")
  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  kept <- unique(as.character(pd$stratum))
  mj <- tapply(m$data$income, as.character(m$data$stratum), mean)
  own <- native_traj(m, kept, grid, ~ (wave | stratum), tweak = function(nd) {
    nd$income <- mj[[as.character(nd$stratum[1])]]
    nd
  })
  glob <- native_traj(m, kept, grid, ~ (wave | stratum), tweak = function(nd) {
    nd$income <- mean(m$data$income)
    nd
  })
  expect_equal(pd$value, own, tolerance = 1e-10)
  # Teeth: the global mean is a materially different answer.
  expect_gt(max(abs(own - glob)), 0.1)
})

test_that("a dimension reachable only through a transformed term warns", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  # factor(gender) stores no raw `gender` column in the model frame, so the grid
  # cannot set it per stratum and it stays at a representative value. Pre-existing
  # behaviour of the model-frame-based grid -- but it must not be silent.
  m <- fit_q(wellbeing ~ wave + factor(gender) + ethnicity + education +
               (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave")
  s <- suppressMessages(summary(m))
  expect_warning(pd <- plot_stratum_trajectories(m, s)$data, "representative value")
  pd <- suppressWarnings(plot_stratum_trajectories(m, s)$data)
  # The other two dimensions DO follow their stratum, so the fix still applies.
  expect_gt(stats::sd(pd$value[pd$time == 0]), 0.3)
  # ... and gender is genuinely frozen: strata differing only in gender then differ
  # by their random effects alone.
  pr <- m$strata_info[, c("stratum", "gender", "ethnicity", "education")]
  pair <- merge(pr, pr, by = c("ethnicity", "education"))
  pair <- pair[as.character(pair$stratum.x) < as.character(pair$stratum.y), ]
  expect_gt(nrow(pair), 0)
  re <- MAIHDA:::maihda_longitudinal_stratum_re(m)
  at0 <- function(j) pd$value[as.character(pd$stratum) == j & pd$time == 0]
  u0 <- function(j) re$coef[[which(as.character(re$stratum) == j)]][1]
  gap <- vapply(seq_len(nrow(pair)), function(k) {
    x <- as.character(pair$stratum.x[k]); y <- as.character(pair$stratum.y[k])
    (at0(x) - at0(y)) - (u0(x) - u0(y))
  }, numeric(1))
  expect_equal(max(abs(gap)), 0, tolerance = 1e-10)
})

# ---- shapes the per-stratum grid must not break -----------------------------

test_that("a formula offset is still re-evaluated on the per-stratum grid", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  m <- fit_q(wellbeing ~ wave + gender + ethnicity + education + offset(0.5 * wave) +
               (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave")
  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  kept <- unique(as.character(pd$stratum))
  # predict.merMod re-evaluates offset(0.5 * wave) itself here (`wave` is in the
  # frame), so the native prediction already carries the time-tracking offset.
  expect_equal(pd$value, native_traj(m, kept, grid, ~ (wave | stratum)),
               tolerance = 1e-10)
  # It tracks time rather than being frozen at mean(0.5 * wave).
  E <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid, strata = kept)
  noff <- unlist(lapply(kept, function(j) {
    rj <- which(as.character(m$data$stratum) == j)[1]
    nd <- m$data[rep(rj, length(grid)), , drop = FALSE]
    nd$wave <- grid
    MAIHDA:::maihda_lme4_fixed_link(m$model, nd, offset = NULL)
  }))
  expect_equal(as.numeric(E), noff + rep(0.5 * grid, length(kept)), tolerance = 1e-10)
  expect_gt(max(abs(as.numeric(E) -
                      (noff + mean(0.5 * maihda_long_data$wave)))), 0.5)
})

test_that("an external offset= is still added as its sample mean, per stratum", {
  skip_on_cran()
  set.seed(9)
  nid <- 160
  base <- do.call(rbind, lapply(seq_len(nid), function(i) {
    data.frame(id = i, g1 = sample(c("F", "M"), 1), g2 = sample(c("A", "B", "C"), 1),
               time = 0:2, logE = log(stats::runif(3, 1, 6)),
               stringsAsFactors = FALSE)
  }))
  base$y <- stats::rpois(nrow(base), exp(-0.2 + 0.15 * base$time + base$logE))
  m <- fit_q(y ~ 1 + g1 + g2 + (1 | g1:g2), data = base, family = "poisson",
             id = "id", time = "time", offset = logE)
  grid <- suppressMessages(summary(m))$longitudinal$time_grid
  st <- as.character(MAIHDA:::maihda_longitudinal_stratum_re(m)$stratum)
  E <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid, strata = st)
  noff <- unlist(lapply(st, function(j) {
    rj <- which(as.character(m$data$stratum) == j)[1]
    nd <- m$data[rep(rj, length(grid)), , drop = FALSE]
    nd$time <- grid
    MAIHDA:::maihda_lme4_fixed_link(m$model, nd, offset = NULL)
  }))
  expect_equal(as.numeric(E), noff + mean(base$logE), tolerance = 1e-10)
  expect_gt(stats::sd(E[1, ]), 0.01)
})

test_that("centered time and a quadratic growth curve stay aligned", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  d <- maihda_long_data
  d$year <- d$wave + 2000
  m <- fit_q(wellbeing ~ year + gender + ethnicity + education +
               (1 | gender:ethnicity:education),
             data = d, id = "id", time = "year", time_degree = 2)
  # Internal centering is active, so the grid rows must carry BOTH the raw grid time
  # and the derived centered column.
  expect_equal(MAIHDA:::maihda_lng_time_center(m$longitudinal_info), 2000)
  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  kept <- unique(as.character(pd$stratum))
  ref <- unlist(lapply(kept, function(j) {
    rj <- which(as.character(m$data$stratum) == j)[1]
    nd <- m$data[rep(rj, length(grid)), , drop = FALSE]
    nd$.maihda_ctime <- grid - 2000
    as.numeric(stats::predict(
      m$model, newdata = nd,
      re.form = ~ (.maihda_ctime + I(.maihda_ctime^2) | stratum),
      allow.new.levels = TRUE))
  }))
  expect_equal(pd$value, ref, tolerance = 1e-9)
})

test_that("the brms engine gets the same per-stratum fixed part", {
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # The nd construction is shared above the engine branch, so brms needs the same
  # per-stratum split -- and the comparison below is an algebraic identity on
  # WHATEVER posterior is drawn, not a distributional claim: posterior_linpred is
  # linear, so its posterior mean equals (posterior-mean fixed part) + (posterior-mean
  # stratum RE), which is exactly what the plot composes. Convergence is therefore
  # irrelevant here and this fixture stays deliberately cheap, unlike the mixing-
  # sensitive brms longitudinal block in test-longitudinal.R.
  set.seed(3)
  nid <- 120
  g <- sample(c("F", "M"), nid, TRUE)
  e <- sample(c("A", "B", "C"), nid, TRUE)
  d <- do.call(rbind, lapply(seq_len(nid), function(i) {
    data.frame(id = i, g = g[i], e = e[i], wave = 0:2, stringsAsFactors = FALSE)
  }))
  b_g <- c(F = 0, M = -0.9)[d$g]
  b_e <- c(A = 0, B = 0.8, C = -1.2)[d$e]
  s_g <- c(F = 0, M = 0.35)[d$g]
  s_e <- c(A = 0, B = -0.30, C = 0.45)[d$e]
  u <- stats::rnorm(nid, 0, 0.5)[d$id]
  d$y <- 5 + b_g + b_e + (0.2 + s_g + s_e) * d$wave + u +
    stats::rnorm(nrow(d), 0, 0.6)

  m <- suppressWarnings(suppressMessages(
    fit_maihda(y ~ wave + g + e + wave:g + wave:e + (1 | g:e), data = d,
               id = "id", time = "wave", engine = "brms",
               chains = 1, iter = 600, warmup = 300, refresh = 0, seed = 7)))
  # suppressWarnings, not a tighter adapt_delta: this fixture is deliberately short
  # and brms' sampler diagnostics (e.g. a divergent transition) are irrelevant to an
  # identity that holds on whatever posterior is drawn -- but they would otherwise
  # surface through summary() as a test warning.
  s <- suppressWarnings(suppressMessages(summary(m)))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  kept <- unique(as.character(pd$stratum))

  dat <- m$data
  nat <- unlist(lapply(kept, function(j) {
    rj <- which(as.character(dat$stratum) == j)[1]
    nd <- dat[rep(rj, length(grid)), , drop = FALSE]
    nd$wave <- grid
    colMeans(brms::posterior_linpred(m$model, newdata = nd,
                                     re_formula = ~ (wave | stratum),
                                     allow_new_levels = TRUE))
  }))
  expect_equal(pd$value, nat, tolerance = 1e-10)
  # The pre-fix construction -- one population curve + each stratum's RE -- collapses
  # the spread on brms harder than on lme4 (SD 0.04 against 1.04 here).
  e0 <- MAIHDA:::maihda_longitudinal_fixed_trajectory(m, grid)
  re <- MAIHDA:::maihda_longitudinal_stratum_re(m)
  old <- unlist(lapply(seq_len(nrow(re)), function(i) {
    co <- re$coef[[i]]
    e0 + vapply(grid, function(t) sum(co * t^(0:(length(co) - 1))), numeric(1))
  }))
  expect_gt(max(abs(old - nat)), 0.5)
  expect_gt(stats::sd(pd$value[pd$time == grid[1]]), 5 *
              stats::sd(old[seq(1, length(old), by = length(grid))]))
})

test_that("stratum_slope = FALSE still matches the native prediction", {
  skip_on_cran()
  data(maihda_long_data, package = "MAIHDA")
  m <- fit_q(wellbeing ~ wave + gender + ethnicity + education +
               (1 | gender:ethnicity:education),
             data = maihda_long_data, id = "id", time = "wave",
             stratum_slope = FALSE)
  s <- suppressMessages(summary(m))
  grid <- s$longitudinal$time_grid
  pd <- plot_stratum_trajectories(m, s)$data
  kept <- unique(as.character(pd$stratum))
  expect_equal(pd$value, native_traj(m, kept, grid, ~ (1 | stratum)),
               tolerance = 1e-10)
})
