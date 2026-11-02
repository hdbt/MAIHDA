# Audit pass 2026-11-02: bootstrap refits of an lme4 fit with a variance parameter just above
# zero.
#
# lme4::refit() starts the optimiser at the fit's own variance parameters. Under nloptwrap
# (lmer's default) and bobyqa (glmer with nAGQ = 0) the first step for one just above its bound
# of zero is a fraction of its distance from the bound, on a criterion that is flat there, so
# every refit stayed at the fit's value: summary(bootstrap = TRUE) of a fit at theta 8.3e-9
# printed the VPC interval [7e-17, 2e-16] where fresh fits of the same draws reach 0.036, and
# beside larger variance parameters one at 1e-4 to 1e-3 stayed where it was in a third of the
# draws. maihda_refit_draw() now starts, under those two optimisers, a random-effect standard
# deviation below 0.01 of the residual one at zero (maihda_refit_start()); the other parameters
# start at the fit, and fits made with other optimisers are refitted as before (nlminb, started
# on the bound, stays there; glmer's Nelder-Mead refits are not held near it).
#
# The fits below are lme4's own, started just above zero on data whose strata means are equal, so
# they stop there on every platform. Random numbers are drawn only inside with_seed(), which puts
# the suite's stream back.

with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv()) else
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  force(expr)
}

quiet <- function(e) suppressWarnings(suppressMessages(e))
criterion_of <- function(m) {
  as.numeric(m@devcomp$cmp[[if (lme4::isLMM(m) && lme4::isREML(m)) "REML" else "dev"]])
}

# Twelve strata of 25 rows (a 2 x 3 x 2 design). Every stratum holds the same 25 residuals in
# another order, so the strata means are equal and the REML stratum variance is zero; `shift`
# adds an effect per stratum (in the order of the 12 cells) and `x` a covariate on 1, 2, 3.
flat_strata <- function(shift = rep(0, 12)) {
  g <- expand.grid(a = c("a1", "a2"), b = c("b1", "b2", "b3"), c = c("c1", "c2"),
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  d <- g[rep(seq_len(nrow(g)), each = 25), ]
  e <- stats::qnorm(stats::ppoints(25))
  cell <- rep(seq_len(nrow(g)), each = 25)
  d$y <- unlist(lapply(seq_len(nrow(g)), function(j) e[(seq_along(e) + 7 * j) %% 25 + 1])) +
    shift[cell]
  d$x <- rep(c(1, 2, 3), length.out = nrow(d))
  rownames(d) <- NULL
  d
}
# Twelve stratum effects with a standard deviation of about 0.5.
stratum_shift <- c(-0.8, 0.6, -0.3, 0.9, 0.1, -0.5, 0.4, -0.9, 0.7, -0.1, 0.3, -0.4)

# The lme4 fit `m` as lme4 would have left it at the variance parameters `theta`: in its
# predictor, where lme4::refit() starts, in its theta slot and in its optimiser's record
# (lme4::getME(, "par")).
fit_at <- function(m, theta) {
  theta <- unname(theta)
  pp <- m@pp$copy()
  pp$setTheta(theta)
  m@pp <- pp
  m@theta <- theta
  m@optinfo$val[seq_along(theta)] <- theta
  m
}

# A fit_maihda() fit whose stratum variance parameter lme4 left just above zero.
near_zero_fit <- function(...) {
  quiet(fit_maihda(y ~ 1 + (1 | a:b:c), data = flat_strata(), start = list(theta = 1e-8), ...))
}
# A fresh fit of `model`'s formula on its own frame with the response `y`.
fresh_fit <- function(model, y, ...) {
  mf <- stats::model.frame(model)
  mf[[1]] <- y
  if (lme4::isLMM(model)) {
    quiet(lme4::lmer(stats::formula(model), data = mf, REML = lme4::isREML(model)))
  } else {
    quiet(lme4::glmer(stats::formula(model), data = mf, family = stats::family(model), ...))
  }
}

test_that("a fit whose variance parameter is just above zero is refitted to a fresh fit's optimum", {
  fit <- near_zero_fit()
  m <- fit$model
  theta <- lme4::getME(m, "theta")
  # The premise: lme4 left the fit just above zero.
  expect_true(theta > 0 && theta <= 1e-8)
  pp_theta <- m@pp$theta

  y_new <- flat_strata(stratum_shift)$y
  r <- quiet(maihda_refit_draw(m, y_new))
  ref <- fresh_fit(m, y_new)
  expect_gt(lme4::getME(ref, "theta"), 0.3)
  expect_lt(abs(criterion_of(r) - criterion_of(ref)), 1e-6)
  expect_lt(abs(lme4::getME(r, "theta") / lme4::getME(ref, "theta") - 1), 1e-4)

  # The fit itself is not touched.
  expect_identical(lme4::getME(m, "theta"), theta)
  expect_identical(m@pp$theta, pp_theta)
})

test_that("summary(bootstrap = TRUE) of such a fit varies the stratum variance as fresh fits do", {
  skip_slow_test()
  fit <- near_zero_fit()
  s <- with_seed(3, quiet(summary(fit, bootstrap = TRUE, n_boot = 40)))

  # The same draws, each fitted afresh.
  sims <- with_seed(3, quiet(maihda_simulate_lme4(fit$model, nsim = 40)))
  ref <- vapply(sims, function(y) {
    f <- fresh_fit(fit$model, y)
    vc <- as.numeric(lme4::VarCorr(f)$stratum)
    c(vpc = vc / (vc + stats::sigma(f)^2), var = vc)
  }, numeric(2))
  ref_vpc <- stats::quantile(ref["vpc", ], c(0.025, 0.975), names = FALSE)
  ref_var <- stats::quantile(ref["var", ], c(0.025, 0.975), names = FALSE)
  expect_gt(ref_vpc[2], 1e-3)

  expect_lt(abs(s$vpc$ci_lower - ref_vpc[1]), 1e-6)
  expect_lt(abs(s$vpc$ci_upper - ref_vpc[2]), 1e-5)
  row <- s$variance_components[s$variance_components$component == "Between-stratum (random)", ]
  expect_lt(abs(row$lower - ref_var[1]), 1e-6)
  expect_lt(abs(row$upper - ref_var[2]), 1e-5)
})

test_that("the PCV bootstrap refits an adjusted fit just above zero as fresh fits do", {
  skip_slow_test()
  d <- flat_strata(stratum_shift)
  # Additive strata: the adjusted model's main effects take the stratum effects whole.
  d$y <- d$y - stratum_shift[rep(seq_len(12), each = 25)] +
    0.6 * (d$a == "a2") + c(0, 0.5, -0.4)[match(d$b, c("b1", "b2", "b3"))] - 0.3 * (d$c == "c2")
  null <- quiet(fit_maihda(y ~ 1 + (1 | a:b:c), data = d))
  adj <- quiet(fit_maihda(y ~ a + b + c + (1 | a:b:c), data = d, start = list(theta = 1e-8)))
  theta <- lme4::getME(adj$model, "theta")
  expect_true(theta > 0 && theta <= 1e-8)

  p <- with_seed(5, quiet(calculate_pcv(null, adj, bootstrap = TRUE, n_boot = 30)))
  # The adjusted fit's draws carry no stratum variance, but fresh fits of them give a positive
  # one about half of the time, so the PCV interval reaches well below 1.
  sims <- with_seed(5, quiet(maihda_simulate_lme4(adj$model, nsim = 30)))
  ref <- vapply(sims, function(y) {
    v1 <- as.numeric(lme4::VarCorr(fresh_fit(null$model, y))$stratum)
    v2 <- as.numeric(lme4::VarCorr(fresh_fit(adj$model, y))$stratum)
    (v1 - v2) / v1
  }, numeric(1))
  ref_ci <- stats::quantile(ref, c(0.025, 0.975), names = FALSE)
  expect_lt(ref_ci[1], 0.99)
  expect_lt(abs(p$ci_lower - ref_ci[1]), 1e-4)
  expect_lt(abs(p$ci_upper - ref_ci[2]), 1e-4)
})

test_that("a variance parameter between 1e-4 and 0.01 beside larger ones is refitted from zero", {
  d <- flat_strata(stratum_shift)
  m <- quiet(lme4::lmer(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b:c), data = d))
  theta <- lme4::getME(m, "theta")
  # A fit that left the dimension `a` at 1e-3 beside the others.
  theta[["a.(Intercept)"]] <- 1e-3
  m <- fit_at(m, theta)

  # A response in which `a` has a clear effect.
  y_new <- d$y + ifelse(d$a == "a2", 1.2, -1.2)
  r <- quiet(maihda_refit_draw(m, y_new))
  ref <- fresh_fit(m, y_new)
  expect_gt(lme4::getME(ref, "theta")[["a.(Intercept)"]], 0.3)
  expect_lt(abs(criterion_of(r) - criterion_of(ref)), 1e-4)
  expect_lt(abs(lme4::getME(r, "theta")[["a.(Intercept)"]] /
                  lme4::getME(ref, "theta")[["a.(Intercept)"]] - 1), 1e-2)
})

test_that("a glmer fit with nAGQ = 0 just above zero is refitted to a fresh fit's optimum", {
  d <- flat_strata()
  # Equal event rates in every stratum: the same 0 / 1 pattern in each.
  d$z <- rep(rep(c(1, 0, 0, 1, 0), 5), 12)
  m <- quiet(lme4::glmer(z ~ 1 + (1 | a:b:c), data = d, family = stats::binomial(),
                         nAGQ = 0, start = list(theta = 1e-8)))
  theta <- lme4::getME(m, "theta")
  expect_true(theta > 0 && theta <= 1e-8)

  # A response with clearly different event rates by stratum.
  rate <- rep(c(0.1, 0.8, 0.3, 0.9, 0.2, 0.6, 0.15, 0.7, 0.4, 0.85, 0.25, 0.5), each = 25)
  z_new <- as.numeric(rep(seq_len(25) / 25, 12) <= rate)
  r <- quiet(maihda_refit_draw(m, z_new))
  ref <- fresh_fit(m, z_new, nAGQ = 0)
  expect_gt(lme4::getME(ref, "theta"), 0.3)
  expect_lt(abs(criterion_of(r) - criterion_of(ref)), 1e-3)
})

test_that("only a variance parameter strictly between 0 and 0.01 starts at zero", {
  d <- flat_strata(stratum_shift)
  m <- quiet(lme4::lmer(y ~ x + (1 + x | a:b:c) + (1 | a), data = d))
  lower <- lme4::getME(m, "lower")
  diag_i <- which(lower == 0)
  off_i <- which(lower != 0)
  expect_length(off_i, 1L)
  base <- lme4::getME(m, "theta")
  # lme4 stops the slope's diagonal entry and the `a` term where the REML criterion is flat at
  # zero: on it on x86_64, 3.5e-5 above it on the arm64 macOS runner. Both are put on it, so that
  # the entries placed below are the only ones strictly between 0 and 0.01.
  base[diag_i[-1]] <- 0
  place <- function(theta) fit_at(m, theta)

  # Below the cut-off, two of them: both started at zero; the other parameters as fitted, and
  # the fit not touched.
  th <- base
  th[diag_i[c(1, 3)]] <- c(0.009, 2e-6)
  placed <- place(th)
  started <- maihda_refit_start(placed)
  want <- th
  want[diag_i[c(1, 3)]] <- 0
  expect_identical(started@pp$theta, unname(want))
  expect_identical(placed@pp$theta, unname(th))
  expect_identical(lme4::getME(placed, "theta"), lme4::getME(started, "theta"))

  # Fits made with an optimiser whose first step is not tied to the bound: as fitted.
  for (opt in c("nlminbwrap", "Nelder_Mead")) {
    other <- placed
    other@optinfo$optimizer <- opt
    expect_identical(maihda_refit_start(other), other, info = opt)
  }

  # At or above the cut-off, at zero, or an unbounded (covariance) parameter near zero: as fitted.
  for (v in c(0.01, 0.011, 0)) {
    th <- base
    th[diag_i[1]] <- v
    placed <- place(th)
    expect_identical(maihda_refit_start(placed), placed)
  }
  th <- base
  th[off_i] <- 1e-9
  placed <- place(th)
  expect_identical(maihda_refit_start(placed), placed)
})

test_that("with precision weights the cut-off is on the residual standard deviation reported", {
  d <- flat_strata(stratum_shift)
  # Weights of 1e4: the residual variance reported is sigma^2 / 1e4, so theta measures the
  # stratum standard deviation against a residual one 100 times smaller than sigma.
  d$w <- 1e4
  m <- quiet(lme4::lmer(y ~ 1 + (1 | a:b:c), data = d, weights = w))
  # A stratum standard deviation of 0.8 of the reported residual one: as fitted.
  inside <- fit_at(m, 0.008)
  expect_identical(maihda_refit_start(inside), inside)
  # 0.005 of it: started at zero.
  near <- fit_at(m, 5e-5)
  expect_identical(maihda_refit_start(near)@pp$theta, 0)
})

test_that("an nlminb fit just above zero is refitted from the fit, as lme4::refit() refits it", {
  m <- quiet(lme4::lmer(y ~ 1 + (1 | a:b:c), data = flat_strata(), start = list(theta = 1e-4),
                        control = lme4::lmerControl(optimizer = "nlminbwrap")))
  expect_identical(m@optinfo$optimizer, "nlminbwrap")
  m <- fit_at(m, 1e-4)
  y_new <- flat_strata(stratum_shift)$y
  ours <- quiet(maihda_refit_draw(m, y_new))
  theirs <- quiet(lme4::refit(m, newresp = y_new))
  expect_identical(lme4::getME(ours, "theta"), lme4::getME(theirs, "theta"))
  expect_identical(criterion_of(ours), criterion_of(theirs))
})

test_that("a block whose optimiser parameters are not theta is left as fitted", {
  skip_if(utils::packageVersion("lme4") < "2.0-0")
  d <- flat_strata(stratum_shift)
  m <- quiet(lme4::lmer(y ~ x + cs(1 + x | a:b:c), data = d))
  # cs() holds a correlation among its parameters, so lme4::getME(, "par") is not theta, and a
  # theta entry of 0.005 is not a variance parameter of its own: the fit is returned as it is.
  theta <- lme4::getME(m, "theta")
  theta[1] <- 0.005
  pp <- m@pp$copy()
  pp$setTheta(unname(theta))
  m@pp <- pp
  m@theta <- unname(theta)
  expect_false(isTRUE(all.equal(unname(lme4::getME(m, "par")), unname(theta))))
  expect_identical(maihda_refit_start(m), m)
})

test_that("a fit without such a parameter is refitted exactly as lme4::refit() refits it", {
  d <- flat_strata(stratum_shift)
  m <- quiet(lme4::lmer(y ~ 1 + (1 | a:b:c), data = d))
  expect_gt(lme4::getME(m, "theta"), 0.01)
  y_new <- flat_strata(rev(stratum_shift))$y
  ours <- quiet(maihda_refit_draw(m, y_new))
  theirs <- quiet(lme4::refit(m, newresp = y_new))
  expect_identical(lme4::getME(ours, "theta"), lme4::getME(theirs, "theta"))
  expect_identical(lme4::fixef(ours), lme4::fixef(theirs))
  expect_identical(criterion_of(ours), criterion_of(theirs))
})
