# Audit pass 2026-10-27: lme4::bootMer() and confint(method = "boot") on the lme4 fits
# fit_maihda() makes.
#
# fit_maihda() binds its `...` values in the formula's environment and names them in the stored
# call (control = .maihda_arg_control). bootMer() evaluates the call's control where it is called,
# so it stopped ("object '.maihda_arg_control' not found") on every fit made with a control and on
# every glmer fit a refit replaced. Per the user, the call gives a control with the user's own
# optimiser settings -- the user's control, and that of a refitted glmer fit's last pass when the
# pass took Nelder-Mead's own first steps -- as lme4 writes a control given in the call, a call to
# lmerControl() or glmerControl(); a control with first steps of the package's own keeps its name,
# since bootMer() would refit each draw from the fit with those steps and stop short, as does a
# control whose optimiser is a function. lme4 2.1-0's bootMer() also looks the control up in the
# formula's environment, so there it finds a control kept by its name and refits each draw with it.
#
# Random numbers are drawn only inside with_seed(), which puts the suite's stream back.

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
deviance_of <- function(m) -2 * as.numeric(stats::logLik(m))
estimates_of <- function(f) c(lme4::fixef(f), lme4::getME(f, "theta"))

# Twelve strata of 25 rows with a covariate on a scale of 1, a normal, a 0 / 1 and a count outcome.
small_data <- function() {
  with_seed(11, {
    s <- rep(seq_len(12), each = 25)
    x <- stats::rnorm(300)
    u <- stats::rnorm(12, 0, 0.5)
    data.frame(y = 1 + u[s] + 0.3 * x + stats::rnorm(300),
               b = stats::rbinom(300, 1, stats::plogis(u[s] + 0.3 * x)), x = x, stratum = factor(s),
               n = stats::rnbinom(300, size = 3, mu = exp(0.5 + u[s] + 0.2 * x)))
  })
}

# The control of a refitted glmer fit's last pass: the user's (glmer's default here) without glmer's
# first stage, and with first steps of 0.1% of each estimate where Nelder-Mead's own fail.
finish_control <- function(steps = FALSE) {
  ctrl <- lme4::glmerControl()
  ctrl$nAGQ0initStep <- FALSE
  if (steps) {
    ctrl$optCtrl$xst <- c(1e-3, 1e-6, 2.5e-4)
    ctrl$optCtrl$xt <- ctrl$optCtrl$xst * 5e-4
  }
  ctrl
}
# The control of a refit by bobyqa with a first step of its own.
bobyqa_control <- function() {
  ctrl <- finish_control()
  ctrl$optimizer[2] <- "bobyqa"
  ctrl$optCtrl <- list(rhobeg = 1e-5, rhoend = 1e-8, maxfun = 1e5)
  ctrl
}

# ---- the call's control ----------------------------------------------------------------------

test_that("a control is written as the call to lmerControl() or glmerControl() that gives it back", {
  env <- new.env()
  user <- as.name(".maihda_arg_control")
  controls <- list(lme4::lmerControl(optimizer = "bobyqa"), lme4::lmerControl(optimizer = NULL),
                   lme4::lmerControl(check.conv.singular = "ignore"), lme4::glmerControl(),
                   lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                   lme4::glmerControl(tolPwrss = 1e-8, check.rankX = "ignore"))
  for (ctl in controls) {
    assign(".maihda_arg_control", ctl, envir = env)
    written <- maihda_lme4_call_control(user, env)
    expect_true(is.call(written))
    expect_identical(eval(written, baseenv()), ctl)
  }
  # Only the settings that differ from the function's defaults are written.
  assign(".maihda_arg_control", lme4::glmerControl(optimizer = "bobyqa"), envir = env)
  expect_identical(maihda_lme4_call_control(user, env),
                   bquote(lme4::glmerControl(optimizer = .(c("bobyqa", "bobyqa")))))
  assign(".maihda_arg_control", lme4::lmerControl(), envir = env)
  expect_identical(maihda_lme4_call_control(user, env), quote(lme4::lmerControl()))
})

test_that("a last pass's control is written only with the user's own optimiser settings", {
  env <- new.env()
  finish <- as.name(".maihda_arg_finish_control")
  # No control of the user's: glmer's default without its first stage.
  assign(".maihda_arg_finish_control", finish_control(), envir = env)
  expect_identical(maihda_lme4_call_control(finish, env),
                   quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  # With first steps of its own, the name is kept.
  assign(".maihda_arg_finish_control", finish_control(steps = TRUE), envir = env)
  expect_identical(maihda_lme4_call_control(finish, env), finish)
  # The user's own settings, an evaluation limit among them.
  user <- lme4::glmerControl(optCtrl = list(maxfun = 2e5), tolPwrss = 1e-8)
  assign(".maihda_arg_control", user, envir = env)
  pass <- user
  pass$nAGQ0initStep <- FALSE
  assign(".maihda_arg_finish_control", pass, envir = env)
  written <- maihda_lme4_call_control(finish, env)
  expect_true(is.call(written))
  expect_identical(eval(written, baseenv()), pass)
  pass$optCtrl$xst <- c(1e-3, 1e-6, 2.5e-4)
  assign(".maihda_arg_finish_control", pass, envir = env)
  expect_identical(maihda_lme4_call_control(finish, env), finish)
})

test_that("a bobyqa refit's control keeps its name, as does a control the call does not give back", {
  env <- new.env()
  for (nm in c(".maihda_arg_refit_control", ".maihda_arg_search_control")) {
    assign(nm, bobyqa_control(), envir = env)
    expect_identical(maihda_lme4_call_control(as.name(nm), env), as.name(nm))
  }
  # Also where the refit took the user's own bobyqa settings.
  own <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  assign(".maihda_arg_control", own, envir = env)
  own$nAGQ0initStep <- FALSE
  assign(".maihda_arg_refit_control", own, envir = env)
  expect_identical(maihda_lme4_call_control(as.name(".maihda_arg_refit_control"), env),
                   as.name(".maihda_arg_refit_control"))
  rm(".maihda_arg_control", envir = env)
  user <- as.name(".maihda_arg_control")
  # Not bound, not a control of lme4's, a setting no argument gives.
  expect_identical(maihda_lme4_call_control(user, env), user)
  assign(".maihda_arg_control", list(optimizer = "bobyqa"), envir = env)
  expect_identical(maihda_lme4_call_control(user, env), user)
  odd <- lme4::glmerControl()
  odd$extra <- 1
  assign(".maihda_arg_control", odd, envir = env)
  expect_identical(maihda_lme4_call_control(user, env), user)
  # An optimiser given as a function.
  for (ctl in list(lme4::lmerControl(optimizer = lme4::nloptwrap),
                   lme4::glmerControl(optimizer = lme4::Nelder_Mead),
                   lme4::glmerControl(optimizer = list("bobyqa", lme4::Nelder_Mead)))) {
    assign(".maihda_arg_control", ctl, envir = env)
    expect_identical(maihda_lme4_call_control(user, env), user)
  }
  # A call without a control, or with one written out, is left as it is.
  expect_null(maihda_lme4_call_control(NULL, env))
  expect_identical(maihda_lme4_call_control(quote(lme4::lmerControl()), env),
                   quote(lme4::lmerControl()))
})

# ---- bootMer() and confint() on the fits ------------------------------------------------------

test_that("bootMer() and confint(method = 'boot') run on lme4 fits made with a control, as on lme4's own", {
  skip_on_cran()
  skip_slow_test()
  d <- small_data()
  lmer_control <- lme4::lmerControl(optimizer = "bobyqa")
  glmer_control <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  m <- quiet(fit_maihda(y ~ x + (1 | stratum), data = d, control = lmer_control))$model
  g <- quiet(fit_maihda(b ~ x + (1 | stratum), data = d, family = "binomial",
                        control = glmer_control))$model
  expect_identical(stats::getCall(m)$control, quote(lme4::lmerControl(optimizer = "bobyqa")))
  expect_identical(eval(stats::getCall(g)$control, baseenv()), glmer_control)
  # lme4's own fits name the control; bootMer() finds it in this frame.
  own_m <- lme4::lmer(y ~ x + (1 | stratum), data = d, control = lmer_control)
  own_g <- quiet(lme4::glmer(b ~ x + (1 | stratum), data = d, family = stats::binomial,
                             control = glmer_control))
  expect_identical(estimates_of(m), estimates_of(own_m))
  expect_identical(estimates_of(g), estimates_of(own_g))
  for (pair in list(list(m, own_m), list(g, own_g))) {
    b <- with_seed(1, quiet(lme4::bootMer(pair[[1]], estimates_of, nsim = 4, seed = 3)))
    b_own <- with_seed(1, quiet(lme4::bootMer(pair[[2]], estimates_of, nsim = 4, seed = 3)))
    expect_identical(unname(b$t), unname(b_own$t))
    expect_false(anyNA(b$t))
  }
  ci <- with_seed(1, quiet(stats::confint(g, method = "boot", nsim = 4, seed = 3, parm = "beta_")))
  expect_identical(dim(ci), c(2L, 2L))
  expect_false(anyNA(ci))
  # update() makes the fit again; summary() prints the control.
  expect_identical(estimates_of(stats::update(m)), estimates_of(m))
  expect_true(any(grepl("^Control: lme4::lmerControl\\(optimizer = \"bobyqa\"\\)$",
                        utils::capture.output(summary(m)))))
})

test_that("growth, negative-binomial, context and maihda() fits give the control written out too", {
  skip_on_cran()
  skip_slow_test()
  lmer_control <- lme4::lmerControl(optimizer = "bobyqa")
  glmer_control <- lme4::glmerControl(optimizer = "bobyqa")
  written_lmer <- quote(lme4::lmerControl(optimizer = "bobyqa"))
  growth <- quiet(fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
                             data = maihda_long_data, id = "id", time = "wave",
                             control = lmer_control))$model
  expect_identical(stats::getCall(growth)$control, written_lmer)
  nb <- quiet(fit_maihda(n ~ x + (1 | stratum), data = small_data(), family = "negbinomial",
                         control = glmer_control))$model
  expect_identical(stats::getCall(nb)$control,
                   bquote(lme4::glmerControl(optimizer = .(c("bobyqa", "bobyqa")))))
  b <- with_seed(1, quiet(lme4::bootMer(nb, estimates_of, nsim = 2, seed = 3)))
  expect_false(anyNA(b$t))
  context <- quiet(fit_maihda(math ~ 1 + (1 | gender:ses), data = maihda_country_data,
                              context = "country", control = lmer_control))$model
  expect_identical(stats::getCall(context)$control, written_lmer)
  two <- quiet(maihda(health_outcome ~ age + (1 | gender:race:education), data = maihda_sim_data,
                      control = lmer_control))
  expect_identical(stats::getCall(two$model$model)$control, written_lmer)
  expect_identical(stats::getCall(two$model_adjusted$model)$control, written_lmer)
})

test_that("a refitted glmer fit's last control is written with the user's steps, not with its own", {
  skip_on_cran()
  skip_slow_test()
  d <- small_data()
  # Stand-ins for maihda_fit_glmer_bound(): the fit as a refit makes it, the call naming the
  # refit's control bound in the fit's environment.
  made_by <- function(name, control) {
    function(fit_call, fit_env, tol = 1e-4) {
      assign(name, control, envir = fit_env)
      fit_call$control <- as.name(name)
      eval(fit_call, fit_env)
    }
  }
  fit_by <- function(name, control) {
    with_mocked_bindings(quiet(fit_maihda(b ~ x + (1 | stratum), data = d, family = "binomial")),
                         maihda_fit_glmer_bound = made_by(name, control))$model
  }
  f <- fit_by(".maihda_arg_finish_control", finish_control())
  expect_identical(stats::getCall(f)$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  b <- with_seed(1, quiet(lme4::bootMer(f, estimates_of, nsim = 3, seed = 3)))
  expect_false(anyNA(b$t))
  expect_identical(estimates_of(stats::update(f)), estimates_of(f))
  for (kept in list(list(".maihda_arg_finish_control", finish_control(steps = TRUE)),
                    list(".maihda_arg_search_control", bobyqa_control()),
                    list(".maihda_arg_refit_control", bobyqa_control()))) {
    s <- fit_by(kept[[1]], kept[[2]])
    expect_identical(stats::getCall(s)$control, as.name(kept[[1]]))
    if (utils::packageVersion("lme4") < "2.1-0") {
      expect_error(with_seed(1, lme4::bootMer(s, estimates_of, nsim = 1, seed = 3)), kept[[1]],
                   fixed = TRUE)
    } else {
      b <- with_seed(1, quiet(lme4::bootMer(s, estimates_of, nsim = 1, seed = 3)))
      draw <- with_seed(1, stats::simulate(s, nsim = 1, seed = 3))[[1]]
      expect_equal(unname(b$t[1, ]),
                   unname(estimates_of(quiet(lme4::refit(s, draw, control = kept[[2]])))))
    }
    expect_identical(estimates_of(stats::update(s)), estimates_of(s))
  }
})

# ---- a live fit ----------------------------------------------------------------------------

# Eighteen strata of 80 rows, each row with a birth year between 1940 and 2000 and a Poisson count
# whose stratum slope deviates around 1970 (the data of test-audit-2026-10-26-refit-default-step.R):
# fit_maihda() refits glmer's fit from its searches and then by Nelder-Mead from there, with
# Nelder-Mead's own first steps.
year_counts <- function(seed) {
  with_seed(seed, {
    lev <- expand.grid(d1 = 1:3, d2 = 1:3, d3 = 1:2)
    add <- stats::rnorm(3, 0, 0.3)[lev$d1] + stats::rnorm(3, 0, 0.3)[lev$d2] +
      stats::rnorm(2, 0, 0.3)[lev$d3]
    v <- stats::rnorm(18, 0, 0.01)
    st <- rep(seq_len(18), each = 80)
    x <- sample(1940:2000, length(st), replace = TRUE)
    eta <- 0.5 + add[st] + (0.01 + v[st] / 2) * (x - 1970)
    data.frame(y = stats::rpois(length(st), exp(eta)), x = x, stratum = factor(st))
  })
}

test_that("bootMer() refits the draws of a fit fit_maihda() refitted by Nelder-Mead", {
  skip_on_cran()
  skip_slow_test()
  m <- quiet(fit_maihda(y ~ x + (1 + x | stratum), data = year_counts(1), family = "poisson"))$model
  expect_identical(m@optinfo$optimizer, "Nelder_Mead")
  expect_identical(stats::getCall(m)$control, quote(lme4::glmerControl(nAGQ0initStep = FALSE)))
  b <- with_seed(1, quiet(lme4::bootMer(m, deviance_of, nsim = 2, seed = 1)))
  expect_false(anyNA(b$t))
  expect_equal(deviance_of(quiet(stats::update(m))), deviance_of(m), tolerance = 1e-8)
})
