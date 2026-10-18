# The realised strata of a brms crossed-dimensions fit: finite-population summaries
# (Gelman 2005, eq. 8) of the effects of the strata and dimension levels in the data.
# In each posterior draw stratum j has the effect zeta_j, the sum of its dimensions'
# level effects and its own (interaction) effect, and with J strata:
#   between_var       sum_j (zeta_j - mean(zeta))^2 / (J - 1)
#   vpc               between_var / (between_var + within-stratum variance)
#   interaction_var   residual sum of squares of the least-squares fit of the zeta_j
#                     on an intercept and the dimensions' levels, over its residual
#                     degrees of freedom (df_interaction)
#   additive_share    1 - interaction_var / between_var; interaction_share the rest
#   per_dim           the variance of each dimension's level effects (divisor
#                     levels - 1)
#   mor               exp(median over the pairs of strata of |zeta_j - zeta_j'|)
# each summarised by the posterior median and the equal-tailed interval. The oracle
# numbers below come from a separate implementation of these definitions that does not
# call the package.

# Evaluate `expr` under a fixed seed and put the caller's RNG state back, so this file
# leaves the random stream the later test files read as it found it (see setup-seed.R).
crs_with_seed <- function(seed, expr) {
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
crs_print <- function(expr) cli::ansi_strip(utils::capture.output(expr))

# Five draws of a 2 x 3 grid of strata, built without random numbers.
crs_draws <- function() {
  cells <- expand.grid(a = c("p", "q"), b = c("x", "y", "z"), stringsAsFactors = FALSE)
  nd <- 5L
  Ua <- matrix(sin(1:10), nd, 2, dimnames = list(NULL, c("p", "q")))
  Ub <- matrix(cos(1:15) / 2, nd, 3, dimnames = list(NULL, c("x", "y", "z")))
  V <- matrix(sin(seq_len(30) * 2.3) / 3, nd, 6)
  Z <- V + Ua[, match(cells$a, colnames(Ua))] + Ub[, match(cells$b, colnames(Ub))]
  list(Z = Z, U = list(a = Ua, b = Ub), V = V, cells = cells)
}

crs_q <- function(x, conf_level = 0.95) {
  a <- 1 - conf_level
  c(stats::median(x), stats::quantile(x, c(a / 2, 1 - a / 2), names = FALSE))
}

# A crossed-dimensions maihda_model shell around a brmsfit stub that carries its own
# posterior draws (an `r_` column marks that level effects were kept) and level-effect
# draws for the mocked reader.
crs_stub <- function(family = "bernoulli", link = "logit", level_draws = TRUE,
                     context = NULL) {
  dr <- crs_draws()
  strata <- paste0("s", seq_len(nrow(dr$cells)))
  dat <- dr$cells
  dat$stratum <- strata
  draws <- data.frame(b_Intercept = rep(0, 5), sd_stratum__Intercept = seq(0.1, 0.5, 0.1),
                      sd_ga__Intercept = seq(0.2, 1, 0.2), sd_gb__Intercept = seq(0.3, 0.7, 0.1),
                      sigma = seq(0.9, 1.3, 0.1))
  if (!is.null(context)) {
    draws$sd_site__Intercept <- seq(0.4, 0.8, 0.1)
  }
  if (level_draws) {
    draws[["r_stratum[s1,Intercept]"]] <- dr$V[, 1]
  }
  V <- dr$V
  colnames(V) <- strata
  levels <- list(ga = dr$U$a, gb = dr$U$b, stratum = V)
  names(dat)[1:2] <- c("ga", "gb")
  fam <- list(family = family, link = link)
  structure(
    list(model = structure(list(formula = y ~ 1, family = fam, draws = draws,
                                levels = levels), class = "brmsfit"),
         engine = "brms", formula = y ~ 1, family = fam, data = dat,
         cc_info = list(dim_groups = c(sex = "ga", region = "gb"),
                        interaction_group = "stratum"),
         context_info = if (!is.null(context)) list(context_vars = context)),
    class = "maihda_model"
  )
}

crs_mock_brms <- function(env = parent.frame()) {
  local_mocked_bindings(
    maihda_posterior_draws_brms = function(model) model$draws,
    maihda_level_draws_brms = function(model, groups) model$levels[groups],
    .env = env
  )
}

test_that("the realised quantities follow their definitions", {
  dr <- crs_draws()
  r <- maihda_cc_realised_from_draws(dr$Z, dr$U, dr$cells, within = pi^2 / 3, mor = TRUE)

  # The separate implementation's values, at full precision.
  expect_equal(c(r$between_var, r$between_var_ci),
               c(0.44656435949565426, 0.092994521902311847, 0.55258524542524146),
               tolerance = 1e-12)
  expect_equal(c(r$vpc, r$vpc_ci),
               c(0.11951623916912939, 0.027145555431101198, 0.14379480870481617),
               tolerance = 1e-12)
  expect_equal(c(r$interaction_var, r$interaction_var_ci),
               c(0.043768542796975382, 0.041630493247242391, 0.044700504295260804),
               tolerance = 1e-12)
  expect_equal(c(r$additive_share, r$additive_share_ci),
               c(0.89984403799278434, 0.43266224971002115, 0.92354778488610967),
               tolerance = 1e-12)
  expect_equal(c(r$interaction_share, r$interaction_share_ci),
               c(0.10015596200721566, 0.076452215113890332, 0.56733775028997913),
               tolerance = 1e-12)
  expect_equal(c(r$per_dim[["a"]], unname(r$per_dim_ci["a", ])),
               c(0.35975405468066696, 0.037254571056467803, 0.67768864183552668),
               tolerance = 1e-12)
  expect_equal(c(r$per_dim[["b"]], unname(r$per_dim_ci["b", ])),
               c(0.098142174164853296, 0.059088052946498479, 0.21568321313586841),
               tolerance = 1e-12)
  expect_equal(c(r$mor, r$mor_ci),
               c(2.3188500506607195, 1.3932691309267047, 2.4927692319903141),
               tolerance = 1e-12)

  # The same, from the definitions written out with base R.
  pairs <- utils::combn(6, 2)
  between <- apply(dr$Z, 1, stats::var)
  # J - 1 is half the mean squared difference between two distinct strata.
  expect_equal(between, apply(dr$Z, 1, function(z) mean((z[pairs[1, ]] - z[pairs[2, ]])^2) / 2))
  inter <- apply(dr$Z, 1, function(z) {
    f <- stats::lm(z ~ a + b, data = cbind(dr$cells, z = z))
    sum(stats::resid(f)^2) / f$df.residual
  })
  expect_equal(c(r$between_var, r$between_var_ci), crs_q(between))
  expect_equal(c(r$vpc, r$vpc_ci), crs_q(between / (between + pi^2 / 3)))
  expect_equal(c(r$interaction_var, r$interaction_var_ci), crs_q(inter))
  expect_equal(c(r$additive_share, r$additive_share_ci), crs_q(1 - inter / between))
  expect_equal(c(r$interaction_share, r$interaction_share_ci), crs_q(inter / between))
  expect_equal(c(r$per_dim[["a"]], unname(r$per_dim_ci["a", ])), crs_q(apply(dr$U$a, 1, stats::var)))
  expect_equal(c(r$per_dim[["b"]], unname(r$per_dim_ci["b", ])), crs_q(apply(dr$U$b, 1, stats::var)))
  expect_equal(c(r$mor, r$mor_ci), crs_q(exp(apply(dr$Z, 1, function(z) {
    stats::median(abs(z[pairs[1, ]] - z[pairs[2, ]]))
  }))))
  expect_identical(r$n_strata, 6L)
  expect_identical(r$df_interaction, 2L)
  expect_identical(r$conf_level, 0.95)
  expect_identical(r$method, "posterior")
  expect_identical(names(r), c(
    "between_var", "between_var_ci", "vpc", "vpc_ci", "additive_share",
    "additive_share_ci", "interaction_share", "interaction_share_ci", "interaction_var",
    "interaction_var_ci", "per_dim", "per_dim_ci", "mor", "mor_ci", "n_strata",
    "df_interaction", "conf_level", "method"))
  expect_identical(dimnames(r$per_dim_ci), list(c("a", "b"), c("lower", "upper")))

  # The interval follows conf_level, and the point estimate can be the posterior mean.
  r80 <- maihda_cc_realised_from_draws(dr$Z, dr$U, dr$cells, conf_level = 0.8)
  expect_equal(c(r80$between_var, r80$between_var_ci), crs_q(between, 0.8))
  rm <- maihda_cc_realised_from_draws(dr$Z, dr$U, dr$cells, point = "mean")
  expect_equal(rm$between_var, 0.37361931307412299, tolerance = 1e-12)
  expect_identical(rm$between_var_ci, r$between_var_ci)
  # A within-stratum variance per draw is used draw by draw.
  w <- c(1, 2, 3, 4, 5)
  rw <- maihda_cc_realised_from_draws(dr$Z, dr$U, dr$cells, within = w)
  expect_equal(c(rw$vpc, rw$vpc_ci), crs_q(between / (between + w)))
})

test_that("strata missing from the grid, no interaction degrees of freedom, few strata", {
  dr <- crs_draws()
  # A 2 x 3 grid with one stratum missing: least squares over the five strata present.
  r <- maihda_cc_realised_from_draws(dr$Z[, -3], dr$U, dr$cells[-3, ])
  expect_identical(r$df_interaction, 1L)
  expect_identical(r$n_strata, 5L)
  expect_equal(c(r$interaction_var, r$interaction_var_ci),
               c(0.025504425680117446, 0.0056137043621405361, 0.071136512237367666),
               tolerance = 1e-12)
  inter <- apply(dr$Z[, -3], 1, function(z) {
    sum(stats::resid(stats::lm(z ~ a + b, data = cbind(dr$cells[-3, ], z = z)))^2)
  })
  expect_equal(r$interaction_var, stats::median(inter))
  # Asked for neither, the block has no VPC and no median odds ratio.
  expect_false(any(c("vpc", "vpc_ci", "mor", "mor_ci") %in% names(r)))

  # Three strata of a 2 x 2 grid: the dimensions fit them exactly, so there is no
  # interaction to estimate; the between-strata variance is still given.
  c3 <- data.frame(a = c("p", "q", "p"), b = c("x", "x", "y"), stringsAsFactors = FALSE)
  r3 <- maihda_cc_realised_from_draws(dr$Z[, 1:3], list(a = dr$U$a), c3, within = 1)
  expect_identical(r3$df_interaction, 0L)
  expect_true(is.na(r3$interaction_var) && all(is.na(r3$interaction_var_ci)))
  expect_true(is.na(r3$additive_share) && is.na(r3$interaction_share))
  expect_equal(r3$between_var, stats::median(apply(dr$Z[, 1:3], 1, stats::var)))
  expect_true(is.finite(r3$vpc))

  # A dimension with one level among the strata adds nothing to the fit and has no
  # variance of its own.
  c1 <- cbind(dr$cells, k = "only")
  r1 <- maihda_cc_realised_from_draws(dr$Z, c(dr$U, list(k = matrix(1, 5, 1))), c1)
  r0 <- maihda_cc_realised_from_draws(dr$Z, dr$U, dr$cells)
  expect_identical(r1$df_interaction, r0$df_interaction)
  expect_equal(r1$interaction_var, r0$interaction_var)
  expect_true(is.na(r1$per_dim[["k"]]))

  # Fewer than two strata: nothing to summarise.
  expect_null(maihda_cc_realised_from_draws(dr$Z[, 1, drop = FALSE], dr$U, dr$cells[1, ]))
})

test_that("the level effects are read from brms's draws by level name", {
  # What brms::ranef(summary = FALSE) returns: draws x levels x coefficients arrays.
  # Levels out of alphabetical order, the intercept not the first coefficient.
  sex <- array(c(1:6, 101:106), dim = c(3, 2, 2),
               dimnames = list(NULL, c("m", "f"), c("age", "Intercept")))
  one <- array(7:9, dim = c(3, 1, 1), dimnames = list(NULL, "only", "Intercept"))
  re <- list(sex = sex, one = one)
  lv <- maihda_level_draws_from_ranef(re, c("one", "sex"))
  expect_identical(names(lv), c("one", "sex"))
  expect_identical(lv$sex, matrix(as.integer(101:106), 3, 2,
                                  dimnames = list(NULL, c("m", "f"))))
  expect_identical(lv$one, matrix(7:9, 3, 1, dimnames = list(NULL, "only")))
  # One draw keeps its matrix shape.
  lv1 <- maihda_level_draws_from_ranef(list(sex = sex[1, , , drop = FALSE]), "sex")
  expect_identical(lv1$sex, matrix(c(101L, 104L), 1, 2, dimnames = list(NULL, c("m", "f"))))
  expect_error(maihda_level_draws_from_ranef(re, "region"), "No 'region' random effects")
  expect_error(maihda_level_draws_from_ranef(list(sex = sex[, , "age", drop = FALSE]), "sex"),
               "no intercept draws")
})

test_that("each stratum's effects are matched by name, whatever order brms keeps", {
  skip_if_not_installed("brms")
  crs_mock_brms()
  m <- crs_stub()
  dr <- crs_draws()
  ref <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)$decomposition$realised
  # The same level effects with every grouping factor's columns in another order.
  lv <- m$model$levels
  m$model$levels <- list(ga = lv$ga[, c("q", "p")], gb = lv$gb[, c("z", "x", "y")],
                         stratum = lv$stratum[, c(4, 6, 1, 5, 3, 2)])
  expect_equal(maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)$decomposition$realised,
               ref)
  # and the fitted data's rows in yet another order
  m$data <- m$data[c(5, 2, 6, 1, 4, 3), ]
  expect_equal(maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)$decomposition$realised,
               ref)
})

test_that("each stratum of the fitted data is read once, with its levels", {
  d <- data.frame(ga = c("p", "p", "q", NA, "q", "p"), gb = c("x", "x", "y", "y", "y", "z"),
                  stratum = c("s1", "s1", "s4", "s4", "s4", "s5"), y = 1:6,
                  stringsAsFactors = FALSE)
  d$gb <- factor(d$gb)
  cells <- maihda_cc_realised_cells(d, c("ga", "gb"), "stratum")
  # Rows repeat a stratum; a row with a missing level is not one; everything is text.
  expect_identical(cells, data.frame(ga = c("p", "q", "p"), gb = c("x", "y", "z"),
                                     stratum = c("s1", "s4", "s5"),
                                     stringsAsFactors = FALSE))
  # A stratum with two levels on a dimension is not a stratum of the crossed design.
  d$ga[2] <- "q"
  expect_error(maihda_cc_realised_cells(d, c("ga", "gb"), "stratum"), "more than one level")
  expect_error(maihda_cc_realised_cells(d, c("ga", "nope"), "stratum"),
               "do not hold the stratum and dimension columns")
  expect_error(maihda_cc_realised_cells(NULL, "ga", "stratum"),
               "do not hold the stratum and dimension columns")
})

test_that("the median odds ratio is taken over every pair, in blocks of draws", {
  dr <- crs_draws()
  pairs <- utils::combn(6, 2)
  brute <- apply(dr$Z, 1, function(z) stats::median(abs(z[pairs[1, ]] - z[pairs[2, ]])))
  # Fifteen pairs: blocks of one, two and all five draws give the same medians.
  for (cells in c(15, 30, 1e6)) {
    expect_equal(maihda_pair_median_draws(dr$Z, max_cells = cells), brute)
  }
  # Two strata: one pair.
  expect_equal(maihda_pair_median_draws(dr$Z[, 1:2]), abs(dr$Z[, 1] - dr$Z[, 2]))

  # Every draw while the work is small, otherwise at least 1000 evenly spaced ones.
  expect_identical(maihda_mor_draw_rows(4000, 8), seq_len(4000))
  rows <- maihda_mor_draw_rows(4000, 500)
  expect_length(rows, 1000L)
  expect_identical(range(rows), c(1L, 4000L))
  expect_false(is.unsorted(rows, strictly = TRUE))
  expect_length(maihda_mor_draw_rows(4000, 125), 2580L)
  expect_identical(maihda_mor_draw_rows(300, 500), seq_len(300))
})

test_that("a brms crossed summary carries the realised strata of its draws", {
  skip_if_not_installed("brms")
  crs_mock_brms()
  m <- crs_stub()
  dr <- crs_draws()
  expect_no_warning(res <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95))
  d <- res$decomposition
  # Placed after the level counts; the elements of today are untouched.
  expect_identical(names(d), c(
    "additive_var", "interaction_var", "between_var", "within_var", "additive_share",
    "interaction_share", "additive_share_ci", "interaction_share_ci", "per_dim",
    "dim_levels", "realised", "bootstrap", "method", "conf_level"))
  m0 <- crs_stub(level_draws = FALSE)
  d0 <- maihda_cc_summary_brms(m0, m0$cc_info, conf_level = 0.95)$decomposition
  expect_identical(d[names(d) != "realised"], d0)
  expect_identical(res$vpc_result, maihda_cc_summary_brms(m0, m0$cc_info, 0.95)$vpc_result)

  # The stratum effects are the stratum's own effect plus its levels' effects, named by
  # the dimensions, and the summary is the one the pure helper gives on them.
  z <- dr$V + dr$U$a[, match(dr$cells$a, c("p", "q"))] +
    dr$U$b[, match(dr$cells$b, c("x", "y", "z"))]
  expect_equal(d$realised, maihda_cc_realised_from_draws(
    z, list(sex = dr$U$a, region = dr$U$b), stats::setNames(dr$cells, c("sex", "region")),
    within = pi^2 / 3, mor = TRUE))
  expect_identical(names(d$realised$per_dim), c("sex", "region"))
  expect_equal(d$realised$between_var, stats::median(apply(z, 1, stats::var)))

  # A credible level reaches the block.
  r90 <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.9)$decomposition$realised
  expect_identical(r90$conf_level, 0.9)
  expect_equal(r90$between_var_ci, crs_q(apply(z, 1, stats::var), 0.9)[2:3])
})

test_that("the VPC line and the median odds ratio follow the family and the context", {
  skip_if_not_installed("brms")
  crs_mock_brms()
  dr <- crs_draws()
  z <- dr$V + dr$U$a[, match(dr$cells$a, c("p", "q"))] +
    dr$U$b[, match(dr$cells$b, c("x", "y", "z"))]
  between <- apply(z, 1, stats::var)
  give <- function(family, link, context = FALSE, within = 1) {
    m <- crs_stub(family, link)
    maihda_cc_realised_brms(m, m$cc_info, m$model$draws, rep(within, 5),
                            conf_level = 0.95, context = context)
  }
  # Latent scales: the within-stratum variance is a constant of the link.
  r <- give("bernoulli", "logit", within = pi^2 / 3)
  expect_equal(r$vpc, stats::median(between / (between + pi^2 / 3)))
  expect_false(is.null(r$mor))
  r <- give("binomial", "logit")
  expect_false(is.null(r$vpc) || is.null(r$mor))
  r <- give("cumulative", "logit")
  expect_false(is.null(r$vpc) || is.null(r$mor))
  r <- give("binomial", "probit")
  expect_false(is.null(r$vpc))
  expect_null(r$mor)
  r <- give("bernoulli", "cloglog")
  expect_false(is.null(r$vpc))
  expect_null(r$mor)
  # A Gaussian fit's residual variance per draw.
  r <- give("gaussian", "identity", within = 1.7)
  expect_equal(r$vpc, stats::median(between / (between + 1.7)))
  expect_null(r$mor)
  # Count families: their within-stratum variance is built from the variance components.
  for (fam in c("poisson", "negbinomial")) {
    r <- give(fam, "log")
    expect_false(any(c("vpc", "vpc_ci", "mor", "mor_ci") %in% names(r)))
    expect_equal(r$between_var, stats::median(between))
  }
  # A context would enter the denominator: no VPC; the rest is the same.
  r_ctx <- give("bernoulli", "logit", context = TRUE, within = pi^2 / 3)
  r_no <- give("bernoulli", "logit", within = pi^2 / 3)
  expect_false(any(c("vpc", "vpc_ci") %in% names(r_ctx)))
  expect_identical(r_ctx[setdiff(names(r_no), c("vpc", "vpc_ci"))],
                   r_no[setdiff(names(r_no), c("vpc", "vpc_ci"))])

  # Through the summary: a contextual brms crossed fit puts the realised variance of the
  # context effects into the VPC's denominator.
  m <- crs_stub(context = "site")
  m$data$site <- rep(c("u", "v"), 3)
  site <- matrix(c(0.3, -0.2, 0.1, 0.4, -0.1, -0.3, 0.2, 0, -0.4, 0.1), 5, 2,
                 dimnames = list(NULL, c("u", "v")))
  m$model$levels$site <- site
  expect_no_warning(d <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)$decomposition)
  ctx <- apply(site, 1, stats::var)
  expect_equal(c(d$realised$vpc, d$realised$vpc_ci),
               crs_q(between / (between + ctx + pi^2 / 3)))
  expect_false(is.null(d$realised$mor))
  # Without the context's level effects there is no VPC, and a warning says why.
  m$model$levels$site <- NULL
  expect_warning(d0 <- maihda_cc_summary_brms(m, m$cc_info, conf_level = 0.95)$decomposition,
                 "realised contextual partition")
  expect_false(is.null(d0$realised))
  expect_null(d0$realised$vpc)
  expect_false(is.null(d0$realised$mor))
})

test_that("no level-effect draws, no block; a failure to read them warns", {
  skip_if_not_installed("brms")
  crs_mock_brms()
  m0 <- crs_stub(level_draws = FALSE)
  expect_no_warning(d0 <- maihda_cc_summary_brms(m0, m0$cc_info, 0.95)$decomposition)
  expect_false("realised" %in% names(d0))
  expect_null(maihda_cc_realised_brms(m0, m0$cc_info, m0$model$draws, rep(1, 5), 0.95))

  # Level effects the posterior should hold but cannot give: the block is dropped with
  # a warning, and everything else is summarised as usual.
  local_mocked_bindings(maihda_level_draws_brms = function(model, groups) {
    stop("no level draws here")
  })
  m <- crs_stub()
  expect_warning(d <- maihda_cc_summary_brms(m, m$cc_info, 0.95)$decomposition,
                 "no level draws here")
  expect_false("realised" %in% names(d))
  expect_identical(d, d0)

  # A stratum or level of the data the posterior lacks is an error the summary turns
  # into that warning.
  local_mocked_bindings(maihda_level_draws_brms = function(model, groups) {
    lv <- model$levels[groups]
    lv$stratum <- lv$stratum[, -1, drop = FALSE]
    lv
  })
  expect_warning(maihda_cc_summary_brms(m, m$cc_info, 0.95), "no effect in the brms posterior")
})

test_that("print() shows the block, and the few-levels note points to it", {
  dr <- crs_draws()
  r <- maihda_cc_realised_from_draws(dr$Z, list(sex = dr$U$a, region = dr$U$b),
                                     stats::setNames(dr$cells, c("sex", "region")),
                                     within = pi^2 / 3, mor = TRUE)
  d <- list(additive_var = 1.5, interaction_var = 0.5, between_var = 2, within_var = 3.29,
            additive_share = 0.75, interaction_share = 0.25,
            additive_share_ci = c(0.5, 0.9), interaction_share_ci = c(0.1, 0.5),
            per_dim = c(sex = 1, region = 0.5), dim_levels = c(sex = 2L, region = 3L),
            realised = r, bootstrap = FALSE, method = "posterior", conf_level = 0.95)
  out <- crs_print(maihda_print_cc_decomposition(d))
  f4 <- function(x) sprintf("%.4f", x)
  pc <- function(x) sprintf("%.1f%%", 100 * x)
  block <- c(
    "Realised strata (the 6 strata in the data):",
    "  (Posterior 95% credible interval)",
    sprintf("  Between-strata variance:  %s [%s, %s]", f4(r$between_var),
            f4(r$between_var_ci[1]), f4(r$between_var_ci[2])),
    sprintf("  VPC:                      %s [%s, %s]", f4(r$vpc), f4(r$vpc_ci[1]),
            f4(r$vpc_ci[2])),
    sprintf("  Interaction variance:     %s [%s, %s]", f4(r$interaction_var),
            f4(r$interaction_var_ci[1]), f4(r$interaction_var_ci[2])),
    sprintf("  Additive share:           %s [%s, %s]", pc(r$additive_share),
            pc(r$additive_share_ci[1]), pc(r$additive_share_ci[2])),
    sprintf("  Interaction share:        %s [%s, %s]", pc(r$interaction_share),
            pc(r$interaction_share_ci[1]), pc(r$interaction_share_ci[2])),
    "  Variance of each dimension's level effects:",
    sprintf("    sex: %s [%s, %s]", f4(r$per_dim[["sex"]]), f4(r$per_dim_ci["sex", 1]),
            f4(r$per_dim_ci["sex", 2])),
    sprintf("    region: %s [%s, %s]", f4(r$per_dim[["region"]]),
            f4(r$per_dim_ci["region", 1]), f4(r$per_dim_ci["region", 2])),
    sprintf("  Median odds ratio:        %.3f [%.3f, %.3f]", r$mor, r$mor_ci[1], r$mor_ci[2]),
    "  Note: finite-population summaries of the effects of these strata and levels",
    "  (Gelman 2005); they depend little on the prior when each stratum has enough",
    "  data. See ?maihda.",
    "")
  expect_identical(utils::tail(out, length(block)), block)
  # The few-levels note above it ends by pointing to it.
  note <- maihda_cc_few_levels_note(c(sex = 2L, region = 3L), pointer = TRUE)
  i <- match(note[1], out)
  expect_identical(out[i + seq_along(note) - 1L], note)
  expect_identical(out[i + length(note)], "")
  expect_identical(out[i + length(note) + 1L], block[1])
  expect_true(grepl(paste("; the realised-strata block below depends little on it when",
                          "each stratum has enough data. See ?maihda."),
                    paste(trimws(note), collapse = " "), fixed = TRUE))

  # Without the block, the note keeps its own ending and nothing follows it.
  d0 <- d
  d0$realised <- NULL
  out0 <- crs_print(maihda_print_cc_decomposition(d0))
  expect_identical(out0, c(out[seq_len(i - 1L)],
                           maihda_cc_few_levels_note(c(sex = 2L, region = 3L)), ""))
  expect_false(any(grepl("realised", out0, fixed = TRUE)))

  # A block without a VPC or a median odds ratio leaves their lines out; one without
  # interaction degrees of freedom says so; a share below zero prints as it is.
  r2 <- r
  r2[c("vpc", "vpc_ci", "mor", "mor_ci")] <- NULL
  r2$additive_share <- -0.123
  out2 <- crs_print(maihda_print_cc_realised(r2))
  expect_false(any(grepl("^  (VPC|Median odds ratio):", out2)))
  expect_true("  Additive share:           -12.3% [43.3%, 92.4%]" %in% out2)
  r3 <- r
  r3$df_interaction <- 0L
  out3 <- crs_print(maihda_print_cc_realised(r3))
  expect_false(any(grepl("^  (Interaction variance|Additive share|Interaction share):",
                         out3)))
  expect_true(all(c("  Interaction variance and shares: NA (the strata in the data leave no",
                    "  degrees of freedom for an interaction)") %in% out3))
})

test_that("print() of the analysis and of its summary show the block", {
  g <- expand.grid(d1 = c("a", "b"), d2 = c("p", "q", "r"), stringsAsFactors = FALSE)
  dg <- g[rep(seq_len(nrow(g)), each = 20), ]
  rownames(dg) <- NULL
  dg$y <- rep(c(0, 0.5, 0.3, 0.9, 0.6, 1.4), each = 20) + sin(seq_len(nrow(dg)) * 12.9898)
  a <- crs_with_seed(1, suppressWarnings(suppressMessages(
    maihda(y ~ 1 + (1 | d1:d2), data = dg, decomposition = "crossed-dimensions"))))
  out_a <- crs_print(print(a))
  out_s <- crs_print(print(a$summary))
  # An lme4 decomposition has no block.
  expect_false("realised" %in% names(a$decomposition))
  expect_false(any(grepl("Realised strata", c(out_a, out_s), fixed = TRUE)))

  # The same analysis with the decomposition of a brms fit that carries the block:
  # both print methods show it below the decomposition and the note, and nothing else
  # in either print moves.
  dr <- crs_draws()
  r <- maihda_cc_realised_from_draws(dr$Z, list(d1 = dr$U$a, d2 = dr$U$b),
                                     stats::setNames(dr$cells, c("d1", "d2")),
                                     within = 1)
  b <- a
  b$decomposition$method <- "posterior"
  b$decomposition$dim_levels <- c(d1 = 2L, d2 = 3L)
  b$decomposition$realised <- r
  b$summary$decomposition <- b$decomposition
  note <- maihda_cc_few_levels_note(c(d1 = 2L, d2 = 3L), pointer = TRUE)
  block <- crs_print(maihda_print_cc_realised(r))
  for (pair in list(list(crs_print(print(b)), out_a),
                    list(crs_print(print(b$summary)), out_s))) {
    out_b <- pair[[1]]
    i <- match(note[1], out_b)
    expect_identical(out_b[i + seq_along(note) - 1L], note)
    k <- i + length(note) + 1L
    expect_identical(out_b[k + seq_along(block) - 1L], block)
    expect_identical(out_b[-c(i + seq_along(note) - 1L, k + seq_along(block) - 1L)],
                     pair[[2]])
  }
})

test_that("the help page defines the realised strata and cites the source", {
  skip_on_cran()
  man <- testthat::test_path("..", "..", "man")
  # Present when the suite runs from the package tree; R CMD check copies only tests/.
  skip_if_not(dir.exists(man), "man/ is not available from this test run")
  rd <- gsub("[[:space:]]+", " ",
             paste(readLines(file.path(man, "maihda.Rd"), warn = FALSE), collapse = " "))
  expect_true(grepl("\\section{Realised strata}{", rd, fixed = TRUE))
  expect_true(grepl("10.1214/009053604000001048", rd, fixed = TRUE))
  expect_true(grepl("(J - 1)", rd, fixed = TRUE))
})

test_that("a real brms fit: the block from the fit's own level effects", {
  # Compiles a Stan model, so it is OPT-IN: set MAIHDA_TEST_BRMS=true to run it.
  skip_on_cran()
  skip_if(Sys.getenv("MAIHDA_TEST_BRMS") != "true",
          "brms Stan tests are opt-in; set MAIHDA_TEST_BRMS=true to run them")
  skip_if_not_installed("brms")

  # 2 x 3 x 4 strata (one left out) of 40 binomial trials each, additive effects plus a
  # small stratum effect. brms draws its sampler seed from R's stream, so the fit sits
  # inside the seeded block too.
  a <- crs_with_seed(20261014, {
    g <- expand.grid(d1 = c("a", "b"), d2 = c("c", "d", "e"), d3 = c("f", "g", "h", "i"),
                     stringsAsFactors = FALSE)
    eta <- -0.9 + c(0, 0.5)[match(g$d1, c("a", "b"))] +
      c(0, 0.3, 0.6)[match(g$d2, c("c", "d", "e"))] +
      c(0, 0.2, 0.4, 0.7)[match(g$d3, c("f", "g", "h", "i"))] + stats::rnorm(24, 0, 0.2)
    g$n <- 40L
    g$s <- stats::rbinom(24, 40, stats::plogis(eta))
    # One stratum missing, so that a level read under another level's name changes the
    # strata's effects rather than permuting a full grid.
    g <- g[-1, ]
    suppressWarnings(suppressMessages(
      maihda(s | trials(n) ~ 1 + (1 | d1:d2:d3), data = g, engine = "brms",
             family = "binomial", decomposition = "crossed-dimensions",
             chains = 2, iter = 1000, refresh = 0)))
  })
  r <- a$decomposition$realised
  expect_false(is.null(r))
  expect_identical(a$summary$decomposition$realised, r)
  expect_identical(r$n_strata, 23L)
  expect_identical(r$df_interaction, 16L)

  # The stratum effects assembled here from brms's own level-effect draws.
  re <- brms::ranef(a$model$model, summary = FALSE)
  dat <- unique(a$model$data[, c("d1", "d2", "d3", "stratum")])
  z <- re$stratum[, as.character(dat$stratum), "Intercept"]
  for (k in c("d1", "d2", "d3")) {
    z <- z + re[[k]][, as.character(dat[[k]]), "Intercept"]
  }
  between <- apply(z, 1, stats::var)
  inter <- apply(z, 1, function(v) {
    sum(stats::resid(stats::lm(v ~ d1 + d2 + d3, data = cbind(dat, v = v)))^2) / 16
  })
  pairs <- utils::combn(23, 2)
  mor <- exp(apply(z, 1, function(v) stats::median(abs(v[pairs[1, ]] - v[pairs[2, ]]))))
  expect_equal(c(r$between_var, r$between_var_ci), crs_q(between), tolerance = 1e-10)
  expect_equal(c(r$vpc, r$vpc_ci), crs_q(between / (between + pi^2 / 3)), tolerance = 1e-10)
  expect_equal(c(r$interaction_var, r$interaction_var_ci), crs_q(inter), tolerance = 1e-10)
  expect_equal(c(r$additive_share, r$additive_share_ci), crs_q(1 - inter / between),
               tolerance = 1e-10)
  expect_equal(c(r$mor, r$mor_ci), crs_q(mor), tolerance = 1e-10)
  expect_equal(r$per_dim[["d2"]], stats::median(apply(re$d2[, , "Intercept"], 1, stats::var)),
               tolerance = 1e-10)

  # Both print methods show the block, under a note that points to it.
  note <- maihda_cc_few_levels_note(c(d1 = 2L, d2 = 3L, d3 = 4L), pointer = TRUE)
  for (out in list(crs_print(print(a)), crs_print(print(a$summary)))) {
    i <- match(note[1], out)
    expect_identical(out[i + seq_along(note) - 1L], note)
    expect_identical(out[i + length(note) + 1L], "Realised strata (the 23 strata in the data):")
  }
})
