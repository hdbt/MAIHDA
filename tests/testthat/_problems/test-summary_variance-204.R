# Extracted from test-summary_variance.R:204

# test -------------------------------------------------------------------------
bars <- reformulas::findbars(y ~ x + (1 | g) + (t | id))
d <- data.frame(g = c("a", "b", "a"), id = c(1, 1, 2), t = c(0, 2, 1))
blocks <- list(
    matrix(0.5, 1, 1, dimnames = list("(Intercept)", "(Intercept)")),
    matrix(c(1, 0.2, 0.2, 0.3), 2, 2,
           dimnames = list(c("(Intercept)", "t"), c("(Intercept)", "t")))
  )
v <- MAIHDA:::maihda_latent_re_variance_rows(bars, blocks, d)
expect_equal(v, 0.5 + 1 + 2 * d$t * 0.2 + d$t^2 * 0.3)
