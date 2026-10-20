# data-raw/maihda_sparse_data.R
# Simulate `maihda_sparse_data`, the data of vignette("bayesian_sparse_maihda"): a
# cross-sectional intersectional dataset in which many strata hold only a handful of
# individuals, with a known intersectional interaction.
#
#   * 4 dimensions -> 36 intersectional strata (gender x ethnicity x education x age)
#   * deliberately skewed stratum sizes (median ~6, several singletons)
#   * stratum effects built from an additive part (one effect per level of each
#     dimension) and an interaction part, each standardised over the 36 strata and
#     weighted for an interaction share of 40% of a between-stratum variance of 0.35
#     on a Gaussian outcome `y` and of 1.5 on the latent scale of a binary outcome
#     `event`, the two outcomes sharing one intersectional structure.
#
# The two parts are drawn independently, but over 36 strata they come out correlated,
# and part of the interaction part lies along the dimensions' main effects, so the
# drawn effects do not carry those targets exactly. attr(maihda_sparse_data, "truth")
# holds what they carry: the intercept, the 36 stratum effects, their variance, its
# split into additive and interaction parts by mean squares (as the realised-strata
# block of a brms crossed-dimensions summary splits it), the interaction share and the
# VPC. R/maihda_sparse_data.R quotes these values.
#
# Base R only. Run from the package root:  Rscript data-raw/maihda_sparse_data.R

set.seed(2026)

## ---- intersectional dimensions: 2 x 3 x 2 x 3 = 36 strata ----
gender    <- c("Women", "Men")
ethnicity <- c("White", "Black", "Asian")
education <- c("Low", "High")
age_group <- c("Young", "Mid", "Older")

grid <- expand.grid(gender = gender, ethnicity = ethnicity,
                    education = education, age_group = age_group,
                    stringsAsFactors = FALSE)
S <- nrow(grid)
stratum_key <- apply(grid, 1, paste, collapse = ":")

## ---- generative stratum components (standardised to unit variance) ----
# Additive part: a sum of independent dimension main effects, the part the
# crossed-dimensions model's per-dimension random intercepts describe.
dim_fx <- function(levels, sd) stats::setNames(stats::rnorm(length(levels), 0, sd), levels)
b_gender <- dim_fx(gender, 0.30); b_eth <- dim_fx(ethnicity, 0.30)
b_edu    <- dim_fx(education, 0.30); b_age <- dim_fx(age_group, 0.30)
raw_add  <- b_gender[grid$gender] + b_eth[grid$ethnicity] +
            b_edu[grid$education] + b_age[grid$age_group]
# Interaction part: a per-stratum effect drawn independently of the additive part.
raw_int  <- stats::rnorm(S, 0, 0.45)

std <- function(x) (x - mean(x)) / stats::sd(x)
add0 <- std(raw_add)   # var 1 across the 36 strata
int0 <- std(raw_int)   # var 1; its correlation with add0 is -0.12 in this draw

## ---- place the components on each outcome's scale, aiming at the targets ----
target_share   <- 0.40   # interaction's share of the between-stratum variance
# Gaussian: target between-stratum variance 0.35 (VPC 0.26 with residual SD 1).
between_g <- 0.35
s_add <- sqrt(between_g * (1 - target_share))
s_int <- sqrt(between_g * target_share)
sigma_e <- 1.0
# Binary (latent logit scale): target between-stratum variance 1.5 (latent VPC 0.31).
between_b <- 1.5
a_add <- sqrt(between_b * (1 - target_share))
a_int <- sqrt(between_b * target_share)
beta0 <- 0.0             # ~50% marginal prevalence

g_stratum <- 50 + s_add * add0 + s_int * int0       # Gaussian stratum means
b_eta     <- beta0 + a_add * add0 + a_int * int0     # binary latent stratum effect
names(g_stratum) <- names(b_eta) <- stratum_key

## ---- skewed allocation of N individuals -> sparse cells ----
N <- 240L
p <- stats::runif(S, 0.25, 1)^1.5
p <- p / sum(p)
who <- sample(seq_len(S), N, replace = TRUE, prob = p)

maihda_sparse_data <- data.frame(
  gender    = factor(grid$gender[who],    levels = gender),
  ethnicity = factor(grid$ethnicity[who], levels = ethnicity),
  education = factor(grid$education[who],  levels = education),
  age_group = factor(grid$age_group[who], levels = age_group),
  y         = g_stratum[who] + stats::rnorm(N, 0, sigma_e),
  event     = factor(ifelse(stats::runif(N) < stats::plogis(b_eta[who]), "Yes", "No"),
                     levels = c("No", "Yes")),
  stringsAsFactors = FALSE
)
rownames(maihda_sparse_data) <- NULL

## ---- record what the drawn stratum effects carry ----
# The interaction variance is the residual mean square of the effects' least-squares
# fit on the four dimensions' levels (36 - 7 = 29 degrees of freedom), the between
# variance the variance of the 36 effects (35); the additive part is the difference.
X <- stats::model.matrix(~ gender + ethnicity + education + age_group, grid)
drawn_truth <- function(intercept, effects, residual_var) {
  between <- stats::var(effects)
  interaction <- sum(stats::lm.fit(X, effects)$residuals^2) / (S - ncol(X))
  list(additive_var = between - interaction, interaction_var = interaction,
       between_var = between, residual_var = residual_var,
       interaction_share = interaction / between,
       vpc = between / (between + residual_var),
       intercept = intercept, stratum_effects = effects)
}
truth <- list(
  target_interaction_share = target_share,
  gaussian = drawn_truth(50, g_stratum - 50, sigma_e^2),
  binary_latent = drawn_truth(beta0, b_eta - beta0, pi^2 / 3),
  n_strata = S, n = N
)
attr(maihda_sparse_data, "truth") <- truth

## ---- report (paste these numbers into the roxygen doc) ----
cells <- table(factor(apply(maihda_sparse_data[1:4], 1, paste, collapse = ":"),
                      levels = stratum_key))
cat(sprintf("N=%d  strata occupied=%d/%d  median cell=%g  min=%g  singletons=%d  empty=%d  below 5=%d\n",
            N, sum(cells > 0), S, stats::median(as.numeric(cells)),
            min(as.numeric(cells)), sum(cells == 1), sum(cells == 0), sum(cells < 5)))
cat(sprintf("Gaussian: between=%.3f resid=%.3f VPC=%.3f interaction share=%.1f%%\n",
            truth$gaussian$between_var, sigma_e^2, truth$gaussian$vpc,
            100 * truth$gaussian$interaction_share))
cat(sprintf("Binary (latent): between=%.3f VPC=%.3f interaction share=%.1f%%  prevalence=%.2f\n",
            truth$binary_latent$between_var, truth$binary_latent$vpc,
            100 * truth$binary_latent$interaction_share,
            mean(maihda_sparse_data$event == "Yes")))

usethis::use_data(maihda_sparse_data, overwrite = TRUE)
