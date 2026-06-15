# data-raw/maihda_sparse_data.R
# Simulate `maihda_sparse_data`: a cross-sectional intersectional dataset built to
# showcase **Bayesian (brms) MAIHDA for sparse intersections**.
#
# The point of the dataset is the regime where frequentist ML (lme4) breaks down:
# many intersectional strata, MANY of them with only a handful of individuals (some
# singletons). In that regime the maximum-likelihood estimate of the *interaction*
# between-stratum variance is pinned at the boundary (a singular fit) and carries no
# uncertainty -- so the additive-vs-interaction split is both unstable and falsely
# precise. Weakly-informative priors (engine = "brms") regularise the variance off
# the boundary and return a calibrated credible interval.
#
# To make that demonstrable, the data carry a KNOWN, non-trivial interaction:
#   * 4 dimensions -> 36 intersectional strata (gender x ethnicity x education x age)
#   * deliberately skewed stratum sizes (median ~6, several singletons)
#   * a true interaction share of 40% of the between-stratum variance, on BOTH a
#     Gaussian outcome `y` and the latent scale of a binary outcome `event`, sharing
#     one intersectional structure.
# The realised true variances are stored as attr(maihda_sparse_data, "truth") and
# documented in R/maihda_sparse_data.R so the vignette can claim *recovery*, not just
# report numbers.
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
# Additive part: a sum of independent dimension main effects -> recoverable by the
# crossed-dimensions model's per-dimension random intercepts.
dim_fx <- function(levels, sd) stats::setNames(stats::rnorm(length(levels), 0, sd), levels)
b_gender <- dim_fx(gender, 0.30); b_eth <- dim_fx(ethnicity, 0.30)
b_edu    <- dim_fx(education, 0.30); b_age <- dim_fx(age_group, 0.30)
raw_add  <- b_gender[grid$gender] + b_eth[grid$ethnicity] +
            b_edu[grid$education] + b_age[grid$age_group]
# Interaction part: a genuine per-stratum departure from additivity (independent of
# the additive component).
raw_int  <- stats::rnorm(S, 0, 0.45)

std <- function(x) (x - mean(x)) / stats::sd(x)
add0 <- std(raw_add)   # var 1 across the 36 strata
int0 <- std(raw_int)   # var 1, ~uncorrelated with add0

## ---- place the components on each outcome's scale to hit exact targets ----
target_share   <- 0.40   # interaction's share of the between-stratum variance
# Gaussian: between-stratum variance 0.35 (VPC ~= 0.26 with residual SD 1).
between_g <- 0.35
s_add <- sqrt(between_g * (1 - target_share))
s_int <- sqrt(between_g * target_share)
sigma_e <- 1.0
# Binary (latent logit scale): between-stratum variance 1.5 (latent VPC ~= 0.31).
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

## ---- record the ground truth (exact, because we built it) ----
truth <- list(
  target_interaction_share = target_share,
  gaussian = list(additive_var = s_add^2, interaction_var = s_int^2,
                  between_var = between_g, residual_var = sigma_e^2,
                  interaction_share = s_int^2 / between_g,
                  vpc = between_g / (between_g + sigma_e^2)),
  binary_latent = list(additive_var = a_add^2, interaction_var = a_int^2,
                       between_var = between_b, residual_var = pi^2 / 3,
                       interaction_share = a_int^2 / between_b,
                       vpc = between_b / (between_b + pi^2 / 3)),
  n_strata = S, n = N
)
attr(maihda_sparse_data, "truth") <- truth

## ---- report (paste these numbers into the roxygen doc) ----
cells <- table(factor(apply(maihda_sparse_data[1:4], 1, paste, collapse = ":"),
                      levels = stratum_key))
cat(sprintf("N=%d  strata occupied=%d/%d  median cell=%g  min=%g  singletons=%d  empty=%d\n",
            N, sum(cells > 0), S, stats::median(as.numeric(cells)),
            min(as.numeric(cells)), sum(cells == 1), sum(cells == 0)))
cat(sprintf("Gaussian: between=%.3f resid=%.3f VPC=%.3f interaction share=%.1f%%\n",
            between_g, sigma_e^2, truth$gaussian$vpc, 100 * truth$gaussian$interaction_share))
cat(sprintf("Binary (latent): between=%.3f VPC=%.3f interaction share=%.1f%%  prevalence=%.2f\n",
            between_b, truth$binary_latent$vpc, 100 * truth$binary_latent$interaction_share,
            mean(maihda_sparse_data$event == "Yes")))

usethis::use_data(maihda_sparse_data, overwrite = TRUE)
