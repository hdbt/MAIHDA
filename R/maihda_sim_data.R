#' Simulated Health Data for MAIHDA Use
#'
#' A simulated dataset for demonstrating Multilevel Analysis of Individual
#' Heterogeneity and Discriminatory Accuracy (MAIHDA).
#'
#' @format A data frame with 500 rows and 6 variables:
#' \describe{
#'   \item{id}{Unique participant identifier.}
#'   \item{gender}{Gender of the participant.}
#'   \item{race}{Simulated race/ethnicity category.}
#'   \item{education}{Educational attainment level.}
#'   \item{age}{Age in years, a continuous covariate.}
#'   \item{health_outcome}{A continuous simulated health outcome.}
#' }
#'
#' @source Simulated for the purpose of the MAIHDA package.
#'
#' @examples
#' data(maihda_sim_data)
#' strata_result <- make_strata(maihda_sim_data, vars = c("gender", "race", "education"))
"maihda_sim_data"
