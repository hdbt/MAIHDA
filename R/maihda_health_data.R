#' NHANES Health Data Subset for MAIHDA Use
#'
#' A pedagogical subset of the National Health and Nutrition Examination Survey
#' (NHANES) dataset, serving as a real-world example for Multilevel Analysis
#' of Individual Heterogeneity and Discriminatory Accuracy (MAIHDA). Contains
#' selected records demonstrating intersectional demographic health inequalities.
#'
#' @format A data frame with 3,000 rows and 7 variables:
#' \describe{
#'   \item{BMI}{Body Mass Index (kg/m^2), a continuous outcome variable.}
#'   \item{Obese}{Factor indicating obesity status (No/Yes).}
#'   \item{Age}{Age in years at screening, a continuous covariate.}
#'   \item{Gender}{Gender of the participant (male/female).}
#'   \item{Race}{Self-reported race/ethnicity.}
#'   \item{Education}{Educational attainment level.}
#'   \item{Poverty}{Poverty to income ratio, a continuous covariate. Some values
#'     may be missing.}
#' }
#'
#' @source
#' Derived from the \code{NHANES} R package. Original data collected by the
#' Centers for Disease Control and Prevention (CDC).
#'
#' @note
#' This is a teaching/illustration dataset only. It is a non-random subsample and
#' does \strong{not} carry the NHANES survey weights or complex sampling design,
#' so results are \strong{not} survey-representative and should not be used for
#' substantive population inference. (For your own survey data, the package
#' supports design-weighted MAIHDA via the \code{sampling_weights} argument of
#' \code{\link{fit_maihda}} / \code{\link{maihda}}.)
#'
#' @examples
#' data(maihda_health_data)
#'
#' # Example usage:
#' # strata_result <- make_strata(maihda_health_data, vars = c("Gender", "Race", "Education"))
#' # model <- fit_maihda(BMI ~ Age + (1 | stratum), data = strata_result$data)
"maihda_health_data"
