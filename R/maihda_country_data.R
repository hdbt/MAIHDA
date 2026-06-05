#' Cross-National Educational Achievement Data for MAIHDA
#'
#' A cross-national dataset for demonstrating how Multilevel Analysis of
#' Individual Heterogeneity and Discriminatory Accuracy (MAIHDA) can be used to
#' compare intersectional inequality \emph{across} a higher-level grouping
#' variable (here, country) with \code{\link{compare_maihda_groups}} and
#' \code{\link{maihda}}. Each row is a 15-year-old student; the intersectional
#' strata are formed by \code{gender} and socioeconomic status (\code{ses}), and
#' the outcome is the PISA mathematics score.
#'
#' Intersectional inequality (the between-stratum share of variance, VPC/ICC) in
#' mathematics achievement differs across the six countries, which is what makes
#' the dataset a useful showcase for the group-comparison workflow.
#'
#' @format A data frame with 3,600 rows (600 students in each of 6 countries) and
#'   7 variables:
#' \describe{
#'   \item{country}{Factor; one of Finland, Germany, United Kingdom, Italy, Japan,
#'     Mexico. The higher-level grouping variable.}
#'   \item{gender}{Factor; student gender (female/male). A stratum dimension.}
#'   \item{ses}{Factor; socioeconomic status as global tertiles (Low/Medium/High)
#'     of \code{escs}, computed on the pooled sample so a band means the same in
#'     every country. A stratum dimension.}
#'   \item{escs}{Numeric; the PISA index of economic, social and cultural status
#'     (the continuous measure underlying \code{ses}).}
#'   \item{math}{Numeric; PISA mathematics score (first plausible value). The
#'     primary outcome.}
#'   \item{reading}{Numeric; PISA reading score (first plausible value).}
#'   \item{low_math}{Factor; "Yes" if \code{math} is below 420 (PISA proficiency
#'     Level 2 baseline), else "No". A binary outcome for logistic examples.}
#' }
#'
#' @details
#' The intersectional strata are \code{gender:ses} (2 x 3 = 6 strata). A canonical
#' MAIHDA "null" model is \code{math ~ 1 + (1 | gender:ses)}; comparing its VPC
#' across countries quantifies how much joint gender-by-class inequality in
#' achievement varies between countries.
#'
#' @source
#' Derived from the OECD Programme for International Student Assessment (PISA)
#' 2018 student questionnaire data (OECD (2019), \emph{PISA 2018 Database}),
#' accessed and cleaned via the \pkg{learningtower} R package (MIT licensed),
#' \url{https://CRAN.R-project.org/package=learningtower}. A balanced random
#' subsample of 600 complete-case students per country was taken (seed 2026). The
#' data preparation script is in \code{data-raw/maihda_country_data.R}.
#'
#' @note
#' This is a teaching/illustration dataset only. It uses a single PISA plausible
#' value for each score and the analysis functions ignore the PISA survey weights
#' and complex sampling design, so results are \strong{not} survey-representative
#' and should not be used for substantive cross-national inference.
#'
#' @examples
#' \donttest{
#' data(maihda_country_data)
#'
#' # Compare intersectional (gender x SES) inequality across countries
#' analysis <- maihda(
#'   math ~ 1 + (1 | gender:ses),
#'   data = maihda_country_data,
#'   group = "country"
#' )
#' analysis
#' plot(analysis, type = "group_vpc")
#' }
"maihda_country_data"
