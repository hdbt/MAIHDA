#' Simulated Longitudinal Data for MAIHDA
#'
#' A simulated long-format (one row per person-occasion) panel for demonstrating
#' \strong{longitudinal / growth-curve MAIHDA} (\code{fit_maihda(id =, time =)} and
#' \code{maihda(decomposition = "longitudinal")}). 600 individuals are each
#' measured over five waves on a continuous wellbeing score, within 12
#' intersectional strata defined by gender x ethnicity x education. The
#' between-stratum \emph{trajectory} differences are constructed to be mostly
#' additive (each dimension's main effect on the baseline level and on the rate of
#' change) with one genuine multiplicative interaction, so the longitudinal PCV --
#' a high but sub-1 \code{PCV_slope} -- is demonstrable.
#'
#' @format A data frame with 3000 rows (600 persons x 5 waves) and 8 variables:
#' \describe{
#'   \item{id}{Person identifier (level 2); repeated across waves.}
#'   \item{wave}{Measurement occasion, 0 to 4 (the numeric time variable).}
#'   \item{gender}{Gender (\code{Women}/\code{Men}); a stratum dimension.}
#'   \item{ethnicity}{Ethnicity (\code{EthA}/\code{EthB}/\code{EthC}); a stratum dimension.}
#'   \item{education}{Education (\code{Low}/\code{High}); a stratum dimension.}
#'   \item{age}{Baseline age in years, a time-invariant covariate.}
#'   \item{wellbeing}{Continuous wellbeing outcome (the growth-curve response).}
#'   \item{low_wellbeing}{Binary companion outcome (1 = bottom 40\% of wellbeing),
#'     for the logistic longitudinal path.}
#' }
#'
#' @source Simulated for the purpose of the MAIHDA package. The growth structure
#'   follows the longitudinal MAIHDA of Bell, Evans, Holman & Leckie (2024)
#'   \doi{10.1016/j.socscimed.2024.116955}.
#'
#' @examples
#' data(maihda_long_data)
#' \donttest{
#' # Time-varying VPC from a 3-level growth model:
#' m <- fit_maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
#'                 data = maihda_long_data, id = "id", time = "wave")
#' summary(m)
#'
#' # Additive-vs-multiplicative PCV (null vs adjusted growth model):
#' a <- maihda(wellbeing ~ wave + (1 | gender:ethnicity:education),
#'             data = maihda_long_data, id = "id", time = "wave",
#'             decomposition = "longitudinal")
#' a$pcv
#' }
"maihda_long_data"
