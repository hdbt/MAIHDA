#' Simulated Cross-National Youth Data for MAIHDA
#'
#' A simulated, cross-sectional sample of 16--29-year-olds in four European
#' countries, for demonstrating MAIHDA on the questions youth-research institutes
#' study: the school-to-work transition (\strong{NEET} status -- Not in Employment,
#' Education or Training) and youth \strong{subjective well-being}, at the
#' intersection of \strong{gender, migration background, and social origin}
#' (parental education). It supports the two-model decomposition
#' (\code{\link{maihda}}), the binary discriminatory-accuracy workflow
#' (\code{family = "binomial"}), and the cross-group comparison
#' (\code{maihda(group = "country")}).
#'
#' The 3 stratum dimensions form 2 x 2 x 3 = 12 intersectional strata, evaluated
#' within each of 4 countries (the higher-level grouping variable). Two outcomes
#' are generated from one shared latent disadvantage, so they correlate as in
#' reality. The between-stratum differences are mostly additive (the dimensions'
#' main effects) with a genuine intersectional interaction -- the migration penalty
#' is steeper at low social origin, and again for young women -- built orthogonal to
#' the main effects so it survives adjustment and is recoverable (a sub-1 PCV). The
#' intensity of intersectional inequality is amplified in the higher-NEET countries,
#' so the VPC/ICC genuinely differs across the grouping variable.
#'
#' @format A data frame with 3000 rows (750 per country) and 8 variables:
#' \describe{
#'   \item{id}{Person identifier.}
#'   \item{country}{Country (\code{Luxembourg}/\code{Germany}/\code{France}/\code{Italy});
#'     the higher-level grouping variable for \code{maihda(group = "country")}.}
#'   \item{gender}{Gender (\code{Women}/\code{Men}); a stratum dimension.}
#'   \item{migration}{Migration background (\code{None}/\code{Migration background});
#'     a stratum dimension.}
#'   \item{parental_edu}{Parental education / social origin
#'     (\code{Low}/\code{Medium}/\code{High}); a stratum dimension.}
#'   \item{age}{Age in years (16--29); an individual-level covariate.}
#'   \item{neet}{NEET status (\code{No}/\code{Yes}); the binary outcome.}
#'   \item{wellbeing}{Subjective well-being (a life-satisfaction index, higher is
#'     better); the continuous outcome.}
#' }
#'
#' @details The exact generative parameters are stored as
#'   \code{attr(maihda_youth_data, "truth")}: the additive main effects, the
#'   standardised interaction per stratum (\code{interaction_by_stratum}), the
#'   between-stratum variances and VPCs (at amplitude 1), the per-country amplitude
#'   multipliers, and the calibration targets. The split a fitted model
#'   \emph{recovers} (the PCV) differs from the construction target by outcome,
#'   because the logit link attenuates the interaction for the binary outcome and
#'   pooling countries of different amplitude inflates it for the Gaussian one --
#'   read the PCV the model reports rather than the nominal target.
#'
#' @source Simulated for the MAIHDA package; \strong{not real microdata}. The design
#'   (age band, dimensions, and outcomes) follows established youth surveys -- the
#'   Youth Survey Luxembourg
#'   (\url{https://www.youth-in-luxembourg.lu/project/youth-survey-luxembourg-ysl/}),
#'   the German Youth Institute's AID:A
#'   (\url{https://www.dji.de/en/about-us/projects/projekte/aida.html}), and the
#'   youth sub-sample of the European Social Survey -- and the NEET prevalences and
#'   life-satisfaction levels are calibrated to Eurostat figures for 15--29-year-olds
#'   (\url{https://ec.europa.eu/eurostat/web/education-and-training/data}). The real
#'   surveys require data-use agreements and cannot be redistributed in a package, so
#'   the bundled data are synthetic, following the same convention as the package's
#'   other simulated datasets (\code{\link{maihda_sim_data}},
#'   \code{\link{maihda_long_data}}, \code{\link{maihda_sparse_data}}).
#'
#' @seealso \code{\link{maihda}}, \code{\link{maihda_discriminatory_accuracy}},
#'   \code{\link{maihda_interactions}}
#'
#' @examples
#' data(maihda_youth_data)
#' \donttest{
#' # Well-being: the two-model VPC/ICC + PCV (additive vs intersectional):
#' wb <- maihda(wellbeing ~ gender + migration + parental_edu +
#'                (1 | gender:migration:parental_edu), data = maihda_youth_data)
#' wb
#'
#' # NEET: the binary discriminatory-accuracy workflow (AUC / MOR ride along):
#' neet <- maihda(neet ~ gender + migration + parental_edu +
#'                  (1 | gender:migration:parental_edu),
#'                data = maihda_youth_data, family = "binomial")
#' neet
#'
#' # Compare intersectional inequality across countries:
#' by_country <- maihda(neet ~ gender + migration + parental_edu +
#'                        (1 | gender:migration:parental_edu),
#'                      data = maihda_youth_data, group = "country",
#'                      family = "binomial")
#' by_country
#' }
"maihda_youth_data"
