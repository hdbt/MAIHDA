\name{MAIHDA-package}
\alias{MAIHDA-package}
\alias{MAIHDA}
\docType{package}
\title{
Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy
}
\description{
The MAIHDA package provides tools for conducting Multilevel Analysis of Individual
Heterogeneity and Discriminatory Accuracy. This approach is useful for examining
intersectional inequalities in health and other outcomes.
}
\section{Main Functions}{
\itemize{
  \item \code{\link{make_strata}}: Create intersectional strata from multiple variables
  \item \code{\link{fit_maihda}}: Fit multilevel models using lme4 or brms
  \item \code{\link{summary_maihda}}: Summarize models with variance partition and stratum estimates
  \item \code{\link{predict_maihda}}: Make predictions at individual or stratum level
  \item \code{\link{plot_maihda}}: Visualize results with various plot types
  \item \code{\link{compare_maihda}}: Compare models with bootstrap confidence intervals
}
}
\author{
MAIHDA Package Authors
}
\references{
Evans, C. R., Williams, D. R., Onnela, J. P., & Subramanian, S. V. (2018). 
A multilevel approach to modeling health inequalities at the intersection of 
multiple social identities. Social Science & Medicine, 203, 64-73.

Merlo, J. (2018). Multilevel analysis of individual heterogeneity and 
discriminatory accuracy (MAIHDA) within an intersectional framework. 
Social Science & Medicine, 203, 74-80.
}
\keyword{package}
