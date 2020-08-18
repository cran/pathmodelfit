# Package documentation

#' @aliases pathmodelfit-package
#' @import lavaan
#' @references
#' Hancock, G. R., & Mueller, R. O. (2011). The reliability paradox in assessing structural relations within covariance structure models. Educational and Psychological Measurement, 71(2), 306-324.
#'
#'McNeish, D., & Hancock, G. R. (2018). The effect of measurement quality on targeted structural model fit indices: A comment on Lance, Beck, Fan, and Carter (2016). Psychological Methods, 23(1), 184–190. https://doi.org/10.1037/met0000157
#'
#' O'Boyle, E. H., Jr., & Williams, L. J. (2011). Decomposing model fit: Measurement vs. theory in organizational research using latent variables. Journal of Applied Psychology, 96(1), 1–12. https://doi.org/10.1037/a0020539
#'
#' Williams, L. J., & O’Boyle, E. H. (2011). The myth of global fit indices and alternatives for assessing latent variable relations. Organizational Research Methods, 14, 350-369.
#'
#' Williams, L. J., O’Boyle, E. H., & Yu, J. (2020). Condition 9 and 10 tests of model confirmation: A review of James, Mulaik, and Brett (1982) and contemporary alternatives. Organizational Research Methods, 23, 1, 6-29.
#'
#' @examples
#'
#' library(lavaan)
#'
#' model4 <- '
#' Ldrrew =~ LdrrewI1 + LdrrewI2 + LdrrewI3
#' Jobcom =~ JobcomI1 + JobcomI2 + JobcomI3
#' Jobsat =~ JobsatI1 + JobsatI2 + JobsatI3
#' Orgcom =~ OrgcomI1 + OrgcomI2 + OrgcomI3
#' Jobsat ~ Ldrrew + Jobcom
#' Orgcom ~ Jobsat'
#'
#' data(mediationVC)
#'
#' fit <- sem(model4, sample.cov = mediationVC, sample.nobs = 232)
#' pathmodelfit(fit)
#'
"_PACKAGE"
