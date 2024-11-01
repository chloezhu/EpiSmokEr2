#' Coefficients of 450K model
#'
#' Coefficents of EpiSmokEr 450K model from Bollepalli et al., 2019.
#'
#' \itemize{
#'   \item coef.450k.ls: Coefficients of intercept, sex term, and 121 CpGs for current smoker (cs), former smoker (fs), and never smoker (ns).
#' }
#'
#' @docType data
#' @keywords coefficients
#' @name coef_450k
#' @usage data("coef_450k")
#' @format A list of 3 vectors, each contains 123 coefficients
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL

#' Coefficients of EPIC model
#'
#' Coefficents of EpiSmokEr2 EPIC model
#'
#' \itemize{
#'   \item coef.epic.ls: Coefficients of intercept, sex term, and 511 CpGs for current smoker (cs), former smoker (fs), and never smoker (ns).
#' }
#'
#' @docType data
#' @keywords coefficients
#' @name coef_EPIC
#' @usage data("coef_EPIC")
#' @format A list of 3 vectors, each contains 513 coefficients
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL

#' Quantiles of 450K data
#'
#' Quantiles of intensity values based on each probe type and color channel from training data (DILGOM)
#'
#' \itemize{
#'   \item quantile.450k.ls: Quantiles of intensity values based on each probe type and color channel
#' }
#'
#' @docType data
#' @keywords quantiles
#' @name quan_450k
#' @usage data("quan_450k")
#' @format A list of 6 vectors
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL

#' Quantiles of EPIC data
#'
#' Quantiles of intensity values based on each probe type and color channel from training data (YFS)
#'
#' \itemize{
#'   \item quantile.epic.ls: Quantiles of intensity values based on each probe type and color channel
#' }
#'
#' @docType data
#' @keywords quantiles
#' @name quan_epic
#' @usage data("quan_epic")
#' @format A list of 6 vectors
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL

#' Probe list for 450K data
#'
#' Type I (two color channel) and Type II probe names of 450K training data (DILGOM)
#'
#' \itemize{
#'   \item probe.450k.ls: Type I (two color channel) and Type II probe names of 450K training data (DILGOM)
#' }
#'
#' @docType data
#' @keywords probes
#' @name probe_450k
#' @usage data("probe_450k")
#' @format A list of 3 vectors
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL

#' Probe list for EPIC data
#'
#' Type I (two color channel) and Type II probe names of EPIC training data (YFS)
#'
#' \itemize{
#'   \item probe.epic.ls: Type I (two color channel) and Type II probe names of EPIC training data (YFS)
#' }
#'
#' @docType data
#' @keywords probes
#' @name probe_epic
#' @usage data("probe_epic")
#' @format A list of 3 vectors
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL

#' Test data
#'
#' Test EPIC beta value matrix of 511 CpGs and 10 samples and the sex information.
#'
#' \itemize{
#'   \item beta.m: Simulated beta value matrix of 511 CpGs and 10 samples
#'   \item sex.v: Sex information (0=Female, 1=Male)
#' }
#'
#' @docType data
#' @keywords testdata
#' @name TestData
#' @usage data("TestData")
#' @format A matrix and a vector
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
NULL
