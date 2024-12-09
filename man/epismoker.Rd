% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epismoker.R
\name{epismoker}
\alias{epismoker}
\title{Epigenetic Smoking status Estimator}
\usage{
epismoker(data.m, sex.v, array)
}
\arguments{
\item{data.m}{Pre-processed beta value matrix with CpGs as rows and samples as columns. Missing values is imputed as 0.5.}

\item{sex.v}{A vector indicating the sex, should be marked as 0 and 1 representing female and male respectively. Missing "sex" information is imputed as 0.5.}

\item{array}{Array type for the input data. "450K" for Illumina HumanMethylation450 data and "EPIC" for Illumina MethylationEPIC (or EPIC v2.0) data.}
}
\value{
A data frame containing following columns for each sample:
\item{result.df}{}
\item{CS}{the probability of being current smoker}
\item{FS}{the probability of being former smoker}
\item{NS}{the probability of being never smoker}
\item{Class}{the classification of smoking status (CS, FS, NS)}
\item{CvsN}{the classification of smoking status (CS, NS)}
}
\description{
A function to estimate the smoking status based on the Illumina 450K/EPIC/EPICv2 methylation profiles generated from whole blood.
}
\details{
This function utilizes the LASSO model, trained on YFS data, to assess the smoking status of each sample in the input data. It computes the probability of each sample being classified as a current smoker, former smoker, or never smoker, and assigns each sample to the category with the highest probability.
}
\examples{
#pred.df<-epismoker(beta.m,sex.v,array='EPIC')


}
\references{
Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
\emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
Epigenomics (2019) 11(13):1469-1486.
doi:\href{https://doi.org/10.2217/epi-2019-0206}{
10.2217/epi-2019-0206}.
}
\author{
Tianyu Zhu, Teodóra Faragó, Sailalitha Bollepalli
}
