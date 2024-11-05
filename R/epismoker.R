#' @name epismoker
#'
#' @aliases epismoker
#'
#' @title Epigenetic Smoking status Estimator
#'
#' @description
#' A function to estimate the smoking status based on the Illumina 450K/EPIC/EPICv2 methylation profiles generated from whole blood.
#'
#' @usage epismoker(data.m, sex.v, array)
#'
#' @param data.m
#' Pre-processed beta value matrix with CpGs as rows and samples as columns. Missing values is imputed as 0.5.
#'
#' @param sex.v
#' A vector indicating the sex, should be marked as 0 and 1 representing female and male respectively. Missing "sex" information is imputed as 0.5.
#'
#' @param array
#' Array type for the input data. "450K" for Illumina HumanMethylation450 data and "EPIC" for Illumina MethylationEPIC (or EPIC v2.0) data.
#'
#' @details
#' This function utilizes the LASSO model, trained on YFS data, to assess the smoking status of each sample in the input data. It computes the probability of each sample being classified as a current smoker, former smoker, or never smoker, and assigns each sample to the category with the highest probability.
#'
#' @return
#' A data frame containing following columns for each sample:
#' \item{result.df}{}
#' \item{CS}{the probability of being current smoker}
#' \item{FS}{the probability of being former smoker}
#' \item{NS}{the probability of being never smoker}
#' \item{Class}{the classification of smoking status (CS, FS, NS)}
#' \item{CvsN}{the classification of smoking status (CS, NS)}
#'
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
#' @author Tianyu Zhu, Teodóra Faragó, Sailalitha Bollepalli
#'
#' @examples
#' #pred.df<-epismoker(beta.m,sex.v,array='EPIC')
#'
#'
#' @export
#'
#'

epismoker<-function(data.m,sex.v,array){

  if (array=='450K') {
    coef.ls<-coef.450k.ls
    sex.v[is.na(sex.v)]<-0.5 # impute NA as 0.5 for 450k model
  }

  if (array=='EPIC') {
    coef.ls<-coef.epic.ls
    sex.v<-ifelse(sex.v==0,1,2) # female=1, male=2 for EPIC model
    sex.v[is.na(sex.v)]<-1.5 # impute NA as 1.5 for EPIC model
  }

  odds.ls<-list()
  common.v <- intersect(rownames(data.m),names(coef.ls[[1]]))
  print(paste0('Number of represented CpGs (max=',length(coef.ls[[1]])-2,')=',length(common.v)))
  if (length(common.v)<460) {print('>10% missing CpGs detected, recommend using QN method')}
  for (i in 1:3){
    rep.beta.m<-data.m[match(names(coef.ls[[i]]),rownames(data.m)),]
    rep.beta.m[1,]<-1 # Intercept term
    rep.beta.m[2,]<-sex.v
    if (sum(is.na(rep.beta.m))!=0){rep.beta.m[is.na(rep.beta.m)]<-0.5} # impute NAs with 0.5
    #if (sum(is.na(rowSums(rep.beta.m)))!=0){rep.beta.m [is.na(rowSums(rep.beta.m)),] <- mean.v[is.na(rowSums(rep.beta.m))]} # impute NAs with population mean
    multi.beta.m<-rep.beta.m * coef.ls[[i]]
    logodds.v<-colSums(multi.beta.m,na.rm = T)
    odds.ls[[i]]<-exp(logodds.v)
  }
  sumodds.v<-odds.ls[[1]]+odds.ls[[2]]+odds.ls[[3]]
  results.df<-data.frame('CS'=odds.ls[[1]]/sumodds.v,
                         'FS'=odds.ls[[2]]/sumodds.v,
                         'NS'=odds.ls[[3]]/sumodds.v)
  results.df$Class<-c('CS','FS','NS')[apply(results.df,1,which.max)]
  results.df$CvsN<-ifelse(results.df$CS>results.df$NS,'CS','NS')
  return(results.df)
}
