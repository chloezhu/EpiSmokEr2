#' @name normalizeData
#'
#' @aliases normalizeData
#'
#' @title Color and channel specific quantile normalization
#'
#' @description
#' Performs quantile normaliation on signal intensity values for each color channel and probe type for Illumina 450K or EPIC data.
#'
#' @usage normalizeData(RGset, normMethod=c("QN","ILM","SQN","ALL"), array = c("450K","EPIC"))
#'
#' @param RGset
#' RGChannelSet object from minfi.
#'
#' @param normMethod
#' The method for quantile normalization, can be chosen from "QN", "ILM", "SQN", "ALL". "QN" is quantile normalization for each color channel and probe type with training data as reference. "ILM" and "SQN" are from minfi. "ALL" means using all 3 methods.
#'
#' @param array
#' Array type for the input CpGs. "450K" for Illumina HumanMethylation450 data and "EPIC" for Illumina MethylationEPIC data.
#'
#' @details
#' This function does quantile normalization to input data.
#'
#' @return \item{dataset_QN}{Beta value matrix after quantile normalization. Available if normMethod = "QN" or "ALL".}
#'
#' @return \item{dataset_ILM}{Beta value matrix after ILM normalization. Available if normMethod = "ILM" or "ALL".}
#'
#' @return \item{dataset_SQN}{Beta value matrix after SQN normalization. Available if normMethod = "SQN" or "ALL".}
#'
#' @references
#' Bollepalli S, Korhonen T, Kaprio J, Ander S, Ollikainen M.
#' \emph{EpiSmokEr: a robust classifier to determine smoking status from DNA methylation data.}
#' Epigenomics (2019) 11(13):1469-1486.
#' doi:\href{https://doi.org/10.2217/epi-2019-0206}{
#' 10.2217/epi-2019-0206}.
#'
#' @author Sailalitha Bollepalli, Tianyu Zhu
#'
#' @examples
#' #dataset <- normalizeData(RGset = RGset, normMethod="QN", array = "EPIC")
#'
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylationEPICmanifest
#'
#' @export
#'

normalizeData <- function (RGset=RGset, normMethod=c("QN", "ILM", "SQN", "ALL"), array = c("450K", "EPIC")){

  message("Loading RGset object.")
  RGset <- RGset
  message(sprintf("Dataset has %s samples.",dim(RGset)[2]))

if(normMethod=="QN")
{
  TypeII.Name <- minfi::getProbeInfo(RGset, type = "II")$Name
  TypeII.Green <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
  TypeII.Red <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
  rownames(TypeII.Red) <- TypeII.Name
  colnames(TypeII.Red) <- sampleNames(RGset)
  rownames(TypeII.Green) <- TypeII.Name
  colnames(TypeII.Green) <- sampleNames(RGset)
  TypeI.Green.Name <- minfi::getProbeInfo(RGset, type = "I-Green")$Name
  TypeI.Green.M <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "I-Green")$AddressB,]
  rownames(TypeI.Green.M) <- TypeI.Green.Name
  colnames(TypeI.Green.M) <- sampleNames(RGset)
  TypeI.Green.U <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "I-Green")$AddressA,]
  rownames(TypeI.Green.U) <- TypeI.Green.Name
  colnames(TypeI.Green.U) <- sampleNames(RGset)
  TypeI.Red.Name <- minfi::getProbeInfo(RGset, type = "I-Red")$Name
  TypeI.Red.M <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "I-Red")$AddressB,]
  rownames(TypeI.Red.M) <- TypeI.Red.Name
  colnames(TypeI.Red.M) <- sampleNames(RGset)
  TypeI.Red.U <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "I-Red")$AddressA,]
  rownames(TypeI.Red.U) <- TypeI.Red.Name
  colnames(TypeI.Red.U) <- sampleNames(RGset)

  # Subsetting rownames
  if (array=="450K"){
    probe.ls<-probe.450k.ls
    quantiles.ls<-quantile.450k.ls
  }

  if (array=='EPIC'){
    probe.ls<-probe.epic.ls
    quantiles.ls<-quantile.epic.ls
  }

  int.ls<-list('TypeI.Green.M'=TypeI.Green.M[probe.ls[[1]],],
               'TypeI.Green.U' = TypeI.Green.U[probe.ls[[1]],],
               'TypeI.Red.M' = TypeI.Red.M[probe.ls[[2]],],
               'TypeI.Red.U' = TypeI.Red.U[probe.ls[[2]],],
               'TypeII.Green' = TypeII.Green[probe.ls[[3]],],
               'TypeII.Red' = TypeII.Red[probe.ls[[3]],])

  norm.ls<-list()
  #Quantile normalise TypeII probes using quantiles from Reference data
  for (i in 1:length(int.ls)){
    # get normalized intensity matrix
    rank.m<-apply(int.ls[[i]],2,rank,ties.method='random')
    norm.ls[[i]]<-apply(rank.m,2,function(x){quantiles.ls[[i]][x]})
  }
  QN_beta_TypeII <- norm.ls[[5]]/(norm.ls[[5]]+norm.ls[[6]]+100)
  rownames(QN_beta_TypeII) <- rownames(int.ls[[5]])
  QN_beta_GreenTypeI <- norm.ls[[1]]/(norm.ls[[1]]+norm.ls[[2]]+100)
  rownames(QN_beta_GreenTypeI)<- rownames(int.ls[[1]])
  QN_beta_RedTypeI <- norm.ls[[3]]/(norm.ls[[3]]+norm.ls[[4]]+100)
  rownames(QN_beta_RedTypeI) <- rownames(int.ls[[3]])
  dataset_QN <- rbind(QN_beta_TypeII, QN_beta_GreenTypeI, QN_beta_RedTypeI)

  return(dataset_QN)

}else if(normMethod == "ILM") # Illumina Normalisation
{
  message("Performing Illumina Normalisation.")
  MSet.illumina <- minfi::preprocessIllumina(RGset, bg.correct = FALSE, normalize = "controls")
  dataset_ILM <- minfi::getBeta(MSet.illumina)
  return(dataset_ILM)

}else if(normMethod == "SQN") # Stratified Quantile Normalisation (Touleimat & Tost)
{
  message("Performing Subset Quantile Normalisation.")
  gset.quantile <- minfi::preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                             quantileNormalize = TRUE, stratified = TRUE,mergeManifest = FALSE, sex = NULL)
  dataset_SQN<- minfi::getBeta(gset.quantile)
  return(dataset_SQN)
}else if(normMethod=="ALL")
{
  # QN
  TypeII.Name <- minfi::getProbeInfo(RGset, type = "II")$Name
  TypeII.Green <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
  TypeII.Red <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "II")$AddressA,]
  rownames(TypeII.Red) <- TypeII.Name
  colnames(TypeII.Red) <- sampleNames(RGset)
  rownames(TypeII.Green) <- TypeII.Name
  colnames(TypeII.Green) <- sampleNames(RGset)
  TypeI.Green.Name <- minfi::getProbeInfo(RGset, type = "I-Green")$Name
  TypeI.Green.M <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "I-Green")$AddressB,]
  rownames(TypeI.Green.M) <- TypeI.Green.Name
  colnames(TypeI.Green.M) <- sampleNames(RGset)
  TypeI.Green.U <- minfi::getGreen(RGset)[minfi::getProbeInfo(RGset, type = "I-Green")$AddressA,]
  rownames(TypeI.Green.U) <- TypeI.Green.Name
  colnames(TypeI.Green.U) <- sampleNames(RGset)
  TypeI.Red.Name <- minfi::getProbeInfo(RGset, type = "I-Red")$Name
  TypeI.Red.M <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "I-Red")$AddressB,]
  rownames(TypeI.Red.M) <- TypeI.Red.Name
  colnames(TypeI.Red.M) <- sampleNames(RGset)
  TypeI.Red.U <- minfi::getRed(RGset)[minfi::getProbeInfo(RGset, type = "I-Red")$AddressA,]
  rownames(TypeI.Red.U) <- TypeI.Red.Name
  colnames(TypeI.Red.U) <- sampleNames(RGset)
  # Subsetting rownames
  if (array=="450K"){
    probe.ls<-probe.450k.ls
    quantile.ls<-quantile.450k.ls}

  if (array=='EPIC'){
    probe.ls<-probe.epic.ls
    quantile.ls<-quantile.epic.ls}

  int.ls<-list('TypeI.Green.M'=TypeI.Green.M[probe.ls[[1]],],
               'TypeI.Green.U' = TypeI.Green.U[probe.ls[[1]],],
               'TypeI.Red.M' = TypeI.Red.M[probe.ls[[2]],],
               'TypeI.Red.U' = TypeI.Red.U[probe.ls[[2]],],
               'TypeII.Green' = TypeII.Green[probe.ls[[3]],],
               'TypeII.Red' = TypeII.Red[probe.ls[[3]],])

  norm.ls<-list()
  #Quantile normalise TypeII probes using quantiles from Reference data
  for (i in 1:length(int.ls)){
    # get normalized intensity matrix
    rank.m<-apply(int.ls[[i]],2,rank,ties.method='random')
    norm.ls[[i]]<-apply(rank.m,2,function(x){quantiles.ls[[i]][x]})
  }
  QN_beta_TypeII <- norm.ls[[5]]/(norm.ls[[5]]+norm.ls[[6]]+100)
  rownames(QN_beta_TypeII) <- rownames(int.ls[[5]])
  QN_beta_GreenTypeI <- norm.ls[[1]]/(norm.ls[[1]]+norm.ls[[2]]+100)
  rownames(QN_beta_GreenTypeI)<- rownames(int.ls[[1]])
  QN_beta_RedTypeI <- norm.ls[[3]]/(norm.ls[[3]]+norm.ls[[4]]+100)
  rownames(QN_beta_RedTypeI) <- rownames(int.ls[[3]])
  dataset_QN <- rbind(QN_beta_TypeII, QN_beta_GreenTypeI, QN_beta_RedTypeI)

  ###################################################################################################
  # SQN
  message("Performing Subset Quantile Normalisation.")
  gset.quantile <- minfi::preprocessQuantile(RGset,quantileNormalize = TRUE, stratified = TRUE,mergeManifest = FALSE, sex = NULL, verbose=TRUE)
  dataset_SQN<- minfi::getBeta(gset.quantile)
  ###################################################################################################
  # ILM
  message("Performing Illumina Normalisation.")
  MSet.illumina <- minfi::preprocessIllumina(RGset, bg.correct = FALSE, normalize = "controls")
  dataset_ILM <- minfi::getBeta(MSet.illumina)
  return(list(dataset_QN=dataset_QN, dataset_ILM=dataset_ILM, dataset_SQN=dataset_SQN))

}

}


