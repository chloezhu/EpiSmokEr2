#' @title
#' Reads IDAT files
#'
#' @description
#' loadData function reads IDAT files and saves as a minfi RGSet object. Details can be found in the help page of *read.metharray.sheet* and *read.metharray.exp* function from *minfi* R package.
#'
#' @param idatPath
#' Requires a path to the folder with the IDAT files and sample sheet included.
#'
#' @return
#' A RGChannelSet object
#'
#' @examples
#' #rawdata <- loadData(idatPath)
#'
#'
#' @import minfi
#' @import IlluminaHumanMethylation450kmanifest
#' @import IlluminaHumanMethylationEPICmanifest
#' @export
#'
loadData <- function (idatPath){
  message("Loading idat files.")
  baseDir <- file.path(idatPath)
  targets <- minfi::read.metharray.sheet(baseDir)
  RGset <- minfi::read.metharray.exp(base = baseDir, recursive = TRUE)
  message(sprintf("Dataset has %s samples.",dim(RGset)[2]))
  return(RGset)
}
