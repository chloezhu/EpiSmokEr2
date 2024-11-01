## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message = FALSE----------------------------------------------------------
library(EpiSmokEr2)
data('TestData')

## -----------------------------------------------------------------------------
pred.df<-epismoker(data.m = beta.m, sex.v = sex.v, array = 'EPIC');pred.df

## ----message = FALSE----------------------------------------------------------
rawdata <- loadData(idatPath = '~/Desktop/204957740136')

## ----message = FALSE----------------------------------------------------------
dataset_QN <- normalizeData(RGset=rawdata, normMethod = "QN",array = 'EPIC') 

## -----------------------------------------------------------------------------
pred.df <- epismoker(data.m = dataset_QN, sex.v = c(0,0,1,0,1,1,0,1), array = 'EPIC');pred.df

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

