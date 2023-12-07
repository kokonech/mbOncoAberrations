##### This script provides the input as required for pyABC
rm(list=ls())

#source("Settings.R")
rdata.directory <- "~/pyABC/Medulloblastoma/"
##### load functions

library(NBevolution)

##### observed data

load(paste0(rdata.directory, "Input_data_MB_initiation.RData"))

mySumStatData <- list(P.MRCA=P.MRCA, P.ECA=P.ECA)


model <- "contraction"

source("~/pyABC/Medulloblastoma/ABC_fit.R")

