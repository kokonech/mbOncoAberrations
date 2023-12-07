###### Input data to fit G3/G4 MB initiation
library(openxlsx)
library(moments)

source("./RScripts//Settings.R")
load(paste0(rdata.directory, "MRCA_timing.RData"))


g3g4.tumors <- rownames(sample.information[sample.information$mnp11 %in% c("MB, G3", "MB, G4"),])

## collect the cumulative distribution of MRCA densities
P.MRCA = data.frame(Density=mutation.time.mrca[g3g4.tumors,]$Mean)
rownames(P.MRCA) <- g3g4.tumors
P.MRCA <- P.MRCA[order(P.MRCA$Density),,drop=F]
P.MRCA <- P.MRCA[!is.na(P.MRCA$Density),,drop=F]
P.MRCA$P <- seq(1, nrow(P.MRCA))/nrow(P.MRCA)
## lower and upper bounds
P.MRCA$P.upper = sapply(P.MRCA$Density, function(x){
  sum(mutation.time.mrca[g3g4.tumors,]$Min <= x, na.rm = T)
})/nrow(P.MRCA)
P.MRCA$P.lower = sapply(P.MRCA$Density, function(x){
  sum(mutation.time.mrca[g3g4.tumors,]$Max <= x, na.rm = T)
})/nrow(P.MRCA)

## collect the cumulative distribution of ECA densities
P.ECA = data.frame(Density=mutation.time.eca[g3g4.tumors,]$Mean)
rownames(P.ECA) <- g3g4.tumors
P.ECA <- P.ECA[order(P.ECA$Density),,drop=F]
P.ECA <- P.ECA[!is.na(P.ECA$Density),,drop=F]
P.ECA$P <- seq(1, nrow(P.ECA))/nrow(P.ECA)
## lower and upper bounds
P.ECA$P.upper = sapply(P.ECA$Density, function(x){
  sum(mutation.time.eca[rownames(P.ECA),]$Min <= x, na.rm = T)
})/nrow(P.ECA)
P.ECA$P.lower = sapply(P.ECA$Density, function(x){
  sum(mutation.time.eca[rownames(P.ECA),]$Max <= x, na.rm = T)
})/nrow(P.ECA)


save(P.MRCA, P.ECA, 
     file=paste0(rdata.directory, "Input_data_MB_initiation.RData"))

