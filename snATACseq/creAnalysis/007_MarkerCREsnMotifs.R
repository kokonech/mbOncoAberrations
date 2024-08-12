## Marker peak analysis per sample by clonal clusters
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
  library(scran)
  library(scater)
})

args = commandArgs(trailingOnly=TRUE) #sample id

addArchRThreads(threads =6)
DIR_IN <- "/path/to/output/ATACana/"
DIR_OUT <- "/path/to/output/ATACana/CREres/"

DIR_OTH <- "/path/to/other/files/"

setwd(DIR_IN)
### load and subset ArchR project
proj_com <- loadArchRProject("comATAC/")

ann <- proj_com$Clonal_clusters
k <- grep("_Normal",ann)
annx <- ann[k]
ann[k] <- "Normal"
proj_com$Clonal_clusters <- ann

## subsetting to sample but using all normals
keep1 <- proj_com$Sample==args
keep2 <- ann=="Normal"
table(keep1,keep2)
keep <- keep1 | keep2
table(keep)
proj_com <- proj_com[keep,]
table(proj_com$Clonal_clusters)

###-----
# ## Marker peaks: ArchR implementation 
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_com, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clonal_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
peaks <-data.frame(rowData(markersPeaks))
pn <- row.names(peaks) <- paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)
AUC <- assay(markersPeaks,"AUC")
FDR <- assay(markersPeaks,"FDR")
row.names(AUC) <- row.names(FDR) <- row.names(peaks)

##This gets the same peaks as obtained by plotMarkerHeatmap
mps <- c()
clones <- colnames(AUC)
for(x in clones){
  keep1 <- AUC[,x]>=0.52
  keep2 <- FDR[,x]<0.01
  del <- pn[keep1 & keep2]
  rg <- peaks[del,c("seqnames","start","end")]
  write.table(rg, file = paste0("CREres/",args,"/",x,"_CREs.bed"),
              row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
  rm(rg)
  rx <- rep(x,length(del))
  del <- cbind(del,rx)
  mps <- rbind(mps,del)
  rm(keep1,keep2,del,rx)
}
mps <- data.frame(mps)
colnames(mps) <- c("Peaks","Cluster")
write.table(mps, file = paste0("CREres/",args,"/MarkerPeaks.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)
q()
