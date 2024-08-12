## This script filter cells to cell present in RNA and add annotation by KO
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

addArchRThreads(threads =6)
#Next, I load the reference genome
addArchRGenome("hg38")

DIR_OUT <- "/path/to/output/ATACana/"
DIR_OTH <- "/path/to/other/files/"

setwd(DIR_OUT)
##load subsetted ArchR project------- 
proj <- loadArchRProject("comATAC")
cells <- getCellNames(proj)
plot.data <- data.frame(getCellColData(proj))
load(paste0(DIR_OTH,"clonalann.RData")) ##ann data.frame object, clonal annotation data
ann <- ann[cells,]

proj$Clonal_clusters <- 
  plot.data$Clonal_clusters <-ann$V2

#Creating PseudoBulk Replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  sampleRatio = 0.8, 
  returnGroups = F,
  force = T,
  groupBy = "Sample" , 
  minCells = 2000,
  maxCells = 5000,
  minReplicates = 2,
  maxReplicates = 5,
  maxFragments = 100 * 10^6,
  useLabels= TRUE ## groups are already defined per sample
)
print("Saving project after group coverage!")
saveArchRProject(proj)
###### Peak Calling ##########

#Define path to MACS2
pathToMacs2 <- findMacs2()

#Peak Calling
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Sample",
  maxPeaks = 150000,
  pathToMacs2 = pathToMacs2,
  reproducibility = "1", ##here there only 2 replicates per sample
  excludeChr = c("chrY", "chrMT"), 
  method = "q",
  cutOff = 0.01, 
  extendSummits = 250,
  force = T
)

#Adding Peak matrix
proj <- addPeakMatrix(proj,
                          ceiling = 5,
                          binarize = F,
                          force = T)
print("Saving project after adding peak matrix!")

saveArchRProject(ArchRProj = proj, outputDirectory = "comATAC")
q()



