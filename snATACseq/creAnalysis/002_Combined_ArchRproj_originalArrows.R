## This script is used for creating a combined project from all the arrow files
## for integrated analysis of all samples together.
## Also perform Iterative LSI and clustering. Iterative LSI is required for gene expression integration.
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

addArchRThreads(threads =6)
#Next, I load the reference genome
addArchRGenome("hg38")
Sample <- "comATAC" ##project name of combined archr project

DIR_IN <- "/path/to/output/ATACana/Arrows/"
DIR_OUT <- "/path/to/output/ATACana/comATAC/"
DIR_OTH <- "/path/to/other/files/"

setwd(DIR_OUT)

#Define Arrow files
sams <- readLines(paste0(DIR_OTH,"libids_sel.txt")) #IDs of samples

ArrowFiles <- paste0(sams,".arrow")
ArrowFiles <- paste0(DIR_IN,ArrowFiles)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = Sample,
  copyArrows = TRUE
)

### filter for cells from combined metadata
load(paste0(DIR_OTH,"plotdataATAC.RData")) ##loads plot.data data.frameobject, these are ATAC cells with clonal annotations
cn <- getCellNames(proj)
keep <- cn %in% row.names(plot.data)
table(keep)
cn <- cn[keep]
proj <- proj[cn,]

#saveArchRproj in between

print("Saving ArchR project:")
saveArchRProject(proj)

######################## LSI ################################
#calculate iterative LSI
proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI_int",
                        iterations=5,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4, 0.8),
                          sampleCells = 20000, ##this causes all the cells to be used as all samples are less than 29k
                          n.start = 10
                        ),
                        varFeatures = 100000,
                        dimsToUse = 1:100,
                        totalFeatures = 500000,
                        seed = 1,
                        LSIMethod = 1,
                        scaleDims = FALSE,
                        corCutOff = 0.75,
                        excludeChr = c("chrX", "chrY", "chrMT"),
                        binarize = T,
                        force = T)

#saveArchRproj in between
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample)

### post LSI ------
##iidentifying clusters to perform second QC pass with high resolution=2
print("Adding clusters:")
proj <- addClusters(input = proj,
                    name = "Clusters_int",
                    reducedDims = "IterativeLSI_int",
                    method = "Seurat",
                    force = T,
                    resolution=1,
                    corCutOff = 0.75,
                    scaleDims = FALSE,
                    seed = 1)


# UMAp calculation
proj <- addUMAP(ArchRProj = proj,
                name = "UMAP_int",
                reducedDims = "IterativeLSI_int",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
#saveArchRproj
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample)

print("Finished Combining Arrow Files to one proj")
q()
