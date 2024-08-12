## This script used RNA-Seq data to identify Peak to Gene links

#P2G links ArchR

suppressPackageStartupMessages({
  library(ArchR)
 library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(scran)
  library(scater)
  library(Matrix)
})

addArchRThreads(threads =8)
addArchRGenome("hg38")

DIR_OUT <- "/path/to/output/ATACana/"
DIR_OTH <- "/path/to/other/files/"

setwd(DIR_OUT)

set.seed(456)
### To add peak to gene links first I need to add RNA expression data
## For cells with true-multiomic data

#Load ArchR project
proj <- loadArchRProject("comATAC/")
cells <- getCellNames(proj)
##load Gene expression matrix and subset for valid ATAC cells
load(paste0(DIR_OTH,"plotdata_G34MB.RData")) ## this loads plot.data object, this is metadata from RNA cells
load(paste0(DIR_OTH,"lgcountsG34MB.RData")) ## this loads mdt dgCMatrix, this is log-noralized RNA counts 
##converting RNA cells to ATAC cells id
plot.data$RNA_cells <- rna_cells <- row.names(plot.data)
rna_cells <- paste0(gsub("_","#",rna_cells),"-1")

keep <- rna_cells %in% cells
table(keep)
plot.data <- plot.data[keep,]
mdt <- mdt[,row.names(plot.data)]

colnames(mdt) <-row.names(plot.data) <- rna_cells[keep]

## subset ArchR project
subsetArchRProject(ArchRProj = proj,
                   cells = row.names(plot.data),
                   outputDirectory = "comATAC2")
rm(proj)
## load
proj <- loadArchRProject("comATAC2/")

#save the combined matrix as SE
gr_obj <- rtracklayer::import(paste0(DIR_OTH,"gencdH38p13r37CR_genesN.filtered.bed"))
SE <- SummarizedExperiment(assays = list(logcounts=mdt), rowRanges = gr_obj)
rm(mdt)

## add gene expression data
proj <- addGeneExpressionMatrix(input = proj,
                                seRNA = SE)


#add Peak tokeep#add Peak to gene linkskeep#add Peak to gene links
proj<- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "IterativeLSI_int",
  useMatrix = "GeneExpressionMatrix"
)

#save
saveArchRProject(proj, outputDirectory = "comATAC2/")

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE #FALSE, df is returned
)

p2g_granges <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE #TRUE, granges returned
)

write.table(p2g,"comATAC/Peak2Gene_df.txt",col.names = TRUE,sep="\t")
write.table(p2g_granges,"comATAC/Peak2Gene_Granges.bed",col.names = T,sep = "\t")

#### use from github - add peak and gene names
p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "_", end(.))})[p2geneDF$idxATAC]
write.table(p2geneDF,"comATAC/Peak2Gene_df_withNames.txt",sep = "\t",col.names = TRUE,row.names = F)







