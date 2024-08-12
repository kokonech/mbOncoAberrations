##These scripts creates arrow file for each of the sample
## I have used my custome genome annotation derived from Gencode hg38 p13 r37
suppressPackageStartupMessages({
  library(ArchR)
  library(tibble)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

set.seed(456)
addArchRThreads(threads =1) 
#define the Sample ID

args = commandArgs(trailingOnly=TRUE) #sample id

##Dirs
DIR_OUT <- "/path/to/output/ATACana/Arrows/"
DIR_OTH <- "/path/to/other/files/"
DIR_FQ <- "/path/to/cellranger_arc/output/"
#Next, I load the reference genome
addArchRGenome("hg38")
load(paste0(DIR_OTH,"GENCDH38p13r37_ann_4_archr.RData"))

#define the input file
inputFiles <- paste0(DIR_FQ,args[1],"/outs/atac_fragments.tsv.gz") 

#generate ArrowFiles
setwd(DIR_OUT)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = args[1],
  geneAnnotation=new_gen_ann,
  minTSS = 3,  
  minFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = "001_barcode_qc",
  promoterRegion = c(2000,100),
  force = T
)
warning()
print("Generated Arrow File ")


q()

