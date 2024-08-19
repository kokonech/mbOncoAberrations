library(Seurat)
library(Signac)

# This script plots InferCnv dervided results in UMAP per sample 
# also it supports meta-cells back assignment
# It depends on the location of the results 

## INPUT DATA

# InferCNV
cnvDir = "/PATH/TO/cnvProfiling/globalInferCNV/"

# snRNA per sample 
snDir = "/PATH/TO/perSampleV2/"
resName = paste0(sId,"_meta_clusters_n3_hg38")

# snATAC per sample
#snDir = "/PATH/TO/atacPerSample/"
#resName = paste0(sId,"_ATAC_meta_clusters_n3")

# input ID
args =  commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
}


# analysis

cnvPath = paste0(cnvDir, resName,"/infercnv.observation_groupings.txt")
cnvInfo <- read.delim(cnvPath,sep=" ")
cnvInfo$metacell_mapping = rownames(cnvInfo)

mb <- readRDS(paste0(snDir,sId,"_obj.RDS"))

# snRNA
metaInfo <- read.delim(paste0(snDir,sId,"_metacells_mapping.tsv"))
# snATAC
#metaInfo <- read.delim(paste0(snDir,sId,"_ATAC_metacells_mapping.tsv"))


metaInfo$cellId = rownames(metaInfo)


mergedDf <- merge(metaInfo, cnvInfo, by = "metacell_mapping")
rownames(mergedDf) <- mergedDf$cellId

mb$cnvCluster = "Unclear"
mb@meta.data[ rownames(mergedDf), "cnvCluster" ] = paste0("cl",mergedDf$Dendrogram.Group)

pdf(paste0(cnvDir,resName,"/",resName,".UMAP.pdf"),width=8,height=6)
DimPlot(mb, reduction = "umap", label = F,group.by="cnvCluster")
dev.off()


