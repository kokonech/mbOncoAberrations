## R 4.2.2
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(rprojroot)

resDir <- Sys.getenv("ONCO_AB_RESDIR")
if (nchar(resDir) > 0) {
    setwd(resDir)
} 
print(paste("Result dir:",resDir))

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Error! Input data is is not provided. Required format : <SID> <PATH_TO_INPUT>", call.=FALSE)
} else {
  sId = args[1]
  inPath = args[2]
}

print("Loading data...")

inData <- Read10X(inPath)
annData <- NULL

TARG="cnvBlock" # could be adjusted based on annotation
targColumn = "seurat_clusters" # default

if (length(args) > 2) {
  annData <- read.delim(args[3])
  if (! (TARG %in% colnames(annData) ) ) {
    stop(paste("Error! Required annotation column is not available:", TARG))
  }
  print(paste("Using target annoitation column:",TARG))
  targColumn = TARG
}

## process


print("Processing data...")

mb <- CreateSeuratObject(inData, project=sId, min.cells = 5, min.features = 300, assay = "RNA")
print("Loaded input")
print(mb)

mb  <- PercentageFeatureSet(mb, pattern = "^MT-", col.name = "percent.mt")

if (!is.null(annData)) {
  #print(head(annData))
  cIds = intersect(rownames(annData), rownames(mb@meta.data))
  if (length(cIds) == 0) {
      print("Error: no overlapping cells IDs in annotation")
  }
  print("Adjust for custom annoitation")
  mb <- mb[,cIds]
  mb@meta.data <- cbind(mb@meta.data, annData[cIds,targColumn])
  colnames(mb@meta.data)[ncol(mb@meta.data)] <- targColumn
} else {
  print("QC filtering")
  # minimum QC filtering, assuming 10X v3 protocol
  geneLim <- 4000
  transLim <- 8000
  mb <- subset(mb, subset = nFeature_RNA < geneLim  & nCount_RNA < transLim )
}

#head(mb@meta.data)
print("After adjustment...")
print(mb)

# standard processing

# default n feature = 2500
mb <- SCTransform(mb, vars.to.regress = "percent.mt", 
                    variable.features.n = 2500, verbose = FALSE)


# default = 20
ndim = 20
mb <- RunPCA(mb, verbose = FALSE)
mb <- FindNeighbors(mb, reduction = "pca", dims = 1:ndim)
# default 0.3
mb <- FindClusters(mb, resolution = 0.3)
mb <- RunUMAP(mb, reduction = "pca", dims = 1:ndim)

pdf(paste0(sId,".clusters_UMAP.pdf"), width=8, height = 6)
DimPlot(mb, reduction = "umap", label = TRUE)
if (nchar(targColumn) > 0) {
  DimPlot(mb, reduction = "umap", label = TRUE,group.by = targColumn)
}
dev.off()

markers <- FindAllMarkers(mb, only.pos = TRUE)
fName = paste0(sId,".DEG_table.txt")
write.table(markers, fName, sep="\t", quote=F)

# check target genes expression

spGenes <- c("MYC", "MYCN", "PRDM6", "SNCAIP")
spGenes <- spGenes[ spGenes %in% rownames(mb)]

pdf(paste0(sId,".RB_targets.pdf"), width=6, height = 4)
for (gName in spGenes) {
    print(gName)
    print(FeaturePlot(mb, features = gName, reduction = "umap",raster=T))
}

dev.off()

# save object
saveRDS(mb, paste0(sId,"_obj.RDS" ))

# save result for InferCNV
gz1 <- gzfile(paste0(sId,"_raw_counts.txt.gz"), "w")
rawCounts <- as.matrix(mb@assays$RNA@counts)
write.table(rawCounts,gz1,sep="\t",quote=F)
close(gz1)

# based on the cluster
annTable <- mb@meta.data
# minor fix to avoid mix numbers vs factor
if (nchar(targColumn) == 0) {  
  annTable2 <- annTable[,"seurat_clusters",drop=F]
  annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
} else {
  annTable2 <- annTable[,targColumn,drop=F]
}
write.table(annTable2, paste0(sId,"_cnv_ann.txt") ,col.names = F,sep="\t",quote=F)



