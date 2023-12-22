# R 4.2.2
setwd("/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial")
library(Seurat)
library(dplyr)
library(RImageJROI)

# main cohort details : sample ids, image sizes
fullInfo <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/MB_RB_runs_info.280623.txt")

args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
  sId <- args[2]
  runId <- args[1]
}


sInfo <- fullInfo[ fullInfo$Cohort == runId, ]
rownames(sInfo) <- sInfo$RB_ID


resId <- sInfo[sId,]$TUM_ID
iHeight <- sInfo[sId,]$Height
print(paste(runId,sId,resId,iHeight))

countsFile <- paste0("cellposeRes/",runId,"/",sId,"_counts.csv")
resDir = "cellposeRes/seuratAnalysisR3v2/" # full cell borders





print(paste("Input counts data:",countsFile))
cellData <- read.csv(countsFile)


rownames(cellData) <- paste0("Cell",cellData$cell)

print("Prepare expression...")
exprData <- t(cellData[,6:ncol(cellData)])

# filter the FP values
exprData <- exprData[!startsWith(rownames(exprData),"FP."),]


# centroids
if (FALSE) {
print("Prepare centroids...")
cellCoords <- cellData[,c("centroid.y", "centroid.x")]
colnames(cellCoords) <- c("x","y")
cellCoords$cell <- rownames(cellCoords)
cellCoords$x <- iHeight - cellCoords$x 
}

# full coords
#if (FALSE) {
require(data.table)
print("Prepare full coords")
roiPath <- paste0("cellposeRes/",runId,"/",sId,"_ROIs.zip")
cell_roi <- read.ijzip(roiPath, names = TRUE, list.files = FALSE, verbose = FALSE)
cell_roi <- lapply(cell_roi, `[[`, "coords")
cell_roi <- lapply(cell_roi, as.data.table)
cell_roi <- rbindlist(cell_roi, idcol = "cell")

cell_roi[, x2 := y] # reorienting the coordinates
cell_roi[, y := x]
cell_roi[, x := iHeight - x2]
cell_roi[, x2 := NULL]
cell_roi <- as.data.frame(cell_roi)
# fix cell names
cell_roi$cell <- paste0("Cell", cell_roi$cell)
#}

print("Running analysis...")
mb <- CreateSeuratObject(exprData,project=sId, assay = "Spatial")
mb@meta.data$sample <- resId
mb@meta.data$nArea_Spatial <- cellData$area 
print(mb)

mb@images <- list()

if (FALSE) {
mb@images["cen"] <- CreateFOV(
  coords = cellCoords,
  type = "centroids",
  nsides = 0L,
  radius = 1L,
  theta = 0L,
  molecules = NULL, 
  assay = "Spatial",
  key = NULL,
  name = NULL)
}

#if (FALSE) {
mb@images["seg"] <- CreateFOV(
  coords = cell_roi,
  type = "segmentation",
  #nsides = 0L,
 #radius = 1L,
  #theta = 0L,
  #molecules = NULL,
  assay = "Spatial",
  #key = NULL,
  #name = NULL)
)
#}

pdf(paste0(resDir, resId,"_",sId,".QC_VlnPlot.pdf"), width=6, height = 4)
VlnPlot(mb, features = c("nCount_Spatial","nFeature_Spatial","nArea_Spatial"), pt.size = 0.1)
dev.off()

# collect and check stats
#quit(save = "no",status=0)

mb <- subset(mb, subset = nCount_Spatial > 5) # filter control, default was 4
mb <- subset(mb, subset = nCount_Spatial < 300)
mb <- subset(mb, subset = nArea_Spatial < 10000) 



pdf(paste0(resDir, resId,"_",sId,".QC_VlnPlot.after_filtering.pdf"), width=6, height = 4)
VlnPlot(mb, features = c("nCount_Spatial","nFeature_Spatial", "nArea_Spatial"), pt.size = 0.1)
dev.off()

mb <- SCTransform(mb, assay = "Spatial", verbose = FALSE)

mb <- RunPCA(mb, assay = "SCT", verbose = FALSE,npcs = 30, approx=FALSE)
#ElbowPlot(mb)
ndim = 20
mb <- FindNeighbors(mb, reduction = "pca", dims = 1:ndim)
mb <- FindClusters(mb, resolution = 0.3) # default 0.25
mb <- RunUMAP(mb, reduction = "pca", dims = 1:ndim)


pdf(paste0(resDir, resId,"_",sId,".UMAP.pdf"), width=6, height = 4)
DimPlot(mb, reduction = "umap", label = TRUE)
dev.off()

require(scales)
cTypes <- levels(mb@active.ident)
selCols = hue_pal()(length(cTypes))


fName <- paste0(resDir, resId, "_",sId,".spatial_clusters.pdf")
pdf(fName,width=8,height = 8)


# K.O. standard
print(ImageDimPlot(mb, cols=selCols, border.color= NA))

# K.O. possible to remove it, but result looks a bit weird 
#print(ImageDimPlot(mb,cols = "parade", border.size=0, border.color=NA, fov = "seg"))

for (i in 1:length(cTypes)) {
  cType = cTypes[i]
  print(cType)
  print(ImageDimPlot(mb,  cols = selCols[i] , border.color = NA,
                     cells = WhichCells(mb, idents = cType)))
}

dev.off()



gNames <- rownames(mb)

pdf(paste0(resDir, resId,"_",sId,".UMAP.gene_info.pdf"), width=6, height = 4)
for (gName in gNames) {
  print(gName)
  print(FeaturePlot(mb, features = gName, reduction = "umap"))
}
dev.off()


markers <- FindAllMarkers(mb, only.pos = TRUE)
fName = paste0(resDir, resId,"_",sId,".DEG_table.txt")
write.table(markers, fName, sep="\t", quote=F)

markers2 <- markers %>% group_by(cluster) 

pdf(paste0(resDir, resId,"_",sId,".heatmap.pdf"),width=16,height = 14)
DoHeatmap(mb,features = markers2$gene)
dev.off()

markers3 <- markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)

subDir = paste0(resDir, resId, "_",sId,"_spatial_DEGs/")
if (!file.exists(subDir)){
  dir.create(subDir)
}

# NOTE: control the visualization

#for (gName in markers3$gene) {
for (gName in gNames) {
  fName <- paste0(subDir, resId, "_",sId,".spatial_", gName,".png")
  png(fName,width=800,height = 600)
  print(ImageFeaturePlot(mb,features = gName, #size=1 )) # for circles
                            border.color = NA ))  # for precise borders
  dev.off()
}

saveRDS(mb,file=paste0(resDir, resId, "_",sId, ".result.rds"))
# mb <- readRDS(paste0(resDir, resId, "_",sId, ".result.rds"))


print("Done!")




