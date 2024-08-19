library(Seurat)
library(dplyr)

# NOTE: this code depends on server envirnoment with data locations

resDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/snCombRes/"

# text file with prperties of samples
sInfo <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/MB_MYC_MYCN_PRDM6_cohort.210723.txt")
sIds <- sInfo$Sample
rownames(sInfo) <- sIds


# requires up to 160 GB RAM
# reading from per sample result
grpName = "MB_G34_somatic"
ob.list <- list()
for (sId in sIds) {
    print(sId)
    ob.list[[sId]] <- readRDS(paste0("perSampleV2/",sId,"_obj.RDS"))
}


mb.merged <- merge(ob.list[[1]],
                y = ob.list[2: length(ob.list)],
                project = grpName)


# no batch effect adjustment

DefaultAssay(mb.merged) <- "RNA"
mb.merged[['SCT']] <- NULL

ob.list <- NULL
gc()

summary(rownames(sInfo) %in% mb.merged$Sample)
require(tidyverse)
mb.merged$subgroup <- map_chr(mb.merged$Sample, function(x) sInfo[x,"Subgroup"] )
mb.merged$group <- map_chr(mb.merged$Sample, function(x) sInfo[x,"Group"] )


# INPUT GENERATION FROM BELOW IS REQUIRED !!!!


# analysis

mb.merged <- NormalizeData(mb.merged)
all.genes <- rownames(mb.merged)
mb.merged <- ScaleData( mb.merged, all.genes )
n = 2500 # default
mb.merged <- FindVariableFeatures(mb.merged, selection.method = "vst", nfeatures=n)

mb.merged <- RunPCA(mb.merged, features = VariableFeatures(object = mb.merged))

png(paste0(resDir,grpName,"_PcElbow.png"), width=600, height=400)
ElbowPlot(mb.merged, ndims = 20)
dev.off()
print("Find clusters...")

ndim = 20
mb.merged <- FindNeighbors(mb.merged, dims = 1:ndim)
mb.merged <- FindClusters(mb.merged, resolution = 0.5)
mb.merged <- RunUMAP(mb.merged, reduction = "pca", dims = 1:ndim)


colInfo <- read.csv("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/subtypeG34.csv")
colorSubgroups <- colInfo$Color
names(colorSubgroups) <- colInfo$Subtype
colorGroups <- c("G3"="goldenrod1", "G4"="darkgreen")


pdf(paste0(resDir,grpName,"_UMAP.clusters.pdf"),width = 8, height = 6)
DimPlot(object = mb.merged, pt.size=1, label=T,reduction = 'umap')
DimPlot(mb.merged, reduction = "umap",group.by = "group",
       cols = colorGroups)
DimPlot(object = mb.merged, label=F,reduction = 'umap', group.by="subgroup", 
            cols = colorSubgroups)
DimPlot(object = mb.merged, label=T,reduction = 'umap', group.by="Sample")
dev.off()




gInfo <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/spatial_gene_selection_MB.160322.txt")

gNames <- gInfo$Gene
gNames <- gNames[gNames %in% rownames(mb.merged)]


pdf(paste0(resDir, grpName,".UMAP.gene_expr.pdf"), width=8, height = 6)

for (gName in gNames) {
    print(gName)
    print(FeaturePlot(mb.merged, features = gName, reduction = "umap"))
}
dev.off()

refCellTypes = names(summary(as.factor(mb.merged@meta.data$RF_finlab)))

umapCoords <- mb.merged@reductions$umap@cell.embeddings
maxX <- max(umapCoords[,1])
minX <- min(umapCoords[,1])
maxY <- max(umapCoords[,2])
minY <- min(umapCoords[,2])
require(ggplot2)

pdf(paste0(resDir,grpName,"_UMAP.RF_finlab.pdf"),width = 8, height = 6)
for (ct in refCellTypes) {
    print(ct)
    targCells  = Cells(mb.merged)[mb.merged@meta.data$RF_finlab == ct]
    if (length(targCells) == 1) {
        next
    }
    print(DimPlot(mb.merged, reduction = "umap",group.by = "RF_finlab", cells=targCells) +xlim(minX,maxX)+ylim(minY,maxY) )
}
dev.off()

saveRDS(mb.merged, paste0(resDir,grpName,"_object.Rdata"))

