library(Seurat)
library(dplyr)
setwd("/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/cellposeRes")

resDir = "combAnalysis/"

# used for supplementary table
dataDir = "seuratAnalysisR3v2/"


extractPseudobulkMatrix <- function(rawTable, annData, clId ) {
    resTable <- NULL
    clusters <- sort(unique(as.character(annData[,clId] )) )
    print(clusters)
    targCells <- as.character(annData[,clId])
    for (k in clusters) {
        cNames <- rownames(annData)[ (targCells %in% k)  ]
        print(paste("cluster:",k,"num cells:",length(cNames)))
        vals <- rowSums(rawTable[,cNames])
        resTable <- cbind(resTable,vals)
    }
    colnames(resTable) <- paste0(clusters)
    resTable
}

geneLengths = read.table("/b06x-isilon/b06x-m/mbCSF/annotation/hg19/hg19_genes.length", header=1)



convToRpkm <- function(clusterMatrix, geneLengths) {
    ns=colSums(clusterMatrix)
    rpkmTable=NULL
    # issue with different genes
    commonGenes <- intersect(rownames(clusterMatrix) , rownames(geneLengths))
    geneLengths <- geneLengths[commonGenes, ,drop=FALSE]
    clusterMatrix <- clusterMatrix[commonGenes,]

    geneLengths <- geneLengths[row.names(clusterMatrix), ,drop=FALSE]
    for(i in 1:ncol(clusterMatrix)){
        rpkmTable = cbind(rpkmTable,
        round((clusterMatrix[,i]*10^9)/(geneLengths$Length*ns[i]), digits=3))
    }
    colnames(rpkmTable) = colnames(clusterMatrix)
    rownames(rpkmTable) = rownames(clusterMatrix)
    rpkmTable
}
 

# select samples
fullInfo <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/MB_RB_runs_info.280623.txt", comment.char="#")
rownames(fullInfo) <- paste0(fullInfo$TUM_ID,"_",fullInfo$RB_ID)

# specific MYC/MYCN/PRDM6 cohort
sAnn <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/MB_MYC_MYCN_PRDM6_cohort.210723.txt")
sAnn$Sample <- gsub("I023_015_P",  "I023_015_it1_prim",sAnn$Sample )
rownames(sAnn) <- sAnn$Sample
grpName = "MB_G34_somatic"
sIds <- rownames(sAnn)


# adjusted IDS
targInfo <- fullInfo[fullInfo$TUM_ID %in% sIds,]


ob.list <- list()
for (sId in rownames(targInfo)) {
    sPath = paste0(dataDir,sId,".result.rds")
    print(sId)
    ob.list[[sId]] <- readRDS(sPath)
}


mb.comb <- merge(ob.list[[1]],
                y = ob.list[2: length(ob.list)],
                #add.cell.ids = sIds,
                project = grpName)

# clean preivous result
mb.comb$SCT_snn_res.0.3 <- NULL
mb.comb$seurat_clusters <- NULL

# avoid unique RB run genes
tagMtx <- mb.comb@assays$Spatial@data[ inGenes, ]
tCells = colSums(tagMtx) != 0  # some cells lose expression
annInfo <- mb.comb@meta.data[tCells,]
mb.comb <-  CreateSeuratObject(tagMtx[,tCells],project=grpName, assay="Spatial")
mb.comb$orig.ident <- annInfo$orig.ident
mb.comb$sample <- annInfo$sample




### no batch effect adjustment

resName = paste0(grpName,"_merged")

mb.merged <- SCTransform(mb.comb, assay = "Spatial", verbose = FALSE,ncells=25000)

# add additional annotation subgroups and groups
summary(rownames(sAnn) %in% mb.merged$sample)

require(tidyverse)
#library(purr)
mb.merged$subgroup <- map_chr(mb.merged$sample, function(x) sAnn[x,"Subgroup"] )
mb.merged$group <- map_chr(mb.merged$sample, function(x) sAnn[x,"Group"] ) 

colInfo <- read.csv("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/subtypeG34.csv")
colorSubgroups <- colInfo$Color
names(colorSubgroups) <- colInfo$Subtype
colorGroups <- c("G3"="goldenrod1", "G4"="darkgreen")


## main analysis

mb.merged <- RunPCA(mb.merged, assay = "SCT", verbose = FALSE,npcs = 30, approx=FALSE)
# default 20
ndim = 30
mb.merged <- FindNeighbors(mb.merged, reduction = "pca", dims = 1:ndim)
mb.merged <- FindClusters(mb.merged, resolution = 0.4)
mb.merged <- RunUMAP(mb.merged, reduction = "pca", dims = 1:ndim)



pdf(paste0(resDir, resName,".UMAP.pdf"), width=6, height = 4)
DimPlot(mb.merged, reduction = "umap", label = TRUE)
DimPlot(mb.merged, reduction = "umap",group.by = "sample")
dev.off()


pdf(paste0(resDir, resName,".UMAP_samples.pdf"), width=6, height = 4)

DimPlot(mb.merged, reduction = "umap",group.by = "group", 
       cols = colorGroups)

DimPlot(mb.merged, reduction = "umap",group.by = "subgroup", 
       cols = colorSubgroups)

dev.off()

# genes inspection

gNames <- rownames(mb.merged)

pdf(paste0(resDir, resName,".UMAP.gene_expr.pdf"), width=6, height = 4)

for (gName in gNames) {
    print(gName)
    print(FeaturePlot(mb.merged, features = gName, reduction = "umap", raster = TRUE))
}
dev.off()

markers <- FindAllMarkers(mb.merged, only.pos = TRUE)
fName = paste0(resDir, resName,".DEG_table.txt")
write.table(markers, fName, sep="\t", quote=F)





pdf(paste0(resDir, resName,".gene_vln_profiles_per.pdf"), width=12, height = 4)
for (gName in gNames) {
    print(gName)
    print(VlnPlot(mb.merged, features = gName, combine = TRUE,pt.size=0))
}
dev.off()


pdf(paste0(resDir, resName,".gene_vln_profiles_per_sample.pdf"), width=12, height = 4)
for (gName in gNames) {
    print(gName)
    print(VlnPlot(mb.merged, features = gName, group.by="sample", combine = TRUE,pt.size=0))
}
dev.off()


saveRDS(mb.merged, paste0(resDir,resName,".object.rds"))


### Suppl Table : QC

# resuires input targInfo  :  19 samples


resInfo <- NULL
for (sId in rownames(targInfo)) {
     sPath = paste0(dataDir,sId,".result.rds")
     print(sId)
     sObj <- readRDS(sPath)
     qcInfo <- data.frame( numCells = ncol(sObj), nCount = median(sObj$nCount_Spatial), nGenes = median(sObj$nFeature_Spatial) )
     resInfo <- rbind(resInfo, qcInfo)
}


fullResInfo <- cbind(targInfo, resInfo)


resName = paste0(resDir, "MB_spatial.MYC_MYCN_PRDM6_cohort.QC_table.txt" )
write.table(fullResInfo, resName, quote=F, sep="\t")




