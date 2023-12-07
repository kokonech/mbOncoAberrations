# R 4.2.2
library(Seurat)
library(SingleCellExperiment)
library(dplyr)

setwd("/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/perSampleV2")

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
}

print("Loading data...")

load(annData) # plot.data
load(inputData) # com.sce

targAnn <- plot.data[ plot.data$Sample %in% sIds, c( "Batch", "Sample", "SVM_prilabKD", "Nor") ]


targAnn <- targAnn[targAnn$Nor != "ND",]

targMtx <- round( counts(com.sce)[,rownames(targAnn)] )
targMtx <- targMtx[ rowSums(targMtx) >0, ]
sel2 <- colSums(targMtx) != 0 # avoid zero columns
targMtx <- targMtx[, sel2]
targAnn <- targAnn[sel2,]


## process


print("Processing data...")

mb <- CreateSeuratObject(targMtx, project=sId, assay = "RNA")
mb@meta.data <- cbind(mb@meta.data, targAnn)

mb  <- PercentageFeatureSet(mb, pattern = "^MT-", col.name = "percent.mt")
# SCT -> required for FindAnchors 
# default n feature = 2500
mb <- SCTransform(mb, vars.to.regress = "percent.mt", 
                    variable.features.n = 2500, verbose = FALSE)

# default = 20
ndim = 20
mb <- RunPCA(mb, verbose = FALSE)#, features = VariableFeatures(mb))
mb <- FindNeighbors(mb, reduction = "pca", dims = 1:ndim)
# default 0.3
mb <- FindClusters(mb, resolution = 0.3)
mb <- RunUMAP(mb, reduction = "pca", dims = 1:ndim)

pdf(paste0(sId,".clusters_UMAP.pdf"), width=8, height = 6)
DimPlot(mb, reduction = "umap", label = TRUE)
# latest udpate
DimPlot(mb, reduction = "umap", label = F,group.by="Nor")
# pairs
DimPlot(mb, reduction = "umap", label = F,group.by="Sample")
dev.off()

pdf(paste0(sId,".split_clusters_UMAP.pdf"), width=12, height = 10)
DimPlot(mb, reduction = "umap", label = F,split.by = 'seurat_clusters',ncol=5)
dev.off()


markers <- FindAllMarkers(mb, only.pos = TRUE)
fName = paste0(sId,".DEG_table.txt")
write.table(markers, fName, sep="\t", quote=F)


top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

pdf(paste0(sId,".DEG_heatmap.pdf"),width=16,height = 14)
DoHeatmap(mb,features = top10$gene)
dev.off()

gList <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/spatial_gene_selection_MB.160322.txt")

spGenes <- gList$Gene[gList$Gene %in% rownames(mb)]

pdf(paste0(sId,".RB_targets.pdf"), width=6, height = 4)
for (gName in spGenes) {
    print(gName)
    print(FeaturePlot(mb, features = gName, reduction = "umap",raster=T))
}

dev.off()




saveRDS(mb, paste0(sId,"_obj.RDS" ))

# for InferCNV
gz1 <- gzfile(paste0(sId,"_raw_counts.txt.gz"), "w")
rawCounts <- as.matrix(mb@assays$RNA@counts)
write.table(rawCounts,gz1,sep="\t",quote=F)
close(gz1)

annTable <- mb@meta.data
annTable2 <- annTable[,"seurat_clusters",drop=F]
annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
write.table(annTable2, paste0(sId,"_cnv_ann.txt") ,col.names = F,sep="\t",quote=F)

# only for combintations
annTable2 <- annTable[,"Sample",drop=F]
write.table(annTable2, paste0(sId,"_sample_ann.txt") ,col.names = F,sep="\t",quote=F)


