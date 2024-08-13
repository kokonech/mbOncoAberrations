# requirements
# for Seurat: module load R/4.2.2-foss-2022a

library(Seurat)
library(dplyr)
library(config)
options(stringsAsFactors=FALSE)


##### MAIN ANALYSIS

args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
    #dataPath = args[2]
}

print(paste("Sample ID:", sId))

dataDir = "/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/rna/Result/"


timeStr = format(Sys.time(), "%Y%m%d")

# default
dataPath = paste0(dataDir,sId,"/counts/counts_",sId,"_full.gn.txt.gz")
resDir = paste0(dataDir,sId,"/",sId,"_S4_res.", timeStr, "/")

print("Read initial data...")
mols <- read.delim(dataPath) 

if (!file.exists(resDir)){
    dir.create(resDir)
}


ndim = 20

print("Start analysis...")

cb <- CreateSeuratObject(mols, min.cells = 5, min.features = 300, project = sId)


cb[["percent.mito"]] <- PercentageFeatureSet(cb, pattern = "^MT-")


print("QC report...")
pdf(paste0(resDir, sId,"_premrna_VlnPlot.pdf"), width=10, height=5)
VlnPlot(cb, c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()


print("Normalization...")

cb <- NormalizeData(cb, normalization.method = "LogNormalize", scale.factor = 10000)

cb <- FindVariableFeatures(cb, selection.method = "vst", nfeatures = 1000)

top10 <- head(VariableFeatures(cb), 10)
plot1 <- VariableFeaturePlot(cb)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(paste0(resDir, sId,"_top10_variable.pdf"), width=12, height=6 )
CombinePlots(plots = list(plot1, plot2))
dev.off()


all.genes <- rownames(cb)
cb <- ScaleData( cb, features = all.genes )


print("Dimensioional reduction...")
cb <- RunPCA(cb, features = VariableFeatures(object = cb))


png(paste0(resDir,sId,"_PcElbow.png"), width=600, height=400)
ElbowPlot(cb, ndims = 20)
dev.off()


print("Find clusters...")
ndim = 20
cb <- FindNeighbors(cb, dims = 1:ndim)
cb <- FindClusters(cb, resolution = 0.8)
table(Idents(cb))


print("tSNE...")
cb <- RunTSNE(cb, dims.use = 1:ndim)
pdf(paste0(resDir,sId,"_tSNE.pdf"))
TSNEPlot(cb, pt.size = 1, label = T)
dev.off()

print("UMAP...")
cb <- RunUMAP(object = cb, dims = 1:ndim)
pdf(paste0(resDir,sId,"_UMAP.pdf"),width = 8, height = 6)
DimPlot(object = cb, pt.size=1, reduction = 'umap')
dev.off()

print("Save result...")
saveRDS(cb, paste0(resDir,sId,"_object.Rdata"))

# cb <- readRDS(paste0(resDir,sId,"_object.Rdata")) 

gNames <- c("EOMES", "LMX1A", "MYC", "DDX1", "MYCN", "ELP4", "IMMP1L", "PAX6", "SNCAIP", "PRDM6", "MKI67", "TOP2A", "FOXP2", "CD44", "ITGAM", "TMEM119", "GPC5", "IGFBP7","OLIG2","AQP4", "TNC","PTPRC" )

gNames <- gNames[gNames %in% rownames(cb)]

pdf(paste0(resDir, sId,"_UMAP.gene_info.pdf"), width=8, height=6)
for (gName in gNames) {
  print(gName)
  print(FeaturePlot(cb, features = gName, reduction = "umap"))
}
dev.off()

print("Find DEGs...")
cb.markers <- FindAllMarkers(cb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(cb.markers, paste0(resDir, sId,"_clusters_marker_genes.txt") , sep = "\t", quote = FALSE)


require(SingleR)
require(SingleCellExperiment)

input_counts <- cb@assays$RNA@data
input_target <- SingleCellExperiment(assays = list(counts = input_counts))

rId = "CB_fetal"
cbRef = "/b06x-isilon/b06x-m/mbCSF/results/humanCbData/extFetal/Cerebellum_cell_types_pseudobulk.gn.log2.txt"


logCbCtrl <- read.delim(cbRef)
cb_ref <- SummarizedExperiment( assays=list(logcounts = logCbCtrl) )
colData(cb_ref)$celltypes <-  colnames(logCbCtrl)

pred.cb <- SingleR(test=input_target, ref = cb_ref , assay.type.test=1,
    labels = cb_ref$celltypes )

cb@meta.data$ref_cb <- pred.cb$pruned.labels

pdf(paste0(  resDir,sId,"_UMAP.SingleR_",rId,".pdf"), width=8, height = 6)
DimPlot(cb, reduction = "umap",group.by = "ref_cb")
dev.off()

resName2 = paste0( resDir,sId,".comparison_to_", rId ,".SingleR.tsv")
write.table(pred.cb, resName2, sep="\t", quote=F)


print("Done!")




