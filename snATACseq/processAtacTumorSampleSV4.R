# requirements
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)


options(stringsAsFactors=FALSE)


##### MAIN ANALYSIS


resDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/atacPerSample/"


args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
}


print(paste("Sample ID:", sId))
sInfo = read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/atacSeq/MB_targ_cohort.auto.txt",header=F)
colnames(sInfo) <- c("sId","mId")
rownames(sInfo) <- sInfo$sId

mId = sInfo[sId,"mId"]
adjId = gsub("MSM","MSN", mId)

# for full path
dataPath = paste0("/omics/odcf/analysis/OE0290_projects/mb_10x/ATACana/snMULTI/",pathId, "/outs/")

print(paste("Input data:", dataPath))

if (!(file.exists(dataPath))) {
    print("Error! Data not found, exiting...")
    quit(save = "no",status=1)   
}

print("Read initial data...")
countsPath = paste0(dataPath,"/filtered_feature_bc_matrix.h5")
counts <- Read10X_h5(countsPath)

fragpath = paste0(dataPath, "/atac_fragments.tsv.gz")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))


# create ATAC assay and add it to the object
chrom_assay  <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

mb <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  project = sId,
)


print("Start analysis...")

mb <- NucleosomeSignal(mb)
mb <- TSSEnrichment(mb)

resName = paste0(resDir, sId,"_VlnPlot_QC.pdf")
pdf(resName, width=8,height=6)
VlnPlot(
  object = mb,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0
)
dev.off()


# requires inspection
# geneLim <- 6000
# transLim <- 20000

# filter out low quality cells
mb <- subset(
  x = mb,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

pdf(paste0(resDir, sId,"_VlnPlot_QC.after_filter.pdf"), width=8, height=6)
VlnPlot(mb, features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),     pt.size=0)
dev.off()


print("Normalization...")

mb <- RunTFIDF(mb)
# q0 - all features
# q75  - top 25% : result is similar
mb <- FindTopFeatures(mb, min.cutoff = 'q0')
mb <- RunSVD(mb)

fName = paste0(resDir,sId, "_seq_depth_cor.pdf")
pdf(fName,width=6,height=4)
DepthCor(mb)
dev.off()

print("Dimensioional reduction...")

ndim = 30 # default 30

mb <- RunUMAP(object = mb, reduction = 'lsi', dims = 2:ndim)
mb <- FindNeighbors(object = mb, reduction = 'lsi', dims = 2:ndim)
mb <- FindClusters(object = mb, verbose = FALSE, algorithm = 3)

pdf(paste0(resDir,sId,"_UMAP.pdf"),width = 8, height = 6)
DimPlot(object = mb, pt.size=1, label=T)
dev.off()


print("Check snRNA-seq annoitation...")

snResDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/perSampleV2/"
blocksAnn = paste0(snResDir,sId,".CNV_blocks_ann.txt")

if (file.exists(blocksAnn)) {

    cnvAnn <- read.delim(blocksAnn)
    rownames(cnvAnn) <- gsub(paste0(adjId,"_"), "", paste0(rownames(cnvAnn),"-1"))
    summary(rownames(mb@meta.data) %in% rownames(cnvAnn))
    targAnn <- cnvAnn[ intersect(rownames(mb@meta.data) , rownames(cnvAnn)), ]

    mb$refAnn <- "Unclear"
    mb@meta.data[ rownames(targAnn) , ]$refAnn <- targAnn$cnvBlock

    pdf(paste0(resDir,sId,"_UMAP.ref_CNV_blocks.pdf"),width = 8, height = 6)
    DimPlot(object = mb, group.by = "refAnn")
    dev.off()
}




print("Save signal...")

saveRDS(mb, paste0(resDir,sId,"_obj.RDS" ))
# mb <- readRDS(paste0(resDir,sId,"_obj.RDS"))

# for InferCNV
gz1 <- gzfile(paste0(resDir,sId,"_ATAC_raw_counts.txt.gz"), "w")
rawCounts <- as.matrix(mb@assays$ATAC@counts)

rawCounts <- rawCounts[ grep("chr",rownames(rawCounts)),]
write.table(rawCounts,gz1,sep="\t",quote=F)
close(gz1)

annTable <- mb@meta.data
annTable2 <- annTable[,"seurat_clusters",drop=F]
annTable2$seurat_clusters <- paste0("cl",annTable2$seurat_clusters)
write.table(annTable2, paste0(resDir,sId,"_ATAC_cnv_ann.txt") ,col.names = F,sep="\t",quote=F)

peakRegions <- GRanges(sub("-",":",rownames(rawCounts), fixed=TRUE))
seqlevels(peakRegions) <- paste0("chr",c(1:22,"X","Y")) # fix order
peakRegions <- sort(peakRegions)


peakDf <- data.frame(peakRegions)[,1:3]
rownames(peakDf) <- paste0(  peakDf[,1],"-",peakDf[,2],"-",peakDf[,3] )
peakDf$seqnames <- gsub("chr","", peakDf$seqnames)
summary(rownames(rawCounts) %in% rownames(peakDf))
write.table(peakDf, paste0(resDir,sId,"_ATAC_cnv_ref.txt") ,col.names = F,sep="\t",quote=F)


# particular gene check : approach #1 signal strength
print("Inspect promoters...")

gene.activities <- GeneActivity(mb)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
mb[['RNA']] <- CreateAssayObject(counts = gene.activities)
mb <- NormalizeData(
  object = mb,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(mb$nCount_RNA)
)

DefaultAssay(mb) <- 'RNA'


gList <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/spatial_gene_selection_MB.160322.txt")

spGenes <- gList$Gene[gList$Gene %in% rownames(mb)]

pdf(paste0(resDir,sId,"_UMAP.target_genes.pdf"),width = 8, height = 6)

for (gName in spGenes) {
    print(gName)
    print( FeaturePlot( mb, features = gName, pt.size = 0.1,
        max.cutoff = 'q95',raster =T) 
    )
}

dev.off()

