require(Seurat)
require(dplyr)
require(stringr)

# input

args =  commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("Input data is is not provided.", call.=FALSE)
} else {
    mainId  = args[1]
}

# trick
blocks = strsplit(mainId, "_\\s*(?=[^_]+$)", perl=TRUE)
tId = blocks[[1]][1]
sId = blocks[[1]][2]


print("load sample  RB  matrix")

# latest result on the same imsages for all samples
rbResDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/cellposeRes/seuratAnalysisR3/"

# results for precise cell border
#rbResDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/cellposeRes/seuratAnalysisR3v2/"

rbObj = readRDS(paste0(rbResDir, tId, "_", sId, ".result.rds"))

print("load refrence single cell matrix")

snRNAseqDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/perSampleV2/"
adj = "" # default
if (tId == "MB165") {
    adj ="_v2" # fix for MB165
}

tId2 = paste0(tId,adj)
targPath = paste0(snRNAseqDir, tId2,"_obj.RDS")
if (!(file.exists(targPath)) ) {
    print("No reference found!")
    stop()
}
targRefObj <- readRDS(targPath)

refId = "Tumor_prolif_CNV" # manual CNV segments with prolif
cnvAnn <- read.delim(paste0(snRNAseqDir, tId2,".CNV_prolif_blocks_ann.txt"))


if (FALSE) {

refId = "Tumor_CNV" # manual CNV segments
cnvAnn <- read.delim(paste0(snRNAseqDir, tId2,".CNV_blocks_ann.txt"))
#}

refId = "Tumor_detailed_CNV" # manual CNV segments with prolif, prog, diff
cnvAnn <- read.delim(paste0(snRNAseqDir, tId2,".CNV_ABC_blocks_ann.txt"))

}

summary(rownames(cnvAnn) == rownames(targRefObj@meta.data))
targRefObj@meta.data$cnvBlock <- cnvAnn$cnvBlock


# run TransferAnchor

anchors <- FindTransferAnchors(reference = targRefObj, query = rbObj, normalization.method = "SCT",approx=F,dims=1:20,reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors,
    refdata = targRefObj$cnvBlock,
    weight.reduction = rbObj[["pca"]], dims = 1:20)

# save result
resName = paste0( rbResDir,mainId,".anchor_predict.ref_",refId ,".txt")
write.table(predictions, resName, sep="\t", quote=F)

# predictions <- read.delim(resName)


rbObj$predict <- predictions$predicted.id

require(ggplot2)
refCellTypes = summary(as.factor(rbObj@meta.data$predict))
print(refCellTypes)
refCellTypes <- refCellTypes[ refCellTypes > 1] # fix issue

refPreds = colnames(predictions)[2:(ncol(predictions)-1)]

resName2=paste0( rbResDir,mainId,".UMAP_anchor_predict.ref_",refId ,".pdf")
pdf(resName2, width=8, height = 6)
print(DimPlot(rbObj, reduction = "umap",group.by = "predict" ))
for (ct in refPreds) {
    print(ct)
    rbObj@meta.data$pred.score = predictions[ , ct ]
    print(FeaturePlot(rbObj, reduction = "umap", feature = "pred.score" ) + labs(title=ct))
}
dev.off()



# fix order
rbObj@meta.data$predict <- as.factor(rbObj@meta.data$predict)

fName <- paste0(rbResDir, mainId,".spatial_anchor.ref_",refId,".pdf")
pdf(fName, width = 8, height = 6)

print(ImageDimPlot(rbObj, group.by="predict", 
        border.color = NA)) # for the detailed loci
        #size=1) ) # for points

for (ct in names(refCellTypes)) {
    print(ct)

    targCells  = Cells(rbObj)[rbObj@meta.data$predict == ct]
    print( ImageDimPlot(rbObj, group.by="predict",
           border.color = NA, # for detailed cell border
           #size=1,  # for points
           cells = targCells ) )   
}

dev.off()

# PNG more fitting for precise cell border

subDir = paste0(rbResDir, mainId,"_anchor_figures/")
if (!file.exists(subDir)){
  dir.create(subDir)
}

fName <- paste0(subDir, mainId,".spatial_anchor.ref_",refId,".png")
png(fName, width = 1000, height = 800)
print(ImageDimPlot(rbObj, group.by="predict", 
                    #cols=c("darkorchid4",  "green4"), # PRDM6  
                    #dark.background = F,
                    #size = 1)) # for points
                    border.color = NA )) # for precise cell border
dev.off()


for (ct in names(refCellTypes)) {
    print(ct)
    fName <- paste0(subDir,mainId,".spatial_anchor.ref_", refId,"_", ct,".png")
    png(fName,width=800,height = 600)
    targCells  = Cells(rbObj)[rbObj@meta.data$predict == ct]
    print( ImageDimPlot(rbObj, group.by="predict", 
               #size = 1, # For points
                border.color = NA, # for precise cell border
                cells = targCells ) )   
    dev.off()    
}


print("DONE!")




