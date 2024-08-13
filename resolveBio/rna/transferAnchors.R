# R 4.1.1
require(Seurat)
require(dplyr)
require(stringr)

# input

args =  commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("Input data is is not provided.", call.=FALSE)
} else {
    inPath  = args[1]
}


resDir = dirname(inPath)
sId =  gsub("_object.Rdata", "", basename(inPath))

rObj = readRDS(inPath)

print("load refrence single cell matrix")



snRNAseqDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/perSampleV2/"
adjId = sId  # default

targPath = paste0(snRNAseqDir, adjId,"_obj.RDS")
if (!(file.exists(targPath)) ) {
    print("No reference found!")
    stop()
}
targRefObj <- readRDS(targPath)

refId = "Tumor_CNV" # manual CNV segments
cnvAnn <- read.delim(paste0(snRNAseqDir, adjId,".CNV_blocks_ann.txt"))



summary(rownames(cnvAnn) == rownames(targRefObj@meta.data))
targRefObj@meta.data$cnvBlock <- cnvAnn$cnvBlock


print("run analysis")

# run TransferAnchor

anchors <- FindTransferAnchors(reference = targRefObj, query = rObj, normalization.method = "SCT",approx=F,dims=1:20,reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors,
    refdata = targRefObj$cnvBlock,
    weight.reduction = rObj[["pca"]], dims = 1:20)

# save result
resName = paste0( resDir,"/",sId,".anchor_predict.ref_",refId ,".txt")
write.table(predictions, resName, sep="\t", quote=F)


# plot result

rObj$predict <- predictions$predicted.id

require(ggplot2)
refCellTypes = summary(as.factor(rObj@meta.data$predict))
print(refCellTypes)
refCellTypes <- refCellTypes[ refCellTypes > 1] # fix issue

refPreds = colnames(predictions)[2:(ncol(predictions)-1)]

resName2=paste0( resDir,"/", sId,".UMAP_anchor_predict.ref_",refId ,".pdf")
pdf(resName2, width=8, height = 6)
print(DimPlot(rObj, reduction = "umap", pt.size = 2, group.by = "predict" ))
for (ct in refPreds) {
    print(ct)
    rObj@meta.data$pred.score = predictions[ , ct ]
    print(FeaturePlot(rObj, reduction = "umap", feature = "pred.score" ) + labs(title=ct))
}
dev.off()


print("DONE!")

