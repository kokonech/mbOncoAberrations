# NOTE: this code depends on server envirnoment with data locations

library(Seurat)
library(Signac)
library(dplyr)

options(stringsAsFactors=FALSE)
setwd("/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial")

resDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/atacMerged/"

sInfo <- read.delim("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/MB_MYC_MYCN_PRDM6_cohort.210723.txt")
sIds <- sInfo$Sample
rownames(sInfo) <- sIds

# requires up to 150 (???) GB RAM
grpName = "MB_G34_somatic"

ob.list <- list()
ids.list <- c()
for (sId in sIds) {
    print(sId)
    inFile = paste0("atacPerSample/",sId,"_obj.RDS") 
    if ( file.exists(inFile) ) {
        mb <- readRDS(inFile)
        #ob.list[[sId]] <- mb
        gene.activities <- GeneActivity(mb)
        mbObj <- CreateAssayObject(counts = gene.activities)
        ids.list <- c(ids.list,rep(sId,ncol(ob.list[[sId]])) )
    } else {
        print("Input does not exist, skipping...")
    }
}


mb.merged <- merge(ob.list[[1]],
                y = ob.list[2: length(ob.list)],
                project = grpName)

mb.merged$Sample = ids.list

summary(rownames(sInfo) %in% mb.merged$Sample)
require(tidyverse)
mb.merged$subgroup <- map_chr(mb.merged$Sample, function(x) sInfo[x,"Subgroup"] )
mb.merged$group <- map_chr(mb.merged$Sample, function(x) sInfo[x,"Group"] )



mb.merged <- RunTFIDF(mb.merged)
# q0 - all features
# q75  - top 25% : result is similar
mb.merged <- FindTopFeatures(mb.merged, min.cutoff = 'q0')
mb.merged <- RunSVD(mb.merged)


ndim = 30 # default 30

fName = paste0(resDir,grpName, "_seq_depth_cor.pdf")
pdf(fName,width=6,height=4)
DepthCor(mb.merged)
dev.off()

mb.merged <- RunUMAP(object = mb.merged, reduction = 'lsi', dims = 2:ndim)
mb.merged <- FindNeighbors(object = mb.merged, reduction = 'lsi', dims = 2:ndim)
mb.merged <- FindClusters(object = mb.merged, verbose = FALSE, algorithm = 3)



colorSubgroups <- c("I"="purple","II"="violetred", "III"="sienna",
                 "IV"="yellow", "V"="yellowgreen","VI"="cyan",
                 "VII"="skyblue", "VIII"= "darkgreen")
names(colorSubgroups) <- paste0("G34_", names(colorSubgroups))
colorGroups <- c("G3"="yellow2", "G4"="green")


colInfo <- read.csv("/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/subtypeG34.csv")
colorSubgroups <- colInfo$Color
names(colorSubgroups) <- colInfo$Subtype
colorGroups <- c("G3"="goldenrod1", "G4"="darkgreen")

pdf(paste0(resDir,grpName,"_UMAP.clusters.pdf"),width = 8, height = 6)
DimPlot(object = mb.merged, pt.size=1, label=T,reduction = 'umap')
# v1
#DimPlot(object = mb.merged, pt.size=1, reduction = 'umap', group.by="orig.ident")
# v2
DimPlot(mb.merged, reduction = "umap",group.by = "group",
       cols = colorGroups)
DimPlot(object = mb.merged, label=F,reduction = 'umap', group.by="subgroup",
            cols = colorSubgroups)
DimPlot(object = mb.merged, label=T,reduction = 'umap', group.by="Sample")
dev.off()


saveRDS(mb.merged, paste0(resDir,grpName,"_object.Rdata"))
# mb.merged <- readRDS( paste0(resDir,grpName,"_object.Rdata") )


