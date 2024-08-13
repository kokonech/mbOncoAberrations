# R 4.2.2
library(ComplexHeatmap)
library(gplots)
library(stringr)
col16 <- colorpanel(16,"green", "grey", "red")


args =  commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  stop("Input data is is not provided.", call.=FALSE)
} else {
    sId = args[1]
}

# 100 Kbp bin
binSize = "100000"
inDir = paste0("/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/dna/Result/BJ-DNA-QC/",sId,"_100Kbp")

## load annotation from RNA

rnaDir = "/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/rna/Result/"
annInfo <- read.delim(paste0(rnaDir, sId, "/",sId,"_S4_res/",sId,".anchor_predict.ref_Tumor_CNV.txt"))

## load input

cellDirs = list.files(inDir)

first = TRUE
regions <- NULL
segNorm <- NULL

for (cId in cellDirs) {

    print(cId)
    dirPath = paste0(inDir,"/",cId)
    targFile = paste0("ginkgo_res.binsize_",binSize,".RDS" )
    found_files <- list.files(dirPath, pattern = targFile , recursive = TRUE, full.names = TRUE)
    fObj <- readRDS(found_files[1])
    sn <- fObj$SegCopy
    if (first) {
        regions <- paste0(sn$CHR,":",sn$START,"-",sn$END)
        first = FALSE
    }
    #print(max(sn$final))
    segNorm <- cbind(segNorm, sn$final )
      

}


rownames(segNorm) <- regions
colnames(segNorm) <- cellDirs


segNorm <- t(segNorm)
chrInfo <- factor(gsub("chr","",sn$CHR), levels = c(1:22,"X","Y"))
rownames(segNorm) <- str_split_fixed(rownames(segNorm),"_",3)[,1]


segNormAdj <- segNorm

# adjust signal and exclude outliers

threshold <- mean(segNorm) + 2*sd(segNorm)
segNormAdj[segNorm > threshold] <- round(threshold)

cutLim <- 3 # cells with low signal
summary(rowMeans(segNorm) > cutLim)
segNormAdj <- segNormAdj[ rowMeans(segNorm) > cutLim, ]


# remove X,Y
sel <- grepl("chrY|chrX",colnames(segNormAdj))
segNormAdj <- segNormAdj[ ,!sel]
chrInfo <- chrInfo[!sel]

hc <- hclust(dist(segNormAdj), method = "ward.D2")


resDir = "/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/dna/Result/BJ-DNA-QC/CNV/"

row_ha = rowAnnotation(Type =  annInfo[ rownames(segNormAdj), "predicted.id" ] )
fName = paste0(resDir,sId,"_Norm_CNV.bin_",binSize,".no_XY.RNA_ann.cheatmap.pdf")

pdf(fName,width=15,height = 10)

ht1 <- Heatmap(as.matrix(segNormAdj), name=sId,
        cluster_columns = FALSE,
        cluster_rows =  hc,#TRUE,
         show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_title_side = "bottom",
        clustering_method_rows = "ward.D2",
        column_title_gp = gpar(fontsize = 8),
      column_split = chrInfo
)
ht1 + row_ha
dev.off()


## Target region

sel <- grepl("chr2:",colnames(segNormAdj)) # MYCN
#sel <- grepl("chr8:",colnames(segNormAdj)) # MYC
segNormTarg <- segNormAdj[,sel]

# use clustering from main
fName = paste0(resDir,sId,"_Norm_CNV.bin_",binSize,".targ_ref_hc.cheatmap.pdf")


pdf(fName,width=12,height = 6)

ht1 <- Heatmap(as.matrix(segNormTarg), name=sId,
        cluster_columns = FALSE,
        cluster_rows = hc,
         show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_title_gp = gpar(fontsize = 8)
)
ht1
dev.off()

fName2 = paste0(resDir, sId,"_Norm_CNV.bin_",binSize,".targ.tsv")
write.table(segNormTarg, fName2, sep="\t",row.names=F,quote=F)

