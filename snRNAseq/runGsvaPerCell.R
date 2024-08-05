# R/4.2.2-foss-2022a
library(Seurat)
library(GSVA)


args =  commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("Input data is is not provided.", call.=FALSE)
} else {
    inFile  = args[1] # rds seurat file
}

snResDir = paste0(dirname(inFile),"/")
sId = gsub("_obj.RDS","",basename(inFile))
print(paste("Processing:",sId))
mb <- readRDS(inFile)


print("Load reference")

scriptDir = Sys.getenv("SRC")
print(paste("SRC dir:",scriptDir))
refId = "MB_programs_Hov"
refData <- paste0(scriptDir,"/ann/MB_Hovestadt_groups.symbols.gmt")

# load gene list 
gDf <- read.delim(refData, header = F,row.names = 1)
gDf$V2 <- NULL

gList <- list() 
for (i in 1:nrow(gDf)) {
  # avoid empty
  genes <- as.character(gDf[i,])
  genes <- genes[ genes != ""]
  gList[[i]]   <- genes
}

# Normalized
targData <- mb@assays$SCT@data
gsvascore <- gsva(targData, gList, method="ssgsea", parallel.sz = 4)

rownames(gsvascore) <- rownames(gDf)
mb@meta.data <- cbind(mb@meta.data,t(data.frame(gsvascore)))

print("Plot result...")

fName = paste0(snResDir, sId,"_UMAP.gsva_scores.",refId,".pdf")
pdf(fName, width=8,height=6)

for ( fName in rownames(gDf) )  {
    print(  FeaturePlot(mb, reduction = "umap",features  = fName,cols = c("grey","red")) )
}

dev.off()


fName2 = paste0(snResDir, sId,"_gsva_scores.",refId,".txt")
write.table(mb@meta.data[, rownames(gDf)] , fName2, sep="\t",quote=F)



print("Done!")




