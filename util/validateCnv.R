library(GenomicRanges)
library(stringr)

# This code describes how validation of InferCNV results was performed
# by comparing to bulk CNV result and computing correlation 

# Input is the main InferCNV output matrix per sample

# As reference it requires bulk methylation arrays CNV output from
# from conumee (https://bioconductor.org/packages/release/bioc/html/conumee.html)


# NOTE: all path are server specific

refDir = "/b06x-isi/g703/g-i/IlmnArrayDB/CNanalysis6/"
resDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/cnvProfiling/bulkVerif/"
 

# gene coords for RNA
geneRegDf <- read.delim("/omics/odcf/analysis/OE0290_projects/Ependymoma/annotations/cell10x/refdata-cellranger-GRCh38-3.0.0/genes/genes_ref.InferCnv.txt", header=0)
colnames(geneRegDf) <- c("id","chr","start","end")
geneRegions <- GRanges(geneRegDf[,2:4])
names(geneRegions) <- geneRegDf[,1]

# for atac signal coords
atacDir = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/atacPerSample/"

getAtacCoords <- function(sId) {
    regDf <- read.delim(paste0(atacDir, sId, "_ATAC_cnv_ref.txt"),header=0)
    colnames(regDf) <- c("id","chr","start","end")
    atacRegions <- GRanges(regDf[,2:4])
    names(atacRegions) <- regDf[,1]
    
    return (atacRegions)    
}


# compute correlation function
computeCor <- function(refId, resName, fullRegions) {
    # load meth CNV
   refPath <- paste0(refDir,substr(refId,1,1),"XX/",str_split(refId,"_")[[1]][1],"/",refId,".bins.igv")
    refInfo <- read.delim(refPath)
    tiles <- GRanges(refInfo[,1:3])
    seqlevels(tiles) <- gsub("chr","", seqlevels(tiles))

    # load single cell CNV
    inDir = paste0("/b06x-isilon/b06x-m/mbCSF/results/humanTumor/cnvProfiling/globalInferCNV/",resName)
    cnvMtx <- read.table(paste0(inDir,"/infercnv.observations.txt"))
    targVals = rowMeans(cnvMtx)
    targRegions <- fullRegions[names(targVals)]
    ovlps <- findOverlaps(targRegions, tiles)

    binnedVals = c()
    for (i in 1:length(tiles)) {
        gidx <- queryHits(ovlps)[ subjectHits(ovlps) == i]
        binnedVals <- c(binnedVals, mean(targVals[gidx]) )
    }
    rgSel <- !(is.na(binnedVals))
    # compure cor
    res <- cor.test( binnedVals[rgSel], refInfo[rgSel, 5])# ,method = "spearman" ) 
    fName = paste0(resDir,resName,".",refId,".CNV_cor.pdf")
    pdf(fName,width=6,height=6)
    mTitle = sprintf("CNV cor = %.3f",res$estimate) 
    plot(binnedVals[rgSel], refInfo[rgSel,5], pch=20, col="black",
        xlab= paste("Ref:",refId) , ylab=paste("Targ:",resName), 
             main = mTitle )
    dev.off()

    return(res)
}


# TEST example
#resName = "MB272_meta_clones_hg38"
# resName = "MB272_ATAC_meta_clusters_n3"
# geneRegions <- getAtacCoords("MB272")
#refId <- "9761749012_R04C01" # MB272 
#res <- computeCor(refId, resName, geneRegions) # test


## ICGC MB MYC/MYCN/PRDM6
inFile = "/b06x-isilon/b06x-m/mbCSF/scripts/tumorAnalysis/spatial/MB_MYC_MYCN_PRDM6_cohort.210723.with_meth.txt"

inDf <- read.delim(inFile)
# IDs fixing
inDf$Sample <- gsub("MB165","MB165_v2",inDf$Sample) # for RNA
inDf$Sample <- gsub("N4059", "NB4059", inDf$Sample) # for RNA
inDf <- inDf[ inDf$Sample != "MB91" ,] # for ATAC



resDf <- NULL
n = nrow(inDf)
for (i in 11:n) {
    resVals = c()
    # RNA
    #tId = paste0(inDf[i,1],"_meta_clusters_n3_hg38") 
    # ATAC
    tId = paste0(inDf[i,1],"_ATAC_meta_clusters_n3") 
    geneRegions <- getAtacCoords(inDf[i,1])    

    for (j in 1:n ) {
        mId = inDf[j,"MethId"]
        print(paste(tId,":",mId))
        r <- computeCor(mId, tId,geneRegions) 
        print(r$estimate)
        resVals <- c(resVals,r$estimate) 
    }
    resDf <- rbind(resDf, resVals)
}


colnames(resDf) <- paste0(inDf$Sample,"_meth")

# RNA
rownames(resDf) <- paste0(inDf$Sample,"_snRNA")
write.table(resDf, paste0(resDir,"comparison_meth_snRNA.250424.txt"), sep="\t",quote=F)

# ATAC
rownames(resDf) <- paste0(inDf$Sample,"_snATAC")
write.table(resDf, paste0(resDir,"comparison_meth_snATAC.250424.txt"), sep="\t",quote=F)



library(RColorBrewer)
library(gplots)

fName = paste0(resDir,"comparison_meth_snRNA_heatmap.250424.pdf")
fName = paste0(resDir,"comparison_meth_snATAC_heatmap.260424.pdf")

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 32)

pdf(fName, width=8,height=8)
heat.out <- heatmap.2(as.matrix(resDf ),
                      main = "Comparison: snRNAseq to meth CNV", # heat map title
                      notecol="black",      
                      density.info="none",  
                      trace="none",         
                      scale="row",    
                      keysize=1,
                      #margins =c(5,5),     # widens margins around plot
                      dendrogram="none",
                      Rowv = F, Colv = F,
                      col = my_palette,
                      #lhei=c(3,7),
                      cexRow = 0.5, cexCol = 0.5
                      #lwid=c(3,9),
                     
)
dev.off()


