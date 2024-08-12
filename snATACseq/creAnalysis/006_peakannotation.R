## This script annotate peaks to best matching protein coding gene for GO analysis
suppressPackageStartupMessages({
  library(rGREAT)
  library(ggplot2)
  library(gridExtra)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
})
## In this hybrid approach I will combine GREAT approach with
## peak2gene links based on correlation
DIR_OUT <- "/path/to/output/ATACana/"
DIR_OTH <- "/path/to/other/files/"

setwd(DIR_OUT)
# ###################GREAT Annotation ###################
###GREAT annotation ------
##loading peak matrix
load("comATAC/Peaks_robust.RData")

bed <- read.delim("comATAC/G34_allPeaks_robust.fil.sorted.bed", header = FALSE) #robust peaks after removing blacklisted regions
sel <- paste0(bed$V1,":",bed$V2,"-",bed$V3)
head(sel)
##filter to only robust peaks and removing blacklisted region
PS$Final <- "No"
keep <- row.names(PS) %in% sel
table(keep)
PS[keep,"Final"] <- "Yes"
PS <- PS[PS$Final=="Yes",]
###converting to GRanges
peaks <- row.names(PS)
table(duplicated(peaks))##all unique!!
PS_gr <- makeGRangesFromDataFrame(PS,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames"),
                         start.field="start",
                         end.field="end",
                         strand.field="strand",
                         starts.in.df.are.0based=TRUE)
##dummy geneset to avoid GO analysis
gs <- list(What=c("ENSG00000225880","ENSG00000288531","ENSG00000230368","ENSG00000272438","ENSG00000223764","ENSG00000212907"))
res <- great(PS_gr,gs,"Gencode_v37")

plotRegionGeneAssociations(res)

R2G <- getRegionGeneAssociations(res)
R2G <- data.frame(R2G)
row.names(R2G) <- paste0(R2G$seqnames,":",R2G$start-1,"-",R2G$end) ##-1 in Peaks
head(R2G)

table(row.names(R2G) %in% row.names(PS)) ## 5 regions not present
##adding gene names
##loading ENS ids and getting gene name associated with them
ens <- read.delim(paste0(DIR_OTH,"gencdH38p13r37CR_genes.txt"),header = FALSE)
colnames(ens) <- c("Gene.ID","Gene_name")
ens <- ens %>% separate(Gene.ID,c("Gene.ID"))
head(ens)
ens <- ens[!duplicated(ens$Gene.ID),]
row.names(ens) <- ens$Gene.ID

peaks <- row.names(R2G)
R2G$Gene.ID <- R2G$Genes <- NA

for(x in peaks){
  gx <- unlist(R2G[x,"annotated_genes"])
  gene_names <- ens[gx,"Gene_name"]
  gene_names <- paste0(gene_names,collapse = ",")
  gx <- paste0(gx,collapse = ",")
  R2G[x,"Gene.ID"] <- gx
  R2G[x,"Genes"] <- gene_names
  rm(gx, gene_names)
}
rm(x)
head(R2G)

save(PS, PS_gr, R2G,file="OUT/Peaks_robust_GREAT.RData") ##association GREAT
###ChIPPeakannotation-----------

dd<- GenomicFeatures::makeTxDbFromGFF(paste0(DIR_OTH,"gencdH38p13r37CR.gtf"), format = "gtf",
                     dataSource="Gencodep13r37",
                     organism = "Homo sapiens")


out <- ChIPpeakAnno::genomicElementDistribution(PS_gr,
                                  TxDb = dd,
                                  #TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                  promoterRegion=c(upstream=2000, downstream=500),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-2000, -1000, -500, 0, 500),
                                    labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                               "upstream <500b", "TSS - 500b"),
                                    colors = c("#FFE5CC", "#FFCA99",
                                               "#FFAD65", "#FF8E32")))
out2 <- ChIPseeker::annotatePeak(PS_gr, tssRegion=c(-500, 500),
                             TxDb=dd)
ANN <- out$peaks
ANN2 <- data.frame(out2@anno)
row.names(ANN2) <- names(out2@anno)
save(PS, PS_gr, R2G,ANN,ANN2, file="comATAC/Peaks_robust_GREAT_annotated.RData") ##association GREAT
# q()
# #################### Peak to Gene links added GREAT annotation ###############
#load("comATAC/Peaks_robust_GREAT_annotated.RData")

P2G <- read.delim("comATAC/Peak2Gene_df_withNames.txt")
table(is.na(P2G$Correlation))
P2G[is.na(P2G$Correlation), "Correlation"] <- 0
##filter for positive corr
P2G <- P2G[P2G$Correlation>0.1,]

head(P2G)

peaks <- row.names(PS)
peaks2 <- row.names(R2G)

##getting genes associated with each peak
R2G$Genes_cor <- NA
R2G$Genes_int <- NA
R2G$MaxCor <- NA

for(x in peaks2){
  y <- gsub("-","_",x)
  keep <- P2G$peakName==y
  if(any(keep)){
    ##max Correlation and then a window around it
    del <- P2G[keep,]
    co <- max(del[,"Correlation"])
    R2G[x,"MaxCor"] <- co
    co <- co - 0.05
    genes <- del[del$Correlation>=co,"geneName"]
    if(length(genes)>0){
      R2G[x,"Genes_cor"] <- paste(genes,collapse=",")
      #intersect correlated to adjacent
      gx <- R2G[x,"Genes"]
      gx <- unlist(strsplit(gx,","))
      gx2 <- gx[gx %in% genes]
      if(length(gx2)>0){
        R2G[x,"Genes_int"] <- paste(gx2,collapse=",")
      }
      rm(y,keep,genes,gx,gx2)
    }
    rm(del,co)
  }
  
}
rm(x)
##repalcing NA values in Genes_int with best correlated or closes gene
keep1 <- is.na(R2G$Genes_cor) 
table(keep1)
R2G[keep1,"Genes_int"] <- R2G[keep1,"Genes"]
keep2 <- !is.na(R2G$Genes_cor) & is.na(R2G$Genes_int)
table(keep2)
R2G[keep2,"Genes_int"] <- R2G[keep2,"Genes_cor"]

save(PS, PS_gr, R2G,ANN,ANN2, file="comATAC/Peaks_robust_GREAT_annotated_cor.RData") ##association GREAT, 
q()
