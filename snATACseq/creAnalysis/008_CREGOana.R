## Marker peak analysis per sample by clonal clusters
suppressPackageStartupMessages({
  library(enrichR)
  library(tidyr)
  library(dplyr)
  library(msigdbr)
  library(clusterProfiler)
})

args = commandArgs(trailingOnly=TRUE)
DIR_IN <- "/path/to/output/ATACana/"
DIR_OUT <- "/path/to/output/ATACana/CREres/"

DIR_OTH <- "/path/to/other/files/"

setwd(DIR_IN)
##load marker CREs-
Markers <- read.delim(paste0("CREres/",args,"/MarkerPeaks.txt"))

head(Markers)
## get a list of CREs associated with subclones
Marker_list <- list()
clones <- unique(Markers$Cluster)
for(x in clones){
  Marker_list[[x]] <- Markers[Markers$Cluster==x,"Peaks"]
}
rm(x)
lengths(Marker_list)

## get genes associated with regions from peak 2 Gene links
load("comATAC/Peaks_robust_GREAT_annotated_cor.RData") ##generated in 006_peakannotation.R
rm(ANN,ANN2,PS_gr,PS)
##using intera
GSEA <- list()
CREs2Gene <- c()
for(x in clones){
  del <- Marker_list[[x]]
  del <-del [del %in% row.names(R2G)]
  del_cre <- R2G[del,c("seqnames","start","end","Genes_int")]
  del_cre$Clone <- x
  CREs2Gene <- rbind(CREs2Gene,del_cre)
  gene_sel <- R2G[del,"Genes_int"]
  GSEA[[x]] <- unique(unlist(strsplit(gene_sel,",")))
  rm(del,gene_sel,del_cre)
}
rm(x)
lengths(GSEA)

##EnrichR
# terms <- c()
# for(x in clones){
#   genes <- GSEA[[x]]
#   ks <- enrichr(genes = genes,databases = "Chromosome_Location")
#   ks <- ks$Chromosome_Location
#   ks <- ks[order(ks$P.value),]
#   ks <- ks[1:20,c("Term","P.value","Genes")]
#   ks$Cluster <- x
#   terms <- rbind(terms,ks)
#   rm(genes,ks)
# }
#rm(x)

### MSigDB
pid_gene_sets = msigdbr(species = "Homo sapiens", category = "C1") #positional geneset
head(pid_gene_sets)
msigdbr_t2g = pid_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

clones <- names(GSEA)
terms <- c()
for(x in clones){
  ks <- enricher(gene = GSEA[[x]], TERM2GENE = msigdbr_t2g)
  res <- ks@result
  keep <- res$p.adjust<=0.1
  if(any(keep)){
    res <- res[keep,c("ID","GeneRatio","p.adjust","geneID")]
    row.names(res) <- NULL
    res$Cluster <- x
    terms <- rbind(terms,res)
    rm(ks,keep,res)
  }
  
}
rm(x)


save(Markers, Marker_list, GSEA,terms,CREs2Gene, file =paste0("CREres/",args,"/ChrEnrch.RData"))

q()

