##############################################################################################################################################
## Script to extract functional mutations and likely driver mutations from the detection tumor cohort

##############################################################################################################################################
## Load settings

source("./RScripts/Settings.R")
library(ComplexHeatmap)
load("./RScripts/RData/Purity_ploidy.RData")
##############################################################################################################################################
## For each of the tumors, extract coding mutations; structural variations, high-level amplifications and deletions

functional.mutations <- data.frame(CHROM=c(), POS=c(), ID=c(), REF=c(), ALT=c(), QUAL=c(), FIlTER=c(), INFO=c(), FORMAT=c(), FORMAT_INFO=c(),
                                   ANNOVAR_FUNCTION=c(), GENE=c(), EXONIC_CLASSIFICATION=c(), READS_REF=c(), READS_ALT=c(), 
                                   CoverageRatio = c(), BAF = c(), Purity = c(), Ploidy = c(), SAMPLE=c(), AA_change=c())

## store to generate a table later on
translocations.for.output <- c()
amplifications.for.output <- c()
deletions.for.output <- c()

translocations <- data.frame(gene1=c(), gene2=c(), Sample=c(), svtype=c())
amplifications <- data.frame(gene=c(), Sample=c())
deletions <- data.frame(gene=c(), Sample=c())

## iterate through the tumors
for(i in tumors){
  
  print(i)
  
  ## Read in combined estimates of ploidy and purity:
  aceseq <- list.files(paste0(data.directory,  i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  purity. <- purity[i]
  ploidy. <- ploidy[i]
  
  ## Read in the mutation file
  files <- list.files(paste0(data.directory,  i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  if(is.na(files)){
    print(i)
    next
  }
  
  mutations <- read.vcf(paste0(data.directory,  i, "/", snv.directory, "/", files))
  mutations$vcf <- mutations$vcf[mutations$vcf$ANNOVAR_FUNCTION %in% c("exonic", "upstream", "splicing"),]
  colnames(mutations$vcf)[10] <- "FORMAT_INFO"
  mutations$vcf <- mutations$vcf[,c(1:10, 16, 17, 18, 19)]
  
  mutations$vcf$AA_change <- Extract.info.from.vcf(mutations, info="AA_change")

  counts <- Extract.info.from.vcf(mutations, info="readcounts")
  
  copy.number.info. <- read.delim(file=paste0(data.directory,  i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios/BAF/genotype/TCN for the mutations of interest
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info.)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios. <- cnv.info.per.mutation$coverage.ratio
  bafs. <- cnv.info.per.mutation$baf
  A <- cnv.info.per.mutation$A
  copy.number. <- cnv.info.per.mutation$tcn
  
  mutations <- mutations$vcf
  mutations$READS_REF <- counts[,1]
  mutations$READS_ALT <- counts[,2]
  mutations$VAF <- counts[,2]/rowSums(counts)
  mutations$CoverageRatio <- coverage.ratios.
  mutations$BAF <- bafs.
  mutations$Purity <- purity.
  mutations$Ploidy <- ploidy.
  mutations$TCN <- copy.number.
  mutations$A <- A
  mutations$SAMPLE <- i
  mutations$ANNOVAR_TRANSCRIPTS <- NULL
  
  mutations$Clonality <- apply(mutations, 1, function(x){
    
    CN <- round((as.numeric(x["CoverageRatio"])*ploidy.*purity.+2*(1-purity.)*(as.numeric(x["CoverageRatio"])-1))/(purity.))
    
    if(is.na(x["CoverageRatio"])|is.na(x["VAF"])|as.numeric(x["CoverageRatio"])==0|round(as.numeric(x["CoverageRatio"])*ploidy.)==0){return("n.d.")}
    depth <- as.numeric(x["READS_REF"]) + as.numeric(x["READS_ALT"])
    if(as.numeric(x["VAF"]) < qbinom(p = 0.05, size = depth,
                                     prob =  1*purity./(purity.*CN + 2*(1-purity.)))/depth){
      "SC"
    }else{  
      if(CN > 2 & as.numeric(x["VAF"]) < qbinom(p = 0.05, size = depth,
                                                 prob =  2*purity./(purity.*CN + 2*(1-purity.)))/depth){
        if(CN > 3){
          "LC"
        }else if(CN==3){
          "C"
        }
      }else if(CN > 2){
        "EC"
      }else{
        "C"
      }
      
    }
  })
  
  functional.mutations <- rbind(functional.mutations, mutations)
  
  
  ## same for indels
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory,  i, "/", indel.directory, "/"), pattern="somatic_indels_conf_8_to_10")[1]
  if(is.na(files)){
    print(i)
    next
  }
  
  mutations <- read.vcf(paste0(data.directory,  i, "/", indel.directory, "/", files))
  
  # for some cases, Control and Tumor columns were swapped. Identify tumor column:
  
  # 1. try column # 10
  
  test <- mutations$vcf
  colnames(test)[10] <- "FORMAT_INFO"
  counts.1 <- Extract.info.from.vcf(test, info="readcounts", type = "indel")
  
  # 2. try column # 11
  
  test <- mutations$vcf
  colnames(test)[11] <- "FORMAT_INFO"
  counts.2 <- Extract.info.from.vcf(test, info="readcounts", type = "indel")
  
  tumor.col <- which.max(c(mean(counts.1[,2]), mean(counts.2[,2])))
  
  if(i!="MB99"){
    mutations$vcf <- as.data.frame(mutations$vcf)[,c(1:9, 9 + tumor.col, 17, 18, 19, 20)]
  }else{
    mutations$vcf <- as.data.frame(mutations$vcf)[,c(1:9, 9 + tumor.col, 14, 15, 16, 17)]
  }
  rm(test)
  rm(counts.1)
  rm(counts.2)
  
  
  colnames(mutations$vcf)[10] <- "FORMAT_INFO"
  mutations$vcf$AA_change <- Extract.info.from.vcf(mutations, info="AA_change", type = "indel")
  
  if(nrow(mutations$vcf)==0){next}
  
  counts <- Extract.info.from.vcf(mutations, info="readcounts", type = "indel")
  
  copy.number.info. <- read.delim(file=paste0(data.directory,  i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios for the mutations of interest
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info.)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios. <- cnv.info.per.mutation$coverage.ratio
  bafs. <- cnv.info.per.mutation$baf
  A <- cnv.info.per.mutation$A
  copy.number. <- cnv.info.per.mutation$tcn
  
  mutations <- mutations$vcf
  mutations$READS_REF <- counts[,1]
  mutations$READS_ALT <- counts[,2]
  mutations$VAF <- counts[,2]/rowSums(counts)
  mutations$CoverageRatio <- coverage.ratios.
  mutations$BAF <- bafs.
  mutations$Purity <- purity.
  mutations$Ploidy <- ploidy.
  mutations$TCN <- copy.number.
  mutations$A <- A
  mutations$SAMPLE <- i
  mutations$ANNOVAR_TRANSCRIPTS <- NULL
  
  ## classify mutations as clonal or subclonal based on VAF
  
  mutations$Clonality <- apply(mutations, 1, function(x){
    
    CN <- round((as.numeric(x["CoverageRatio"])*ploidy.*purity.+2*(1-purity.)*(as.numeric(x["CoverageRatio"])-1))/(purity.))
    
    if(is.na(x["CoverageRatio"])|is.na(x["VAF"])|as.numeric(x["CoverageRatio"])==0|round(as.numeric(x["CoverageRatio"])*ploidy.)==0){return("n.d.")}
    depth <- as.numeric(x["READS_REF"]) + as.numeric(x["READS_ALT"])
    if(as.numeric(x["VAF"]) < qbinom(p = 0.05, size = depth,
                                     prob =  1*purity./(purity.*CN + 2*(1-purity.)))/depth){
      "SC"
    }else{  
      "C"
    }
  })
  
  
  functional.mutations <- rbind(functional.mutations, mutations)
  
  ## look up structural variants 
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory,  i, "/", sv.directory, "/"), pattern="filtered_somatic_minEventScore3.tsv")[1]
  if(is.na(files)){
    warning(paste0("no SV info for sample", i))
    print(paste0("no SV info for sample", i))
    next
  }
  
  mutations <- read.vcf(paste0(data.directory,  i, "/", sv.directory, "/", files))
  mutations <- mutations$vcf
  mutations$genepos1 <- sapply(mutations$gene1, function(x){strsplit(x, split="_")[[1]][2]})
  mutations$genepos2 <- sapply(mutations$gene2, function(x){strsplit(x, split="_")[[1]][2]})
  
  mutations$gene1 <- sapply(mutations$gene1, function(x){strsplit(x, split="_")[[1]][1]})
  mutations$gene2 <- sapply(mutations$gene2, function(x){strsplit(x, split="_")[[1]][1]})
  
  ## Take translocations between 2 genes or within 1 gene, but not just intronic deletions/amplifications
  mutations <- mutations[mutations$gene1 != mutations$gene2 | mutations$genepos1 != mutations$genepos2,]
  ## Require an event score of >= 5
  mutations <- mutations[(mutations$gene1 %in% driver.genes | mutations$gene2 %in% driver.genes) &
                           as.numeric(mutations$eventScore) >= 5,,drop=F]
  if(nrow(mutations)>0){
    mutations$Sample <- i
    translocations <- rbind(translocations, mutations[,c("gene1", "gene2", "Sample", "svtype")])
    translocations.for.output <- rbind(translocations.for.output, mutations)
  }
  
  ### Look for high-level amplifications and deletions
  
  amp <- unlist(apply(driver.genes.with.genetic.pos, 1, function(x){
    tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                     (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                     (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
    
    if(nrow(tmp)>0){
      ## Require at least 10 copies
      if(mean(as.numeric(tmp$tcnMean))>=10 ){
        return(x[1])
      }else{
        return(c())
      }
    }else{
      return(c())
    }
  }))
  

  if(length(amp)>0){
    amplifications <- rbind(amplifications, data.frame(gene=amp, Sample=i))
    
    amp.for.output <- do.call("rbind", apply(driver.genes.with.genetic.pos, 1, function(x){
      tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                       (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                       (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
      
      if(nrow(tmp)>0){
        ## Require at least 10 copies
        if(mean(as.numeric(tmp$tcnMean))>=10 ){
          tmp$Gene <- x[1]
          tmp <- tmp[tmp$tcnMean>=10,]
          return(tmp)
        }else{
          return()
        }
      }else{
        return()
      }
    }))
    
    amp.for.output$Sample <- i
    amplifications.for.output <- rbind(amplifications.for.output, amp.for.output)
  }
  
  ## Homozygous deletions; require <0.9 copy numbers
  del <- unlist(apply(driver.genes.with.genetic.pos, 1, function(x){
    tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                     (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                     (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
    
    if(nrow(tmp)>0){
      if(mean(as.numeric(tmp$tcnMean))<=0.9){
        return(x[1])
      }else{
        return(c())
      }
    }else{
      return(c())
    }
  }))
  
 
  
  if(length(del)>0){
    deletions <- rbind(deletions, data.frame(gene=del, Sample=i))
    
    del.for.output <- do.call("rbind", apply(driver.genes.with.genetic.pos, 1, function(x){
      tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                       (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                       (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
      
      if(nrow(tmp)>0){
        if(mean(as.numeric(tmp$tcnMean))<=0.9){
          tmp$Gene <- x[1]
          tmp <- tmp[tmp$tcnMean<=0.9,]
          
          return(tmp)
        }else{
          return()
        }
      }else{
        return()
      }
    }))
    del.for.output$Sample <- i
    deletions.for.output <- rbind(deletions.for.output, del.for.output)
  }
  
  
}

translocations <- unique(translocations)
translocations$SV <- "SV"

## filter non-promoter TERT mutations
functional.mutations <- functional.mutations[-which(functional.mutations$GENE=="TERT" &
                                                      !functional.mutations$POS %in% c(1295228, 1295250)),]

## keep splice-site mutations, nonsynonymous SNVs and indels

to.keep <- functional.mutations$ANNOVAR_FUNCTION=="splicing" | 
(  functional.mutations$ANNOVAR_FUNCTION=="exonic" &
  functional.mutations$EXONIC_CLASSIFICATION != "synonymous SNV" )|
  functional.mutations$EXONIC_CLASSIFICATION =="ncRNA_exonic" |
  (functional.mutations$ANNOVAR_FUNCTION=="upstream" &
     functional.mutations$GENE=="TERT" & functional.mutations$POS %in% c(1295228, 1295250))


filtered.functional.mutations <- functional.mutations[to.keep,]

write.table(filtered.functional.mutations[filtered.functional.mutations$GENE %in% driver.genes,], 
            file=paste0(meta.data, "Driver_mutations_cohort.tsv"), sep="\t", quote=F, row.names = F)


translocations <- unique(translocations)
translocations$SV <- "SV"


## simplify by adding splice sites to Exonic classification
filtered.functional.mutations$EXONIC_CLASSIFICATION[is.na(filtered.functional.mutations$EXONIC_CLASSIFICATION)] <- filtered.functional.mutations$ANNOVAR_FUNCTION[is.na(filtered.functional.mutations$EXONIC_CLASSIFICATION)]


##############################################################################################################################################
## transform the mutation info into a mutation matrix, where rows correspond to genes/chromosomes and columns to samples
mat <- matrix("", nrow=length(unique(driver.genes)), ncol=length(tumors),
              dimnames=list(c(unique(driver.genes)), tumors))

translocations <- translocations[,c("gene1", "gene2", "Sample", "SV")]
translocations <- unique(translocations)

for(i in 1:nrow(mat)){
  
  goi <- rownames(mat)[i]
  tmp <- as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==goi,]$SAMPLE)
  tmp <- intersect(tmp, tumors)
  
  for(k in tmp){
    mat[i,k] <- paste(filtered.functional.mutations[filtered.functional.mutations$GENE==goi &
                                                  filtered.functional.mutations$SAMPLE ==k,]$EXONIC_CLASSIFICATION, collapse = ";")
    
  }
    
    
    
    mat[i,intersect(tumors, as.character(deletions[deletions$gene==goi,]$Sample))] <- paste(mat[i,intersect(tumors, as.character(deletions[deletions$gene==goi,]$Sample))], "DEL", sep=";")
    mat[i,intersect(tumors, as.character(amplifications[amplifications$gene==goi,]$Sample))] <- paste(mat[goi,intersect(tumors, as.character(amplifications[amplifications$gene==goi,]$Sample))], "AMP", sep=";")
    mat[i,as.character(translocations[(translocations$gene1==goi | translocations$gene2==goi) &
                                        translocations$Sample %in% tumors,]$Sample)] <- 
      paste(mat[goi,as.character(translocations[(translocations$gene1==i | translocations$gene2==goi) & translocations$Sample %in% tumors,]$Sample)], translocations[(translocations$gene1==goi | translocations$gene2==goi)&
                                                                                                                                                              translocations$Sample %in% tumors,]$SV, sep=";")
}


## Manual adjustments:

## Take information from other arrays
mat["MYC",rownames(sample.information[sample.information$MYC==1,])] <- paste(mat["MYC",rownames(sample.information[sample.information$MYC==1,])], "AMP", sep=";")
mat["MYCN",rownames(sample.information[sample.information$MYCN==1,])] <- paste(mat["MYCN",rownames(sample.information[sample.information$MYCN==1,])], "AMP", sep=";")
mat["SNCAIP",rownames(sample.information[sample.information$PRDM6==1,])] <- paste(mat["SNCAIP",rownames(sample.information[sample.information$PRDM6==1,])], "AMP", sep=";")


##############################################################################################################################################
## Oncoprint

column_annotation <- data.frame(Group = sample.information[tumors,]$mnp11,
                                Subgroup = sample.information[tumors,]$mnp12)
rownames(column_annotation) <- tumors


col = c("black", "black", "firebrick", c(brewer.pal(n=8, name="Paired")), c("green", "#B2DF8A"))
names(col) <- c("Gain", "Loss", "nonsynonymous SNV", "splicing", "stopgain", "frameshift deletion", "nonframeshift deletion",
                "DUP", "AMP", "DEL", "SV", "upstream", "frameshift insertion")

#palette_check(col2hcl(col), plot = TRUE)


alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = "white", col = NA))
  },
  "Gain" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["Gain"], col = NA))
  },
  "Loss" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["Loss"], col = NA))
  },
  "AMP" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  "nonsynonymous SNV" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = col["nonsynonymous SNV"], col = NA))
  },
  "upstream" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = col["upstream"], col = NA))
  },
  "DEL" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["DEL"], col = NA))
  },
  "SV" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = col["SV"], col = NA))
  },
  "nonframeshift deletion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["nonframeshift deletion"], col = NA))
  },
  "splicing" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["splicing"], col = NA))
  },
  "stopgain" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["stopgain"], col = NA))
  },
  "frameshift deletion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["frameshift deletion"], col = NA))
  },
  "frameshift insertion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["frameshift insertion"], col = NA))
  },
  "DUP" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["DUP"], col = NA))
  }
)

column_title = "OncoPrint for Medulloblastoma"


annotation_colors <- list(Group = cbp1[1:length(unique(sample.information$mnp11))],
                          Subgroup = hcl.colors(n=length(unique(sample.information$mnp12))))

names(annotation_colors$Group) <- unique(sample.information$mnp11)
names(annotation_colors$Subgroup) <- unique(sample.information$mnp12)

pdf(paste0(output.directory, "Oncoprint.pdf"), width=14, height=9, useDingbats = F)

oncoPrint(mat,
          bottom_annotation = HeatmapAnnotation(df=column_annotation, col=annotation_colors),
          alter_fun = alter_fun, col = col, row_names_gp = gpar(fontsize = 8),
          remove_empty_columns = F, remove_empty_rows = TRUE, show_column_names = T,
          column_title = column_title, show_pct = F)



dev.off()


##############################################################################################################################################
## We plot a second oncoprint that is not concerned with the type of mutations but whether they are clonal or subclonal
mat.clonal <- matrix("", nrow=length(unique(driver.genes)), ncol=length(tumors),
              dimnames=list(c(unique(driver.genes)), tumors))


for(i in 1:nrow(mat.clonal)){
  
  goi <- rownames(mat.clonal)[i]
  tmp <- as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==goi,]$SAMPLE)
  tmp <- intersect(tmp, tumors)
  
  for(k in tmp){
    mat.clonal[i,k] <- paste(filtered.functional.mutations[filtered.functional.mutations$GENE==goi &
                                                      filtered.functional.mutations$SAMPLE ==k,]$Clonality, collapse = ";")
    
  }
  
  mat.clonal[i,intersect(tumors, as.character(deletions[deletions$gene==goi,]$Sample))] <- paste(mat.clonal[i,intersect(tumors, as.character(deletions[deletions$gene==goi,]$Sample))], "n.d.", sep=";")
  mat.clonal[i,intersect(tumors, as.character(amplifications[amplifications$gene==goi,]$Sample))] <- paste(mat.clonal[goi,intersect(tumors, as.character(amplifications[amplifications$gene==goi,]$Sample))], "n.d.", sep=";")
  mat.clonal[i,as.character(translocations[(translocations$gene1==goi | translocations$gene2==goi) &
                                      translocations$Sample %in% tumors,]$Sample)] <- 
    paste(mat.clonal[goi,as.character(translocations[(translocations$gene1==i | translocations$gene2==goi) & translocations$Sample %in% tumors,]$Sample)], "n.d.", sep=";")
}


## Manual adjustments:

## Take information from other arrays
mat.clonal["MYC",rownames(sample.information[sample.information$MYC==1 & sample.information$Clonality_of_drivers%in% "clonal",])] <- "clonal"
mat.clonal["MYC",rownames(sample.information[sample.information$MYC==1 & sample.information$Clonality_of_drivers %in% "subclonal",])] <- "subclonal"
mat.clonal["MYC",rownames(sample.information[sample.information$MYC==1 & is.na(sample.information$Clonality_of_drivers),])] <-  "n.d."

mat.clonal["MYCN",rownames(sample.information[sample.information$MYCN==1 & sample.information$Clonality_of_drivers%in% "clonal",])] <- "clonal"
mat.clonal["MYCN",rownames(sample.information[sample.information$MYCN==1 & sample.information$Clonality_of_drivers %in% "subclonal",])] <- "subclonal"
mat.clonal["MYCN",rownames(sample.information[sample.information$MYCN==1 & is.na(sample.information$Clonality_of_drivers),])] <- "n.d."

mat.clonal["SNCAIP",rownames(sample.information[sample.information$PRDM6==1 & sample.information$Clonality_of_drivers%in% "clonal",])] <- "clonal"
mat.clonal["SNCAIP",rownames(sample.information[sample.information$PRDM6==1 & sample.information$Clonality_of_drivers %in% "subclonal",])] <- "subclonal"
mat.clonal["SNCAIP",rownames(sample.information[sample.information$PRDM6==1 & is.na(sample.information$Clonality_of_drivers),])] <- "n.d."



##############################################################################################################################################
## Oncoprint

column_annotation <- data.frame(Group = sample.information[tumors,]$mnp11,
                                Subgroup = sample.information[tumors,]$mnp12)
rownames(column_annotation) <- tumors


col_clonal = c(time.colors, "cornflowerblue", "orange", "grey", "orange", "cornflowerblue")
names(col_clonal) <- c("LC", "EC", "C", "SC", "n.d.", "subclonal", "clonal")

#palette_check(col2hcl(col), plot = TRUE)


alter_fun_clonal = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = "white", col = NA))
  },
  "EC" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["EC"], col = NA))
  },
  "LC" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["LC"], col = NA))
  },
  "n.d." = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["n.d."], col = NA))
  },
  "C" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["C"], col = NA))
  },
  "SC" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["SC"], col = NA))
  },
  "subclonal" = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col["subclonal"], col = NA))
  },
  "clonal" = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col["clonal"], col = NA))
  }
)


column_title = "OncoPrint for Medulloblastoma"


annotation_colors <- list(Group = cbp1[1:length(unique(sample.information$mnp11))],
                          Subgroup = hcl.colors(n=length(unique(sample.information$mnp12))))

names(annotation_colors$Group) <- unique(sample.information$mnp11)
names(annotation_colors$Subgroup) <- unique(sample.information$mnp12)

pdf(paste0(output.directory, "Oncoprint_clonal.pdf"), width=14, height=9, useDingbats = F)

oncoPrint(mat.clonal,
          bottom_annotation = HeatmapAnnotation(df=column_annotation, col=annotation_colors),
          alter_fun = alter_fun_clonal, col = col_clonal, row_names_gp = gpar(fontsize = 8),
          remove_empty_columns = F, remove_empty_rows = TRUE, show_column_names = T,
          column_title = column_title, show_pct = F)



dev.off()


##############################################################################################################################################
## Oncoprint for G3/G4 tumors only
# only G3/G4 tumors
source("./RScripts/CNV_timing_overview.R")

pdf(paste0(output.directory, "Oncoprint_G3G4.pdf"), width=14, height=9, useDingbats = F)

mat.g3g4 <- rbind(cnv.mat.g3g4, mat[,colnames(cnv.mat.g3g4)])

col <- c(col, time.colors, "ECA/MRCA" = "violet", clonal = "cornflowerblue", "n.d." = "grey", subclonal="orange")

alter_fun$ECA <- function(x, y, w, h) {
  grid.rect(x, y, w, h, 
            gp = gpar(fill = col["ECA"], col = NA))
}
alter_fun$MRCA <- function(x, y, w, h) {
  grid.rect(x, y, w, h, 
            gp = gpar(fill = col["MRCA"], col = NA))
}
alter_fun$`ECA/MRCA` <- function(x, y, w, h) {
  grid.rect(x, y, w, h, 
            gp = gpar(fill = col["ECA/MRCA"], col = NA))
}
alter_fun$`n.d.` <- function(x, y, w, h) {
  grid.rect(x, y, w, h, 
            gp = gpar(fill = col["n.d."], col = NA))
}
alter_fun$subclonal <- function(x, y, w, h) {
  grid.rect(x, y, w, h, 
            gp = gpar(fill = col["subclonal"], col = NA))
}
alter_fun$clonal <- function(x, y, w, h) {
  grid.rect(x, y, w, h, 
            gp = gpar(fill = col["clonal"], col = NA))
}

oncoPrint(mat.g3g4,
          bottom_annotation = HeatmapAnnotation(df=column_annotation[g3g4.tumors,], col=annotation_colors),
          alter_fun = alter_fun, col = col, row_names_gp = gpar(fontsize = 8),
          remove_empty_columns = F, remove_empty_rows = TRUE, show_column_names = T,
          column_title = column_title, show_pct = F)


# add annotation ECA yes/no
load("./RScripts/RData/MRCA_timing.RData")

column_annotation$ECA <- !is.na(mutation.time.eca[colnames(mat),]$Mean)
annotation_colors$ECA <- c("TRUE"="black", "FALSE"="white")
column_annotation$MRCA <- mutation.time.mrca[colnames(mat),]$Mean/3.3/10^3

oncoPrint(mat.g3g4,
          bottom_annotation = HeatmapAnnotation(df=column_annotation[column_annotation$Group %in% c("MB, G3", "MB, G4"),], col=annotation_colors),
          alter_fun = alter_fun, col = col, row_names_gp = gpar(fontsize = 8),
          row_split = ifelse(rownames(mat.g3g4) %in% rownames(cnv.mat.g3g4), "CNV", "small mutation"),
          remove_empty_columns = F, remove_empty_rows = TRUE, show_column_names = T,
          column_title = column_title, show_pct = F, column_split = column_annotation[column_annotation$Group %in% c("MB, G3", "MB, G4"),"ECA"])

dev.off()


##############################################################################################################################################
## Also for the G3/G4 tumors, we plot an oncoprint with the clonality information
mat.g3g4.clonal <- rbind(cnv.mat.g3g4, mat.clonal[,colnames(cnv.mat.g3g4)])

column_annotation <- data.frame(Subgroup = sample.information[g3g4.tumors,]$mnp12)
rownames(column_annotation) <- g3g4.tumors
column_annotation$ECA <- !is.na(mutation.time.eca[colnames(mat.g3g4.clonal),]$Mean)
annotation_colors$ECA <- c("TRUE"="black", "FALSE"="white")
column_annotation$MRCA <- mutation.time.mrca[colnames(mat.g3g4.clonal),]$Mean/3.3/10^3

col_clonal = c(time.colors, "cornflowerblue", "orange", "grey", "orange", "cornflowerblue")
names(col_clonal) <- c("LC", "EC", "C", "SC", "n.d.", "subclonal", "clonal")
col_clonal <- c(col_clonal, time.colors, "ECA/MRCA" = "violet")

#palette_check(col2hcl(col), plot = TRUE)


alter_fun_clonal = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = "white", col = NA))
  },
  "ECA" = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col_clonal["ECA"], col = NA))
  },
  "MRCA" = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col_clonal["MRCA"], col = NA))
  },
  `ECA/MRCA` = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col_clonal["ECA/MRCA"], col = NA))
  },
  "subclonal" = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col_clonal["subclonal"], col = NA))
  },
  "clonal" = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = col_clonal["clonal"], col = NA))
  },
  "n.d." = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["n.d."], col = NA))
  },
  "LC" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["LC"], col = NA))
  },
  "EC" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["EC"], col = NA))
  },
  "C" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["C"], col = NA))
  },
  "SC" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col_clonal["SC"], col = NA))
  }
)


annotation_colors <- list(Subgroup = g34.subgroup.colors)

pdf(paste0(output.directory, "Oncoprint_clonal_G34.pdf"), width=14, height=9, useDingbats = F)

oncoPrint(mat.g3g4.clonal,
          bottom_annotation = HeatmapAnnotation(df=column_annotation, col=annotation_colors),
          alter_fun = alter_fun_clonal, col = col_clonal, row_names_gp = gpar(fontsize = 8),
          row_split = ifelse(rownames(mat.g3g4) %in% rownames(cnv.mat.g3g4), "CNV", "small mutation"),
          remove_empty_columns = F, remove_empty_rows = TRUE, show_column_names = T, show_pct = F, 
          column_split = column_annotation[,"ECA"])

dev.off()

##############################################################################################################################################
## Summary statistics: number of clonal events per tumor; stratified by small mutation or CNV

to.plot <- rbind(data.frame(Clonal = apply(mat.g3g4.clonal[rownames(cnv.mat.g3g4),], 2, function(x){
  sum((grepl("clonal", x) & !grepl("subclonal", x)) | grepl("MRCA", x) | grepl("ECA", x) )}), 
  Subclonal = apply(mat.g3g4.clonal[rownames(cnv.mat.g3g4),], 2, function(x){
    sum(grepl("subclonal", x))}),Type = "CNV",
  NoData = apply(mat.g3g4.clonal[rownames(cnv.mat.g3g4),], 2, function(x){
    sum(grepl("n.d.", x))}),
  ID = colnames(mat.g3g4.clonal)),
  data.frame(Clonal = apply(mat.g3g4.clonal[setdiff(rownames(mat.g3g4.clonal),rownames(cnv.mat.g3g4)),], 
                            2, function(x){
    sum((!grepl("SC", x) & grepl("C", x)) | (!grepl("subclonal", x) & grepl("clonal", x)))}), 
    Subclonal = apply(mat.g3g4.clonal[setdiff(rownames(mat.g3g4.clonal),rownames(cnv.mat.g3g4)),], 
                   2, function(x){
                     sum(grepl("SC", x) | grepl("subclonal", x))}), Type = "small mutation",
    NoData = apply(mat.g3g4.clonal[setdiff(rownames(mat.g3g4.clonal),rownames(cnv.mat.g3g4)),], 
                      2, function(x){
                        sum((!grepl("C", x) & grepl("n.d.", x)) | (!grepl("clonal", x) & grepl("n.d.", x)))}),
    ID = colnames(mat.g3g4.clonal)))

to.plot <- reshape2::melt(to.plot, id.vars=c("Type", "ID"))
  
pdf(paste0(output.directory, "Clonal_drivers_G3G4.pdf"), width=5, height=4, useDingbats = F)

ggplot(to.plot, aes(x=variable, y=value, col=Type)) + geom_boxplot() + 
  scale_color_manual(values = c(CNV = "orange", `small mutation` = "violet")) +
  scale_y_continuous("Count per tumor") + stat_compare_means(paired = T, method = "wilcox.test", size=2)

# do stratified by ECA yes/no
to.plot$ECA <- ifelse(to.plot$ID %in% rownames(mutation.time.eca[!is.na(mutation.time.eca$Mean),]),
                      "ECA", "no ECA")

ggplot(to.plot, aes(x=variable, y=value, col=Type)) + geom_boxplot() + 
  scale_color_manual(values = c(CNV = "orange", `small mutation` = "violet")) +
  scale_y_continuous("Count per tumor") + stat_compare_means(paired = T, method = "wilcox.test", size=2) +
  facet_wrap(~ECA, nrow=2, scales = "free") 

# number of tumors with 1, 2, 3, ... clonal CNV

to.plot.2 <- rbind(data.frame(Count = seq(0, 10), N_tumors = sapply(seq(0, 10), function(x){
  sum(to.plot[to.plot$variable=="Clonal" & to.plot$Type=="CNV",]$value == x)}),
  Type = "CNV"
), data.frame(Count = seq(0, 10), N_tumors = sapply(seq(0, 10), function(x){
  sum(to.plot[to.plot$variable=="Clonal" & to.plot$Type=="small mutation",]$value == x)}),
  Type = "small mutation"
))

ggplot(to.plot.2, aes(x=Count, y=N_tumors, col=Type)) + geom_point() + geom_line() +
  scale_color_manual(values = c(CNV = "orange", `small mutation` = "violet")) +
  scale_x_continuous("Number of clonal drivers") + scale_y_continuous("Number of tumors")

dev.off()
