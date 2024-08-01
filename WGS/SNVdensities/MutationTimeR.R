##############################################################################################################################################
## source the settings
source("./RScripts/Settings.R")
library(MutationTimeR)
library(VariantAnnotation)


##############################################################################################################################################
## load purity/ploidy per case
load(paste0(rdata.directory, "Purity_ploidy.RData"))

# iterate through all cases and compute mutation densities with MutationTimeR
for(i in tumors){
  
  print(i)
  purity. <- purity[i]
  
  if(!dir.exists(paste(data.directory, i, "MutationTimeR", sep="/"))){
    dir.create(paste(data.directory, i, "MutationTimeR", sep="/"))
  }
  
  vcf.file <- list.files(paste(data.directory, i, snv.directory, sep = "/"), pattern="somatic_snvs_conf_8_to_10", full.names = T)[1]
  if(!file.exists(vcf.file)){next}
  
  vcf <- readVcf(vcf.file) # Point mutations, needs `geno` entries `AD` and `DP` or `info` columns t_alt_count t_ref_count.
  
  info(vcf)$t_ref_count <- sapply(info(vcf)$DP4, function(x){
    sum(x[c(1,2)])
  })
  info(vcf)$t_alt_count <- sapply(info(vcf)$DP4, function(x){
    sum(x[c(3,4)])
  })
  
  coverage.file <- list.files(paste0(data.directory, i, "/ACEseq"), pattern="comb_pro_extra")[1]
  coverage <- read.delim(paste0(data.directory, i, "/ACEseq/",coverage.file))
  coverage <- coverage[!is.na(coverage$A) & !is.na(coverage$B) & coverage$A!="sub" & coverage$B!="sub",]
  coverage$A <- as.numeric(coverage$A)
  coverage$B <- as.numeric(coverage$B)
  coverage <- coverage[coverage$A!=0 ,]
  ## merge consecutive segments with equal copy number
  coverage.smoothed <- data.frame()
  for(j in unique(coverage$X.chromosome)){
    tmp <- coverage[coverage$X.chromosome==j,]
    output.tmp <- tmp[1,]
    k <- 2
    while(k <= nrow(tmp)){
      if(tmp[k,"A"]==output.tmp[nrow(output.tmp),"A"] &
         tmp[k,"B"]==output.tmp[nrow(output.tmp), "B"]){
        output.tmp[nrow(output.tmp),"end"] <- tmp[k,"end"]
      }else{
        output.tmp <- rbind(output.tmp, tmp[k,])
      }
      k <- k+1
    }
    coverage.smoothed <- rbind(coverage.smoothed, output.tmp)
  }
  coverage <- coverage.smoothed
  coverage$start <- coverage$start+1
  
  bb <- GRanges(seqnames = coverage$X.chromosome, IRanges(start=coverage$start, end=coverage$end),
                major_cn = coverage$A, minor_cn = coverage$B, clonal_frequency = purity.)
  
  mt <- mutationTime(vcf, bb,  n.boot=10)
  
  vcf <- addMutTime(vcf, mt$V)
  mcols(bb) <- cbind(mcols(bb),mt$T)
  
  ## extract early clonal mutations per copy
  early.clonals.per.segment <- countOverlaps(bb, vcf[!is.na(mt$V$CLS) & 
                                                       mt$V$CLS=="clonal [early]",])
  ## extract mutations that are clonal on 1 copy per segment
  single.copy.clonals.per.segment <- countOverlaps(bb, vcf[!is.na(mt$V$CLS) & 
                                                       mt$V$CLS %in% c("clonal [late]", "clonal [NA]"),])
  
  subclonals.per.segment <- countOverlaps(bb, vcf[!is.na(mt$V$CLS) & 
                                                             mt$V$CLS %in% c("subclonal"),])
  
  mcols(bb) <- cbind(mcols(bb), EarlyClonals=early.clonals.per.segment)
  mcols(bb) <- cbind(mcols(bb), SingleCopyClonals=single.copy.clonals.per.segment)
  
  ## save the vcf including the clonality state
  writeVcf(vcf, filename=paste0(data.directory, i, "/MutationTimeR/", i, "_clonality.vcf"))
  ## also write the bb output
  
  df <- data.frame(seqnames=seqnames(bb),
                   starts=start(bb)-1,
                   ends=end(bb),
                   names=c(rep(".", length(bb))),
                   major_cn = elementMetadata(bb)$major_cn,
                   minor_cn = elementMetadata(bb)$minor_cn,
                   clonal_frequency = elementMetadata(bb)$clonal_frequency,
                   type = elementMetadata(bb)$type,
                   time = elementMetadata(bb)$time,
                   time.lo = elementMetadata(bb)$time.lo,
                   time.up = elementMetadata(bb)$time.up,
                   time.2nd = elementMetadata(bb)$time.2nd,
                   time.2nd.lo = elementMetadata(bb)$time.2nd.lo,
                   time.2nd.up = elementMetadata(bb)$time.2nd.up,
                   time.star = elementMetadata(bb)$time.star,
                   n.snv_mnv = elementMetadata(bb)$n.snv_mnv,
                   n.early.clonals = elementMetadata(bb)$EarlyClonals,
                   n.single.copy.clonals = elementMetadata(bb)$SingleCopyClonals)
  
  write.table(df, file=paste0(data.directory, i, "/MutationTimeR/", i, "_MutationTimeR_w_clonality.bed"), quote=F,
              sep="\t", row.names=F, col.names=T)
  
  pdf(paste0(data.directory, i, "/MutationTimeR/", i, "_MutationTimeR.pdf"))
  p1 <- plotSample(vcf,bb)
  print(p)
  dev.off()
  save(p1, file = paste0(data.directory, i, "/MutationTimeR/", i, "_MutationTimeRplot.rds"))
  
  bb <- bb[width(bb)>=10^7,]
  bb <- bb[!seqnames(bb) %in% c("X", "Y") & (bb$major_cn + bb$minor_cn) <=4,]
  bb <- bb[order(as.numeric(as.character(seqnames(bb)))),]
  if(length(bb)>0){
    to.plot <- data.frame(Time = bb$time, Time.lo = bb$time.lo, Time.up = bb$time.up,
                          SegNo = 1:length(bb), SegId = paste(seqnames(bb), paste0(bb$major_cn + bb$minor_cn, "N"),
                                                              c("I", "II", "III", "IV")[bb$major_cn], sep="_"))
    to.plot <- to.plot[!is.na(to.plot$Time) ,]
    to.plot$SegId <- factor(to.plot$SegId, levels=unique(to.plot$SegId))
    
  }else{
    to.plot <- data.frame()
  }
  
  if(nrow(to.plot)>0){
    to.plot$SegNo <- 1:nrow(to.plot)
    
    p2 <- ggplot(to.plot, aes(x = Time, xmin = Time.lo, xmax = Time.up, y = SegNo, col = SegId)) + 
      geom_point() + geom_errorbarh(height = 0) + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), 
                                                        axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
      geom_vline(xintercept = 1, col = time.colors["MRCA"], linetype = 2) + 
      scale_x_continuous(name = "Mutation time relative to MRCA", limits = c(0, 2)) +
      guides(col=guide_legend(ncol=2))
    
    write.table(to.plot, file = paste("./mbOncoAberApp/data/WGS/MutationTimeR/Timeline_MutationTimeR_", i, ".txt", sep=""), quote = F, row.names = F, sep="\t")
    
  }else{
    p2 <- ggplot(NULL) + scale_x_continuous(name = "Mutation time relative to MRCA", limits = c(0,1))
  }
  
  save(p1, p2, file = paste0(data.directory, i, "/MutationTimeR/", i, "_MutationTimeRplot.rds"))
  
                                                          
                                                          
}


##############################################################################################################################################
## Compare the mutation densities between our method and MutationTimeR segment-wise
load(paste0(rdata.directory, "MRCA_timing.RData"))

p1 <- list()
p2 <- list()
cor_coeff <- c()
mean_dev <- c()
percent_no_overlap <- c()

for(i in tumors){
  
  lachesis <- mutation.time[[i]]
  mutationtimer <- read.delim(paste0(data.directory, i, "/MutationTimeR/", i, "_MutationTimeR_w_clonality.bed"))
  mutationtimer$copy_number <- apply(mutationtimer, 1, function(x){
    cn <- as.numeric(x["major_cn"]) + as.numeric(x["minor_cn"])
    if(cn == 4){
      "tetrasomic"
    }else if(cn == 3){
      "trisomic"
    }else if(cn == 2){
      "disomic"
    }else if(cn == 1){
      "monosomic"
    }else{
      "out_of_range"
    }
  })
  mutationtimer$B_allele <- apply(mutationtimer, 1, function(x){
    a <- as.numeric(x["major_cn"])
    if(a == 1){
      "I"
    }else if(a == 2){
      "II"
    }else if(a == 3){
      "III"
    }else if(a == 4){
      "IV"
    }else{
      "out_of_range"
    }
  })
  mutationtimer$Segment <- paste("chr", mutationtimer$seqnames, mutationtimer$copy_number,
                                 mutationtimer$B_allele, sep="_")
  
  mutationtimer <- mutationtimer[mutationtimer$copy_number!="out_of_range" &
                                   (mutationtimer$ends - mutationtimer$starts) >= 10^7,]
  
  to.plot <- data.frame(Segment = intersect(lachesis$Segment, mutationtimer$Segment))
  to.plot$Lachesis_min <- sapply(to.plot$Segment, function(x){
    lachesis[lachesis$Segment==x,]$Min/mrca.eca[[i]]$mutation.time.mrca.upper
  })
  to.plot$Lachesis_max <- sapply(to.plot$Segment, function(x){
    lachesis[lachesis$Segment==x,]$Max/mrca.eca[[i]]$mutation.time.mrca.lower
  })
  to.plot$Lachesis_mean <- sapply(to.plot$Segment, function(x){
    lachesis[lachesis$Segment==x,]$Mean/mrca.eca[[i]]$mutation.time.mrca
  })
  to.plot$MutationTimer_min <- sapply(to.plot$Segment, function(x){
    mean(mutationtimer[mutationtimer$Segment==x,]$time.lo)
  })
  to.plot$MutationTimer_max <- sapply(to.plot$Segment, function(x){
    mean(mutationtimer[mutationtimer$Segment==x,]$time.up)
  })
  to.plot$MutationTimer_mean <- sapply(to.plot$Segment, function(x){
    mean(mutationtimer[mutationtimer$Segment==x,]$time)
  })

  to.plot <- to.plot[!is.na(to.plot$MutationTimer_mean),,drop=F]

  
  if(nrow(to.plot)==0){next}
  
  cor_coeff[i] <- cor(to.plot$Lachesis_mean, to.plot$MutationTimer_mean)
  mean_dev[i] <- mean(to.plot$Lachesis_mean - to.plot$MutationTimer_mean)
  percent_no_overlap[i] <- sum((to.plot$Lachesis_max < to.plot$MutationTimer_min) |
                                 to.plot$MutationTimer_max < to.plot$Lachesis_min)/nrow(to.plot) * 100
  
  
  p1[[i]] <- ggplot(to.plot, aes(x = Lachesis_mean, xmin = Lachesis_min, xmax = Lachesis_max,
                      y = MutationTimer_mean, ymin = MutationTimer_min, ymax = MutationTimer_max)) +
    geom_errorbar(width = 0) + geom_errorbarh(height = 0) + geom_point() + 
    geom_abline(slope = 1, intercept = 0, linetype = 2) + theme(aspect.ratio = 1) +
    ggtitle(paste(i, "Pearson's r =", round(cor_coeff[i] , digits = 2))) +
    scale_x_continuous(limits = c(0, max(c(to.plot$Lachesis_max, to.plot$MutationTimer_max))))+
    scale_y_continuous(limits = c(0, max(c(to.plot$Lachesis_max, to.plot$MutationTimer_max))))
  
  p2[[i]] <- ggplot(to.plot, aes(x = Lachesis_mean, xmin = Lachesis_min, xmax = Lachesis_max,
                      y = MutationTimer_mean, ymin = MutationTimer_min, ymax = MutationTimer_max)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) + theme(aspect.ratio = 1) +
    ggtitle(paste(i, "Pearson's r =", round(cor_coeff[i] , digits = 2)))+
    scale_x_continuous(limits = c(0, max(c(to.plot$Lachesis_mean, to.plot$MutationTimer_mean))))+
    scale_y_continuous(limits = c(0, max(c(to.plot$Lachesis_mean, to.plot$MutationTimer_mean))))

}

save(p1, p2, cor_coeff, mean_dev, percent_no_overlap, file = paste(rdata.directory, "Compare_to_MutationTimeR.RData"))

hist(cor_coeff, breaks=20)
abline(v=median(cor_coeff, na.rm=T), col = "red")

pdf(paste(output.directory, "Compare_densities_to_MutationTimeR.pdf"), width = 3, height = 2.5)

ggplot(data.frame(Percent = 100-percent_no_overlap,
                  Sample = names(percent_no_overlap)), aes(x = Percent)) + 
  geom_histogram() + scale_x_continuous("Overlapping mutation densities (%)") +
  scale_y_continuous("Number of tumors")

dev.off()

