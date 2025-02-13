##############################################################################################################################################
### Analyze gains and losses across the cohort
##############################################################################################################################################
library(GenomicRanges)
library(ggbio)
## homozygous-deletions due to artefacts need to be filtered
homdel.artifacts <- read.delim("refdata/artifact.homoDels.potentialArtifacts.txt", sep="\t", header = T)

resolution <- 10^6
##############################################################################################################################################
## collect for each genomic position the number of tumors harboring a gain or loss at this position

# for plotting
gains <- GRanges()
losses <- gains

# for enrichment test
cnv.summary <- data.frame(Chr = c(), Arm = c(), CNV = c(), ID = c(), CN = c(), Clonality = c())

for(i in tumors){
  ## ACEseq results files are named ...comb_pro...
  aceseq <- list.files(paste0("./data/",  i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  
  purity. <- Extract.purity.ploidy.from.ACEseq(aceseq)$purity
  ploidy. <- Extract.purity.ploidy.from.ACEseq(aceseq)$ploidy
  
  ploidy. <- round(as.numeric(ploidy.))
  
  aceseq <-  read.delim(paste0("./data/",  i, "/", cnv.directory, "/", aceseq),  sep="\t", stringsAsFactors = F)
  
  aceseq$X.chromosome <- replace(aceseq$X.chromosome, aceseq$X.chromosome=="X", 23)
  aceseq$X.chromosome <- replace(aceseq$X.chromosome, aceseq$X.chromosome=="Y", 24)
  aceseq$X.chromosome <- as.numeric(aceseq$X.chromosome)
  aceseq <- aceseq[order(aceseq$X.chromosome),]
  
  
  ## filter artefacts
  for(j in 1:nrow(homdel.artifacts[1:12,])){
    tmp <- which(aceseq$X.chromosome==homdel.artifacts[j,1] & aceseq$start >= homdel.artifacts[j,2] & aceseq$end <= homdel.artifacts[j,3])
    if(length(tmp)>0){
      aceseq <- aceseq[-tmp,,drop=F]
    }
  }
  
  ## store the copy number info; choose 10kb resolution
  aceseq <-  GRanges(seqnames=aceseq[,1],
                     ranges = IRanges(ceiling(aceseq[,2]/resolution), floor(aceseq[,3]/resolution)),
                     TCN = aceseq[,"TCN"],
                     tumor = i,
                     gain.loss=aceseq[,ncol(aceseq)],
                     cov.ratio = aceseq[,"tcnMeanRaw"])
  aceseq <- aceseq[width(aceseq) > 1,]

  if(any(max(coverage(aceseq)) > 1)){
    stop("DUPLICATE alarm!!!!")
  }
  
  gains.this.tumor <- aceseq[(aceseq$gain.loss %in% c("gain", "LOHgain", "DUP", "DUP;LOH", "TCNneutral;DUP", "TCNneutral;LOH")) |
    (!is.na(aceseq$gain.loss) & aceseq$cov.ratio > 1.1 & aceseq$TCN %in% "sub"),]
  gains <- c(gains, gains.this.tumor)
  
  losses.this.tumor <-  aceseq[(aceseq$gain.loss %in% c("homozygousDel", "loss", "LOH", "cnLOH", "DEL", "HomoDel", "DEL;LOH")) |
                                 (is.na(aceseq$gain.loss) & aceseq$cov.ratio < 0.9 & aceseq$TCN %in%"sub"),]
  losses <- c(losses, losses.this.tumor)
  
  if(length(gains.this.tumor) > 0){
    cnv.summary <- rbind(cnv.summary,
                         data.frame(Chr = unique(seqnames(gains.this.tumor)), 
                                    Arm = sapply(unique(seqnames(gains.this.tumor)), function(j){
                                      begin <- min(start(gains.this.tumor[seqnames(gains.this.tumor)==j,]))
                                      ending <- max(end(gains.this.tumor[seqnames(gains.this.tumor)==j,]))
                                      total.length <- sum(width(gains.this.tumor[seqnames(gains.this.tumor)==j,]))
                                      if(total.length < 0.25 * centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution){
                                        "small level"
                                      }else if(begin <= centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution &
                                               ending >= centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution){
                                        "both arms"
                                      }else if (begin <= centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution){
                                        "p arm"
                                      }else{
                                        "q arm"
                                      }
                                    }), CNV = c("gain"), ID = i, 
                                    CN = sapply(unique(seqnames(gains.this.tumor)), function(j){
                                      mean(as.numeric(gains.this.tumor[seqnames(gains.this.tumor)==j,]$TCN), na.rm=T)
                                    }),
                                    Seglength = sapply(unique(seqnames(gains.this.tumor)), function(j){
                                      total.length <- sum(width(gains.this.tumor[seqnames(gains.this.tumor)==j,]))
                                      total.length * resolution
                                    }),
                                    Clonality = sapply(unique(seqnames(gains.this.tumor)), function(j){
                                      ifelse(any(gains.this.tumor[seqnames(gains.this.tumor)==j,]$TCN != "sub"), "clonal", "subclonal")
                                    })))
  }
  if(length(losses.this.tumor) > 0){
    cnv.summary <- rbind(cnv.summary,
                         data.frame(Chr = unique(seqnames(losses.this.tumor)), 
                                    Arm = sapply(unique(seqnames(losses.this.tumor)), function(j){
                                      begin <- min(start(losses.this.tumor[seqnames(losses.this.tumor)==j,]))
                                      ending <- max(end(losses.this.tumor[seqnames(losses.this.tumor)==j,]))
                                      total.length <- sum(width(losses.this.tumor[seqnames(losses.this.tumor)==j,]))
                                      if(total.length < 0.25 * centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution){
                                        "small level"
                                      }else if(begin <= centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution &
                                               ending >= centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution){
                                        "both arms"
                                      }else if (begin <= centromers[centromers$chrom==paste0("chr", j),]$chromStart/resolution){
                                        "p arm"
                                      }else{
                                        "q arm"
                                      }
                                    }), CNV = c("loss"), ID = i, 
                                    CN = sapply(unique(seqnames(losses.this.tumor)), function(j){
                                      mean(as.numeric(losses.this.tumor[seqnames(losses.this.tumor)==j,]$TCN), na.rm=T)
                                    }),
                                    Seglength = sapply(unique(seqnames(losses.this.tumor)), function(j){
                                      total.length <- sum(width(losses.this.tumor[seqnames(losses.this.tumor)==j,]))
                                      total.length * resolution
                                    }),
                                    Clonality = sapply(unique(seqnames(losses.this.tumor)), function(j){
                                      ifelse(any(losses.this.tumor[seqnames(losses.this.tumor)==j,]$TCN != "sub"), "clonal", "subclonal")
                                    })))
  }
 
}

save(cnv.summary, file="./RScripts/RData/CNV_per_tumor_res_1Mb.RData")

# for the significance test we only care for arm-level or whole chromosome gains
cnv.summary <- cnv.summary[cnv.summary$Arm != "small level",]

gains$y<- 1
gains$Type <- "Gain"
losses$y <- -1
losses$Type <- "Loss"

data(ideoCyto, package = "biovizBase") ## data for ideogram
## use the lengths of the chromosomes from reference data
seqlengths(gains)[1:24] <- as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr",c(1:22, "X", "Y"))])/resolution


## do separate plots for all tumors and per subtype
p <- list()
for(subtype in c("ALL", "MB, SHH CHL AD" ,"MB, SHH INF" , "MB, WNT", "G3/G4")){
  
  if(subtype=="ALL"){
    subset <- sample.information$MRCA_ID
  }else if(subtype == "G3/G4"){
    subset <- rownames(sample.information[sample.information$mnp11 %in% c("MB, G3", "MB, G4"),])
  }else{
    subset <- sample.information$MRCA_ID[sample.information$mnp11 %in% subtype]
  }
  
  ## sum up the number of gains of individual tumors
  gains.cov <- coverage(gains[gains$tumor %in% subset,])
  start <- c()
  end <- c()
  ## extract the positions (start end of each gain)
  for(i in 1:length(gains.cov)){
    position <- unlist(cumsum(unlist(runLength(gains.cov[[i]]))))
    if(length(position)==0){next}
    start <- c(start,1,position[-length(position)])
    end <- c(end,position)
  }
  
  ## convert to GRanges object
  gain.granges <- data.frame(Chr=rep(names(gains.cov), sapply(gains.cov, function(x) length(runLength(x)))), Strand = rep("*", length(unlist(runLength(gains.cov)))),
                             Start = resolution*start, End = resolution*end,
                             Coverage = as.numeric(unlist(runValue(gains.cov))))
  
  gain.granges <- GRanges(seqnames = gain.granges$Chr, ranges = IRanges(start = gain.granges$Start, end = gain.granges$End), strand = 
                            gain.granges$Strand, Coverage = gain.granges$Coverage, Type = "gain")
  
  gain.granges <- keepSeqlevels(gain.granges, as.character(1:24)[as.character(1:24) %in% seqlevels(gain.granges)])
  ## use reference lengths of chromosomes
  seqlevels(gain.granges) <-  as.character(c(1:24))
  seqlengths(gain.granges) <- as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr", c(1:22, "X", "Y"))])

  ## same for losses
  losses.cov <- coverage(losses[losses$tumor %in% subset,])
  start <- c()
  end <- c()
  for(i in 1:length(losses.cov)){
    position <- unlist(cumsum(unlist(runLength(losses.cov[[i]]))))
    if(length(position)==0){next}
    start <- c(start,1,position[-length(position)])
    end <- c(end,position)
  }
  
  
  loss.granges <- data.frame(Chr=rep(names(losses.cov), sapply(losses.cov, function(x) length(runLength(x)))), Strand = rep("*", length(unlist(runLength(losses.cov)))),
                             Start = resolution*start, End = resolution*end,
                             Coverage = as.numeric(unlist(runValue(losses.cov))))
  
  loss.granges <- GRanges(seqnames = loss.granges$Chr, ranges = IRanges(start = loss.granges$Start, end = loss.granges$End), strand = 
                            loss.granges$Strand, Coverage = -loss.granges$Coverage, Type = "loss")
  
  loss.granges <- keepSeqlevels(loss.granges, as.character(1:24)[as.character(1:24) %in% seqlevels(loss.granges)])
  
  seqlengths(loss.granges)<- as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr", c(1:22, "X", "Y"))[as.character(1:24) %in% seqlevels(loss.granges)]])
  
  
  # are particular chromosomes gained more often than expected?
  gain.probs <- p.adjust(sapply(unique(seqnames(gain.granges)), function(j){
    
    pbinom(sum(as.character(cnv.summary$Chr) == as.character(j) & cnv.summary$CNV=="gain" & cnv.summary$ID %in% subset), 
            size = sum(cnv.summary$CNV=="gain" & cnv.summary$ID %in% subset),
           prob = 1/24, lower.tail = F)
    
    # cont.mat <- matrix(c(sum(cnv.summary$Chr == j & cnv.summary$CNV=="gain" & cnv.summary$ID %in% subset), 
    #                      sum(cnv.summary$Chr != j & cnv.summary$CNV=="gain" & cnv.summary$ID %in% subset),
    #                      1, 23),
    #                    nrow = 2, byrow = T)
    # 
    # fisher.test(cont.mat, alternative = "greater")$p.value
   # pbinom(table(seqnames(gain.granges))[j], size = length(gain.granges), 
    #       prob = seqlengths(gain.granges)[j]/sum(seqlengths(gain.granges)), lower.tail = F)
  }))
  names(gain.probs) <- unique(seqnames(gain.granges))
  lost.probs <- p.adjust(sapply(unique(seqnames(loss.granges)), function(j){
    
    pbinom(sum(as.character(cnv.summary$Chr) == as.character(j) & cnv.summary$CNV=="loss" & cnv.summary$ID %in% subset), 
           size = sum(cnv.summary$CNV=="loss" & cnv.summary$ID %in% subset),
           prob = 1/24, lower.tail = F)
    
    #pbinom(table(seqnames(loss.granges))[j], size = length(loss.granges), 
    #       prob = seqlengths(loss.granges)[j]/sum(seqlengths(loss.granges)), lower.tail = F)
  }))
  names(lost.probs) <- unique(seqnames(loss.granges))
  # highlight significant gains and losses
  
  if(any(lost.probs < 0.05) | any(gain.probs < 0.05)){
    if(any(gain.probs < 0.05)){
      highlights <- unique(GRanges(c(names(gain.probs[gain.probs < 0.05])),
                                   IRanges(rep(1, sum(gain.probs < 0.05) ),
                                           width = as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr",c(
                                             names(gain.probs[gain.probs < 0.05])
                                           ))]))))
      
      overlaps <- findOverlaps(query = gain.granges, subject = highlights)
      mcols(gain.granges)$sig <- F
      mcols(gain.granges)$sig[queryHits(overlaps)] <- T
      
    }else{
      mcols(gain.granges)$sig <- F
    }
    
    if(any(lost.probs < 0.05)){
      highlights <- unique(GRanges(c(names(lost.probs[lost.probs < 0.05])),
                                   IRanges(rep(1, sum(lost.probs < 0.05)),
                                           width = as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr",c(
                                             names(lost.probs[lost.probs < 0.05])
                                           ))]))))
      
      overlaps <- findOverlaps(query = loss.granges, subject = highlights)
      mcols(loss.granges)$sig <- F
      mcols(loss.granges)$sig[queryHits(overlaps)] <- T
      
    }else{
      mcols(loss.granges)$sig <- F
      
    }
    
    p[[subtype]] <- plotGrandLinear(c(gain.granges, loss.granges), coord="genome", geom="bar", 
                                    ymax=c(gain.granges, loss.granges)$Coverage/length(subset)*100, 
                                    aes(y=Coverage/length(subset)*100, fill = sig, col=sig)) +  
      theme(panel.background = element_rect(fill = "white", color="black"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + geom_hline(yintercept = 0) +
      scale_color_manual(values = c("TRUE" = "#DA1F28", "FALSE"= "grey")) +scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE"= "grey"))+
      scale_y_continuous("% Tumors") + ggtitle(subtype)
    
  }else{
    p[[subtype]] <- plotGrandLinear(c(gain.granges, loss.granges), coord="genome", geom="bar", 
                                    ymax=c(gain.granges, loss.granges)$Coverage/length(subset)*100, 
                                    aes(y=Coverage/length(subset)*100), col="grey", fill = "grey") +  
      theme(panel.background = element_rect(fill = "white", color="black"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + geom_hline(yintercept = 0) +
      scale_y_continuous("% Tumors") + ggtitle(subtype)
  }


 
  
}

pdf("./Plots/CNVs_per_cohort.pdf", width = 8, height = 3)
p[[1]]
p[[2]]
p[[3]]
p[[4]]
p[[5]]
dev.off()
