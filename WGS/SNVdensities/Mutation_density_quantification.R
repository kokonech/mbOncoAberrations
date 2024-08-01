##############################################################################################################################################
## Quantify the densitiy of amplified and non-amplified clonal mutations
##############################################################################################################################################
#source(paste0(custom.script.directory, "Adjust_purity.R"))
load(paste0(rdata.directory, "/Purity_ploidy.RData"))
library(ggbeeswarm)
# no info on false positive clonals, hence set to 0.
clonal.mutations.false.positives <- c(0,0)
clonal.mutations.all <- c(1,1)
##############################################################################################################################################
## Count the number of clonal mutations stratified by copy number and autosome per tumor

clonal.mutations.all.tumors <- list()

for(i in tumors){

  print(i)
  
  ## Find ACEseq file:
  aceseq <- list.files(paste0(data.directory, i, "/ACEseq"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  aceseq <- paste0(data.directory, "/", i, "/", cnv.directory, "/", aceseq)

  ## Find mutation file
  mutations <- list.files(paste0(data.directory, "/", i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations <- paste0(data.directory, "/", i, "/", snv.directory, "/", mutations)

  clonal.mutations.all.tumors[[i]] <- count.clonal.mutations(aceseq, mutations, chromosomes = c(1:22), purity.=purity[i], ploidy.=ploidy[i])

}

B.allele.indicator <- clonal.mutations.all.tumors[[1]]$B.allele.indicator
copy.number.indicator <- clonal.mutations.all.tumors[[1]]$ copy.number.indicator

save(clonal.mutations.all.tumors, 
     file=paste0(rdata.directory, "Clonal_mutations_different_CNs.RData"))

load(paste0(rdata.directory, "Clonal_mutations_different_CNs.RData"))


##############################################################################################################################################
## Compute the normalized mutational density at each genomic fragment

## We model mutation accumulation as a Poisson process. The number of mutations acquired on a piece of DNA depends on the per-base mutation rate (at this time) and the length of the piece
## using a chisquare-approximation, we obtain confidence intervals for the mutation time.

## store results of all tumors in a list
mutation.time <- list()

for(i in tumors){

  mutation.time[[i]] <- data.frame()

  ## iterate through the chromosomes
  for(j in 1:22){

    ## Get normalized mutation counts per copy and segment on this chromosome
    tmp.mut.count <- Mutation.time.converter(clonal.mutations.all.tumors[[i]]$clonal.mutation.matrix[,j])

    tmp.genome.length <- clonal.mutations.all.tumors[[i]]$segment.length.matrix[,j]
    ## Restrict analysis to fragments > 10^7 bp
    tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>10^7)]
    tmp.genome.length <- tmp.genome.length[tmp.genome.length>10^7]

    if(length(tmp.mut.count)==0){next}

    ## Convert to mutations per haploid genome and store output

    mutation.time[[i]] = rbind(mutation.time[[i]],
                          data.frame(Mean =  tmp.mut.count*3.3*10^9/tmp.genome.length,
                                     Min = 0.5*qchisq(0.025, tmp.mut.count*2)/tmp.genome.length*3.3*10^9,
                                     Max = 0.5*qchisq(0.975, (tmp.mut.count*2+2))/tmp.genome.length*3.3*10^9,
                                     Segment = paste("chr", j, names(tmp.mut.count), sep="_")))

  }

  mutation.time[[i]]$Segment <- factor(mutation.time[[i]]$Segment, levels = mutation.time[[i]]$Segment[order(mutation.time[[i]]$Mean)])

  ## Visualize the mutation density per segment
  p <- ggplot(mutation.time[[i]],
              aes(x=Segment, y=Mean, ymin=Min, ymax=Max)) + geom_col() + geom_errorbar(aes(x=Segment, ymin=Min, ymax=Max)) +
    scale_y_continuous(limits=c(0, max(mutation.time[[i]]$Max)), name = "# Mutations per haploid genome") +
    theme(axis.text.x = element_text(angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

  print(p)

}

##############################################################################################################################################
## Estimate the mutational density at ECA and MRCA testing with negative binomial distributions

##############################################################################################################################################
#### Obtain estimates for the mutational density at the tumor's MRCA and ECA and get the segments conforming to either time point with the function
#### MRCA.ECA.quantification

mutation.time.mrca <- c()
earliest.mutation.time <- c()
mutation.time.eca <- c()

mrca.eca <- list()

for(i in tumors){
  print(i)

  ## build input matrix

  mrca.eca[[i]] <- MRCA.ECA.quantification(clonal.mutations.all.tumors[[i]]$clonal.mutation.matrix,
                                      clonal.mutations.all.tumors[[i]]$segment.length.matrix, )

  mutation.time.mrca <- rbind(mutation.time.mrca, data.frame(Mean = mrca.eca[[i]]$mutation.time.mrca,
                                                             Min = mrca.eca[[i]]$mutation.time.mrca.lower,
                                                             Max = mrca.eca[[i]]$mutation.time.mrca.upper,
                                                             Sample = i))

  if(length(mrca.eca[[i]]$earliest.mutation.time)>0){
    earliest.mutation.time <- rbind(earliest.mutation.time, data.frame(Mean = mrca.eca[[i]]$earliest.mutation.time,
                                                               Min = mrca.eca[[i]]$earliest.mutation.time.lower,
                                                               Max = mrca.eca[[i]]$earliest.mutation.time.upper,
                                                               Sample = i))
  }

  mutation.time.eca <- rbind(mutation.time.eca, data.frame(Mean = mrca.eca[[i]]$mutation.time.eca,
                                                             Min = mrca.eca[[i]]$mutation.time.eca.lower,
                                                             Max = mrca.eca[[i]]$mutation.time.eca.upper,
                                                             Sample = i))

}

rownames(mutation.time.mrca) <- mutation.time.mrca$Sample
rownames(mutation.time.eca) <- mutation.time.eca$Sample
rownames(earliest.mutation.time) <- earliest.mutation.time$Sample

save(mutation.time.mrca, mutation.time.eca, earliest.mutation.time, mrca.eca, mutation.time,
     file=paste0(rdata.directory, "MRCA_timing.RData"))


median(mutation.time.mrca$Mean)/3.3/10^3
median(mutation.time.eca$Mean, na.rm=T)/3.3/10^3
quantile(mutation.time.mrca$Mean)/3.3/10^3
quantile(mutation.time.mrca$Mean[!is.na(mutation.time.eca$Mean)])/3.3/10^3
quantile(mutation.time.eca$Mean, na.rm=T)/3.3/10^3

mutation.time.eca[earliest.mutation.time$Sample,]$Mean <- earliest.mutation.time$Mean
mutation.time.eca[earliest.mutation.time$Sample,]$Min <- earliest.mutation.time$Min
mutation.time.eca[earliest.mutation.time$Sample,]$Max <- earliest.mutation.time$Max

quantile(mutation.time.eca$Mean, na.rm=T)/3.3/10^3

# compute some stats about the ECA

pdf(paste0(output.directory, "ECA_stats.pdf"), width=2, height=2.5, useDingbats = F)
## How many group3/4 tumors show evidence for an ECA?
nrow(mutation.time.eca[!is.na(mutation.time.eca$Mean) & rownames(mutation.time.eca) %in% group34.tumors,])
## How many CNVs map to ECA?
tmp <- unlist(lapply(mrca.eca[group34.tumors], function(x){length(x$gains.uniquely.mapped.to.eca)})) 
mean(tmp[tmp!=0])
range(tmp[tmp!=0])

ggplot(data.frame(N_early = tmp[tmp!=0]), aes(x = "G3/4", y = N_early)) + geom_boxplot() + geom_beeswarm() +
  scale_x_discrete(name = "") + scale_y_continuous(name = "Number of early CNVs") + expand_limits(y=0)

## How many group3/4 tumors show evidence for an ECA but not all CNVs agree to ECA or MRCA?
tmp <- unlist(lapply(mrca.eca[group34.tumors], function(x){(length(x$gains.not.maping.to.eca.or.mrca) +
                                                              length(x$gains.at.earliest.time)+
                                                              length(x$gains.not.mapping.to.earliest.time))/
    (length(x$gains.not.maping.to.eca.or.mrca)+
       length(x$gains.uniquely.mapped.to.eca) +
       length(x$gains.at.earliest.time)+
       length(x$gains.not.mapping.to.earliest.time))}))
tmp <- tmp[!is.na(mutation.time.eca[names(tmp),"Mean"])]
table(tmp)

ggplot(data.frame(Percent_in_ECA = (1-tmp)*100), aes(x = "G3/4", y = Percent_in_ECA)) + geom_boxplot() + geom_beeswarm() +
  scale_x_discrete(name = "") + scale_y_continuous(name = "% early CNVs in ECA") + expand_limits(y=0)

dev.off()

##############################################################################################################################################
## Plot the mutation density distribution

pdf(paste0(output.directory, "MRCA_density.pdf"), width=3.5, height=3.5, useDingbats = F)

max.x <- max(mutation.time.mrca$Max/3.3/10^3)
max.x <- 0.6

# MRCA
to.plot <- cbind(mutation.time.mrca, data.frame(Event="MRCA", PID = rownames(mutation.time.mrca)))
# ECA
to.plot <- rbind(to.plot, cbind(mutation.time.eca, data.frame( Event="ECA", PID = rownames(mutation.time.eca))))

to.plot$Group <- sample.information[to.plot$PID,]$mnp11
to.plot$Subgroup <- sample.information[to.plot$PID,]$mnp12
to.plot$MYC <- sample.information[to.plot$PID,]$MYC
to.plot$MYC[is.na(to.plot$MYC)] <- 0
to.plot$MYCN <- sample.information[to.plot$PID,]$MYCN
to.plot$MYCN[is.na(to.plot$MYCN)] <- 0
to.plot$PRDM6 <- sample.information[to.plot$PID,]$PRDM6
to.plot$PRDM6[is.na(to.plot$PRDM6)] <- 0
to.plot$TERT_prom_mutation <- sample.information[to.plot$PID,]$TERT_prom_mutation
to.plot$TERT_prom_mutation[is.na(to.plot$TERT_prom_mutation)] <- "n.d."
to.plot$Focus <- sample.information[to.plot$PID,]$TargetCohort

# plot Mutation densities at ECA/MRCA for G3/G4 tumors
to.plot. <- rbind(cbind(to.plot[!is.na(to.plot$Group) & to.plot$Group%in%c("MB, G3", "MB, G4") & to.plot$Event=="MRCA",], ecdf.stats(to.plot[!is.na(to.plot$Group) & to.plot$Group%in%c("MB, G3", "MB, G4") & to.plot$Event=="MRCA",], mean = "Mean", lower = "Min", upper = "Max")),
                  cbind(to.plot[!is.na(to.plot$Group) & to.plot$Group%in%c("MB, G3", "MB, G4") & to.plot$Event=="ECA" & !is.na(to.plot$Mean),], ecdf.stats(to.plot[!is.na(to.plot$Group) & to.plot$Group%in%c("MB, G3", "MB, G4") & to.plot$Event=="ECA" & !is.na(to.plot$Mean),], mean = "Mean", lower = "Min", upper = "Max")))

p <- ggplot(to.plot., aes(x=Mean/3.3/10^3, col=Event, fill=Event, group=Event)) + 
  geom_stepribbon(aes(ymin=ecdf.upper, ymax=ecdf.lower, col=NULL), alpha=0.5) + stat_ecdf() + 
  scale_color_manual(values=time.colors) + scale_fill_manual(values=time.colors) + 
  scale_x_continuous(name="#SSNVs/Mb") + scale_y_continuous(name="Fraction of tumors")  +
  ggtitle("MB, G3/G4")

print(p)

# compare MRCA in cases with and without MYC amplification or MYCN amplification or PRDM6 expression


to.plot. <- rbind(cbind(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")])>=1 & to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="MRCA",], 
                        ecdf.stats(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")])>=1 & to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="MRCA",], mean = "Mean", lower = "Min", upper = "Max")),
                  cbind(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")])>=1 & to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="ECA",], 
                        ecdf.stats(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")])>=1 &  to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="ECA",], mean = "Mean", lower = "Min", upper = "Max")))

p1 <- ggplot(to.plot., aes(x=Mean/3.3/10^3, col=Event, fill=Event, group=Event)) + 
  geom_stepribbon(aes(ymin=ecdf.upper, ymax=ecdf.lower, col=NULL), alpha=0.5) + stat_ecdf() + 
  scale_color_manual(values=time.colors) + scale_fill_manual(values=time.colors) + 
  scale_x_continuous(name="#SSNVs/Mb", limits=c(0, max.x)) + scale_y_continuous(name="Fraction of tumors")  +
  ggtitle("MYC/MYCN/PRDM6" ) 

to.plot. <- rbind(cbind(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")]==0)==3 & to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="MRCA",], 
                        ecdf.stats(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")]==0)==3 & to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="MRCA",], mean = "Mean", lower = "Min", upper = "Max")),
                  cbind(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")]==0)==3 & to.plot$Group%in% c("MB, G3", "MB, G4") & to.plot$Event=="ECA",], 
                        ecdf.stats(to.plot[rowSums(to.plot[,c("MYC", "MYCN", "PRDM6")]==0)==3 &  to.plot$Group %in% c("MB, G3", "MB, G4") & to.plot$Event=="ECA",], mean = "Mean", lower = "Min", upper = "Max")))

p2 <- ggplot(to.plot., aes(x=Mean/3.3/10^3, col=Event, fill=Event, group=Event)) + 
  geom_stepribbon(aes(ymin=ecdf.upper, ymax=ecdf.lower, col=NULL), alpha=0.5) + stat_ecdf() + 
  scale_color_manual(values=time.colors) + scale_fill_manual(values=time.colors) + 
  scale_x_continuous(name="#SSNVs/Mb", limits = c(0, max.x)) + scale_y_continuous(name="Fraction of tumors")  +
  ggtitle("No driver in MYC/MYCN/PRDM6"  ) 

print(ggpubr::ggarrange(p1, p2, nrow=2))

## Summary statistics of ECA/MRCA at driver

to.plot. <- to.plot[to.plot$Group %in% c("MB, G3", "MB, G4")  & to.plot$Event=="ECA",]
to.plot.$Driver <- ifelse(rowSums(to.plot.[,c("MYC", "MYCN", "PRDM6")]==0)==3, F, T)

p1 <- ggplot(to.plot., aes(x=Driver, y=Mean/3.3/10^3)) + scale_y_continuous(name="#SSNVs/Mb at ECA") +  
  geom_boxplot() + geom_beeswarm() + scale_x_discrete() + ggtitle("ECA") + theme(aspect.ratio = 1) + 
  stat_compare_means(comparisons = list(c("FALSE","TRUE")), size=1)

wilcox.test(to.plot.[to.plot.$Driver==T,]$Mean, to.plot.[to.plot.$Driver==F,]$Mean)

to.plot. <- to.plot[to.plot$Group %in% c("MB, G3", "MB, G4")  & to.plot$Event=="MRCA",]
to.plot.$Driver <- ifelse(rowSums(to.plot.[,c("MYC", "MYCN", "PRDM6")]==0)==3, F, T)

p2 <- ggplot(to.plot., aes(x=Driver, y=Mean/3.3/10^3)) + scale_y_continuous(name="#SSNVs/Mb at MRCA")+
  geom_boxplot() + geom_beeswarm() + scale_x_discrete() + ggtitle("MRCA") + theme(aspect.ratio = 1) + 
  stat_compare_means(comparisons = list(c("FALSE","TRUE")), size=1)

wilcox.test(to.plot.[to.plot.$Driver==T,]$Mean, to.plot.[to.plot.$Driver==F,]$Mean)

print(ggpubr::ggarrange(p1, p2, ncol=2))


# only mrca against group
to.plot <- cbind(mutation.time.mrca, data.frame(Event="MRCA", PID = rownames(mutation.time.mrca)))
to.plot$Group <- sample.information[to.plot$PID,]$mnp11
to.plot$ecdf.lower <- NA
to.plot$ecdf.mean <- NA
to.plot$ecdf.upper <- NA

for(type in unique(to.plot$Group)){
  if(is.na(type)){next}
  to.plot[!is.na(to.plot$Group) & to.plot$Group==type,c("ecdf.lower",
                                                        "ecdf.mean",
                                                        "ecdf.upper")] <- ecdf.stats(to.plot[!is.na(to.plot$Group) & to.plot$Group==type, ], mean = "Mean", lower = "Min", upper = "Max")
  
}


p<- ggplot(to.plot, aes(x=Mean/3.3/10^3, col=Group, y=ecdf.mean,
                         fill=Group, group=Group)) + 
  geom_stepribbon(data = to.plot, aes(ymin=ecdf.upper, ymax=ecdf.lower, col=NULL), alpha=0.5) + 
  stat_ecdf(data = to.plot) + 
  scale_color_manual(values=group.colors) + scale_fill_manual(values=group.colors) + theme(legend.position = "bottom") +
  scale_x_continuous(name="#SSNVs/Mb at MRCA") + scale_y_continuous(name="Fraction of tumors") + coord_cartesian(xlim=c(0, max.x))
print(p)


# only mrca against subgroup in G3/G4 tumors
to.plot <- cbind(mutation.time.mrca, data.frame(Event="MRCA", PID = rownames(mutation.time.mrca)))
to.plot$Group <- sample.information[to.plot$PID,]$mnp11
to.plot$Subgroup <- sample.information[to.plot$PID,]$mnp12
to.plot <- to.plot[to.plot$Group %in% c("MB, G3", "MB, G4"),]
to.plot$ecdf.lower <- NA
to.plot$ecdf.mean <- NA
to.plot$ecdf.upper <- NA

for(type in unique(to.plot$Subgroup)){
  if(is.na(type)){next}
  to.plot[!is.na(to.plot$Subgroup) & to.plot$Subgroup==type,c("ecdf.lower",
                                                        "ecdf.mean",
                                                        "ecdf.upper")] <- ecdf.stats(to.plot[!is.na(to.plot$Subgroup) & to.plot$Subgroup==type, ], mean = "Mean", lower = "Min", upper = "Max")
  
}


p<- ggplot(to.plot, aes(x=Mean/3.3/10^3, col=Subgroup, y=ecdf.mean,
                        fill=Subgroup, group=Subgroup)) + 
  geom_stepribbon(data = to.plot, aes(ymin=ecdf.upper, ymax=ecdf.lower, col=NULL), alpha=0.5) + 
  stat_ecdf(data = to.plot) + 
  scale_color_manual(values=g34.subgroup.colors) + scale_fill_manual(values=g34.subgroup.colors) + theme(legend.position = "bottom") +
  scale_x_continuous(name="#SSNVs/Mb at MRCA") + scale_y_continuous(name="Fraction of tumors") 
print(p)

# average density per subgroup
mean(to.plot$Mean)/3.3/10^3
sapply(unique(to.plot$Subgroup), function(x){
  mean(to.plot[to.plot$Subgroup==x,]$Mean)/3.3/10^3
})

# are there significant differences between subgroups?
combinations <- as.data.frame(combn(unique(to.plot$Subgroup), 2), stringsAsFactors=FALSE)

sort(p.adjust(apply(combinations, 2, function(x){
  res <- wilcox.test(to.plot[to.plot$Subgroup==x[1],]$Mean/3.3/10^3,
              to.plot[to.plot$Subgroup==x[2],]$Mean/3.3/10^3)
  res$p.value
})))


dev.off()

##############################################################################################################################################
## MRCA density vs age


## all tumors
pdf(paste0(output.directory, "MRCA_density_vs_age.pdf"), width=3.5, height=3.5, useDingbats = F)

to.plot <- cbind(mutation.time.mrca, data.frame(Event="MRCA", PID = rownames(mutation.time.mrca)), 
                 OS = sample.information[rownames(mutation.time.mrca),]$OS, 
                 Age = sample.information[rownames(mutation.time.mrca),]$Age)
to.plot$Group <- sample.information[to.plot$PID,]$mnp11

ggplot(to.plot, aes(x=Max/3.3/10^3, as.numeric(Age), col=Group)) + geom_point() +
  scale_color_manual(values=group.colors) + theme(aspect.ratio = 1)

cor.test(to.plot$Max, to.plot$Age)

dev.off()


## only G3/G4

pdf(paste0(output.directory, "MRCA_density_vs_age_G34.pdf"), width=3.5, height=3.5, useDingbats = F)

to.plot <- cbind(mutation.time.mrca, data.frame(Event="MRCA", PID = rownames(mutation.time.mrca)), 
                 OS = sample.information[rownames(mutation.time.mrca),]$OS, 
                 Age = sample.information[rownames(mutation.time.mrca),]$Age)
to.plot$Group <- sample.information[to.plot$PID,]$mnp11
to.plot <- to.plot[to.plot$Group %in% c("MB, G3", "MB, G4"),]
to.plot$Subgroup <- sample.information[to.plot$PID,]$mnp12

cor.test(to.plot$Max, to.plot$Age, method="spearman")

ggplot(to.plot[to.plot$Mean < 3.3*10^3,], aes(x=Max/3.3/10^3, as.numeric(Age), col=Subgroup)) + geom_point() +
  scale_color_manual(values=g34.subgroup.colors) + theme(aspect.ratio = 1)

dev.off()


##############################################################################################################################################
### Extract for each tumor the VAF distribution for each ploidy state. Exclude sex chromosomes

vafs.all.tumors <- list()
genome.size.all.tumors <- list()

load(paste0(rdata.directory, "Purity_ploidy.RData"))

purities.all.tumors <- purity
ploidies.all.tumors <- ploidy

tert.p <- c()

for(i in tumors){

  print(i)

  ## Read in copy number information and extract ploidy/purity, as before
  aceseq <- list.files(paste0(data.directory, i, "/ACEseq"), pattern="comb_pro_extra")[1]

  ## read in the mutation file
  files <- list.files(paste0(data.directory, "/", i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations <- read.vcf(paste0(data.directory, "/", i, "/", snv.directory, "/", files))
  mutations$vcf <- mutations$vcf[!mutations$vcf$CHROM %in% c("X", "Y"), ]
  
  if(any(c(1295228, 1295250) %in% mutations$vcf$POS)){
    print("++++")
    print(i)
    tert.p <- c(tert.p, i)
  }

  purity <- purities.all.tumors[i]
  ploidy <- ploidies.all.tumors[i]

  copy.number.info <- read.delim(file=paste0(data.directory, "/", i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  copy.number.info <- copy.number.info[!copy.number.info$X.chromosome %in% c("X", "Y"),]
  ## obtain the coverage ratios for the mutations of interest

  ## Extract copy number info for each mutation
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios <- cnv.info.per.mutation$coverage.ratio
  bafs <- cnv.info.per.mutation$baf
  genotype <- cnv.info.per.mutation$genotype
  tcn <- cnv.info.per.mutation$tcn


  ## Extract readcounts of reference and alternative bases
  readcounts <- Extract.info.from.vcf(mutations, info="readcounts")


    #######################################################################
    ## Iterate through all copy number states

  vafs.this.tumor <- list()
  genome.size.this.tumor <- list()
  p <- list()

  ## Plot separately for each copy number
    for(k in unique(copy.number.indicator)){

      expected.coverage.ratio <- (k*purity + (1-purity)*2)/(ploidy*purity+(1-purity)*2)

      readcounts. <- readcounts[((coverage.ratios>(expected.coverage.ratio-0.1) & coverage.ratios<(expected.coverage.ratio+0.1) & !is.na(coverage.ratios)) |
                                   (tcn ==k & !is.na(tcn))) ,,drop=F]


      vafs.this.tumor[[k]] <- readcounts.


      genome.size <- sum(as.numeric(copy.number.info[(copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1)) |
                                                       (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ,]$end)-
                           as.numeric(copy.number.info[(copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1)) |
                                                         (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))),]$start))

      genome.size.this.tumor[[k]] <- genome.size


      for(l in B.allele.indicator[which(copy.number.indicator==k)]){


        readcounts. <- readcounts[((coverage.ratios>(expected.coverage.ratio-0.1) & coverage.ratios<(expected.coverage.ratio+0.1) & !is.na(coverage.ratios)) |
                                     (tcn ==k & !is.na(tcn))) &
                                    (((bafs < (max(l/k, 1-l/k)+0.05) & bafs > (max(l/k, 1-l/k)-0.05)) & !is.na(bafs) | (is.na(bafs) & l==k/2)) |
                                       (genotype==paste(k-l, l, sep=":") & !is.na(genotype)) |
                                       (genotype==paste(l, k-l, sep=":") & !is.na(genotype))),,drop=F]

        if(nrow(readcounts.)==0){next}



        prob.clonal <- l*purity/(purity*k + (1-purity)*2)
        monosomic.prob.clonal <- purity/(purity*k + (1-purity)*2)


      }
    }

  vafs.all.tumors[[i]] <- vafs.this.tumor
  genome.size.all.tumors[[i]] <- genome.size.this.tumor
  save(vafs.all.tumors, genome.size.all.tumors, file=paste0(rdata.directory, "Vafs_all_tumors.RData"))
  
  ## plot a VAF histogram illustrating the estimated density of clonal SNVs stratified by copy number between 1 and 4
  
  p <- list()
  for(CN in 1:4){
    if(length(vafs.all.tumors[[i]][[CN]])==0){
      p[[length(p)+1]] <- ggplot() + theme(axis.line=element_blank())
      next}
    to.plot <-  data.frame(VAF=vafs.all.tumors[[i]][[CN]][,2]/rowSums(vafs.all.tumors[[i]][[CN]]),
                           CN=CN)
    
    avg.depth <- round(mean(rowSums(vafs.all.tumors[[i]][[CN]])))
    
    p[[length(p)+1]] <- ggplot(to.plot, aes(x=VAF)) + geom_histogram(binwidth=1/avg.depth) + 
      scale_x_continuous(limits=c(0,1), name="Variant allele frequency") + ggtitle(paste0("CN = ", CN)) +
      geom_vline(xintercept = seq(1,CN)*purity/(CN*purity + 2*(1-purity)), linetype=2) +
      scale_y_continuous(name="# SNVs")
    
    ## estimated peak sizes from Binomial mixture model:
    peak.sizes <- c()
    for(j in 1:CN){
      peak.sizes <- c(peak.sizes, sum(clonal.mutations.all.tumors[[i]]$clonal.mutation.matrix[clonal.mutations.all.tumors[[i]]$B.allele.indicator==j & 
                                                                                                clonal.mutations.all.tumors[[i]]$copy.number.indicator==CN]))
      
    }
    
    for(B in 1:CN){
      color.this.peak <- ifelse(B==1, time.colors["MRCA"], time.colors["ECA"])
      
      ## the density distribution for the VAF counts follows a binomial distribution 
      
      p[[length(p)]] <-  p[[length(p)]] +
        geom_line(data=data.frame(x=(0:avg.depth)/avg.depth,
                                  p=dbinom(0:avg.depth, 
                                           size = avg.depth, 
                                           prob = purity*(1:CN)[B]/(CN*purity + 2*(1-purity)))*
                                    peak.sizes[B]),
                  aes(x=x, y=p), col=color.this.peak, size=1)
    }
    
    
  }
  
  pdf(paste0(output.directory, "VAF_histogram_w_clonality.pdf"), width=3.5, height=3.5, useDingbats = F)
  print(p)
  dev.off()
  
}

save(vafs.all.tumors, genome.size.all.tumors, file=paste0(rdata.directory, "Vafs_all_tumors.RData"))
