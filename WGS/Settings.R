## General settings, such as libraries, samples taken, color codes, etc.
##############################################################################################################################################
library("RColorBrewer")
library(ggplot2); theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=6, color="black"),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                  axis.ticks = element_line(color = "black", size=0.75),
                                  axis.text = element_text(size=6, color="black")))
library(bedr)
library(openxlsx)
library(NBevolution)
library(ggpubr)
library(ggsci)
library(pammtools)
library(tidyverse)
library(mixtools)
##############################################################################################################################################
## set directories
data.directory <- "./data/"
snv.directory <- "SNVs"
cnv.directory <- "ACEseq/"
indel.directory <- "Indels/"
sv.directory <- "SVs/"
custom.script.directory <- "./"
output.directory <- "./Plots/"
meta.data <- "./Meta_data/"
rdata.directory <- "./RScripts/RData/"

if(!dir.exists(rdata.directory)){
  dir.create(rdata.directory)
}

##############################################################################################################################################
## Define samples to use, meta information and colors
sample.information<- read.xlsx(paste0(meta.data, "Extended Data Table 6.xlsx"), sheet=1)
rownames(sample.information) <- sample.information$MRCA_ID
tumors <- rownames(sample.information)

##############################################################################################################################################
## Color palette:
## early clonal, late clonal, subclonal

time.colors <- c("#4FB12B", "#176A02")
names(time.colors) <- c("MRCA", "ECA")

group.colors <- c(`MB, G3` = "goldenrod1", `MB, G4` = "darkgreen", `MB, SHH CHL AD` = "firebrick",
                  `MB, WNT` = "blue", `MB, SHH INF` = "red")
g34.subgroup.colors <- c(MB_G34_I = "#cd89d9", MB_G34_II = "#A34D7C", MB_G34_III = "#C4A908", 
                     MB_G34_IV = "#FFFF00",
                     MB_G34_V = "#ADCA02", MB_G34_VI = "#89B395", MB_G34_VII = "#9BC2E6", MB_G34_VIII = "#127134")

##############################################################################################################################################
## Driver genes 

candidates.intogen <- read.delim("refdata/2023-05-31_IntOGen-Drivers/Compendium_Cancer_Genes.tsv")
candidates.intogen <- candidates.intogen$SYMBOL[candidates.intogen$CANCER_TYPE %in% c("MBL")]

driver.genes <- c("SNCAIP", "TERT", "MYC", "PRDM6", candidates.intogen)

## genetic positions
gene.positions <- read.delim("refdata/gencode_v19_gene_pos.txt", header=F, row.names = 1)
driver.genes.with.genetic.pos <- data.frame(GENE=unique(driver.genes), gene.positions[unique(driver.genes),])
colnames(driver.genes.with.genetic.pos) <- c("GENE", "CHROM", "START", "END")

## for TERT we want the upstream region
driver.genes.with.genetic.pos["TERT","END"] <- 1297200
##############################################################################################################################################
## Tumors for which we re-estimate the purity
tumor.to.adjust.purity <- c("MB140", "MB188")

##############################################################################################################################################
## Compute lower and upper bounds of an ecdf statistics

# data: the data frame, mean, lower, upper: column indices giving mean, and lower and upper CI of the measurements
ecdf.stats <- function(data, mean, lower, upper){
  if(nrow(data)==0){
    return(NULL)
  }
  l <- ecdf(data[,lower])
  u <- ecdf(data[,upper])
  m <- ecdf(data[,mean])
  data.frame(ecdf.lower = l(data[,mean]),
             ecdf.mean = m(data[,mean]),
             ecdf.upper = u(data[,mean]))
}

# example:
# data <- data.frame(Mean = seq(1,10), Lower = seq(1,10)-1, Upper=seq(1,10)+1)
# ecdf.stats(data, "Mean", "Lower", "Upper")

##############################################################################################################################################
## Centromer information
centromers <- read.delim("./refdata/Centromers.tsv")
centromers <- centromers[centromers$type=="centromere",]
centromers$chrom[centromers$chrom=="chrX"] <- "chr23"
centromers$chrom[centromers$chrom=="chrY"] <- "chr24"

