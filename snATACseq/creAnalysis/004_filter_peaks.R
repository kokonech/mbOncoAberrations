## This is script for filtering peak matrix generated in 003_Peakcalling.R
## The aim is the filter out peaks not detected in at least three samples based on
## cut-off activity

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(Matrix)
})

addArchRThreads(threads =6)
addArchRGenome("hg38")

DIR_OUT <- "/path/to/output/ATACana/"

setwd(DIR_OUT)
#
## Obtain Peak matrix and Peak set-----
## load combined ArchR project
proj <- loadArchRProject("comATAC/")
getAvailableMatrices(proj)

PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
PeakSet <- getPeakSet(proj)
save(PeakMatrix, file = "comATAC/combinedPeakMatrix.RData")
save(PeakSet, file = "comATAC/combinedPeakSet.RData")
#q()
#
########### filtering peaks###########

## Processing Peak matrix and PeakSet ---------
#load("OUT/combinedPeakMatrix.RData") ##peak matrix
#load("OUT/combinedPeakSet.RData") ##peakset

region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)
dim(region_names)
PS <- data.frame(PeakSet) 
row.names(PS) <- paste0(PS$seqnames,":",PS$start,"-",PS$end)
head(PS)
##adding cluster info
cl <- PS$GroupReplicate
cl <-gsub("._."," ",cl,fixed = TRUE)
cl <- data.frame(X=cl)

cl2 <- cl %>% separate(X,c("A","B"), sep=" ")
head(cl2)
PS$Cluster <- cl2$A
head(PS)
rm(cl,cl2)
table(row.names(PS)==row.names(region_names))
###all peaks exist in peakset, and are in same order as in the PeakMatrix.!!!

##obtaining pseudobulk actvitiy---------
plot.data <- data.frame(getCellColData(proj))
dim(plot.data)
clust <- sort(unique(PS$Cluster))
dim(PeakMatrix)
table(colnames(PeakMatrix) %in% row.names(plot.data)) ##all exists!
## rearranging the cell order for plotdata as thts easier
pd <- plot.data[colnames(PeakMatrix),]
table(colnames(PeakMatrix)==row.names(pd)) ##fixed!

##Ioannis's lapply function takes a lot of memory, so using for loop
peak_clust_f <- c()
for(id in clust){
  cell_i <- which(pd$Sample==id)
  PM <- PeakMatrix[,cell_i]
  PM <- assay(PM)
  del  <- rowSums(PM > 0)/length(cell_i)
  peak_clust_f <- cbind(peak_clust_f,del)
  rm(del,cell_i,PM)
}
rm(id)

colnames(peak_clust_f) <- clust
row.names(peak_clust_f) <- row.names(region_names)
peak_clust_f[1:10,1:4]

save(PeakMatrix, PS, pd,peak_clust_f, file = "comATAC/peak_freq_sample.RData")
#q()
## filtering for peak activity by pseudobulk ----------------
#load("comATAC/peak_freq_sample.RData")

PS$max_freq <- apply(peak_clust_f, 1, max)
sum(PS$max_freq >= 0.01)
sum(PS$max_freq >= 0.03)
sum(PS$max_freq >= 0.05)

ggplot(PS, aes(log10(max_freq))) +
  geom_histogram(bins = 30) +
  theme_classic()
p <- ggplot(PS, aes(x=log10(max_freq), y=log10(score))) +
  geom_hex() +
  theme_classic()
tiff("OUT/Peaks_maxFreq_byPeakScore.tiff", unit = "in",  width=5, height=5, res = 300)
print(p)
dev.off()
p <- ggplot(PS, aes(x=as.factor(Reproducibility), y=log10(max_freq), fill=as.factor(Reproducibility))) +
  geom_violin() +
  geom_boxplot(notch = T, width=0.1) +
  scale_fill_brewer(palette = "Greens", name="") +
  geom_hline(yintercept = c(log10(0.03),log10(0.05)), color=c("red","orange"), lty="dashed") +
  xlab("Number of samples supporting peak") +
  theme_classic()

tiff("Peaks_maxFreq_byReproducibility.tiff", unit = "in",  width=6, height=4, res = 300)
print(p)
dev.off()
## cut-off of 5% seems better 

cutoff <- 0.05

PS$robust <- PS$max_freq >= cutoff 

save(PS,  file = "comATAC/Peaks_robust.RData")
##write robust and normal peaks into bed file
bed <- PS[,c("seqnames",  "start",    "end")]
head(bed)
write.table(bed, file="comATAC/G34_allPeaks.bed", row.names = FALSE, col.names = FALSE,quote = FALSE, sep="\t")
bed <- PS[PS$robust,c("seqnames",  "start",    "end")]
head(bed)
write.table(bed, file="comATAC/G34_allPeaks_robust.bed", row.names = FALSE,col.names = FALSE, quote = FALSE, sep="\t")
##these peaks are sorted later

q()
##this generated G34_allPeaks_robust.bed was sorted and filtered to remove blacklisted regions
#sort -k 1,1 -k2,2n G34_allPeaks_robust.bed > G34_allPeaks_robust.sorted.bed
#bedtools subtract -A -a G34_allPeaks_robust.sorted.bed -b hg38-blacklist.v2.bed > G34_allPeaks_robust.fil.sorted.bed
