library(gplots)
library(ggplot2)

# This script is adjusted for SComatic results output
# main goal to visualize per subclone mutations and inspect VAF

# SComatic output

inDir = "PATH/HERE"
resDir = "PATH/HERE"

sId = "targetSample"

# required for results figure name
mainId = paste0(sId, "_ATAC_deep")

# somatic
inPath = paste0(inDir,sId,"_deep.calling.step2.somatic.tsv") 

mutInfo <- read.delim(inPath)
rownames(mutInfo) <- paste0(mutInfo$X.CHROM,":",mutInfo$Start)

vafPath = paste0("~/work/spatialMB/mb_mutations/hg38BulkRes/snvs_ICGC_",sId,"_somatic_snvs_conf_8_to_10.hg38.vaf.txt")

vafInfo <- read.delim(vafPath)
rownames(vafInfo) <- vafInfo$Loci


# clean normal
somaticInfo <- somaticInfo[ somaticInfo[,"Normal"] != 1, ]
somaticInfo <- somaticInfo[ ,colnames(somaticInfo) != "Normal"]

hc <- hclust(dist(somaticInfo))



require(ComplexHeatmap)
require(circlize)

row_ha = rowAnnotation(#refVAF = vafInfo[ rownames(somaticInfo), 2],
                       snVAF = snVafInfo[ rownames(somaticInfo)],
                       col = list(  #refVAF = colorRamp2(c(0, 1), c("white", "red4")),
                                    snVAF = colorRamp2(c(0, 1), c("white", "lightgreen")) 
                                    )
                       )

resName2 = paste0(resDir,mainId,".mut_heatmap.pdf")
pdf(resName2, width=4, height = 6)  

ht1 <- Heatmap(somaticInfo, name = "Mut Heatmap", 
              #rect_gp = gpar(type = "none"),
            #   right_annoitation = row_ha,
                col = c("grey","black"),
               cluster_rows = T,
                show_row_names = F, 
                cluster_columns = F,
               row_names_gp = gpar(fontsize = 10)
                )
ht1 + row_ha
dev.off()


# correlation plot
vafDf <- data.frame( bulkVAF = vafInfo[ rownames(somaticInfo), 2],
      snVAF = snVafInfo[ rownames(somaticInfo)])
     
corVal = cor(vafDf[,1],vafDf[,2])

resName1 = paste0(resDir,mainId,".vaf_variance.pdf")

pdf(resName1,width=8,height = 6)
plot(vafDf,  pch=19,
     col="black",
     xlab = "Bulk VAF", ylab = "snATAC VAF", 
     main = sprintf("%s cor = %.4f",sId,corVal))
dev.off()

# boxplot clusters

clusters <- cutree(hc,4) 
summary(as.factor(clusters))


summary(names(clusters) == rownames(somaticInfo))
bxDf <- data.frame(Cluster=paste0("cl",clusters), 
                    #VAF=vafInfo[ rownames(somaticInfo), 2] ) # v1
                    VAF=snVafInfo[ rownames(somaticInfo)] ) # v2
    

resName3 = paste0(resDir,mainId,".vaf_boxplot.pdf")
require(ggplot2)
pdf(resName3, width=8,height=6)
p <- ggplot(data = bxDf, aes(x=Cluster, y=VAF,fill="red4")) +
    geom_boxplot() +   
    theme_bw() +
    theme(legend.position="none") 
print(p)
dev.off()





