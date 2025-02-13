# Collect the timing information of chromosomal gains of interest, for G3/G4 tumors only
load("./RScripts/RData/MRCA_timing.RData")

g3g4.tumors <- rownames(sample.information[sample.information$mnp11 %in% c("MB, G3", "MB, G4"),])

load("./RScripts/RData/CNV_per_tumor_res_1Mb.RData")

## For the CNVs of interest (Gain of 4, 7, 12, 17q, 1q and 18; Loss of Chr 8, 11, 5q, 10) (Doussouki et al. and own analysis)
cnv.mat.g3g4 <- matrix("", nrow = 14, ncol = length(g3g4.tumors),
                       dimnames = list(c("Gain whole 4", "Gain whole 7", "Gain q-arm 7", "Gain whole 12", "Gain q-arm 17", 
                                         "Gain q-arm 1", "Gain whole 18", "Loss whole 8", "Loss whole 11", "Gain whole 17",
                                         "Loss q-arm 5", "Loss whole 10", "Loss q-arm 10", "Loss p-arm 17"),
                                       g3g4.tumors))

for(i in 1:nrow(cnv.mat.g3g4)){
  type = gsub("\\ .*", "", rownames(cnv.mat.g3g4)[i])
  arm = strsplit(rownames(cnv.mat.g3g4)[i], split=" ")[[1]][2]
  if(arm == "whole"){
    arm <- "both arms"
  }else if(arm == "q-arm"){
    arm <- "q arm"
  }else if(arm == "p-arm"){
    arm <- "p arm"
  }
  chr = strsplit(rownames(cnv.mat.g3g4)[i], split=" ")[[1]][3]
  
  tmp <- cnv.summary[cnv.summary$Chr == chr & cnv.summary$CNV == tolower(type) &
                       cnv.summary$Arm == arm & cnv.summary$ID %in% g3g4.tumors, ]
  if(nrow(tmp)==0){next}
  
  timing <- tmp$Clonality
  # losses cannot be timed --> all n.d.
  
  # if(type == "Gain"){
  timing[(is.na(timing) | timing == "clonal") & tmp$CN < 5 & tmp$CN > 1] <- sapply(tmp$ID[(is.na(timing) | timing=="clonal") & tmp$CN < 5 & tmp$CN > 1], function(id){
    # does the gain map to MRCA?
    if(any(grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.at.mrca)) & 
       any(grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.at.mrca.conforming.eca))){
      res <- "ECA/MRCA"
    }else if(any(grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.uniquely.mapped.to.eca))){
      res <- "ECA"
    }else if(any(grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.at.mrca))){
      res <- "MRCA"
    }else if(any(grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.at.earliest.time))){
      res <- "< ECA"
    }else if(any(grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.not.maping.to.eca.or.mrca))){
      segment <- mrca.eca[[id]]$gains.not.maping.to.eca.or.mrca[grepl(pattern = paste0("_",chr), mrca.eca[[id]]$gains.not.maping.to.eca.or.mrca)]
      if(any(mutation.time[[id]][mutation.time[[id]]$Segment %in% segment,]$Min > mrca.eca[[id]]$mutation.time.mrca/3.3/10^3 &
             mutation.time[[id]][mutation.time[[id]]$Segment %in% segment,]$Min > mrca.eca[[id]]$mutation.time.mrca.upper/3.3/10^3)){
        res <- "> MRCA"
      }else if(any(mutation.time[[id]][mutation.time[[id]]$Segment %in% segment,]$Mean > mrca.eca[[id]]$mutation.time.eca/3.3/10^3 &
                   mutation.time[[id]][mutation.time[[id]]$Segment %in% segment,]$Mean < mrca.eca[[id]]$mutation.time.mrca.upper/3.3/10^3)){
        res <- "> ECA, < MRCA"
      }else{
        res <- "clonal"
      }
    }else{
      res <- "clonal"
    }
    return(res)
  })
  # }else{
  #   timing <- rep("n.d.", length(tmp$ID))
  # }
  cnv.mat.g3g4[i,tmp$ID] <- timing
}


# Also plot for each of these CNVs the absolute mutation timing relative to the MRCA

p <- list()
col_clonal <- c( time.colors, "ECA/MRCA" = "violet", "-" = "grey")


for(i in sort(rownames(cnv.mat.g3g4))){
  type = gsub("\\ .*", "", i)
  arm = strsplit(i, split=" ")[[1]][2]
  if(arm == "whole"){
    arm <- "both arms"
  }else if(arm == "q-arm"){
    arm <- "q arm"
  }else if(arm == "p-arm"){
    arm <- "p arm"
  }
  chr = strsplit(i, split=" ")[[1]][3]
  
  tmp <- cnv.summary[cnv.summary$Chr == chr & cnv.summary$CNV == tolower(type) &
                       cnv.summary$Arm == arm & cnv.summary$ID %in% g3g4.tumors, ]
  if(nrow(tmp)==0){next}
  
  # now, get for each gain the mutation timing with lower and upper bounds
  
  tmp$SNV_density <- apply(tmp, 1, function(x){
    if(is.na(x["CN"])){return(NA)}
    if(x["CN"] > 4){return(NA)}
    tmp2 <- mutation.time[[x["ID"]]]
    tmp2 <- tmp2[grepl(paste("chr", x["Chr"], sep="_"), as.character(tmp2$Segment)),,drop=F]
    tmp2 <- tmp2[sapply(as.character(tmp2$Segment), function(y){
      strsplit(y, split="_")[[1]][4] != "I"
    }),,drop=F]
    mean(tmp2$Mean/3.3/10^3)
  })
  
  tmp$SNV_density_lower <- apply(tmp, 1, function(x){
    if(is.na(x["CN"])){return(NA)}
    if(x["CN"] > 4){return(NA)}
    tmp2 <- mutation.time[[x["ID"]]]
    tmp2 <- tmp2[grepl(paste("chr", x["Chr"], sep="_"), as.character(tmp2$Segment)),,drop=F]
    tmp2 <- tmp2[sapply(as.character(tmp2$Segment), function(y){
      strsplit(y, split="_")[[1]][4] != "I"
    }),,drop=F]
    mean(tmp2$Min/3.3/10^3)
  })
  
  tmp$SNV_density_upper <- apply(tmp, 1, function(x){
    if(is.na(x["CN"])){return(NA)}
    if(x["CN"] > 4){return(NA)}
    tmp2 <- mutation.time[[x["ID"]]]
    tmp2 <- tmp2[grepl(paste("chr", x["Chr"], sep="_"), as.character(tmp2$Segment)),,drop=F]
    tmp2 <- tmp2[sapply(as.character(tmp2$Segment), function(y){
      strsplit(y, split="_")[[1]][4] != "I"
    }),,drop=F]
    mean(tmp2$Max/3.3/10^3)
  })
  
  tmp$MRCA <- sapply(tmp$ID, function(x){mrca.eca[[x]]$mutation.time.mrca/3.3/10^3})
  tmp$MRCA_lower <- sapply(tmp$ID, function(x){mrca.eca[[x]]$mutation.time.mrca.lower/3.3/10^3})
  tmp$MRCA_upper <- sapply(tmp$ID, function(x){mrca.eca[[x]]$mutation.time.mrca.upper/3.3/10^3})
  
  tmp$Assignment <- apply(tmp, 1, function(x){
    if(is.na(x["CN"])){return(NA)}
    if(x["CN"] > 4){return(NA)}
    tmp2 <- mrca.eca[[x["ID"]]]$p.values.mrca
    tmp2 <- tmp2[grepl(paste("chr", x["Chr"], sep="_"), rownames(tmp2)),,drop=F]
    tmp2 <- tmp2[sapply(rownames(tmp2), function(y){
      strsplit(y, split="_")[[1]][4] != "I"
    }),,drop=F]
    mrca <- any(tmp2$adj.p >= 0.01)
    tmp2 <- mrca.eca[[x["ID"]]]$p.values.eca
    if(length(tmp2)>0){
      tmp2 <- tmp2[grepl(paste("chr", x["Chr"], sep="_"), rownames(tmp2)),,drop=F]
    }else{
      if(mrca){
        return("MRCA")
      }else(return("-"))
    }
    if( nrow(tmp2)>0){
      tmp2 <- tmp2[sapply(rownames(tmp2), function(y){
        strsplit(y, split="_")[[1]][4] != "I"
      }),,drop=F]
    }
    if(nrow(tmp2)>0){
      eca <- any(is.na(tmp2$adj.p) | tmp2$adj.p >= 0.01)
      if(eca & mrca){
        "ECA/MRCA"
      }else if(eca){
        "ECA"
      }else if(mrca){
        "MRCA"
      }else("-")
    }else{
      if(mrca){
        "MRCA"
      }else("-")
    }
    
  })
  
  p[[length(p) + 1]] <- ggplot(tmp, aes(x = SNV_density, xmin = SNV_density_lower, xmax = SNV_density_upper,
                                        y = MRCA, ymin = MRCA_lower, ymax = MRCA_upper, col = Assignment)) + geom_point() +
    geom_errorbar(width = 0) + geom_errorbarh(height = 0) + scale_color_manual(values = col_clonal) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) + scale_x_continuous(name = "SNV density at chromosomal gain") +
    scale_y_continuous(name = "SNV density at MRCA") + ggtitle(i) +
    expand_limits(x = 0, y = 0) + theme(aspect.ratio = 1)
  
}


pdf(paste0(output.directory, "/CNV_timing_per_segment.pdf"), width=8, height=5, useDingbats = F)

ggarrange(plotlist = p, nrow = 3, ncol = 5, common.legend = T)

dev.off()