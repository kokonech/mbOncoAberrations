# Collect the timing information of chromosomal gains of interest, for G3/G4 tumors only
load("./RScripts/RData/MRCA_timing.RData")

g3g4.tumors <- rownames(sample.information[sample.information$mnp11 %in% c("MB, G3", "MB, G4"),])

load("./RScripts/RData/CNV_per_tumor_res_1Mb.RData")

## For the CNVs of interest (Gain of 4, 7, 12, 17q, 1q and 18; Loss of Chr 8, 11, 5q, 10) (Doussouki et al. and own analysis)
cnv.mat.g3g4 <- matrix("", nrow = 13, ncol = length(g3g4.tumors),
                       dimnames = list(c("Gain whole 4", "Gain whole 7", "Gain q-arm 7", "Gain whole 12", "Gain q-arm 17", 
                                         "Gain q-arm 1", "Gain whole 18", "Loss whole 8", "Loss whole 11", "Gain whole 17",
                                         "Loss q-arm 5", "Loss whole 10", "Loss q-arm 10"),
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
    timing[(is.na(timing) | timing == "clonal") & tmp$CN < 5] <- sapply(tmp$ID[(is.na(timing) | timing=="clonal") & tmp$CN < 5], function(id){
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
