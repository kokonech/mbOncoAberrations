###################### Infer parameters of a selection model from tumor samples
#i=patient.id # sets the patient ID

load("./RScripts/RData/Purity_ploidy.RData") # load ploidy/purity info from ACEseq
load("./RScripts/RData/Vafs_all_tumors.RData") # load the CN-stratified VAFs


ploidy.all <- ploidy
purity.all <- purity

purity <- purity.all[i]
ploidy <- ploidy.all[i]

############ The fit involves the subclonal tails of the haploid, diploid, triploid and tetraploid fraction


source("./RScripts/Mutation_accumulation_exponential_growth.R")  # source the model functions


## store the genome fraction that is haploid, diploid, triploid and tetraploid and decide whether to fit all or some of them; exclude sex chromosomes


haploid.genome.fraction <- genome.size.all.tumors[[i]][[1]]
diploid.genome.fraction <- genome.size.all.tumors[[i]][[2]]
if(length(genome.size.all.tumors[[i]])>=3){
  triploid.genome.fraction <- genome.size.all.tumors[[i]][[3]]
}else{
  triploid.genome.fraction <- 0
}
if(length(genome.size.all.tumors[[i]])>=4){
  tetraploid.genome.fraction <- genome.size.all.tumors[[i]][[4]]
}else{
  tetraploid.genome.fraction <- 0
}

if(haploid.genome.fraction < 10^8){
  fit.haploid <- F
}else{
  fit.haploid <- T
}
if(diploid.genome.fraction < 10^8){
  fit.diploid <- F
}else{
  fit.diploid <- T
}
if(triploid.genome.fraction < 10^8){
  fit.triploid <- F
}else{
  fit.triploid <- T
}
if(tetraploid.genome.fraction < 10^8){
  fit.tetraploid <- F
}else{
  fit.tetraploid <- T
}

if(length(vafs.all.tumors[[i]][[1]])>1){
  vafs.haploid <- vafs.all.tumors[[i]][[1]][,2]/rowSums(vafs.all.tumors[[i]][[1]])
  vafs.haploid <- vafs.haploid[!is.na(vafs.haploid)]
  depth.haploid <- round(mean(rowSums(vafs.all.tumors[[i]][[1]]), na.rm = T))
  if(length(vafs.haploid)==0){
    fit.haploid <- F
  }
}

if(fit.diploid && length(vafs.all.tumors[[i]][[2]])>1){
  vafs.diploid <- vafs.all.tumors[[i]][[2]][,2]/rowSums(vafs.all.tumors[[i]][[2]])
  vafs.diploid <- vafs.diploid[!is.na(vafs.diploid)]
  
  ## to quantify correctly, we need to cut off higher order peaks and add them to the lower order peaks (otherwise we may underestimate the time prior to
  ## the MRCA as some mutations are "lost" due to amplification)
  depth.diploid <- round(mean(rowSums(vafs.all.tumors[[i]][[2]]), na.rm = T))
  
  homozygous.mutations <- vafs.diploid[vafs.diploid>qbinom(p=0.95, size = depth.diploid, prob = purity/2)/depth.diploid]
  vafs.diploid <- vafs.diploid[vafs.diploid<=qbinom(p=0.95, size = depth.diploid, prob = purity/2)/depth.diploid]
  
  if(length(homozygous.mutations)>0){
    vafs.diploid <- c(vafs.diploid, rep(homozygous.mutations/2, 2))
  }
  
  
}
if(fit.triploid && length(vafs.all.tumors[[i]][[3]])>1){
  vafs.triploid <- vafs.all.tumors[[i]][[3]][,2]/rowSums(vafs.all.tumors[[i]][[3]])
  vafs.triploid <- vafs.triploid[!is.na(vafs.triploid)]
  ## cut off higher order peaks from tri- and tetraploid tumors
  depth.triploid <- round(mean(rowSums(vafs.all.tumors[[i]][[3]]), na.rm = T))
  
  vafs.triploid.two.alleles <- vafs.triploid[vafs.triploid>qbinom(p=0.95, size = depth.triploid, prob = purity/(3*purity + 2*(1-purity)))/depth.triploid &
                                               vafs.triploid<=qbinom(p=0.95, size = depth.triploid, prob = 2*purity/(3*purity + 2*(1-purity)))/depth.triploid]
  vafs.triploid.three.alleles <- vafs.triploid[vafs.triploid>qbinom(p=0.95, size = depth.triploid, prob = 2*purity/(3*purity + 2*(1-purity)))/depth.triploid]
  vafs.triploid <- vafs.triploid[vafs.triploid<=qbinom(p=0.95, size = depth.triploid, prob = purity/(3*purity + 2*(1-purity)))/depth.triploid]
  if(length(vafs.triploid.two.alleles)>0){
    vafs.triploid <- c(vafs.triploid, rep(vafs.triploid.two.alleles*2/3, 2))
  }
  if(length(vafs.triploid.three.alleles)>0){
    vafs.triploid <- c(vafs.triploid, rep(vafs.triploid.three.alleles*1/3, 3))
    
  }
  
}
if(fit.tetraploid && length(vafs.all.tumors[[i]][[4]])>1){
  vafs.tetraploid <- vafs.all.tumors[[i]][[4]][,2]/rowSums(vafs.all.tumors[[i]][[4]])
  vafs.tetraploid <- vafs.tetraploid[!is.na(vafs.tetraploid)]
  ## cut off higher order peaks from tri- and tetraploid tumors
  depth.tetraploid <- round(mean(rowSums(vafs.all.tumors[[i]][[4]]), na.rm = T))
  
  vafs.tetraploid.two.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.95, size = depth.tetraploid, prob = purity/(4*purity + 2*(1-purity)))/depth.tetraploid &
                                                   vafs.tetraploid<=qbinom(p=0.95, size = depth.tetraploid, prob = 2*purity/(4*purity + 2*(1-purity)))/depth.tetraploid]
  vafs.tetraploid.three.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.95, size = depth.tetraploid, prob = 2*purity/(4*purity + 2*(1-purity)))/depth.tetraploid &
                                                     vafs.tetraploid<=qbinom(p=0.95, size = depth.tetraploid, prob = 3*purity/(4*purity + 2*(1-purity)))/depth.tetraploid ]
  vafs.tetraploid.four.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.95, size = depth.tetraploid, prob = 3*purity/(4*purity + 2*(1-purity)))/depth.tetraploid  ]
  
  vafs.tetraploid <- vafs.tetraploid[vafs.tetraploid<=qbinom(p=0.95, size = depth.tetraploid, prob = purity/(4*purity + 2*(1-purity)))/depth.tetraploid]
  if(length(vafs.tetraploid.two.alleles)>0){
    vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.two.alleles/2,2))
  }
  if(length(vafs.tetraploid.three.alleles)>0){
    vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.three.alleles*1/3, 3))
  }
  if(length(vafs.tetraploid.four.alleles)>0){
    vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.four.alleles*1/4, 4))
  }
}


### compute the cumulative distribution functions and do upscaling for the entire genome

  ## haploid fraction
  if(fit.haploid){
    cum.data.haploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.haploid >= x)
    })*3.3*10^9/haploid.genome.fraction
   }else{
    cum.data.haploid <- 0
   }
  if(fit.diploid){
    cum.data.diploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.diploid >= x)
    })*3.3*10^9/diploid.genome.fraction
   }else{
    cum.data.diploid <- 0
   }
  if(fit.triploid){
    cum.data.triploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.triploid >= x)
    })*3.3*10^9/triploid.genome.fraction
   }else{
    cum.data.triploid <- 0
   }
  if(fit.tetraploid){
    cum.data.tetraploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.tetraploid >= x)
    })*3.3*10^9/tetraploid.genome.fraction
   }else{
    cum.data.tetraploid <- 0
   }


mySumStatData <- list(haploid = cum.data.haploid, diploid = cum.data.diploid, triploid = cum.data.triploid, tetraploid = cum.data.tetraploid)

myModel <- function(parms){

## the model functions were sourced from the file Mutation_accumulation_exponential_growth.R
  
  parms$lambda <- 1  
  
  if(selection){
    ## set tau = t - ts; require clone size between 1 and 99% and convert t_s accordingly
    min.tau = log(0.01*10^9)/(parms$lambda - parms$s*parms$delta)
    max.tau = log(0.99*10^9)/(parms$lambda - parms$s*parms$delta)
    
    tau = min.tau + parms$t_s*(max.tau - min.tau)
    
    ## estimate size of the selected subclone at resection
    
    subclonal.size <- exp((parms$lambda - parms$s*parms$delta)*tau) 
    founder.size <- 10^9 - exp((parms$lambda - parms$s*parms$delta)*tau) + exp((parms$lambda - parms$delta)*tau) 
    
    t.max <- log(founder.size)/(parms$lambda - parms$delta)
    
    parms$t_s = t.max - tau
  }else{
    subclonal.size <- 0
  }

  # allow purity adjustment by +-0.05
  parms$purity <- purity + parms$purity_adj
  if(parms$purity>1){
   parms$purity <- 1
  }  
## computation
  if(fit.haploid){
    model <- SimulateVAFDistr(parms, ploidy = 1, depth = depth.haploid, purity = parms$purity, selection = selection)
    cum.sim.haploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.haploid))){
  return(Inf)}
  }else{
    cum.sim.haploid <- 0
  }
  if(fit.diploid){
    model <- SimulateVAFDistr(parms, ploidy = 2, depth = depth.diploid, purity = parms$purity, selection = selection)
    cum.sim.diploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.diploid))){
      return(Inf)}
  }else{
    cum.sim.diploid <- 0
  }
  if(fit.triploid){
    model <- SimulateVAFDistr(parms, ploidy = 3, depth = depth.triploid, purity = parms$purity, selection = selection)
    cum.sim.triploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.triploid))){
     return(Inf)}
  }else{
    cum.sim.triploid <- 0
  }
  if(fit.tetraploid){
    model <- SimulateVAFDistr(parms, ploidy = 4, depth = depth.tetraploid, purity = parms$purity, selection = selection)
    cum.sim.tetraploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.tetraploid))){
     return(Inf)}
  }else{
    cum.sim.tetraploid <- 0
  }


 list(haploid = cum.sim.haploid, diploid = cum.sim.diploid, triploid = cum.sim.triploid, tetraploid = cum.sim.tetraploid, subclonal.size=subclonal.size/10^9)
}


## Distance to data

mySummaryDistance <- function(sumStatSample, sumStatData){

if(selection){  
  ## squared distances as error
sqrt( sum(((sumStatSample$haploid - sumStatData$haploid)*haploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$diploid - sumStatData$diploid)*diploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) + 
    sum(((sumStatSample$triploid - sumStatData$triploid)*triploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$tetraploid - sumStatData$tetraploid)*tetraploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2))/
(length(sumStatSample$haploid )+ length(sumStatSample$diploid) + length(sumStatSample$triploid) + length(sumStatSample$tetraploid) ) + ifelse(sumStatSample$subclonal.size <0.05 | sumStatSample$subclonal.size >0.95, Inf, 0)
}else{
  ## squared distances as error
sqrt( sum(((sumStatSample$haploid - sumStatData$haploid)*haploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$diploid - sumStatData$diploid)*diploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$triploid - sumStatData$triploid)*triploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$tetraploid - sumStatData$tetraploid)*tetraploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2))/
(length(sumStatSample$haploid )+ length(sumStatSample$diploid) + length(sumStatSample$triploid) + length(sumStatSample$tetraploid) ) 

}

}

