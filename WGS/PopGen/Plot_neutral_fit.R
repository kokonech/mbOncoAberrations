##############################################################################################################################################
source("./RScripts/Settings.R")
##############################################################################################################################################
## load data
load(paste0(rdata.directory, "Purity_ploidy.RData"))
load(paste0(rdata.directory, "/Vafs_all_tumors.RData"))
##############################################################################################################################################

selection <- F
tumors.subclonal.evol <- rownames(sample.information)[!is.na(sample.information$Use.for.neutral.fitting) & sample.information$Use.for.neutral.fitting=="yes"]

for(i in tumors.subclonal.evol){
  print(i)
  
  fits <- list.files("Processed_data/Neutral_fits/", pattern=paste0(i, ".csv"), full.names = T)
  if(length(fits)==0){next}
  fits <- read.csv(fits)
  
  source("./RScripts/Mutation_accumulation_fit_pre_clonal_and_clonal.R")
  
  if(fit.haploid){
    sim.haploid <- matrix(0, nrow=100, ncol=length(mySumStatData$haploid))
  }
  if(fit.diploid){
    sim.diploid <- matrix(0, nrow=100, ncol=length(mySumStatData$diploid))
  }
  if(fit.triploid){
    sim.triploid <- matrix(0, nrow=100, ncol=length(mySumStatData$triploid))
  }
  if(fit.tetraploid){
    sim.tetraploid <- matrix(0, nrow=100, ncol=length(mySumStatData$tetraploid))
  }
  
  for(j in 1:100){
    parms <- list(delta=fits$par_delta[j], n_clonal=fits$par_n_clonal[j], mu = fits$par_mu[j],
                  purity_adj = fits$par_purity_adj[j])
    output <- myModel(parms)
    
    if(fit.haploid){
      sim.haploid[j,] <- output$haploid
      max.haploid <- max(apply(sim.haploid, 2, max)*haploid.genome.fraction/(3.3*10^9))
    }else{
      max.haploid <- 0
    }
    if(fit.diploid){
      sim.diploid[j,] <- output$diploid
      max.diploid <- max(apply(sim.diploid, 2, max)*diploid.genome.fraction/(3.3*10^9))
    }else{
      max.diploid <- 0
    }
    if(fit.triploid){
      sim.triploid[j,] <- output$triploid
      max.triploid <- max(apply(sim.triploid, 2, max)*triploid.genome.fraction/(3.3*10^9))
    }else{
      max.triploid <- 0
    }
    if(fit.tetraploid){
      sim.tetraploid[j,] <- output$tetraploid
      max.tetraploid <- max(apply(sim.tetraploid, 2, max)*tetraploid.genome.fraction/(3.3*10^9))
    }else{
      max.tetraploid <- 0
    }
    
  }
  
  y.max <- max(max.haploid, max.diploid, max.triploid, max.tetraploid)
  
  p <- list()
  
  if(fit.haploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$haploid*haploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.haploid, 2, quantile, p=0.025)*haploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.haploid, 2, quantile, p=0.975)*haploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data, 
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                            )) + 
      geom_ribbon(aes(ymin=Mmin, 
                      ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + ggtitle(paste("CN=1, weight=", 
                                                                       round(haploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) + 
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"), list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                     diploid.genome.fraction=diploid.genome.fraction,
                                                                                     triploid.genome.fraction=triploid.genome.fraction,
                                                                                     tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  if(fit.diploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$diploid*diploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.diploid, 2, quantile, p=0.025)*diploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.diploid, 2, quantile, p=0.975)*diploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data,   
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                            )) + 
      geom_ribbon(aes(ymin=Mmin, 
                      ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + ggtitle(paste("CN=2, weight=", 
                                                                       round(diploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) +
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"),
    list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                   diploid.genome.fraction=diploid.genome.fraction,
                                                                                   triploid.genome.fraction=triploid.genome.fraction,
                                                                                   tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  if(fit.triploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$triploid*triploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.triploid, 2, quantile, p=0.025)*triploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.triploid, 2, quantile, p=0.975)*triploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data, 
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                            )) + 
      geom_ribbon(aes(ymin=Mmin, 
                      ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01)  + ggtitle(paste("CN=3, weight=", 
                                                                        round(triploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) + 
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"), list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                     diploid.genome.fraction=diploid.genome.fraction,
                                                                                     triploid.genome.fraction=triploid.genome.fraction,
                                                                                     tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  if(fit.tetraploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$tetraploid*tetraploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.tetraploid, 2, quantile, p=0.025)*tetraploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.tetraploid, 2, quantile, p=0.975)*tetraploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data, 
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                           )) +
      geom_ribbon(aes( ymin=Mmin, 
                       ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + ggtitle(paste("CN=4, weight=", 
                                                                       round(tetraploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) + 
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"),  list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                      diploid.genome.fraction=diploid.genome.fraction,
                                                                                      triploid.genome.fraction=triploid.genome.fraction,
                                                                                      tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  
  pdf(paste0("./Plots/Neutral_evolution/", i, "_fit.pdf"), width=4, height=3)
  print(p)
  dev.off()
  
  save(p, file=paste0("./RScripts/RData/Neutral_evolution/Plot_fit_", i, ".RData"))
}


# parameters


for(i in tumors.subclonal.evol){
  print(i)
  patient.id <- i
  fits <- list.files("Processed_data/Neutral_fits/", pattern=paste0(i, ".csv"), full.names = T)
  if(length(fits)==0){next}
  fits <- read.csv(fits)
  
  parms.of.interest <- c("par_delta", "par_mu", "par_n_clonal")
  parameters.to.plot <- data.frame(parms_i = rep(parms.of.interest, each = length(parms.of.interest)),
                                   parms_j = rep(parms.of.interest, length(parms.of.interest)),
                                   i = rep(1:length(parms.of.interest), each = length(parms.of.interest)),
                                   j = rep(1:length(parms.of.interest), length(parms.of.interest)))
  
  plot.design <- data.frame(Parameter = parms.of.interest,
                            Min = c(0, 0.1, 0),
                            Max = c(0.99, 20, 10000),
                            Label = c("delta/lambda", "mu", "n_clonal"))
  
  
  # fill up
  p <- list()
  
  for(ii in 1:length(parms.of.interest)){
    for(j in 1:length(parms.of.interest)){
      if(ii == j){
        p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits, aes(x=.data[[parms.of.interest[ii]]])) +
          geom_histogram() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                                name = unlist(plot.design[ii, "Label"])) + 
          scale_y_continuous() +
          theme(axis.text.y=element_blank(),  #remove y axis labels
                axis.ticks.y=element_blank(),  #remove y axis ticks
                axis.title.y = element_blank()
          )
      }else if(ii < j){
        p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits, aes(x=.data[[parms.of.interest[ii]]],
                                                                      y=.data[[parms.of.interest[j]]])) +
          geom_point() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                            name = unlist(plot.design[ii, "Label"])) + 
          scale_y_continuous(limits=unlist(plot.design[j,c("Min", "Max")]),
                             name = unlist(plot.design[j, "Label"]))
      }else{
        p[[length(parms.of.interest)*(ii-1) + j]] <- ggplot(fits, aes(x=.data[[parms.of.interest[ii]]],
                                                                      y=.data[[parms.of.interest[j]]])) +
          geom_density2d() + scale_x_continuous(limits=unlist(plot.design[ii,c("Min", "Max")]),
                                                name = unlist(plot.design[ii, "Label"])) + 
          scale_y_continuous(limits=unlist(plot.design[j,c("Min", "Max")]),
                             name = unlist(plot.design[j, "Label"]))
      }
      if(j!=length(parms.of.interest)){
        p[[length(parms.of.interest)*(ii-1) + j]] <-  p[[length(parms.of.interest)*(ii-1) + j]] + theme(axis.text.x=element_blank(),
                                                                                                        axis.ticks.x=element_blank(),
                                                                                                        axis.title.x = element_blank())
      }
      if(ii!=1){
        p[[length(parms.of.interest)*(ii-1) + j]] <- p[[length(parms.of.interest)*(ii-1) + j]] + theme(axis.text.y=element_blank(),
                                                                                                       axis.ticks.y=element_blank(),
                                                                                                       axis.title.y = element_blank())
      }
    }
  }
  
  pdf(paste0("./Plots/Neutral_evolution/", i, "_params.pdf"), width=6, height=6)
  print(egg::ggarrange(plots = p, nrow=length(parms.of.interest), byrow=F)) 
  dev.off()
}
