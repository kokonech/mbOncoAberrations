##############################################################################################################################################
### Learn the dynamics of tumor initiation
##############################################################################################################################################
### Load libraries and settings

source("./RScripts/Settings.R")
library(NBevolution)
library(cdata)
##### observed data

load(paste0(rdata.directory, "Input_data_MB_initiation.RData"))

mySumStatData <- list(P.MRCA=P.MRCA, P.ECA=P.ECA)

load(paste0(rdata.directory, "MRCA_timing.RData"))
load(paste0(rdata.directory,"Clonal_mutations_different_CNs.RData"))
load(paste0(rdata.directory,"Vafs_all_tumors.RData"))

##### take the mutation rate from the fit of tumor initiation
fits <- read.csv("./Processed_data/Dynamic_model/Expansion_decay_continuous_evol.csv")
mutation.rate <- c(mean(fits$par_mu*2), sd(fits$par_mu*2))

load(paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))

##############################################################################################################################################
### The fit was done with pyABC. In order to reproduce it, you need to create the input data by running Input_data_MB.R
### Then run Expansion_decay_continuous_evol.py. You need the files Expansion_decay_continuous_evol.R and Expansion_decay_2_hits_continuous_evol.R

##############################################################################################################################################
## Expansion + decay, ECA

fits <- read.csv("./Processed_data/Dynamic_model/Expansion_decay_continuous_evol.csv")

parameter.samples <- data.frame(N=fits$par_N, delta1=fits$par_delta1, muD1=fits$par_muD1, 
                                muD2=fits$par_muD2, mu=fits$par_mu,
                                delta2=fits$par_delta2, psurv=fits$par_psurv, r=fits$par_r)


model <- "contraction"

source("./RScripts/ABC_fit.R")

## Simulate incidence curves
sim.mrca <- matrix(NA, nrow = 100, ncol = length(mySumStatData$P.MRCA$Density))
sim.eca <- matrix(NA, nrow = 100, ncol = length(mySumStatData$P.ECA$Density))

for(i in 1:100){
  parms <- as.list(parameter.samples[i,])
  
  ## rescale measured mutations into time (gamma-distributed)
  time.points <- sapply(sort(mySumStatData$P.MRCA$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
  
  ## then at each time point, take the probability of having acquired the 2 hits
  P <- sapply(time.points, function(t){
    P.2nd.driver(parms, t, model=model)
  })
  
  sim.mrca[i,] <- P
  
  # ECA
  ## rescale measured mutations into time; MRCA (gamma distributed)
  time.points <- sapply(sort(mySumStatData$P.MRCA[rownames(mySumStatData$P.ECA),]$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
  
  ## and ECA
  time.points.eca <- sapply(sort(mySumStatData$P.ECA$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
  
  ## then sample a time point for the ECA
  x <- runif(length(time.points), 0, 1)
  
  ## Then, choose the associated probability
  tmp <- function(x, parms, t1, t2){
    P.eca.at.t(parms, t1, t2, model=model) - x
  }
  
  t <- apply(rbind(time.points, x), 2, function(x){
    uniroot(tmp, c(0, x[1]), x=x[2], parms=parms, t2 = x[1])$root
  })
  
  incidence <- sapply(sort(time.points.eca), function(x){
    sum(t <= x)
  })/length(time.points.eca)
  
  sim.eca[i,] <- incidence
}

to.plot <- data.frame(x = P.MRCA$Density/3.3/10^3,
                      data = P.MRCA$P*10^-5, 
                      sd = (P.MRCA$P.upper - P.MRCA$P.lower)/2/1.95*10^-5,
                      lower = apply(sim.mrca, 2, quantile, 0.025), 
                      upper = apply(sim.mrca, 2, quantile, 0.975),
                      Event = rep("MRCA", nrow(P.MRCA)))

to.plot.eca <- data.frame(x = P.ECA$Density/3.3/10^3,
                      data = P.ECA$P, 
                      sd = (P.ECA$P.upper - P.ECA$P.lower)/2/1.95,
                      lower = apply(sim.eca, 2, quantile, 0.025), 
                      upper = apply(sim.eca, 2, quantile, 0.975),
                      Event = rep("ECA", nrow(P.ECA)))

pdf(paste0(output.directory, "Model_fit_decay_ECA.pdf"), width=4, height=3, useDingbats = F)

print(ggplot(to.plot.eca, aes(x=x, y = 100*data, ymin = (data-sd)*100, ymax = (data +sd)*100, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x, ymin = 100*lower, ymax = 100*upper), col=NA)  +
        scale_fill_manual(values=time.colors) + scale_color_manual(values=time.colors) + 
        scale_x_continuous(name="#SSNVs/Mb", limits=c(0, max(to.plot.eca$x))) + 
        scale_y_continuous(name = "% Medulloblastoma cases") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(12-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(12-2)*7/2/3.3/10^3, ymin = 0, ymax = 100), fill="grey", 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F)+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38-2)*7/2/3.3/10^3, ymin = 0, ymax = 100), fill="grey", 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F)+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38+52-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38+52-2)*7/2/3.3/10^3, ymin = 0, ymax = 100), fill="grey", 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F))

dev.off()

##############################################################################################################################################
## Expansion + decay, MRCA

pdf(paste0(output.directory, "Model_fit_decay_MRCA.pdf"), width=4, height=3, useDingbats = F)

print(ggplot(to.plot, aes(x=x, y = data, ymin = (data-sd), ymax = (data +sd), col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x, ymin = lower, ymax = upper), col=NA)  +
        scale_fill_manual(values=time.colors) + scale_color_manual(values=time.colors) + 
        scale_x_continuous(name="#SSNVs/Mb") + 
        scale_y_continuous(name = "Medulloblastoma incidence") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(12-2)*7/2/3.3/10^3, 
                                  xmax = estimated.mutation.rate.per.day[3]*(12-2)*7/2/3.3/10^3, ymin = 0, ymax = 1.2*10^-5), fill="grey",
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  inherit.aes = F)+
        geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38-2)*7/2/3.3/10^3, 
                                  xmax = estimated.mutation.rate.per.day[3]*(38-2)*7/2/3.3/10^3, ymin = 0, ymax = 1.2*10^-5), fill="grey",
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  inherit.aes = F)+
        geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38+52*10-2)*7/2/3.3/10^3, 
                                  xmax = estimated.mutation.rate.per.day[3]*(38+52*10-2)*7/2/3.3/10^3, ymin = 0, ymax = 1.2*10^-5), fill="grey",
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  inherit.aes = F))

dev.off()

##############################################################################################################################################
## parameter estimates

## do all correlations
## compute selective advantage from survival probability
fits$par_s <- fits$par_delta2/(1-fits$par_psurv)
parameter.samples$s <- parameter.samples$delta2/(1-parameter.samples$psurv)
## compute geometric mean of mu1 and mu2
parameter.samples$muD1D2 <- log10(sqrt(10^parameter.samples$muD2*10^parameter.samples$muD1))

## parameters to plot
parameters <- c("par_N", "par_delta1", "par_delta2", "par_mu", "par_muD1", "par_muD2", "par_s", "par_r", "muD1D2")

## specify the variables I want to plot
meas_vars <- colnames(parameter.samples)

## a data frame of all combinations of its arguments

controlTable <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = F))

## rename the columns
colnames(controlTable) <- c("x", "y")

## add the key column
controlTable <- cbind(data.frame(par_key = paste(controlTable[[1]], controlTable[[2]]), stringsAsFactors = F), controlTable)

## create the new data frame
to.plot <- rowrecs_to_blocks(parameter.samples, controlTable)

## re-arrange with facet_grid
splt <- strsplit(to.plot$par_key, split=" ", fixed=TRUE)
to.plot$xv <- vapply(splt, function(si) si[[1]], character(1))
to.plot$yv <- vapply(splt, function(si) si[[2]], character(1))

to.plot$xv <- factor(as.character(to.plot$xv), meas_vars)
to.plot$yv <- factor(as.character(to.plot$yv), meas_vars)


## arrange manually

to.plot$xaxis <- F
to.plot$yaxis <- F
to.plot$xaxis[to.plot$yv == to.plot$xv[sqrt(length(unique(to.plot$par_key)))]] <- T
to.plot$yaxis[to.plot$xv==to.plot$xv[1]] <- T
to.plot$topm <- F
to.plot$rightm <- F
to.plot$topm[to.plot$yv == to.plot$xv[1]] <- T
to.plot$rightm[to.plot$xv==to.plot$xv[sqrt(length(unique(to.plot$par_key)))]] <- T

p <- list()

## introduce an artificial top row and right column

for(i in 1:(sqrt(length(unique(to.plot$par_key))))){
  p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
    theme_bw() + theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
  
}

for(i in unique(to.plot$par_key)){
  
  
  tmp <- to.plot[to.plot$par_key==i,]
  
  if(tmp$xv[1]=="psurv" | tmp$yv[1]=="psurv"){next}
  
  if(tmp$xv[1]==tmp$yv[1]){
    p[[length(p)+1]] <- ggplot(tmp, aes(x=x)) + 
      geom_histogram() + scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(legend.position = "none")
    
  }else{
    p[[length(p)+1]] <- ggplot(tmp, aes(x=x, y=y)) + 
      geom_density_2d_filled(col=NA, contour_var = "ndensity", aes( fill = ..level..)) + 
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
  }
  
  ## top-row and right column: adjust margins differently
  if(tmp$rightm[1] & tmp$topm[1]){
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }else if(tmp$rightm[1]){
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }else if(tmp$topm[1]){
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }else{
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }

  if(tmp$xaxis[1]==F){
    p[[length(p)]] <-  p[[length(p)]] +  theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank())
  }
    
  if(tmp$yaxis[1]==F){
    p[[length(p)]] <-  p[[length(p)]] +  theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank())
  }
  
  if(tmp$rightm[1]){
    p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
            legend.position = "none")
    
  }
  
}

pdf(paste0(output.directory, "Parameters_decay_2D.pdf"), width=9, height=9, useDingbats = F)

print(ggpubr::ggarrange(plotlist=p, nrow=10, ncol=10, align="hv"))

dev.off()


##############################################################################################################################################
##### Plot the retraction of neuroblasts and the expansion of the first clone over approximate week of pregnancy

mutation.count <- seq(0, max(P.MRCA$Density))
N.sim <- matrix(0, nrow=1000, ncol=length(mutation.count))

hdi.parameters <- hdi(parameter.samples)

for(i in 1:nrow(parameter.samples)){
  delta1 <- parameter.samples[i, "delta1"]
  delta2 <- parameter.samples[i, "delta2"]
  mu <- parameter.samples[i, "mu"]
  N <- 10^parameter.samples[i, "N"]
  
  
  turning.point <- log(N)/(1-delta1)
  N.sim[i,] <- sapply(mutation.count, function(m){
    ## convert to time
    t <- m/mu
    if(t <=turning.point){
      return(exp((1-delta1)*t))
    }else{
      return(N*exp((1-delta2)*(t-turning.point)))
    }
  })
  
}


N.sim <- N.sim[rowSums(N.sim)!=0,]
test <- apply(N.sim, 2, density)

#N.sim.cells <- data.frame(Mutations=rep(mutation.count, nrow(N.sim)), N=as.vector(t(N.sim)), Sim=rep(1:nrow(N.sim), each=ncol(N.sim)))

N.sim.norm <- t(apply(N.sim, 1, function(x){x/max(x)}))
N.sim.cells <- data.frame(Mutations=mutation.count, ymin=apply(N.sim.norm, 2, quantile, p=0.025), ymax=apply(N.sim.norm, 2, quantile, p=0.975))
N.sim.cells$t <- N.sim.cells$Mutations/3.3/10^3


## Simulate evolution of M1 cells

M1.cells <- function(parms, t){
  delta1 <- parms$delta1
  N <- 10^parms$N
  delta2 <- parms$delta2
  mu1 <- 10^parms$muD1
  r <- parms$r
  lambda1 <- 1
  lambda2 <- 1
  
  turning.point <- log(N)/(1-delta1)
  
  if(t <=  turning.point){
    M1 <- mu1*lambda1*t*exp((lambda1-delta1)*t)
  }else if(t>turning.point){
    M1 <- mu1*lambda1*t*N*exp((lambda2 - delta2*r)*(t-turning.point)) +
      mu1*lambda2/(delta2*(r-1))*N*exp((lambda2-delta2*r)*t)/(exp((lambda2-delta2)*turning.point))*(exp(delta2*(r-1)*t)-
                                                                                                      exp(delta2*(r-1)*turning.point))
  }
  M1
}

M1 <- t(apply(parameter.samples, 1, function(parms){
  parms <- as.list(parms)
  cells <- sapply(mutation.count, function(x){
    M1.cells(parms, x/parms$mu)
  })
}))

M1 <- t(apply(M1,1,function(x){x/max(x)}))
M1 <- apply(M1, 2, quantile, p=c(0.025, 0.975))
M1 <- as.data.frame(t(M1))
M1$t <- mutation.count/3.3/10^3

N.sim.cells$t <- N.sim.cells$Mutations/3.3/10^3


pdf(paste0(output.directory, "First_clone.pdf"), width=4, height=3, useDingbats = F)

p <- ggplot(M1, aes(x=t, ymin=`2.5%`, ymax=`97.5%`)) + geom_ribbon(fill="darkgreen") + 
  geom_ribbon(data=N.sim.cells, aes(x=t, ymin=ymin, ymax=ymax), fill="grey")+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(12-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(12-2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey",
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F)+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38-2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey",
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F)+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38+52*10-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38+52*10-2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey",
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Cell count") +
  scale_x_continuous(name="# SNVs/Mb") 

print(p)

dev.off()

