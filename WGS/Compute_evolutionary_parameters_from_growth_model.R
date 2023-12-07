##############################################################################################################################################
### The fit was done with pyABC. In order to reproduce it, you need to run 
### Neutral_fit_pre_clonal_and_clonal.py, executing Neutral_fit_pre_clonal_and_clonal.R
##############################################################################################################################################
library(ggbeeswarm)

## subset of tumors used to learn the parameters of tumor growth
subset <- sample.information[!is.na(sample.information$Use.for.neutral.fitting) & sample.information$Use.for.neutral.fitting=="yes" ,]

##### take the mutation rate from the fit of tumor initiation
fits <- read.csv("./Processed_data/Dynamic_model/Expansion_decay_continuous_evol.csv")
mutation.rate <- c(mean(fits$par_mu*2), sd(fits$par_mu*2))

## from these, extract estimates of the effective mutation rate and subsequently of the loss rate 
parms <- data.frame(deltas = c(), effective.mutation.rates = c(), division.rate = c(), mutation.rate.per.day=c(),
                    n.clonal = c(), Model = c(), Stat = c(), ID = c())

for(i in rownames(subset)){
  
  fits <- read.csv(paste0("./Processed_data/Neutral_fits/", i, ".csv"))
  
  n.clonal <- c(mean(fits$par_n_clonal), sd(fits$par_n_clonal))
  ## in the growth model, the mutation rate is per 2 cells, but per haploid genome, thus take it as it is (*2/2=1) per cell
  ## The effective mutation rate is per effective division. 
  ## store mean, sd
  effective.mutation.rate <- c(mean(fits$par_mu/(1-fits$par_delta)), sd(fits$par_mu/(1-fits$par_delta)))
  
  ## From this, compute delta; mean and sd
  delta <- 1 - mutation.rate[1]/effective.mutation.rate[1] 
  ## we obtain the error by propagation of uncertainties
  delta[2] <- abs(1/effective.mutation.rate[1] *mutation.rate[2] + mutation.rate[1]/effective.mutation.rate[1]^2 *effective.mutation.rate[2])
  
  mutational.burden.at.mrca <- mutation.time.mrca[i,]$Mean
  ## assume that mutation times are roughly normally distributed. Thus the standard deviation would correspond to 1/3.84 of the 95% CI
  mutational.burden.at.mrca[2] <- (mutation.time.mrca[i,]$Max - mutation.time.mrca[i,]$Min)/(2*1.96)
  age <- subset[i, "Age"]*365
  
  n.generations <- mutational.burden.at.mrca[1]*2/mutation.rate[1] + 9*log(10)/(1-delta[1])
  n.generations[2] <- 2/mutation.rate[1]*mutational.burden.at.mrca[2]*2 +
    mutational.burden.at.mrca[1]*2/mutation.rate[1]^2*mutation.rate[2] + 9*log(10)/(1-delta[1])^2*delta[2]
  
  ## generations until MRCA 
  n.generations.1 <- mutational.burden.at.mrca[1]*2/mutation.rate[1]
  n.generations.1[2] <- 2/mutation.rate[1]*mutational.burden.at.mrca[2]*2+
    mutational.burden.at.mrca[1]*2/mutation.rate[1]^2*mutation.rate[2] 
  n.generations.2 <- 9*log(10)/(1-delta[1])
  n.generations.2[2] <-  9*log(10)/(1-delta[1])^2*delta[2]
  
  ## t.total = age + pregnancy
  t.total <- age + 250
  ## initiation is the number of generations until tumor initiation divided by the total number of generations times the total time
  t.init <- mutational.burden.at.mrca*2/mutation.rate[1]/n.generations[1]*t.total
  t.init[2] <- t.total*(n.generations.1[1]*n.generations.2[1]+n.generations.2[1]*n.generations.1[2])/(n.generations.1[1]+n.generations.2[1])^2
  
  division.rate <- c(n.generations[1]/(age + 250), n.generations[2]/(age + 250))
  
  ## the mutation rate per day is the product of the division rate and the mutation rate
  mut.per.day <-  division.rate[1]*mutation.rate[1]
  mut.per.day[2] <- division.rate[1]*mutation.rate[2] + division.rate[2]*mutation.rate[1]
  
  parms <- rbind(parms,
                 data.frame(deltas = delta[1], effective.mutation.rates = effective.mutation.rate[1], 
                            division.rate = division.rate[1], mutation.rate.per.day=mut.per.day[1],
                            n.clonal = n.clonal[1], Stat = "mean", ID = i),
                 data.frame(deltas = delta[2], effective.mutation.rates = effective.mutation.rate[2], 
                            division.rate = division.rate[2], mutation.rate.per.day=mut.per.day[2],
                            n.clonal = n.clonal[2], Stat = "sd", ID = i))
  
}

parms <- parms[!is.na(parms$mutation.rate.per.day),]
parms$Driver <- ifelse(sample.information[parms$ID,]$MYC==1 | sample.information[parms$ID,]$MYCN==1 |
                         sample.information[parms$ID,]$PRDM6==1, 1, 0)

real.time.parms <- parms
save(real.time.parms, file=paste0(rdata.directory, "Estimated_real_time_parms.RData"))

estimated.mutation.rate.per.day <- c(mean(parms$mutation.rate.per.day[ parms$Stat=="mean"], na.rm = T),
                                     sqrt(sum(parms$mutation.rate.per.day[parms$Stat=="sd"]^2))/
                                       length(parms$mutation.rate.per.day[ parms$Stat=="mean"]))

estimated.mutation.rate.per.day <- c(estimated.mutation.rate.per.day[1] - 1.96*estimated.mutation.rate.per.day[2],
                                     estimated.mutation.rate.per.day[1], 
                                     estimated.mutation.rate.per.day[1] + 1.96*estimated.mutation.rate.per.day[2])

save(estimated.mutation.rate.per.day, file=paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))

