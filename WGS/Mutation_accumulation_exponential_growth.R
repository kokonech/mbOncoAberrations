#########################################################################################################################################
## Compute the accumulation of neutral mutations during an exponentially growing system
#########################################################################################################################################
## Functions

#' Non-critical clone size distribution (exact).
#'
#' @description Exact probability to grow from a clone of size "a" to a clone of size "b" within time "t" according to a non-critical birth-death process.
#'
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t time
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @return The probability that a clone of size a grows to size "b" within "t".
#' @references  Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @examples
#' density.a.b.exact(1, 0, 10, 1, 2)
#' @export

density.a.b.exact <- function(lambda, delta, t, a, b){
  if(a==1){
    if(b==0){
      .alpha(lambda, delta, t)
    }else{
      (1-.alpha(lambda, delta, t))*(1-.beta(lambda, delta, t))*
        .beta(lambda, delta, t)^(b-1)
    }
  }else{
    if(b==0){
      .alpha(lambda, delta, t)^a
    }else{
      sum(sapply(seq(0,min(a,b)), function(j){
        choose(a, j)*choose(a+b-j-1, a-1)*.alpha(lambda, delta, t)^(a-j)*.beta(lambda, delta, t)^(b-j)*(1-.alpha(lambda, delta, t)-
                                                                                                          .beta(lambda, delta, t))^j
      }))
    }
  }
}

.alpha <- function(lambda, delta, t){
  (delta*exp((lambda - delta)*t) - delta)/
    (lambda*exp((lambda - delta)*t) - delta)
}

.beta <- function(lambda, delta, t){
  (lambda*exp((lambda - delta)*t) - lambda)/
    (lambda*exp((lambda-delta)*t) - delta)
}

#' Clone size distribution in a noncritical birth-death process (approximate).
#' @description Probability of a clone of size "a" to grow to size "b" within "t" according to a noncritical linear birth-death process. 
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t time
#' @param a clone size at t=0
#' @param b clone size at t=t
#' @param mode "density" if density distribution is to be returned , "cumulative" if cumulative distribution is to be returned. Defaults to "cumulative"
#' @param approx Approximation to be used. Defaults to "highnumbers"; i.e. the distribution is approximated with a gamma distribution if `a` and `b` are large.
#' @details
#' If `approx="highnumbers"`, the function is approximated with a \eqn{\Gamma}-distribution if `a+b>100` and `mode="density"` or if `a+b>10` and `mode="cumulative"`. The \eqn{\Gamma}-distribution is parametrized with 
#' \eqn{shape = \mu^2/\sigma, scale = \sigma/\mu}, where
#' \eqn{\mu = a e^{(\lambda - \delta)t}, \sigma = a \frac{\lambda + \delta}{\lambda - \delta}e^{(\lambda - \delta)t}(e^{(\lambda - \delta)t}-1)}
#' @references Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @return The probability of growing from size a to size b within t. The Function switches between the exact solution and an approximate solution according to a parametrized gamma distribution.
#' @export
#' 
p.a.b <- function(lambda, delta, t, a, b, mode="cumulative", approx="highnumbers"){
  if(a==0){
    if(mode=="density"){
      if(b==0){
        return(1)
      }else{
        return(0)
      }
    }else{
      return(1)
    }
  }
  if(mode=="density"){
    if(a+b<=100 & approx=="highnumbers"){
      density.a.b.exact(lambda, delta, t, a, b)
    }else{
      ## Approximate with gamma-distribution, parametrized by mean and variance
      mean.g <- a*exp((lambda-delta)*t)
      var.g <- a*(lambda+delta)/(lambda-delta)*exp((lambda-delta)*t)*(exp((lambda-delta)*t)-1)
      dgamma(b, shape=mean.g^2/var.g, scale=var.g/mean.g)
    }
    
  }else if(mode=="cumulative"){
    if(a+b<=10 & approx=="highnumbers"){
      ## cumulative distribution
      sum(sapply(seq(0, b), function(b){
        ## Pab according to Bailey, 1964
        density.a.b.exact(lambda, delta, t, a, b)
      }))
      
    }else{
      ## Approximate with gamma-distribution, parametrized by mean and variance
      mean.g <- a*exp((lambda-delta)*t)
      var.g <- a*(lambda+delta)/(lambda-delta)*exp((lambda-delta)*t)*(exp((lambda-delta)*t)-1)
      pgamma(b, shape=mean.g^2/var.g, scale=var.g/mean.g)
    }
  }
}


#' Mutation accumulation in a growing tissue
#' 
#' @description Expected number of neutral mutations that are present in at least `n.min` cells at `t.end` in an exponentially growing or contracting tissue.
#' @param lambda proliferation rate 
#' @param delta loss rate
#' @param t.end time
#' @param mu mutation rate per cell division
#' @param n.min minimal clone size at `t.end`; can be a value or a vector
#' @param N0 initial population size
#' @param N final population size
#' @param mode if "approx" the sum is approximated by integration. If "exact" the sum is exactly computed for clone sizes between 1 and 10 but beyond that also approximated.
#' @details The expected number of mutations present in at least `n.min` cells is computed as \cr
#' \eqn{M(n_\mathrm{min})=\sum_{n_\mathrm{min}}^N \mu \lambda \int_0^t e^{(\lambda - \delta)(t-t')} P(1,n_{min},t-t')dt'}, which is approximated to \cr
#' \eqn{M(n_\mathrm{min})\approx \mu \lambda \int_0^t e^{(\lambda - \delta)(t-t')} \frac{ P(1,N,t-t') -  P(1,n_\mathrm{min},t-t')}{\log y(t-t')}dt'}, \cr where
#' \eqn{y(t) = \frac{\lambda e^{(\lambda - \delta)t} - \lambda}{\lambda e^{(\lambda - \delta)t}-\delta}}, if `mode=="approx"` or if \eqn{n_\mathrm{min} \le 10}
#' @references  Bailey, NTJ (1964). The elements of stochastic processes with applications to the natural sciences, Wiley (New York).
#' @return The expected number of mutations present in at least `n.min` cells at `t.end` in an exponentially growing tissue. The function assumes that mutations are continuously acquired at a constant rate.
#' @export

mutations.noncritical.bd <- function(lambda, delta, t.end, mu, n.min, N0=1, N=N, mode="approx"){
  
  if(mode=="exact"){
    
    res <- sapply(n.min, function(n.min){
      if(n.min <=10){
        total <- .approx.count(mu, lambda, delta, n.min=11, n.max=100*max(N,N0), t.end, N0) + .exact.count(t, mu, lambda, delta, 10, t.end, N0)
        if(n.min > 1){
          res <- total - .exact.count(t, mu, lambda, delta, n.min-1, t.end)
        }else{
          res <- total 
        }
        return(res)
      }else{
        .approx.count(mu, lambda, delta, n.min=1, n.max=Inf, t.end, N0) - .approx.count(mu, lambda, delta, n.min=1, n.max=n.min, t.end, N0)
      }
      
    })
    return(res)
    
  }
  ## Compute the number of mutations that are present in at least 1 cell and at most n.min cells
  ## The sum necessary in order to compute the cumulative distribution, is here replaced by integration.
  sapply(n.min, function(n.min){
    
    res <- .approx.count(mu, lambda, delta, n.min=1, n.max=100*max(N0, N), t.end, N0) - .approx.count(mu, lambda, delta, n.min=1, n.max=n.min, t.end, N0) 
    return(res)
  })
  
}

## Exact number of mutations present in at least n cells at time t
.exact.count <- function(t, mu, lambda, delta, n, t.end, N0){
  sum(sapply(seq(1,n), function(n){
    integrand <- function(t, mu, lambda, delta, n){
      mu*lambda*N0*exp((lambda - delta)*t)*(density.a.b.exact(lambda, delta, t.end-t, 1, n)) 
      
    }
    ## total number of mutations acquired during exponential growth that survived:
    res <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n=n)$value
    res
  }))
}

## Approximate number of mutations present in at least n.min cells at time t

.approx.count <- function(mu, lambda, delta, n.min, n.max, t.end, N0){
  integrand <- function(t, mu, lambda, delta, n.min, n.max){
    mu*lambda*N0*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, n.max) - 
                                                                               density.a.b.exact(lambda, delta, t.end-t, 1, n.min))
    
  }
  total <- integrate(integrand, lower=0, upper=t.end, mu=mu, lambda=lambda, delta=delta, n.min=n.min, n.max=n.max, rel.tol = .Machine$double.eps^0.1)$value
  total
}


#' Mutation accumulation during exponential expansion with clonal selection. 
#' 
#' @param mu mutation rate per cell division
#' @param N population size
#' @param lambda proliferation rate 
#' @param delta loss rate 
#' @param t.end end point 
#' @param t.s time point at which selective advantage is acquired.
#' @param s selective advantage
#' @param b minimal clone size of interest. Number or vector. 
#' @return This function returns an approximation by first computing the distribution at the transition time within intervals, then averaging the fate of each interval during homeostasis and adding newly acquired mutations in a scenario where a subpopulation is under positive selection. Returns the number of mutations present in at least `b` cells.
#' @export

mutational.burden.selection.expansion=function(mu,lambda,delta,s,t.s,t.end, b){
  if (s<0){
    message("No selective advantage (s<0)")
  }
  
  ## final tissue size
  N <- exp((lambda - delta)*t.end) - exp((lambda - delta)*(t.end - t.s)) + exp((lambda - delta*s)*(t.end - t.s))
  
  ## size of the selected clone at t.end
  sel.size <- exp((lambda - s*delta)*(t.end-t.s))
  
  
  ## Compute the number of mutations that are present in at least 1 cell and at most n.min cells
  ## The sum necessary in order to compute the cumulative distribution, is here replaced by integration.
  mutations.in.selected.clone.prior.t.s <- sapply(b, function(n.min){
    integrand <- function(t, mu, lambda, delta, n){
      p.mut.in.sel <- exp((lambda - delta)*(t.s - t))/exp((lambda - delta)*t.s)
      if(n==0){
        p.mut.in.sel*mu*lambda*exp((lambda - delta)*t)
      }else{
        p.mut.in.sel*mu*lambda*exp((lambda - delta)*t)*((density.a.b.exact(lambda, delta, t.end-t, 1, N*100) -                                                                                                   
                                                           density.a.b.exact(lambda, delta, t.end-t, 1, n))/log(.beta(lambda, delta, t.end-t)))
        
      }
    }
    ## total number of mutations acquired during exponential growth that survived:
    res <- integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=max(0,n.min-sel.size))$value
    return(res)
  })
  
  mutations.not.in.selected.clone.prior.t.s <- sapply(b, function(n.min){
    integrand <- function(t, mu, lambda, delta, n){
      p.mut.in.sel <- exp((lambda - delta)*(t.s - t))/exp((lambda - delta)*t.s)
      (1-p.mut.in.sel)*mu*lambda*exp((lambda - delta)*t)/log(.beta(lambda, delta, t.end-t))*(density.a.b.exact(lambda, delta, t.end-t, 1, 100*N) - 
                                                                                               density.a.b.exact(lambda, delta, t.end-t, 1, n))
      
    }
    
    ## total number of mutations acquired during exponential growth that survived:
    res <- integrate(integrand, lower=0, upper=t.s, mu=mu, lambda=lambda, delta=delta, n=n.min)$value
    return(res)
  })
  
  
  
  mutations.before.t.s <- mutations.in.selected.clone.prior.t.s + mutations.not.in.selected.clone.prior.t.s
  
  
  mutations.from.selected.clone.after.t.s <- sapply(b, function(n.min){
    mutations.noncritical.bd(lambda = lambda, delta = delta *s, t.end = t.end - t.s, 
                             mu = mu, n.min = n.min, N = sel.size)})
  mutations.from.founder.clone.after.t.s <- exp((lambda - delta)*t.s)*sapply(b, function(n.min){
    mutations.noncritical.bd(lambda = lambda, delta = delta, t.end = t.end - t.s, 
                             mu = mu, n.min = n.min, N = exp((lambda - delta)*(t.end - t.s)), mode="exact")})
  
  mutations.before.t.s + mutations.from.selected.clone.after.t.s + mutations.from.founder.clone.after.t.s
  
}

#########################################################################################################################################
## Simulate full cumulative VAF distribution of an exponential expansion with n clonal mutation and subclonal mutation accumulation according to neutral mutation acquisition
## as obtained with next generation sequencing
## parms: a vector consisting of the following elements: 
##  - n_clonal: # clonal mutations per haploid genome
##  - delta: loss rate (relative to division rate)
##  - mu: mutation rate (per cell division)
##  - lambda_s: driver division rate (only required if running in the selection mode)
##  - delta_s: driver loss rate (only required if running the selection mode)
##  - t_s: acquisition time of the driver mutation (only required if running in the selection mode)
## ploidy: the ploidy of the genome fraction we're looking at (1, 2, 3, ...)
## depth: the sequencing depth
## selection: does the model account for a subclone under positive selection? default F
## Note that the computation is re-scaled to the division rate, such that lambda=1!

## previously named as rmodel

SimulateVAFDistr <- function(parms, ploidy, depth, purity=1, selection = F){ 

  ## start with the clonal mutations. For each clonal mutation, sample the sequencing depth from a Poisson distribution
  coverages <- rpois(n=round(parms$n_clonal*ploidy), lambda=depth)
  coverages[coverages==0] = 1
  ## sample the VAFs of the clonal peak according to a binomial distribution. The purity and the ploidy go into the sampling probability
  vafs.clonal <- rbinom(n=round(parms$n_clonal*ploidy), size = coverages, prob = purity/(ploidy*purity + 2*(1-purity)))/coverages

  ## simulate a histogram of bins between 0.05 and 1 with bin size = 0.05; assume a tumor size of 10^9 cells

  parms$lambda <- 1
  
  if(selection){
    parms$delta_s <- parms$delta
    t.max <- uniroot(function(t){
      exp((parms$lambda - parms$delta)*t) - exp((parms$lambda - parms$delta)*(t - parms$t_s)) + exp((parms$lambda - parms$s*parms$delta)*(t - parms$t_s)) - 10^9},
      interval = c(0, 1+9*log(10)/(parms$lambda - parms$delta))
    )$root
    
    vafs.subclonal <- mutational.burden.selection.expansion(mu = parms$mu*ploidy, lambda = parms$lambda, delta = parms$delta, s = parms$s, t.s = parms$t_s, t.end = t.max, b = seq(0.05, 0.99, 0.05)*10^9)
  }else{
    vafs.subclonal <- mutations.noncritical.bd(mu = parms$mu*ploidy, lambda = parms$lambda, delta = parms$delta, t.end = 9*log(10)/(parms$lambda - parms$delta), N = 10^9, n.min = seq(0.05, 0.99, 0.05)*10^9)
  }
    
  # convert to histogram
  vafs.subclonal <- vafs.subclonal - c(vafs.subclonal[-1], 0)
  
  # convert to vafs
  vafs.subclonal <- unlist(apply(rbind(seq(0.05, 0.99, 0.05), vafs.subclonal), 2, function(vaf){
    if(vaf[2] < 1){
      vaf[2] <- sample(x = c(1, 0), size = 1, replace = F, prob = c(vaf[2], 1 - vaf[2]))
    }
    rep(vaf[1], round(vaf[2]))
  }))

  if(any(is.na(vafs.subclonal) | is.infinite(vafs.subclonal))){
    return(NA)
  }
  ## sample the sequencing depth for each VAF from a Poisson distribution
  coverages <- rpois(n=length(vafs.subclonal), lambda=depth)
  coverages[coverages==0] <- 1
  ## Sample the subclonal VAFs from a binomial distribution
  vafs.subclonal <- rbinom(length(vafs.subclonal), size = coverages, prob = vafs.subclonal*purity/(ploidy*purity + 2*(1-purity)))/coverages
 if(any(is.na(vafs.subclonal))){print(purity)} 
  return(c(vafs.clonal, vafs.subclonal))
}


## tests:
# parms = list(n_clonal = 1000, delta = 0.9, mu = 10, delta_s = 0.8, lambda_s = 1, t_s = 30)
# test.sel = SimulateVAFDistr(parms, 2, 90, selection=T)
# hist(test.sel, breaks=100)
# test = SimulateVAFDistr(parms, 2, 90, selection=F)
# hist(test, breaks=100)
# plot(1/seq(0.01, 1, 0.01), sapply(seq(0.01, 1, 0.01), function(x){sum(test>x)}), xlab="1/VAF", ylab="# Mutations")
# lines(1/seq(0.01, 1, 0.01), sapply(seq(0.01, 1, 0.01), function(x){sum(test.sel>x)}))
