dinvgamma <- function(x, shape = 1, rate = 1, scale = 1/rate, log = FALSE) {
    # return( shape^rate / gamma(rate) * exp(-shape/x) * x^(-rate-1) )
    logval = rate*log(shape) - lgamma(rate) - shape/x - (rate+1)*log(x)
    if (log)
	return(logval)
    else
	return(exp(logval))
}

# pinvgamma: cumulative distribution function of the inverse-gamma distribution
pinvgamma <- function(q, shape = 1, rate = 1, scale=1/rate,
			lower.tail = TRUE, log.p = FALSE) {
    return( pgamma(1/q, shape, scale=scale, lower.tail=!lower.tail, log.p=log.p) )
}

rinvgamma <- function(n, shape = 1, rate = 1, scale=1/rate) {
    return( 1 / rgamma(n, shape, scale=scale) )
}



repar <- function(inparv, prior, bound.2.unbound=FALSE){
    ## repar <- function(inparv, R, prior.type, direction){
    ## PURPOSE: performs reparameterisation and computes log prior
    ## -----------------------------------------------------
    ## INPUTS: 
    ##         inparv = parameter vector
    ##         prior = a list containing the name of the prior distribution (dist) and two hyper parameters (p1, p2)
    ## -----------------------------------------------------
    ## RETURNS: 
    ##         outparv =  parameter vector
    ##         lnprior = ln prior 
    ## -----------------------------------------------------
    npar = length(inparv)
    outparv = double(npar)
    lnprior = sum(log(prior$normconst))
    nb = !prior$bounded
    if(!bound.2.unbound){
        ## transformation from from unbounded to bounded
        ret = within(prior,
            {
                outparv[bounded] = (ub[bounded] - lb[bounded]) / (1 + exp(-inparv[bounded])) + lb[bounded]
                lnprior = lnprior + sum(dbeta(outparv[ndb], p1[ndb], p2[ndb], log=TRUE))
                
                ## Gamma
                ndx = nb & ndg
                outparv[ndx] = exp(inparv[ndx])
                lnprior = lnprior + sum(dgamma(outparv[ndg], p1[ndg], p2[ndg], log=TRUE))
                
                ## inverse Gamma
                ndx = nb & ndi
                outparv[ndx] = exp(inparv[ndx])
                lnprior = lnprior + sum(dinvgamma(outparv[ndi], p1[ndi], p2[ndi], log=TRUE))
                
                ## Gaussian
                ndx = nb & ndn
                outparv[ndx] = inparv[ndx]
                lnprior = lnprior +sum( dnorm(outparv[ndn], p1[ndn], p2[ndn], log=TRUE))
            }
            )
    }else{
        ## Transformation from bounded to unbounded (not used for estimation, only for retrieving the likelihood
        ## later)
        ret = within(prior,
            {
                outparv[bounded] = log(inparv[bounded]-lb[bounded]) - log(ub[bounded]-inparv[bounded])
                lnprior = lnprior + sum(dbeta(inparv[ndb], p1[ndb], p2[ndb], log=TRUE))
                
                ## Gamma
                ndx = nb & ndg
                outparv[ndx] = log(inparv[ndx])
                lnprior = lnprior + sum(dgamma(inparv[ndg], p1[ndg], p2[ndg], log=TRUE))
                
                ## inverse Gamma
                ndx = nb & ndi
                outparv[ndx] = log(inparv[ndx])
                lnprior = lnprior + sum(dinvgamma(inparv[ndi], p1[ndi], p2[ndi], log=TRUE))
                
                ## Gaussian
                ndx = nb & ndn
                outparv[ndx] = inparv[ndx]
                lnprior = lnprior +sum( dnorm(inparv[ndn], p1[ndn], p2[ndn], log=TRUE))
            }
            )
    }
    return(list(outparv=ret$outparv, lnprior=ret$lnprior))
}

reject <- function(n, fun, p1, p2, lb, ub, max.samp = 10*n){
    ## rejection samples using fun to generate random numbers
    iter = 0
    i = 1
    a = double(n)
    while(i < n+1){
      if(iter == max.samp) stop('Rejection sampling did not yield enough reasonable values. Increase max.samp')
        test = fun(1, p1, p2)
        if(test>lb && test<ub){
            a[i] = test
            i = i+1
        }
        iter = iter+1
    }
    return(a)
}

start.vals <- function(prior, nstarts){
    dist = prior$dist
    nparm = length(dist)
    p1 = prior$p1
    p2 = prior$p2
    lb = prior$lb
    ub = prior$ub
    bounded = prior$bounded
    out = matrix(NA, nparm, nstarts)
    out1 = matrix(NA, nparm, nstarts)
    for(idx in 1:nparm){
        if(bounded[idx]){
            out[idx,] = switch(dist[idx],
                   beta = reject(nstarts, rbeta, p1[idx], p2[idx], lb[idx], ub[idx]),
                   gamma =  reject(nstarts, rgamma, p1[idx], p2[idx], lb[idx], ub[idx]),
                   igamma = reject(nstarts, rinvgamma, p1[idx], p2[idx], lb[idx], ub[idx]),
                   gaussian = reject(nstarts, rnorm, p1[idx], p2[idx], lb[idx], ub[idx])
                   )
            out1[idx,] = log((out[idx,]-lb[idx])/(ub[idx]-out[idx,]))
        }else{
            switch(dist[idx],
                   beta = {
                       out[idx,] = rbeta(nstarts, p1[idx], p2[idx])
                       out1[idx,] = log(out[idx,]/(1-out[idx,]))
                   },
                   gamma = {
                       out[idx,] = rgamma(nstarts, p1[idx], p2[idx])
                       out1[idx,] = log(out[idx,])
                   },
                   igamma = {
                       out[idx,] = rinvgamma(nstarts, p1[idx], p2[idx])
                       out1[idx,] = log(out[idx,])
                   },
                   gaussian = {
                       out[idx,] = rnorm(nstarts, p1[idx], p2[idx])
                       out1[idx,] = out[idx,]
                   }
                   )
        }
    }
    return(list(bounded=out,unbounded=out1))
}

trans <- function(prior){
    ## takes in name, mean, sd, outputs standard hyperparams for R
    n = prior$dist
    m = prior$m
    sd = prior$sd
    nparm = length(n)
    p1 = double(nparm)
    p2 = double(nparm)
    lb = prior$lb
    ub = prior$ub
    normconst = rep(1, nparm)
    bounded = prior$bounded
    idx = n=='beta'
    for(idx in 1:nparm){
        switch(n[idx],
               beta = {
                   mb = m[idx]
                   vb = sd[idx]^2
                   p1[idx] = (mb^2 * (1-mb) - vb*mb) / vb
                   p2[idx] = p1[idx] * (1-mb) / mb
                   if(bounded[idx]) normconst[idx] = 1 / ( pbeta(ub[idx], p1[idx], p2[idx]) -
                                                              pbeta(lb[idx], p1[idx], p2[idx]) )
               },
               gamma = {
                   mg = m[idx]
                   vg = sd[idx]^2
                   p2[idx] = mg/vg
                   p1[idx] = mg * p2[idx]
                   if(bounded[idx]) normconst[idx] = 1 / ( pgamma(ub[idx], p1[idx], p2[idx]) -
                                                              pgamma(lb[idx], p1[idx], p2[idx]) )
               },
               igamma = {
                   mi = m[idx]
                   vi = sd[idx]^2
                   p1[idx] = mi^2 / vi +2
                   p2[idx] = mi * (p1[idx] - 1)
                   if(bounded[idx]) normconst[idx] = 1/ ( pinvgamma(ub[idx], p1[idx], p2[idx]) -
                                                             pinvgamma(lb[idx], p1[idx], p2[idx]) )
               },
               gaussian = {
                   p1[idx] = m[idx]
                   p2[idx] = sd[idx]
                   if(bounded[idx]) normconst[idx] = 1/ (pnorm(ub[idx], p1[idx], p2[idx]) -
                                                             pnorm(lb[idx], p1[idx], p2[idx]) )
               }
               )
    }
    prior$p1 = p1
    prior$p2 = p2
    prior$normconst = normconst
    return(prior)
}
