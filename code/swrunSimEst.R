## Requires:
## Functions (modelsol.R, reparlogp.R, SimsCodesFort.R)
## Source (priorsetup.R) contains constant definitions for the prior
## Load packages (QZ, FKF)
## Load Estimated Parameters from unpermuted data (StdPvec.Rdata)

## This file is meant to run on the head node using the package BatchJobs

require(BatchExperiments)
load('DSGEs/StdPvec.Rdata')

reg = makeExperimentRegistry(id='SWSimulateEstimate2', packages=c('QZ','FKF'), src.dirs='DSGEs/code/functions',
    src.files=c('DSGEs/code/priorsetup.R'))

generator <- function(static){
    data = generate(nobs=3100, pvec=static$pvec)
    data = data[-(1:1000),]
    return(data)
}
estimator <- function(static, dynamic, nestim){
    out = estim.pred.gen(dynamic, parm2start=static$pvec, prior=static$prior, nestim, npred=1000, maxit=3e4)
}

addProblem(reg, id='Generator', static = list(prior=prior, pvec=pvec), dynamic = generator,seed=12345)
addAlgorithm(reg, id='Estimator',  fun=estimator)

estim.pts = seq(100,1100,by=20)
gen.design = makeDesign('Generator')
estim.design = makeDesign('Estimator', exhaustive = list(nestim=estim.pts))
addExperiments(reg, repls=100, prob.designs=gen.design, algo.designs=estim.design)

chunked = chunk(1:3000, n.chunks=500)
submitJobs(reg,  chunked, resources=list(nodes=1,ppn=1,walltime='24:00:00'))
