## Requires:
## Functions (modelsol.R, reparlogp.R, SimsCodesFort.R)
## Source (priorsetup.R) contains constant definitions for the prior
## Load packages (QZ, FKF)
## Load Estimated Parameters from unpermuted data (StdPvec.Rdata)

## This file is meant to run on the head node using the package BatchJobs

library(batchtools)
pvec = readRDS("update2022/swlong-run-pars.rds")

makeExperimentRegistry(
  file.dir = 'SWSimulateEstimate2022', 
  packages = c('QZ', 'FKF', 'optimr'), 
  source = c(
    list.files("code/functions", full.names = TRUE),
    "code/priorsetup.R"
  )
)

generator <- function(job, data){
    df = generate(nobs = 3100, pvec = data$pvec)
    df = df[, -(1:1000)]
    return(df)
}

estimator <- function(job, data, instance, nestim) {
  estim.pred.gen(
    instance, 
    parm2start = data$pvec, 
    prior = data$prior, 
    nestim, 
    npred = 1000, 
    maxit = 3e4)
}

addProblem('Generator', 
           data = list(prior = prior, pvec = pvec), 
           fun = generator, 
           seed = 12345)
addAlgorithm('Estimator', fun = estimator)

estim.pts = seq(100, 1100, by = 20)
pdes = list(Generator = data.frame())
ades = list(Estimator = data.frame(nestim = estim.pts))
addExperiments(pdes, ades, repls = 100)

ids <- getJobPars()
ids <- tidyr::unnest_wider(ids, algo.pars)
ids <- dplyr::select(ids, job.id, nestim)
ids$chunk <- lpt(ids$nestim, n.chunks = 500) # 
submitJobs(ids, 
           resources = list(
             nodes = 1, ncpus = 1, walltime = '24:00:00',
             memory = "4Gb"
           ))
