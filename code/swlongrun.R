## Requires:
## Functions (modelsol.R, reparlogp.R, SimsCodesFort.R)
## Source (priorsetup.R,dataload.R) contains constant definitions for the prior, loads data (below)
## Load packages (QZ, FKF)
## Load data (SWdataDM.Rdata)

## This file is meant to run on the head node using the package batchtools

library(batchtools)
makeRegistry('SWlong2022', 
             seed = 20190921,
             packages=c('QZ','FKF','optimr'), 
             source = c('code/functions/SimsCodesFort.R', 
                        'code/functions/modelsol.R',
                        'code/functions/reparlogp.R', 
                        'code/priorsetup.R',
                        'code/dataload.R'))


perms = 1:7


    
batchMap(estim.pred.wrap, perms, more.args = list(
  y = y, prior = prior, estim.obs = 1:200, pred.obs = 201:251,
  nstarts = 5, maxit = 1e6)
)

ids = findJobs()

# submitJobs(ids, resources=list(nodes = 1, ncpus=1, walltime='48:00:00'))
