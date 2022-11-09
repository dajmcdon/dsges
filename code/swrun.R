## Requires:
## Functions (modelsol.R, reparlogp.R, SimsCodesFort.R)
## Source (priorsetup.R,dataload.R) contains constant definitions for the prior, loads data (below)
## Load packages (QZ, FKF)
## Load data (SWdataDM.Rdata)

## This file is meant to run on the head node using the package batchtools

library(batchtools)
makeRegistry('SWpermutations2022', 
             seed = 20190921,
             packages=c('QZ','FKF','optimr'), 
             source = c('code/functions/SimsCodesFort.R', 
                        'code/functions/modelsol.R',
                        'code/functions/reparlogp.R', 
                        'code/priorsetup.R',
                        'code/dataload.R'))


perms = perm(7, 7)  ## all 5040 permutations
perms = split(perms, 1:nrow(perms))

    
batchMap(estim.pred.wrap, perms, more.args = list(
  y = y, prior = prior, estim.obs = 1:200, pred.obs = 201:251,
  nstarts = 5, maxit = 3e4)
)

ids = findJobs()

ids$chunk = chunk(ids$job.id, n.chunks = 96)
submitJobs(ids, resources=list(nodes = 1, ppn=1, walltime='48:00:00'))
