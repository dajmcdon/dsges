library(batchtools)
reg = loadRegistry('SWSimulateEstimate-files')
source('DSGEs/code/processClustOutput/SimEstClustfuncs.R')
load('DSGEs/StdPvec.Rdata')
load('DSGEs/data/SWDataDM.Rdata')

SWPerm = c(4,7,5,3,1,2,6)
lnames = rownames(y)[SWPerm]

SummaryStats = reduceResultsExperiments(reg, fun=getClustOutput, parms=pvec)
save(SummaryStats, file='DSGEs/output/SimEst/summaryStats.Rdata')

pred.errs.sc = reduceResultsExperiments(reg, fun=getClustPred, lnames=lnames)
save(pred.errs.sc, file='DSGEs/output/SimEst/predErrs.Rdata')

train.errs.sc = reduceResultsExperiments(reg, fun=getClustTrain, lnames=lnames)
save(train.errs.sc, file='DSGEs/output/SimEst/trainErrs.Rdata')

base.errs.sc = reduceResultsExperiments(reg, fun=getClustBase, lnames=lnames)
save(base.errs.sc, file='DSGEs/output/SimEst/baselineErrs.Rdata')

parmEst = reduceResultsExperiments(reg, fun=getClustParms)
save(parmEst, file = 'DSGEs/output/SimEst/parmEst.Rdata')
