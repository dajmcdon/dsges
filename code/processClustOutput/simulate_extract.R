library(batchtools)
library(tidyverse)

loadRegistry('SWSimulateEstimate2022')
source("code/dataload.R")

SWPerm = c(4,7,5,3,1,2,6)
lnames = rownames(y)[SWPerm]
parmGreek = c(
  'sigma[a]','sigma[b]','sigma[g]','sigma[I]', 'sigma[r]', 
  'sigma[p]','sigma[w]',
  'rho[a]', 'rho[b]', 'rho[g]', 
  'rho[I]', 'rho[r]', 'rho[p]','rho[w]', 'mu[p]',
  'mu[w]', 'phi1', 'sigma[c]', 'h', 'xi[w]', 'sigma[l]', 
  'xi[p]', 'iota[w]',
  'iota[p]', 'Psi', 'Phi', 'r[pi]', 'rho', 'r[y]', 
  'r[Delta][y]', 'bar(pi)', 
  '100(beta^{-1} -1)', 'bar(l)', 'bar(gamma)', 'rho[ga]', 'alpha'
)

allnames = c(
  "lltrain",
  paste("mste", lnames, sep = "_"),
  paste("mspe", lnames, sep = "_"),
  parmGreek,
  "llpred",
  paste(lnames, "sc", sep = "_"),
  paste("base", lnames, sep = "_")
)

pars = getJobPars() %>% 
  unwrap() %>%
  select(job.id, nestim)

res = reduceResultsList()
res = do.call(rbind, res)
colnames(res) = allNames
res = as_tibble(bind_cols(pp, res))

write_rds(res, file = "update2022/simulate_estimate.rds")