library(QZ)
library(FKF)
load('data/SWdataUpdated.Rdata')
source('code/functions/modelsol.R')
source('code/functions/reparlogp.R')
source('code/functions/SimsCodesFort.R')
source('code/priorsetup.R')
load("cluster_output/SWlong.Rdata")

pvec <- res$par_ests
n <- 30000
long <- generate(n, pvec)
filt <- getFilterOutput(long, pvec)
Zt <- filt$ssmats$Zt
GG <- filt$ssmats$GG
Pbar <- filt$Pt[,,n + 1]
Plast <- filt$Pt[,,n]
rm(filt)
(frobenius_norm <- sqrt(sum((Pbar - Plast)^2)))
limit_obs_var <- Zt %*% Pbar %*% t(Zt) + GG
(entropy_rate <- determinant(2 * pi * exp(1) * limit_obs_var)$modulus / 2)

