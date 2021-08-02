library(QZ)
library(FKF)
source('code/functions/modelsol.R')
source('code/functions/reparlogp.R')
source('code/functions/SimsCodesFort.R')
source('code/priorsetup.R')
load("cluster_output/StdPvec.Rdata") # pvec
load("cluster_output/SimEst/parmEst.Rdata")

pvec2 <- unlist(parmEst[nrow(parmEst), -c(1:5)])

n <- 3e6
long <- generate(n, pvec)
ns <- c(1000,2000,3000,5000,10000,15000,20000,30000,100000,150000, 200000, 300000)
truth <- double(length(ns))
other <- double(length(ns))

negll_per_obs <- function(y, simsout) {
  filt = ss.model(y, simsout, output = 'll')
  return(-filt / ncol(y))
}

ss_true <- gensolution(pvec)
ss_est <- gensolution(pvec2)
for (i in seq_along(ns)) {
  print(paste(i," of ", length(ns)))
  truth[i] <- negll_per_obs(long[,seq(ns[i])], ss_true)
  other[i] <- negll_per_obs(long[,seq(ns[i])], ss_est)
}

library(tidyverse)
tib <- tibble(n = ns, true_params = truth, est_params = other) %>% 
  pivot_longer(-n, values_to = "-loglikelihood")
ggplot(tib, aes(n, `-loglikelihood`, color = name)) +
  geom_line() + scale_x_log10(labels = scales::label_number_si()) + 
  geom_point() + theme_bw() +
  scale_color_viridis_d(begin=.25,end=.75) +
  theme(legend.title = element_blank()) +
  scale_y_continuous(trans = scales::pseudo_log_trans())
