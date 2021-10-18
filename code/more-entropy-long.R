library(QZ)
library(FKF)
library(tidyverse)
source('code/functions/modelsol.R')
source('code/functions/reparlogp.R')
source('code/functions/SimsCodesFort.R')
source('code/priorsetup.R')
load("cluster_output/StdPvec.Rdata") # pvec
load("cluster_output/SimEst/parmEst.Rdata")

pvec2 <- unlist(parmEst[nrow(parmEst), -c(1:5)])

n <- 1e6
long <- generate(n, pvec)
log_scale_seq <- function(from = 1, to = 1, length.out = 1) {
  10^seq(log10(from), log10(to), length.out = length.out)
}
ns <- round(log_scale_seq(1000, n, 10), -2)
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
  ii <- seq(ns[i])
  truth[i] <- negll_per_obs(long[,ii], ss_true)
  other[i] <- negll_per_obs(long[,ii], ss_est)
}

tib <- tibble(n = ns, true_params = truth, est_params = other) %>% 
  pivot_longer(-n, values_to = "neg_loglikelihood")
saveRDS(tib, file = "data/entropy-investigation-long.RDS")



# -------------------------------------------------------------------------


# ggplot(tib, aes(n, neg_loglikelihood, color = name)) +
#   geom_line() + 
#   scale_x_log10(labels = scales::label_number_si()) + 
#   geom_point() + 
#   theme_bw() +
#   scale_color_viridis_d(begin=.25,end=.75) +
#   theme(legend.title = element_blank())
# + scale_y_continuous(trans = scales::pseudo_log_trans())
