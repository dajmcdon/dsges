library(tidyverse)
library(QZ)
library(FKF)
source('code/functions/modelsol.R')
source('code/functions/reparlogp.R')
source('code/functions/SimsCodesFort.R')
source('code/priorsetup.R')

simest <- read_rds("update2022/simulate_estimate.rds")
pvec <- read_rds("update2022/swlong-run-pars.rds")


pvec2 <- simest %>% 
  filter(nestim == max(nestim)) %>% 
  select(`sigma[a]`:alpha) %>% 
  summarise(across(everything(), ~ mean(.x))) %>%
  pivot_longer(everything()) %>% 
  pull(value)

rm(simest)

n <- 3e5
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
saveRDS(tib, file = here::here("update2022", "entropy-investigation-long.rds"))
