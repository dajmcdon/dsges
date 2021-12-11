library(tidyverse)
start_year <- 1990
stop_year <- 2021
n_years <- stop_year - start_year + 1
years <- start_year:stop_year
nber <- "https://www2.nber.org/RePEc/nbr/nberwo/"
n_macro <- integer(n_years)
n_dsges <- integer(n_years)
for (i in seq_along(years)) {
  yy <- years[i]
  print(yy)
  u <- url(paste0(nber, "nberwo", yy, ".rdf"))
  f <- readLines(u)
  close(u)
  abst <- f[str_detect(f, "Abstract:")]
  jel <- str_split(f[str_detect(f, "JEL")], ":")
  jel <- sapply(jel, function(x) x[2])
  macro <- str_detect(jel, "E6")
  n_macro[i] <- sum(macro) 
  dsges <- c(grep("dynamic stochastic", abst, ignore.case = TRUE),
             grep("stochastic general", abst, ignore.case = TRUE),
             grep("dynamic general", abst, ignore.case = TRUE),
             grep("general equlibrium", abst, ignore.case = TRUE),
             grep("dynamic equilibrium", abst, ignore.case = TRUE),
             grep("dsge", abst, ignore.case = TRUE))
  n_dsges[i] <- n_distinct(dsges)
}

tib <- tibble(year = years, macro = n_macro, dsges = n_dsges)
ggplot(tib, aes(year, dsges)) + geom_line() + geom_point()
