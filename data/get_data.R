library(fredr)
library(tidyverse)
library(lubridate)


# Download raw series -----------------------------------------------------

# must first obtain a FRED API Key: https://research.stlouisfed.org/docs/api/api_key.html
fredr_set_key("99e75d62216e44f070cf6f934f5f937f")

start_date = mdy('1/1/1956')
end_date = mdy('11/1/2018')
index_date = mdy('7/1/2012') # to match index period for PRS85006023

series_ids = c("GDPC1", "GDPDEF", "PCEC", "FPI", "CE16OV", 
               "FEDFUNDS", "CNP16OV", "PRS85006023", "COMPNFB")

raw = lapply(series_ids, function(x) fredr(
  x,
  observation_start = start_date,
  observation_end = end_date,
  units = 'lin',
  frequency = 'q',
  aggregation_method = 'avg'
))

lapply(raw, function(x) range(x$date)) # check that all series span the date range

raw_df = bind_cols(raw) %>% select_at(vars("date", contains("value"))) 
names(raw_df)[-1] = series_ids


# Transformations ---------------------------------------------------------

# Index population variables to index_date
idx_period = filter(raw_df, date==index_date)
indexed = raw_df %>% mutate(CE16OV_idx = CE16OV/idx_period$CE16OV,
                            CNP16OV_idx = CNP16OV/idx_period$CNP16OV)

# See the SW Dynare documentation, note that some series no longer exist
trans = indexed %>% transmute(
  consumption = log((PCEC/GDPDEF)/CNP16OV_idx)*100,
  investment = log((FPI/GDPDEF)/CNP16OV_idx)*100,
  output = log(GDPC1/CNP16OV_idx)*100,
  hours = log((PRS85006023*CE16OV_idx/100)/CNP16OV_idx)*100, # mistake in docs, use idx instead of raw CE16OV
  inflation = GDPDEF,
  wage = log(COMPNFB/GDPDEF)*100,
  interest_rate = FEDFUNDS/4
)

observed = trans %>% transmute(
  labobs = hours - mean(hours),
  robs = interest_rate,
  pinfobs = log(inflation/lag(inflation))*100,
  dy = output - lag(output),
  dc = consumption - lag(consumption),
  dinve = investment - lag(investment),
  dw = wage - lag(wage)
) %>% filter(!is.na(dc))

y = t(as.matrix(observed))
save(y, file='data/SWdataUpdated.Rdata')
