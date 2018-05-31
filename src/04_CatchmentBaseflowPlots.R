## 03_CatchmentSummarizeBaseflow.R
#' This script is intended to read in CSV files of baseflow for natural catchments
#' (produced with the script 02_CatchmentSelect+Separate.R), calculate some summary
#' statistics, and make some graphs.

source(file.path("src", "paths+packages.R"))

# Read in data ------------------------------------------------------------

# summary data - additional catchment characteristics are in the file catchmentsSummary_GapFill.csv
df.summary <- 
  read.csv(file.path("data", "catchmentSummary_BaseflowStats.csv"), stringsAsFactors=F) %>% 
  left_join(., read.csv(file.path("data", "catchmentSummary_GapFill.csv"), stringsAsFactors=F),
            by="catchment")

df.all <- 
  read.csv(file.path("data", "catchmentSummary_BaseflowStats_AllMethods.csv"), stringsAsFactors=F)

# Filter data and eliminate bad stations ----------------------------------

catchsize.thres <- 10
df.summary <- subset(df.summary, catchsize > catchsize.thres)

df.all <- subset(df.all, catchment %in% df.summary$catchment)

# Derived metrics ---------------------------------------------------------

df.summary$recessionConstantDays_Langbein = -1/log(df.summary$recessionConstant_Langbein)    # eq. 6 from Barlow digital filters tutorial
df.summary$recessionConstantDays_Brutsaert = -1/log(df.summary$recessionConstant_Brutsaert)  # eq. 6 from Barlow digital filters tutorial

# Plots -------------------------------------------------------------------

# distribution of recessionConstant
p.hist.recessonConstantDays <-
  df.summary %>% 
  dplyr::select(catchment, recessionConstant_Langbein, recessionConstant_Brutsaert) %>% 
  melt(id="catchment") %>% 
  ggplot(aes(x=value)) +
  geom_histogram(binwidth=0.05) +
  facet_wrap(~variable)

# distribution of recessionConstantDays
p.hist.recessonConstantDays <-
  df.summary %>% 
  dplyr::select(catchment, recessionConstantDays_Langbein, recessionConstantDays_Brutsaert) %>% 
  melt(id="catchment") %>% 
  ggplot(aes(x=value)) +
  geom_histogram(binwidth=1) +
  facet_wrap(~variable)
sum(df.summary$recessionConstantDays_Langbein >= 30 & df.summary$recessionConstantDays_Langbein <= 60)

# distribution of BFImax
p.hist.BFImax <-
  df.summary %>% 
  ggplot(aes(x=BFImax)) +
  geom_histogram(binwidth=0.05)

# distribution of Q quantiles
p.hist.Q <-
  df.summary %>% 
  dplyr::select(catchment, Q10, Q50, Q90) %>% 
  melt(id="catchment") %>% 
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~variable, scales="free")

# distribution of FDC_slope
p.hist.FDC_slope <-
  df.summary %>% 
  ggplot(aes(x=FDC_slope)) +
  geom_histogram()

# distribution of dQ_dt.slope
p.hist.FDC_slope <-
  df.summary %>% 
  ggplot(aes(x=dQ_dt.slope)) +
  geom_histogram()
