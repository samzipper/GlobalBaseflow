## 03_CatchmentSummarizeBaseflow.R
#' This script is intended to read in CSV files of baseflow for natural catchments
#' (produced with the script 02_CatchmentSelect+Separate.R), calculate some summary
#' statistics, and make some graphs.

source(file.path("src", "paths+packages.R"))

# Read in data ------------------------------------------------------------

# summary data - additional catchment characteristics are in the file catchmentsSummary_GapFill.csv
df.summary <- 
  read.csv(file.path("data", "catchmentSummary_Baseflow.csv"), stringsAsFactors=F)

# Loop through catchments -------------------------------------------------

catchment <- unique(df.summary$catchment)
n.catchment <- length(catchment)
for (cat in 1:n.catchment){
  
  # read in data
  cat.name <- catchment[cat]
  df.cat <- 
    file.path(dir.Q.derived, "Baseflow", paste0(cat.name, "_Daily.csv")) %>% 
    read.csv(., stringsAsFactors=F) %>% 
    transform(date  = ymd(date),
              year  = year(date),
              month = month(date),
              DOY   = yday(date))
  
  # calculate and save monthly means
  df.cat.yr.mo <-
    df.cat %>% 
    group_by(year, month) %>% 
    summarize(Q_mm.d = mean(Q_mm.d),
              HYSEP_fixed = mean(HYSEP_fixed),
              HYSEP_slide = mean(HYSEP_slide),
              HYSEP_local = mean(HYSEP_local),
              UKIH = mean(UKIH),
              BFLOW_1pass = mean(BFLOW_1pass),
              BFLOW_3pass = mean(BFLOW_3pass),
              Eckhardt = mean(Eckhardt))
  path.save <- file.path(dir.Q.derived, "Baseflow", paste0(cat.name, "_Monthly.csv"))
  write.csv(df.cat.yr.mo, path.save, row.names=F)
  
  # calculate long-term average for each day
  df.cat.DOY <-
    df.cat %>% 
    group_by(DOY) %>% 
    summarize(Q_mm.d = mean(Q_mm.d),
              HYSEP_fixed = mean(HYSEP_fixed),
              HYSEP_slide = mean(HYSEP_slide),
              HYSEP_local = mean(HYSEP_local),
              UKIH = mean(UKIH),
              BFLOW_1pass = mean(BFLOW_1pass),
              BFLOW_3pass = mean(BFLOW_3pass),
              Eckhardt = mean(Eckhardt))
  
  # calculate long-term average for each month
  df.cat.mo <-
    df.cat.yr.mo %>% 
    group_by(month) %>% 
    summarize(Q_mm.d = mean(Q_mm.d),
              HYSEP_fixed = mean(HYSEP_fixed),
              HYSEP_slide = mean(HYSEP_slide),
              HYSEP_local = mean(HYSEP_local),
              UKIH = mean(UKIH),
              BFLOW_1pass = mean(BFLOW_1pass),
              BFLOW_3pass = mean(BFLOW_3pass),
              Eckhardt = mean(Eckhardt))
  
  # find maxes and mins
  df.mo.max <-
    df.cat.mo %>% 
    melt(id=c("month")) %>% 
    group_by(variable) %>% 
    filter(value==max(value)) %>% 
    set_colnames(c("maxFlow_mo", "flux", "maxFlow_mo_mm.d")) %>% 
    ungroup()
  
  df.DOY.max <-
    df.cat.DOY %>% 
    melt(id=c("DOY")) %>% 
    group_by(variable) %>% 
    filter(value==max(value)) %>% 
    set_colnames(c("maxFlow_DOY", "flux", "maxFlow_DOY_mm.d")) %>% 
    ungroup()
  
  df.mo.min <-
    df.cat.mo %>% 
    melt(id=c("month")) %>% 
    group_by(variable) %>% 
    filter(value==min(value)) %>% 
    set_colnames(c("minFlow_mo", "flux", "minFlow_mo_mm.d")) %>% 
    ungroup()
  
  df.DOY.min <-
    df.cat.DOY %>% 
    melt(id=c("DOY")) %>% 
    group_by(variable) %>% 
    filter(value==min(value)) %>% 
    set_colnames(c("minFlow_DOY", "flux", "minFlow_DOY_mm.d")) %>% 
    ungroup()
  
  # find slopes and trends
  df.DOY.max.trend <-
    df.cat %>% 
    melt(id=c("date", "DOY", "year", "month"), variable.name="flux") %>% 
    group_by(flux, year) %>% 
    filter(value==max(value)) %>% 
    group_by(flux) %>% 
    do(mod = lm(DOY ~ year, data = .)) %>%
    mutate(maxFlow_DOY_trend = summary(mod)$coeff[2],
           maxFlow_DOY_trend_p = lmp(mod)) %>% 
    ungroup() %>% 
    dplyr::select(-mod)
  
  # combine into single data frame
  df <-
    df.mo.max %>% 
    left_join(df.DOY.max, by="flux") %>% 
    left_join(df.mo.min, by="flux") %>% 
    left_join(df.DOY.min, by="flux") %>% 
    transform(catchment = cat.name,
              units  = "mm_d")
  
  # now: repeat analysis with BFI
  df.cat.BFI <-
    df.cat[,c("date", "year", "month", "DOY")] %>% 
    cbind(df.cat[,c("HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt")]/df.cat$Q_mm.d, by="date")
  
}
