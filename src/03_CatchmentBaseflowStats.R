## 03_CatchmentSummarizeBaseflow.R
#' This script is intended to read in CSV files of baseflow for natural catchments
#' (produced with the script 02_CatchmentSelect+Separate.R), calculate some summary
#' statistics, and make some graphs.

source(file.path("src", "paths+packages.R"))

# Read in data ------------------------------------------------------------

# summary data - additional catchment characteristics are in the file catchmentsSummary_GapFill.csv
df.summary <- 
  read.csv(file.path("data", "catchmentSummary_Baseflow.csv"), stringsAsFactors=F) %>% 
  subset(complete.cases(.)) %>% 
  dplyr::select(catchment, recessionConstant, BFImax) %>% 
  transform(noFlowDays_prc = NaN,  # percent of record with 0 flow
            Q5  = NaN,          # 95th percentile streamflow
            Q10 = NaN,          # 90th percentile streamflow
            Q50 = NaN,          # 50th percentile streamflow
            Q90 = NaN,          # 10th percentile streamflow
            Q95 = NaN,          # 5th percentile streamflow
            FDC_slope = NaN,    # slope of 33rd to 66th percentile of flow-duration curve
            dQ_dt.slope = NaN)  # slope of lower envelope of log(-dQ/dt) vs log(Q) plot

# Loop through catchments -------------------------------------------------

catchment <- unique(df.summary$catchment)
n.catchment <- length(catchment)
start.flag <- T
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
  
  ## fill in Q summary stats
  df.summary$noFlowDays_prc[cat] <- sum(df.cat$Q_mm.d==0)/length(df.cat$Q_mm.d)
  df.summary$Q5[cat]  <- as.numeric(quantile(df.cat$Q_mm.d, 0.95))
  df.summary$Q10[cat] <- as.numeric(quantile(df.cat$Q_mm.d, 0.90))
  df.summary$Q50[cat] <- as.numeric(quantile(df.cat$Q_mm.d, 0.50))
  df.summary$Q90[cat] <- as.numeric(quantile(df.cat$Q_mm.d, 0.10))
  df.summary$Q95[cat] <- as.numeric(quantile(df.cat$Q_mm.d, 0.05))
  
  # slope of FDC - Sawicz et al., 2011, HESS, Eq. 3
  df.summary$FDC_slope[cat] <- 
    (log(as.numeric(quantile(df.cat$Q_mm.d, 0.33))) - 
       log(as.numeric(quantile(df.cat$Q_mm.d, 0.66)))) / 
    (0.66 - 0.33)
  
  # slope of -dQ/dt vs Q relationship
  # calculate lagged difference (dQ/dt) based on before/after point
  dQ_dt <- c(NaN, diff(df.cat$Q_mm.d, lag=2)/2, NaN)
  dQ_dt_left <- c(NaN, diff(df.cat$Q_mm.d))
  
  # screen data for which dQ_dt to calculate recession, based on rules in Brutsaert (2008) WRR Section 3.2
  which_negative <- which(dQ_dt < 0 & dQ_dt_left < 0 & df.cat$Q_mm.d > 0)
  which_positive <- which(dQ_dt >= 0)
  which_positive_with_buffer <- unique(c(which_positive-2, which_positive-1, which_positive,
                                         which_positive+1, which_positive+2, which_positive+3))  # 2 days before and 3 days after a positive or 0 value
  which_positive_with_buffer <- which_positive_with_buffer[which_positive_with_buffer > 0]  # get rid of negative indices; possible because of 2 days before
  which_keep <- which_negative[!(which_negative %in% which_positive_with_buffer)]
  
  # fit quantile regression
  fit.qr <- rq(log10(-dQ_dt[which_keep]) ~ log10(df.cat$Q_mm.d[which_keep]), tau=0.05)
  
  # get slope
  df.summary$dQ_dt.slope[cat] <- coef(fit.qr)[2]
  
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
  
  # calculate annual totals
  df.cat.ann <- 
    df.cat %>% 
    group_by(year) %>% 
    summarize(Q_mm.d = sum(Q_mm.d),
              HYSEP_fixed = sum(HYSEP_fixed),
              HYSEP_slide = sum(HYSEP_slide),
              HYSEP_local = sum(HYSEP_local),
              UKIH = sum(UKIH),
              BFLOW_1pass = sum(BFLOW_1pass),
              BFLOW_3pass = sum(BFLOW_3pass),
              Eckhardt = sum(Eckhardt))
  
  # find maxes and mins
  df.mo.max <-
    df.cat.mo %>% 
    melt(id=c("month")) %>% 
    group_by(variable) %>% 
    filter(value==max(value)) %>% 
    summarize(month = median(month),
              value = median(value)) %>% 
    set_colnames(c("flux", "maxFlow_mo", "maxFlow_mo_mm.d")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.DOY.max <-
    df.cat.DOY %>% 
    melt(id=c("DOY")) %>% 
    group_by(variable) %>% 
    filter(value==max(value)) %>% 
    summarize(DOY = median(DOY),
              value = median(value)) %>% 
    set_colnames(c("flux", "maxFlow_DOY", "maxFlow_DOY_mm.d")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.mo.min <-
    df.cat.mo %>% 
    melt(id=c("month")) %>% 
    group_by(variable) %>% 
    filter(value==min(value)) %>% 
    summarize(month = median(month),
              value = median(value)) %>% 
    set_colnames(c("flux", "minFlow_mo", "minFlow_mo_mm.d")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.DOY.min <-
    df.cat.DOY %>% 
    melt(id=c("DOY")) %>% 
    group_by(variable) %>% 
    filter(value==min(value)) %>% 
    summarize(DOY = median(DOY),
              value = median(value)) %>% 
    set_colnames(c("flux", "minFlow_DOY", "minFlow_DOY_mm.d")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.ann.mean <-
    df.cat.ann %>% 
    melt(id=c("year"), variable.name="flux") %>% 
    group_by(flux) %>% 
    summarize(annFlow_mean_mm.yr = mean(value),
              annFlow_sd_mm.yr = sd(value)) %>% 
    transform(flux = as.character(flux)) %>% 
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
    transform(flux = as.character(flux)) %>% 
    ungroup() %>% 
    dplyr::select(-mod)
  
  df.ann.mean.trend <-
    df.cat.ann %>% 
    melt(id=c("year"), variable.name="flux") %>% 
    group_by(flux) %>% 
    do(mod = lm(value ~ year, data = .)) %>%
    mutate(annFlow_trend_mm.yr = summary(mod)$coeff[2],
           annFlow_trend_p = lmp(mod)) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup() %>% 
    dplyr::select(-mod)
  
  ## now: repeat analysis with BFI
  df.cat.BFI <-
    df.cat[,c("date", "year", "month", "DOY")] %>% 
    cbind(df.cat[,c("HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt")]/df.cat$Q_mm.d)
  df.cat.BFI[df.cat$Q_mm.d==0, c("HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt")] <- 0
  
  df.cat.BFI.yr.mo <-
    df.cat.BFI %>% 
    group_by(year, month) %>% 
    summarize(HYSEP_fixed = mean(HYSEP_fixed),
              HYSEP_slide = mean(HYSEP_slide),
              HYSEP_local = mean(HYSEP_local),
              UKIH = mean(UKIH),
              BFLOW_1pass = mean(BFLOW_1pass),
              BFLOW_3pass = mean(BFLOW_3pass),
              Eckhardt = mean(Eckhardt))
  
  df.cat.BFI.DOY <-
    df.cat.BFI %>% 
    group_by(DOY) %>% 
    summarize(HYSEP_fixed = mean(HYSEP_fixed),
              HYSEP_slide = mean(HYSEP_slide),
              HYSEP_local = mean(HYSEP_local),
              UKIH = mean(UKIH),
              BFLOW_1pass = mean(BFLOW_1pass),
              BFLOW_3pass = mean(BFLOW_3pass),
              Eckhardt = mean(Eckhardt))
  
  df.cat.BFI.mo <-
    df.cat.BFI.yr.mo %>% 
    group_by(month) %>% 
    summarize(HYSEP_fixed = mean(HYSEP_fixed),
              HYSEP_slide = mean(HYSEP_slide),
              HYSEP_local = mean(HYSEP_local),
              UKIH = mean(UKIH),
              BFLOW_1pass = mean(BFLOW_1pass),
              BFLOW_3pass = mean(BFLOW_3pass),
              Eckhardt = mean(Eckhardt))
  
  df.cat.BFI.ann <- 
    df.cat.ann[,"year"] %>% 
    cbind(df.cat.ann[,c("HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt")]/df.cat.ann$Q_mm.d)
  
  df.BFI.mo.max <-
    df.cat.BFI.mo %>% 
    melt(id=c("month")) %>% 
    group_by(variable) %>% 
    filter(value==max(value)) %>% 
    summarize(month = median(month),
              value = median(value)) %>% 
    set_colnames(c("flux", "maxBFI_mo", "maxBFI_mo_BFI")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  # calculate metrics
  df.BFI.DOY.max <-
    df.cat.BFI.DOY %>% 
    melt(id=c("DOY")) %>% 
    group_by(variable) %>% 
    filter(value==max(value)) %>% 
    summarize(DOY = median(DOY),
              value = median(value)) %>% 
    set_colnames(c("flux", "maxBFI_DOY", "maxBFI_DOY_BFI")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.BFI.mo.min <-
    df.cat.BFI.mo %>% 
    melt(id=c("month")) %>% 
    group_by(variable) %>% 
    filter(value==min(value)) %>% 
    summarize(month = median(month),
              value = median(value)) %>% 
    set_colnames(c("flux", "minBFI_mo", "minBFI_mo_BFI")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.BFI.DOY.min <-
    df.cat.BFI.DOY %>% 
    melt(id=c("DOY")) %>% 
    group_by(variable) %>% 
    filter(value==min(value)) %>% 
    summarize(DOY = median(DOY),
              value = median(value)) %>% 
    set_colnames(c("flux", "minBFI_DOY", "minBFI_DOY_BFI")) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.ann.BFI.mean <-
    df.cat.BFI.ann %>% 
    melt(id=c("year"), variable.name="flux") %>% 
    group_by(flux) %>% 
    summarize(annBFI_mean_BFI = mean(value),
              annBFI_sd_BFI = sd(value)) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup()
  
  df.BFI.DOY.max.trend <-
    df.cat.BFI %>% 
    melt(id=c("date", "DOY", "year", "month"), variable.name="flux") %>% 
    group_by(flux, year) %>% 
    filter(value==max(value)) %>% 
    group_by(flux) %>% 
    do(mod = lm(DOY ~ year, data = .)) %>%
    mutate(maxBFI_DOY_trend = summary(mod)$coeff[2],
           maxBFI_DOY_trend_p = lmp(mod)) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup() %>% 
    dplyr::select(-mod)
  
  df.ann.BFI.mean.trend <-
    df.cat.BFI.ann %>% 
    melt(id=c("year"), variable.name="flux") %>% 
    group_by(flux) %>% 
    do(mod = lm(value ~ year, data = .)) %>%
    mutate(annBFI_trend = summary(mod)$coeff[2],
           annBFI_trend_p = lmp(mod)) %>% 
    transform(flux = as.character(flux)) %>% 
    ungroup() %>% 
    dplyr::select(-mod)
  
  ## combine into single data frame
  df <-
    df.mo.max %>% 
    left_join(df.BFI.mo.max, by="flux") %>% 
    left_join(df.DOY.max, by="flux") %>% 
    left_join(df.DOY.max.trend, by="flux") %>% 
    left_join(df.BFI.DOY.max, by="flux") %>% 
    left_join(df.BFI.DOY.max.trend, by="flux") %>% 
    left_join(df.mo.min, by="flux") %>% 
    left_join(df.BFI.mo.min, by="flux") %>% 
    left_join(df.DOY.min, by="flux") %>% 
    left_join(df.BFI.DOY.min, by="flux") %>% 
    left_join(df.ann.mean, by="flux") %>% 
    left_join(df.ann.mean.trend, by="flux") %>% 
    left_join(df.ann.BFI.mean, by="flux") %>% 
    left_join(df.ann.BFI.mean.trend, by="flux") %>% 
    transform(catchment = cat.name)
  
  # round columns
  cols.round <- c("maxFlow_mo_mm.d", "maxBFI_mo_BFI", "maxFlow_DOY_mm.d", "maxFlow_DOY_trend", "maxFlow_DOY_trend_p",
                  "maxBFI_DOY_BFI", "maxBFI_DOY_trend", "maxBFI_DOY_trend_p", "minFlow_mo_mm.d", "minBFI_mo_BFI",
                  "minFlow_DOY_mm.d", "minBFI_DOY_BFI", "annFlow_mean_mm.yr", "annFlow_sd_mm.yr", "annFlow_trend_mm.yr",
                  "annFlow_trend_p", "annBFI_mean_BFI", "annBFI_sd_BFI", "annBFI_trend", "annBFI_trend_p")
  df[,cols.round] <- signif(df[,cols.round], 3)
  
  # save baseflow stats
  path.save <- file.path(dir.Q.derived, "Baseflow", paste0(cat.name, "_BaseflowMethodComparison.csv"))
  write.csv(df, path.save, row.names=F)
  
  if (start.flag){
    df.out <- df
    start.flag <- F
  } else {
    df.out <- rbind(df.out, df)
  }
  
  print(paste0(cat, " of ", n.catchment, " complete, ", Sys.time()))
  
}

# Save output -------------------------------------------------------------

# save to GSAS and git repository
write.csv(df.summary, file.path(dir.Q.derived, "Baseflow", "catchmentSummary_BaseflowStats.csv"), row.names=F)
write.csv(df.summary, file.path("data", "catchmentSummary_BaseflowStats.csv"), row.names=F)

write.csv(df.out, file.path(dir.Q.derived, "Baseflow", "catchmentSummary_BaseflowStats_AllMethods.csv"), row.names=F)
write.csv(df.out, file.path("data", "catchmentSummary_BaseflowStats_AllMethods.csv"), row.names=F)
