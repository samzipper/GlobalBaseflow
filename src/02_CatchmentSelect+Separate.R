## 02_CatchmentSelect+Separate.R
#' This script is intended to read in CSV files for all natural catchments
#' (produced with the script 01_CatchmentGapFill+Stats.R), decide which 
#' catchments for final analysis, and perform baseflow separation on them.

source(file.path("src", "paths+packages.R"))
source(file.path("src", "BaseflowSeparationFunctions.R"))

# Decide catchments to keep ----------------------------------------------------

# read in summary data
df.summary <- 
  read.csv(file.path("data", "catchmentSummary_GapFill.csv"), stringsAsFactors=F) %>% 
  subset(is.finite(yrs.complete))   # only keep sites that had years with no missing data

# # exploratory plots and analysis
# df.summary$yrs.complete.cut <- cut(df.summary$yrs.complete, c(0,10,15,20,500), right=F)
# 
# min(df.summary$yrs.complete)
# max(df.summary$yrs.complete)
# sum(df.summary$yrs.complete>10)
# sum(df.summary$yrs.complete>20)
# 
# ggplot(df.summary, aes(x=yrs.complete)) +
#   geom_histogram(binwidth=1)
# 
# world <- map_data("world")
# ggplot(df.summary, aes(x=lon, y=lat, color=yrs.complete.cut)) +
#   geom_polygon(data = world, aes(x=long, y = lat, group = group), fill=NA, color=col.gray) +
#   geom_point()

# criteria to keep
yrs.complete.min <- 10  # Beck et al. (2013, 2016) WRR used 10 years

# subset summary data
df.summary <- subset(df.summary, yrs.complete >= yrs.complete.min)

cols.hist <- c("catchment", "yrs.complete", "catchsize", "meanelevation", "aridity", "snowfall_fraction", "SoilGrids1km_sand")
p.properties.hist <- 
  df.summary[,cols.hist] %>% 
  melt(id="catchment") %>% 
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~variable, ncol=3, scales="free",
             labeller=as_labeller(c("yrs.complete"="Years Data",
                                    "catchsize"="Catchment Size [km2]",
                                    "meanelevation"="Mean Elevation [m]",
                                    "aridity"="Aridity (=PET/P)",
                                    "snowfall_fraction"="Fraction of Precip as Snow",
                                    "SoilGrids1km_sand"="Soil Sand Content"))) +
  scale_y_continuous(name="# of Catchments") +
  scale_x_continuous(name="Value")
ggsave(file.path("plots", "CatchmentSelect_HistogramProperties.png"),
       p.properties.hist, width=8, height=6, units="in")

df.summary$yrs.complete.cut <- cut(df.summary$yrs.complete, c(10,20,30,500), 
                                   labels=c("10-20", "20-30", ">30"), 
                                   include.lowest=F, right=F)
world <- map_data("world")
p.map.yrs <- 
  ggplot(df.summary, aes(x=lon, y=lat, color=yrs.complete.cut)) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill=NA, color=col.gray) +
  geom_point(shape=21) +
  scale_color_manual(name="Years of Data", values=c("10-20"=col.cat.org, "20-30"=col.cat.grn, ">30"=col.cat.blu)) +
  scale_x_continuous(name=NULL, expand=c(0,0)) +
  scale_y_continuous(name=NULL, expand=c(0,0)) +
  coord_equal() +
  labs(title=paste0(dim(df.summary)[1], " Retained Catchments"), 
       subtitle=paste0(">=", yrs.complete.min, " yrs data & minimal human influence")) +
  theme(legend.position=c(0.01, 0.01),
        legend.justification=c(0, 0))

ggsave(file.path("plots", "CatchmentSelect_Map.png"),
       p.map.yrs, width=8, height=6, units="in")

# Perform baseflow separation on each catchment ---------------------------

# df summary columns to save
df.summary <- subset(df.summary, select=c("catchment", "catchsize"))
df.summary$recessionConstant_Langbein_50  <- NaN
df.summary$recessionConstant_Brutsaert_50 <- NaN
df.summary$recessionConstant_Langbein_90  <- NaN
df.summary$recessionConstant_Brutsaert_90 <- NaN
df.summary$recessionConstant_Langbein_95  <- NaN
df.summary$recessionConstant_Brutsaert_95 <- NaN
df.summary$BFImax_Langbein_50  <- NaN
df.summary$BFImax_Brutsaert_50 <- NaN
df.summary$BFImax_Langbein_90  <- NaN
df.summary$BFImax_Brutsaert_90 <- NaN
df.summary$BFImax_Langbein_95  <- NaN
df.summary$BFImax_Brutsaert_95 <- NaN

n.catchment <- length(df.summary$catchment)
start.flag <- T
for (cat in 1:n.catchment){
    
  # read in data
  cat.name <- df.summary$catchment[cat]
  df.cat <- 
    file.path(dir.Q.derived, "GapFill", paste0(cat.name, ".csv")) %>% 
    read.csv(., stringsAsFactors=F)
  
  # get catchment area
  area_km2 <- df.summary$catchsize[cat]
  
  # unit conversions
  area_mi2 <- area_km2*0.621371*0.621371
  m3s_to_mmd <- 86400*1000/(area_km2*1000*1000)

  # grab flow, normalize to [mm/d]
  df.cat$Q_mm.d <- df.cat$discharge.m3_s*m3s_to_mmd
  
  # make sure no negative Q (very rare)
  df.cat$Q_mm.d[df.cat$Q_mm.d < 0] <- 0
  
  # extract year/month
  df.cat$date <- ymd(df.cat$date)
  df.cat$year <- year(df.cat$date)
  df.cat$month <- month(df.cat$date)
  
  # estimate recession constant and BFImax
  df.summary$recessionConstant_Langbein_5[cat]  <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.05, method="Langbein")
  df.summary$recessionConstant_Brutsaert_5[cat] <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.05, method="Brutsaert")
  df.summary$recessionConstant_Langbein_10[cat]  <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.10, method="Langbein")
  df.summary$recessionConstant_Brutsaert_10[cat] <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.10, method="Brutsaert")
  df.summary$recessionConstant_Langbein_50[cat]  <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.50, method="Langbein")
  df.summary$recessionConstant_Brutsaert_50[cat] <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.50, method="Brutsaert")
  df.summary$recessionConstant_Langbein_90[cat]  <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.90, method="Langbein")
  df.summary$recessionConstant_Brutsaert_90[cat] <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.90, method="Brutsaert")
  df.summary$recessionConstant_Langbein_95[cat]  <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.95, method="Langbein")
  df.summary$recessionConstant_Brutsaert_95[cat] <- baseflow_RecessionConstant(df.cat$Q_mm.d, UB_prc=0.95, method="Brutsaert")
  
  # perform baseflow separations
  df.cat$HYSEP_fixed <- baseflow_HYSEP(Q = df.cat$Q_mm.d, area_mi2 = area_mi2, method="fixed")
  df.cat$HYSEP_slide <- baseflow_HYSEP(Q = df.cat$Q_mm.d, area_mi2 = area_mi2, method="sliding")
  df.cat$HYSEP_local <- baseflow_HYSEP(Q = df.cat$Q_mm.d, area_mi2 = area_mi2, method="local")
  df.cat$UKIH <- baseflow_UKIH(Q = df.cat$Q_mm.d, endrule="B")
  df.cat$BFLOW_1pass <- baseflow_BFLOW(Q = df.cat$Q_mm.d, beta=0.925, passes=1)
  df.cat$BFLOW_3pass <- baseflow_BFLOW(Q = df.cat$Q_mm.d, beta=0.925, passes=3)
  if (is.finite(df.summary$recessionConstant_Langbein_50[cat]) & 
      df.summary$recessionConstant_Langbein_50[cat] < 1){
    
    # calculate BFImax
    df.summary$BFImax_Langbein_50[cat] <- baseflow_BFImax(df.cat$Q_mm.d, 
                                                          k=df.summary$recessionConstant_Langbein_50[cat])
    
    # calculate baseflow
    df.cat$Eckhardt_Langbein_50 <- baseflow_Eckhardt(Q = df.cat$Q_mm.d, 
                                                     BFImax=df.summary$BFImax_Langbein_50[cat], 
                                                     k=df.summary$recessionConstant_Langbein_50[cat])
  } else {
    df.cat$Eckhardt_Langbein_50 <- NaN
  }
  
  if (is.finite(df.summary$recessionConstant_Brutsaert_50[cat]) & 
      df.summary$recessionConstant_Brutsaert_50[cat] < 1){
    
    # calculate BFImax
    df.summary$BFImax_Brutsaert_50[cat] <- baseflow_BFImax(df.cat$Q_mm.d, 
                                                           k=df.summary$recessionConstant_Brutsaert_50[cat])
    
    # calculate baseflow
    df.cat$Eckhardt_Brutsaert_50 <- baseflow_Eckhardt(Q = df.cat$Q_mm.d, 
                                                      BFImax=df.summary$BFImax_Brutsaert_50[cat], 
                                                      k=df.summary$recessionConstant_Brutsaert_50[cat])
  } else {
    df.cat$Eckhardt_Brutsaert_50 <- NaN
  }
  
  if (is.finite(df.summary$recessionConstant_Langbein_90[cat]) & 
      df.summary$recessionConstant_Langbein_90[cat] < 1){
    
    # calculate BFImax
    df.summary$BFImax_Langbein_90[cat] <- baseflow_BFImax(df.cat$Q_mm.d, 
                                                          k=df.summary$recessionConstant_Langbein_90[cat])
    
    # calculate baseflow
    df.cat$Eckhardt_Langbein_90 <- baseflow_Eckhardt(Q = df.cat$Q_mm.d, 
                                                     BFImax=df.summary$BFImax_Langbein_90[cat], 
                                                     k=df.summary$recessionConstant_Langbein_90[cat])
  } else {
    df.cat$Eckhardt_Langbein_90 <- NaN
  }
  
  if (is.finite(df.summary$recessionConstant_Brutsaert_90[cat]) & 
      df.summary$recessionConstant_Brutsaert_90[cat] < 1){
    
    # calculate BFImax
    df.summary$BFImax_Brutsaert_90[cat] <- baseflow_BFImax(df.cat$Q_mm.d, 
                                                           k=df.summary$recessionConstant_Brutsaert_90[cat])
    
    # calculate baseflow
    df.cat$Eckhardt_Brutsaert_90 <- baseflow_Eckhardt(Q = df.cat$Q_mm.d, 
                                                      BFImax=df.summary$BFImax_Brutsaert_90[cat], 
                                                      k=df.summary$recessionConstant_Brutsaert_90[cat])
  } else {
    df.cat$Eckhardt_Brutsaert_90 <- NaN
  }
  
  if (is.finite(df.summary$recessionConstant_Langbein_95[cat]) & 
      df.summary$recessionConstant_Langbein_95[cat] < 1){
    
    # calculate BFImax
    df.summary$BFImax_Langbein_95[cat] <- baseflow_BFImax(df.cat$Q_mm.d, 
                                                          k=df.summary$recessionConstant_Langbein_95[cat])
    
    # calculate baseflow
    df.cat$Eckhardt_Langbein_95 <- baseflow_Eckhardt(Q = df.cat$Q_mm.d, 
                                                     BFImax=df.summary$BFImax_Langbein_95[cat], 
                                                     k=df.summary$recessionConstant_Langbein_95[cat])
  } else {
    df.cat$Eckhardt_Langbein_95 <- NaN
  }
  
  if (is.finite(df.summary$recessionConstant_Brutsaert_95[cat]) & 
      df.summary$recessionConstant_Brutsaert_95[cat] < 1){
    
    # calculate BFImax
    df.summary$BFImax_Brutsaert_95[cat] <- baseflow_BFImax(df.cat$Q_mm.d, 
                                                           k=df.summary$recessionConstant_Brutsaert_95[cat])
    
    # calculate baseflow
    df.cat$Eckhardt_Brutsaert_95 <- baseflow_Eckhardt(Q = df.cat$Q_mm.d, 
                                                      BFImax=df.summary$BFImax_Brutsaert_95[cat], 
                                                      k=df.summary$recessionConstant_Brutsaert_95[cat])
  } else {
    df.cat$Eckhardt_Brutsaert_95 <- NaN
  }
  
  # save gap-filled data as output
  cols.save <- c("date", "Q_mm.d", "HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", 
                 "Eckhardt_Langbein_50", "Eckhardt_Brutsaert_50", 
                 "Eckhardt_Langbein_90", "Eckhardt_Brutsaert_90", 
                 "Eckhardt_Langbein_95", "Eckhardt_Brutsaert_95")
  path.save <- file.path(dir.Q.derived, "Baseflow", paste0(cat.name, "_Daily.csv"))
  write.csv(df.cat[,cols.save], path.save, row.names=F)
  
  ## save some annual statistics to df.summary
  df.ann <-
    df.cat %>% 
    group_by(year) %>% 
    summarize(HYSEP_fixed = sum(HYSEP_fixed),
              HYSEP_slide = sum(HYSEP_slide),
              HYSEP_local = sum(HYSEP_local),
              UKIH = sum(UKIH),
              BFLOW_1pass = sum(BFLOW_1pass),
              BFLOW_3pass = sum(BFLOW_3pass),
              Eckhardt_Langbein_50 = sum(Eckhardt_Langbein_50),
              Eckhardt_Brutsaert_50 = sum(Eckhardt_Brutsaert_50),
              Eckhardt_Langbein_90 = sum(Eckhardt_Langbein_90),
              Eckhardt_Brutsaert_90 = sum(Eckhardt_Brutsaert_90),
              Eckhardt_Langbein_95 = sum(Eckhardt_Langbein_95),
              Eckhardt_Brutsaert_95 = sum(Eckhardt_Brutsaert_95)) %>% 
    melt(id=c("year"), variable.name="method") %>% 
    group_by(method) %>% 
    summarize(bf_mm.y = mean(value)) %>% 
    transform(recessionConstant = c(NaN, NaN, NaN, NaN, NaN, NaN,
                                    df.summary$recessionConstant_Langbein_50[cat],
                                    df.summary$recessionConstant_Brutsaert_50[cat],
                                    df.summary$recessionConstant_Langbein_90[cat],
                                    df.summary$recessionConstant_Brutsaert_90[cat],
                                    df.summary$recessionConstant_Langbein_95[cat],
                                    df.summary$recessionConstant_Brutsaert_95[cat]),
              BFImax = c(NaN, NaN, NaN, NaN, NaN, NaN,
                                    df.summary$BFImax_Langbein_50[cat],
                                    df.summary$BFImax_Brutsaert_50[cat],
                                    df.summary$BFImax_Langbein_90[cat],
                                    df.summary$BFImax_Brutsaert_90[cat],
                                    df.summary$BFImax_Langbein_95[cat],
                                    df.summary$BFImax_Brutsaert_95[cat]))
  
  df.ann$Q_mm.y <- 
    df.cat %>% 
    group_by(year) %>% 
    summarize(Q_mm.y = sum(Q_mm.d)) %>% 
    summarize(mean(Q_mm.y)) %>% 
    as.numeric()
  df.ann$BFI <- df.ann$bf_mm.y/df.ann$Q_mm.y
  df.ann$catchment <- cat.name
  
  # put together new summmary data frame
  if (start.flag){
    df.summary.out <- df.ann
    start.flag <- F
  } else {
    df.summary.out <- rbind(df.summary.out, df.ann)
  }
  
  # status update
  if ((cat %% 5)==0){
    print(paste0(cat, " of ", n.catchment, " complete"))
  }
  
}

# Save summary -------------------------------------------------------------

# round
df.summary.out[,c("bf_mm.y", "recessionConstant", "BFImax", "Q_mm.y", "BFI")] <-
  round(df.summary.out[,c("bf_mm.y", "recessionConstant", "BFImax", "Q_mm.y", "BFI")], 3)

# save to GSAS and git repository
write.csv(df.summary.out, file.path(dir.Q.derived, "Baseflow", "catchmentSummary_Baseflow.csv"), row.names=F)
write.csv(df.summary.out, file.path("data", "catchmentSummary_Baseflow.csv"), row.names=F)
