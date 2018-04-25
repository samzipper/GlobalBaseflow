## 01_CatchmentGapFill+Stats.R
#' This script is intended to read in CSV files for all natural catchments
#' (produced with the script LiberateData.m), calculate some simple statistics
#' which will be added to the table, and trim discharge data to only complete years.

source("src/paths+packages.R")

# Define parameters and script controls -----------------------------------

# what year do CSV files start/end? should be same as LiberateData.m
minyear <- 1950
maxyear <- 2017
dates <- seq(ymd(paste0(minyear, "-01-01")), ymd(paste0(maxyear, "-12-31")),
             by="day")

# what is the maximum number of missing days that should be gap-filled?
max.missing.days <- 7

# Read in summary data ----------------------------------------------------

# only keep the ones defined as 'natural' (CSVs not created for natural==0 catchments)
df.summary <- 
  read.csv(file.path(dir.Q.raw, "summary_21-Mar-2018.csv"), stringsAsFactors=F) %>%
  subset(natural==1)

# get rid of unnecessary columns
cols.drop <- c("biomes",        # all NaN
               "drainagedensity",        # all NaN
               "FTC",                    # all NaN
               "GLWD",                   # all NaN
               "dailyQexists",           # not relevant
               "monthlyQexists",         # not relevant
               "meanQobs_entireperiod",  # not relevant
               "meanQobs_studyperiod",   # not relevant
               "reclength_entireperiod", # not relevant
               "reclength_studyperiod",  # not relevant
               "nonreference",           # all 0
               "soilI"                   # all NaN
               )
df.summary <- df.summary[, !(names(df.summary) %in% cols.drop)]

# set up output columns of interest
df.summary$yrs.complete          <- NaN  # number of years with no missing data
df.summary$yrs.complete.coverage <- NaN  # gaps that were filled within those complete years
df.summary$yrs.complete.first    <- NaN  # first year with no missing data
df.summary$yrs.complete.last     <- NaN  # last year with no missing data

# Cycle through catchments and calculate statistics -----------------------

n.catchment <- length(df.summary$catchment)
for (cat in 1:n.catchment){
  
  cat <- which(df.summary$catchment=="Brazil_14428000")
  
  # read in data
  cat.name <- df.summary$catchment[cat]
  df.cat <- 
    file.path(dir.Q.raw, "DATABASE_WORLD_V3", cat.name, "discharge+met_1950-2017.csv") %>% 
    read.csv(., stringsAsFactors=F)
  
  # add date, year, month columns
  df.cat$date  <- dates
  df.cat$year  <- year(df.cat$date)
  df.cat$month <- month(df.cat$date)
  df.cat$year_mo <- paste(df.cat$year, df.cat$month, sep="_")
  
  # make a column stating whether observation exists for that day
  df.cat$obs <- is.finite(df.cat$discharge.m3_s)
  
  # figure out which years to keep: no months with more than max.missing.days of missing discharge data
  df.cat.mo <-
    df.cat %>% 
    group_by(year, month) %>% 
    summarize(mo.days = mean(days_in_month(date)),
              mo.data.days = sum(is.finite(discharge.m3_s)),
              days.missing = mo.days-mo.data.days)
  
  yrs.keep <- 
    df.cat.mo %>% 
    subset(days.missing <= max.missing.days) %>% 
    group_by(year) %>% 
    summarize(keep = sum(is.finite(mo.data.days))==12) %>% 
    subset(.,keep) %>% 
    .$year

  if (length(yrs.keep)>0){
    
    # subset overall data to only years to keep
    df.cat <- subset(df.cat, year %in% yrs.keep)
    
    # gap fill small gaps using linear interpolation
    gapfill <- "spline"  # "linear" or "spline"
    if (gapfill=="spline"){
      df.cat$discharge.m3_s[!(df.cat$obs)] <- as.numeric(na.spline(df.cat$discharge.m3_s), maxgap=max.missing.days)[!(df.cat$obs)]
    }
    if (gapfill=="linear"){
      df.cat$discharge.m3_s[!(df.cat$obs)] <- as.numeric(na.approx(df.cat$discharge.m3_s), maxgap=max.missing.days)[!(df.cat$obs)]
    }
    
    # fill in summary output
    df.summary$yrs.complete[cat]          <- length(yrs.keep)
    df.summary$yrs.complete.coverage[cat] <- sum(df.cat$obs)/dim(df.cat)[1]
    df.summary$yrs.complete.first[cat]    <- min(yrs.keep)
    df.summary$yrs.complete.last[cat]     <- max(yrs.keep)
    
    # save gap-filled data as output
    cols.save <- c("date", "discharge.m3_s", "obs", "prec.mm_d", "Tmean.C", "Trange.C", "wind.m_s")
    path.save <- file.path(dir.Q.derived, "GapFill", paste0(cat.name, ".csv"))
    write.csv(df.cat[,cols.save], path.save, row.names=F)
    
  }
  
  # status update
  if ((cat %% 10)==0){
    print(paste0(cat, " of ", n.catchment, " complete"))
  }
}

# Save summary -------------------------------------------------------------

# save to GSAS and git repository
write.csv(df.summary, file.path(dir.Q.derived, "GapFill", "catchmentSummary_GapFill.csv"), row.names=F)
write.csv(df.summary, file.path("data", "catchmentSummary_GapFill.csv"), row.names=F)
