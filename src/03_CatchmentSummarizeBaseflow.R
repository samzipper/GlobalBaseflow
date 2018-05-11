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
    read.csv(., stringsAsFactors=F)
  
  # calculate and save monthly means
  df.cat.mo <-
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
  write.csv(df.cat.mo, path.save, row.names=F)
}
