## BaseflowSeparationFunctions.R
#' Script to hold various functions for baseflow separation.

baseflow_HYSEP <- function(Q, area_mi2, method=NULL){
  # R implementation of USGS HYSEP baseflow separation algorithms
  # as described in Sloto & Crouse (1996)
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data)
  #   area_mi2 = area in square miles
  #   method = HYSEP method; options are "fixed", "sliding", "local"
  #
  # Output:
  #   bf = baseflow timeseries, same length as Q
  
  ## package dependencies
  require(caTools)
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  if (is.null(method)){
  }
  
  ## calculate N*
  N <- area_mi2^0.2
  Nstar <- 2*floor(N/2)+1  # round to nearest odd integer
  
  ## calculation depends on method
  if (method=="fixed"){
    
    ## fixed interval of width 2Nstar
    n.ints <- ceiling(length(Q)/(2*Nstar))
    ints <- rep(seq(1,n.ints), each=2*Nstar)[1:length(Q)]
    
    # build data frame
    df <- data.frame(int = ints,
                     Q = Q)
    
    # summarize by interval
    df <- df %>% 
      group_by(int) %>% 
      summarize(bf = min(Q)) %>% 
      left_join(df, ., by="int")
    
    return(df$bf)
    
  } else if (method=="sliding"){
    
    ## sliding interval of width 2Nstar
    bf <- caTools::runmin(Q, 2*Nstar)
    
    return(bf)
    
  } else if (method=="local"){
    
    ## local minimum
    # need to mirror 2Nstar days at the beginning/end to ensure no NaNs
    Q_mirror <- c(Q[(2*Nstar):1], Q, Q[length(Q):(length(Q)-(2*Nstar))])
    
    # local minima are points where Q is equal to the sliding interval minimum
    interval_min <- caTools::runmin(Q_mirror, 2*Nstar)
    i_minima <- which(interval_min==Q_mirror)
    
    # interpolate between values using na.approx
    bf_mirror <- rep(NaN, length(Q_mirror))
    bf_mirror[i_minima] <- Q_mirror[i_minima]
    bf_mirror <- as.numeric(zoo::na.approx(bf_mirror, na.rm=F))
    
    # trim off mirrored portions
    bf <- bf_mirror[(2*Nstar+1):(length(bf_mirror)-(2*Nstar)-1)]
    
    # find any bf>Q and set to Q
    i_tooHigh <- which(bf>Q)
    bf[i_tooHigh] <- Q[i_tooHigh]
    
    return(bf)
    
  } else {
    
    # error
    stop("Wrong or missing method: choose fixed, sliding, or local")
    
  }
  
}

## example: Des Moines River at Fort Dodge
# packages required for sample data/examples
require(dataRetrieval)
require(lubridate)
require(ggplot2)
require(magrittr)
require(reshape2)

# get USGS data for Des Moines River at Fort Dodge, IA
dv <- readNWISdv(siteNumber="05480500",                          # site code (can be a vector of multiple sites)
                 parameterCd="00060",                            # parameter code: "00060" is cubic ft/sec
                 startDate="2013-01-01",endDate="2014-01-01",    # start & end dates (YYYY-MM-DD format)
                 statCd = "00003") %>%                           # statistic code: "00003" is daily mean (default)
  set_colnames(c("agency_cd", "site_no", "Date", "discharge.cfs", "QA.code"))

# area in square miles
area_mi2 <- 4190

# check for missing data
sum(is.na(dv$discharge.cfs))

## perform baseflow separations
dv$HYSEP_fixed <- baseflow_HYSEP(Q = dv$discharge.cfs, area_mi2 = area_mi2, method="fixed")
dv$HYSEP_slide <- baseflow_HYSEP(Q = dv$discharge.cfs, area_mi2 = area_mi2, method="sliding")
dv$HYSEP_local <- baseflow_HYSEP(Q = dv$discharge.cfs, area_mi2 = area_mi2, method="local")

dv.melt <- 
  dv%>% 
  subset(select=c("Date", "discharge.cfs", "HYSEP_fixed", "HYSEP_slide", "HYSEP_local")) %>% 
  melt(id=c("Date", "discharge.cfs"))

p <- 
  ggplot(dv.melt) +
  geom_ribbon(data=dv, aes(x=Date, ymin=0, ymax=discharge.cfs), fill="black") +
  geom_line(aes(x=Date, y=value, color=variable)) +
  scale_y_continuous(name="Discharge [cfs]") +
  scale_x_date(expand=c(0,0)) +
  scale_color_discrete(name="Method") +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        legend.position=c(0.99, 0.99),
        legend.justification=c(1,1))

## calculate BFI
dv.melt %>% 
  group_by(variable) %>% 
  summarize(discharge.sum = sum(discharge.cfs),
            baseflow.sum = sum(value),
            BFI = baseflow.sum/discharge.sum)
