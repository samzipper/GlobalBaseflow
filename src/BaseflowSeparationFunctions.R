## BaseflowSeparationFunctions.R
#' Script to hold various functions for baseflow separation.

baseflow_HYSEP <- function(Q, area_mi2, method=NULL){
  # R implementation of USGS HYSEP baseflow separation algorithms
  # as described in Pettyjohn & Henning (1979) and implemented 
  # in Sloto & Crouse (1996).
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   area_mi2 = area in square miles
  #   method = HYSEP method; options are "fixed", "sliding", "local"
  #
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  if (sum(is.na(Q))>0){
    stop(paste0(sum(is.na(Q)), " missing data points. You need to gap-fill or remove it."))
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
    bf <- rollapply(Q, width=2*Nstar, FUN=min, 
                    align="center", partial=T)
    
    return(bf)
    
  } else if (method=="local"){
    
    ## local minimum
    # need to mirror 2Nstar days at the beginning/end to ensure no NaNs
    Q_mirror <- c(Q[(2*Nstar):1], Q, Q[length(Q):(length(Q)-(2*Nstar))])
    
    # local minima are points where Q is equal to the sliding interval minimum
    interval_min <- rollapply(Q_mirror, width=2*Nstar, FUN=min,
                              align="center", partial=T)
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
    stop("Wrong or missing method. Choose fixed, sliding, or local")
    
  }
  
}

baseflow_UKIH <- function(Q, endrule="NA"){
  # R implementation of UKIH baseflow separation algorithm as
  # described in Piggott et al. (2005)
  #
  # Inputs:
  #   Q = discharge timeseries (no missing data) (any units are OK)
  #   endrule = what to do with endpoints, which will always have NAs?
  #     "NA" (default) = retain NAs
  #     "Q" = use Q of the first/last point
  #     "B" = use bf of the first/last point
  #       
  #
  # Output:
  #   bf = baseflow timeseries, same length and units as Q
  
  ## package dependencies
  require(zoo)
  require(dplyr)
  require(magrittr)
  
  ## fixed interval of width 5
  int_width <- 5
  n.ints <- ceiling(length(Q)/int_width)
  ints <- rep(seq(1,n.ints), each=int_width)[1:length(Q)]
  
  # build data frame
  df <- data.frame(int = ints,
                   day = seq(1,length(ints)),
                   Q = Q)
  
  # summarize by interval
  df <- 
    df %>% 
    group_by(int) %>% 
    summarize(Qmin = min(Q),
              n_int = sum(is.finite(int))) %>% 
    left_join(df, ., by="int") %>% 
    subset(n_int==int_width)
  
  # extract minimum Qmin for each interval; these are
  # candidates to become turning points
  df.mins <- df[df$Q==df$Qmin, ]
  
  # if there are two minima for an interval (e.g. two 
  # days with same Q), choose the earlier one
  df.mins <- df.mins[!duplicated(df.mins$int), ]
  
  ## determine turning points, defined as:
  #    0.9*Qt < min(Qt-1, Qt+1)
  # do this using a weighted rolling min function
  df.mins$iQmin <- rollapply(df.mins$Qmin, width=3, align="center", 
                             fill=NA, FUN=function(z) which.min(z*c(1,0.9,1)))
  df.mins <- subset(df.mins, is.finite(iQmin))
  TP.day <- df.mins$day[df.mins$iQmin==2]
  TP.Qmin <- df.mins$Qmin[df.mins$iQmin==2]
  
  # linearly interpolate to length Q
  bf <- rep(NaN, length(Q))
  bf[TP.day] <- TP.Qmin
  bf <- as.numeric(zoo::na.approx(bf, na.rm=F))
  
  # need to fill in NAs?
  if (endrule=="Q"){
    # start
    bf[1:(TP.day[1]-1)] <- Q[1:(TP.day[1]-1)]
    
    # end
    bf[(TP.day[length(TP.day)]+1):length(Q)] <- 
      Q[(TP.day[length(TP.day)]+1):length(Q)]
    
  } else if (endrule=="B") {
    # start
    bf[1:(TP.day[1]-1)] <- bf[TP.day[1]]
    
    # end
    bf[(TP.day[length(TP.day)]+1):length(Q)] <- 
      bf[TP.day[length(TP.day)]]
  
  } else if (endrule != "NA") {
    
    stop("Invalid endrule")
    
  }
  
  # find any bf>Q and set to Q
  i_tooHigh <- which(bf>Q)
  bf[i_tooHigh] <- Q[i_tooHigh]
  return(bf)
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
dv$UKIH <- baseflow_UKIH(Q = dv$discharge.cfs, endrule="B")

dv.melt <- 
  dv%>% 
  subset(select=c("Date", "discharge.cfs", "HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH")) %>% 
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
p

## calculate BFI
dv.melt %>% 
  group_by(variable) %>% 
  summarize(discharge.sum = sum(discharge.cfs),
            baseflow.sum = sum(value),
            BFI = baseflow.sum/discharge.sum)
