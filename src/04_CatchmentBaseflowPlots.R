## 04_CatchmentBaseflowPlots.R
#' This script is intended to read in CSV files with baseflow summary statistics
#' (produced with the script 03_CatchmentBaseflowStats.R) and make some graphs.

source(file.path("src", "paths+packages.R"))

# prep map data frame
mapWorld <- borders("world", colour="white", fill="#7f7f7f", alpha=0.25) # create a layer of borders

# Read in data ------------------------------------------------------------

# summary data - additional catchment characteristics are in the file catchmentsSummary_GapFill.csv
df.summary <- 
  read.csv(file.path("data", "catchmentSummary_BaseflowStats.csv"), stringsAsFactors=F) %>% 
  left_join(., read.csv(file.path("data", "catchmentSummary_GapFill.csv"), stringsAsFactors=F),
            by="catchment")

df.all <- 
  read.csv(file.path("data", "catchmentSummary_BaseflowStats_AllMethods.csv"), stringsAsFactors=F)

# set factor levels
df.all$flux <- factor(df.all$flux, 
                      levels=c("Q_mm.d", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt", 
                               "HYSEP_fixed", "HYSEP_local", "HYSEP_slide", "UKIH"),
                      labels=c("Discharge [mm/d]", "BFLOW (1 pass)", "BFLOW (3 pass)", 
                               "Eckhardt", "HYSEP (fixed)", "HYSEP (local)", "HYSEP (sliding)", "UKIH"))

# Filter data and eliminate bad stations ----------------------------------

# some catchments have super high runoff - I think due to catchment size error
# eliminate anything with runoff > 5000 mm/yr
df.summary <- subset(df.summary, Q_mm.y < 5000)

df.all <- subset(df.all, catchment %in% unique(df.summary$catchment))

# data frame for site characteristics only
df.info <- 
  df.summary %>% 
  dplyr::select(catchment, Q_mm.y, aridity, catchsize, ETPOT_Hargr, lat, lon) %>% 
  unique()

# Derived metrics ---------------------------------------------------------

# annual range of BFI, bf
df.all$rangeBFI_mo <- df.all$maxBFI_mo_BFI - df.all$minBFI_mo_BFI
df.all$rangeFlow_mo_mm.d <- df.all$maxFlow_mo_mm.d - df.all$minFlow_mo_mm.d

# recession constant, in days
df.summary$recessionConstantDays  <- -1/log(df.summary$recessionConstant)    # eq. 6 from Barlow digital filters tutorial

# summarize differences between different baseflow separation methods
bf_methods <- c("HYSEP_fixed", "HYSEP_slide", "HYSEP_local", "UKIH", "BFLOW_1pass", "BFLOW_3pass", "Eckhardt")
df.methods <- 
  df.all %>% 
  subset(flux %in% bf_methods) %>% 
  group_by(catchment) %>% 
  summarize(BFI_range = (max(BFI) - min(BFI)),
            BFI_sd = sd(BFI),
            BFI_mean = mean(BFI),
            BFI_CV = BFI_sd/BFI_mean,
            bf_range = (max(bf_mm.y) - min(bf_mm.y)),
            bf_sd = sd(bf_mm.y),
            maxFlow_mo_range = min.month.range(maxFlow_mo),
            maxFlow_DOY_range = min.DOY.range(maxFlow_DOY))

# decision tree
df.DecisionTree <-
  df.all %>% 
  subset(flux != "Discharge [mm/d]") %>% 
  left_join(df.summary[,c("catchment", "noFlowDays_prc")], by=c("catchment"))
df.DecisionTree$DecisionTree <- "Unknown"
df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc > 0.1] <- "Int"

df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc < 0.1 &
                               df.DecisionTree$BFI < 0.5 &
                               df.DecisionTree$maxBFI_mo_BFI < 0.5] <- "Per_Qf"

df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc < 0.1 &
                               df.DecisionTree$BFI < 0.5 &
                               df.DecisionTree$maxBFI_mo_BFI >= 0.5 &
                               (df.DecisionTree$maxBFI_mo < 4 | df.DecisionTree$maxBFI_mo > 9)] <- "Per_Qf_highBf_O-M"

df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc < 0.1 &
                               df.DecisionTree$BFI < 0.5 &
                               df.DecisionTree$maxBFI_mo_BFI >= 0.5 &
                               (df.DecisionTree$maxBFI_mo >= 4 | df.DecisionTree$maxBFI_mo <= 9)] <- "Per_Qf_highBf_A-S"

df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc < 0.1 &
                               df.DecisionTree$BFI >= 0.5 &
                               df.DecisionTree$minBFI_mo_BFI >= 0.5] <- "Per_Bf_lowQf"

df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc < 0.1 &
                               df.DecisionTree$BFI >= 0.5 &
                               df.DecisionTree$minBFI_mo_BFI < 0.5 &
                               (df.DecisionTree$minBFI_mo >= 4 | df.DecisionTree$minBFI_mo <= 9)] <- "Per_Bf_highQf_A-S"

df.DecisionTree$DecisionTree[df.DecisionTree$noFlowDays_prc < 0.1 &
                               df.DecisionTree$BFI >= 0.5 &
                               df.DecisionTree$minBFI_mo_BFI < 0.5 &
                               (df.DecisionTree$minBFI_mo < 4 | df.DecisionTree$minBFI_mo > 9)] <- "Per_Bf_highQf_O-M"

table(df.DecisionTree$DecisionTree)

# Plots -------------------------------------------------------------------

# distribution of recessionConstant
df.summary %>% 
  ggplot(aes(x=recessionConstant)) +
  geom_histogram(binwidth=0.01)

# distribution of recessionConstantDays
df.summary %>% 
  ggplot(aes(x=recessionConstantDays)) +
  geom_histogram(binwidth=10) +
  geom_vline(xintercept=45, color="red") +
  annotate("rect", xmin=30, xmax=60, ymin=-Inf, ymax=Inf, fill="red", alpha=0.5)

sum(df.summary$recessionConstantDays >= 30 & df.summary$recessionConstantDays <= 60, na.rm=T)
sum(is.finite(df.summary$recessionConstantDays))

# distribution of BFImax
df.summary %>% 
  ggplot(aes(x=BFImax)) +
  geom_histogram(binwidth=0.05)

# distribution of Q quantiles
df.summary %>% 
  dplyr::select(catchment, Q10, Q50, Q90) %>% 
  melt(id="catchment") %>% 
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~variable)

# distribution of FDC_slope
df.summary %>% 
  ggplot(aes(x=FDC_slope)) +
  geom_histogram()

# distribution of dQ_dt.slope
df.summary %>% 
  ggplot(aes(x=dQ_dt.slope)) +
  geom_histogram()

## comparison of range among different methods
# BFI
df.all %>% 
  ggplot(aes(x=BFI)) +
  geom_density() +
  facet_wrap(~flux, scales="free_y")

df.methods %>% 
  ggplot(aes(x=BFI_range)) +
  geom_histogram(binwidth=0.05)

df.methods %>% 
  ggplot(aes(x=BFI_sd)) +
  geom_histogram()

df.methods %>% 
  ggplot(aes(x=BFI_mean)) +
  geom_histogram()

df.methods %>% 
  ggplot(aes(x=BFI_CV)) +
  geom_histogram()

# bf
df.methods %>% 
  ggplot(aes(x=bf_range)) +
  geom_histogram()

df.methods %>% 
  ggplot(aes(x=bf_sd)) +
  geom_histogram()

# maxFlow_mo
df.methods %>% 
  ggplot(aes(x=maxFlow_mo_range)) +
  geom_histogram(binwidth=1)

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=aridity, y=maxFlow_mo_range)) +
  geom_point() +
  stat_smooth(method="lm")

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=Q_mm.y, y=maxFlow_mo_range)) +
  geom_point() +
  stat_smooth(method="lm")

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  lm(maxFlow_mo_range ~ Q_mm.y, data=.) %>% 
  summary()

## scatterplot matrices
df.all %>% 
  dplyr::select(flux, catchment, maxFlow_mo) %>% 
  dcast(catchment ~ flux, value.var = "maxFlow_mo") %>% 
  ggpairs(columns=levels(df.all$flux))

df.all %>% 
  dplyr::select(flux, catchment, maxFlow_DOY) %>% 
  dcast(catchment ~ flux, value.var = "maxFlow_DOY") %>% 
  ggpairs(columns=levels(df.all$flux))


## maps
df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=Q_mm.y)) +
  mapWorld +
  geom_point(shape=21) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90))

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=BFI_range)) +
  mapWorld +
  geom_point(shape=21) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90))

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=BFI_mean)) +
  mapWorld +
  geom_point(shape=21) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90))

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=BFI_CV)) +
  mapWorld +
  geom_point(shape=21) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90))

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=maxFlow_mo_range)) +
  mapWorld +
  geom_point() +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90))

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, 
             color=cut(maxFlow_mo_range, breaks=c(0,1,2,3,6,12), 
                       include.lowest=T, right=F))) +
  mapWorld +
  geom_point(shape=21) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
  scale_color_discrete(name="Month Range")

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=maxFlow_DOY_range)) +
  mapWorld +
  geom_point() +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90))

df.methods %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, 
             color=cut(maxFlow_DOY_range, breaks=c(0,15,30,60,180,366), 
                       include.lowest=T, right=F))) +
  mapWorld +
  geom_point(shape=21) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
  scale_color_discrete(name="DOY Range")

df.all %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=maxFlow_mo)) +
  mapWorld +
  geom_point(shape=21) +
  facet_wrap(flux~.) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
  scale_color_gradient2(low=col.cat.grn, mid=col.cat.red, high=col.cat.blu, midpoint=6) +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

df.all %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=rangeBFI_mo)) +
  mapWorld +
  geom_point(shape=21) +
  facet_wrap(flux~.) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

df.all %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=rangeFlow_mo_mm.d)) +
  mapWorld +
  geom_point(shape=21) +
  facet_wrap(flux~.) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

df.DecisionTree %>% 
  left_join(df.info, by="catchment") %>% 
  ggplot(aes(x=lon, y=lat, color=DecisionTree)) +
  mapWorld +
  geom_point(shape=21) +
  facet_wrap(flux~.) +
  scale_x_continuous(expand=c(0,0), limits=c(-180, 180)) +
  scale_y_continuous(expand=c(0,0), limits=c(-90, 90)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

