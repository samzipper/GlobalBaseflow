## 02_CatchmentSelect+Separate.R
#' This script is intended to read in CSV files for all natural catchments
#' (produced with the script 01_CatchmentGapFill+Stats.R), decide which 
#' catchments for final analysis, and perform baseflow separation on them.

source("src/paths+packages.R")

# Define parameters and script controls -----------------------------------



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

# make plots of retained catchments
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
ggsave(file.path("plots", "CatchmentSelect_HistographProperties.png"),
       p.properties.hist, width=8, height=6, units="in")
