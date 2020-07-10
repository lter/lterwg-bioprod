#########################################
#' Analyze BIOTIME Biomass Timeseries
#' 
#' 
#' @author Jarrett Byrnes
#' @details 
#' Changelog:
#########################################

library(tidyverse)
library(ggplot2)

###

#### Set some graphical parameters

theme_set(theme_bw(base_size=17))# +
#             theme(panel.grid.major = element_blank(),
#                   panel.grid.minor = element_blank(), 
#                   axis.line = element_line(colour = "black")))

#for a labeller
add_year <- function(x) paste(x,"years", sep=" ")

#### Load data ####
#load(file="../derived_data/biotime_loss_summary_filtered.Rdata")
load(file="../derived_data/biotime_sum_abund.Rdata")


biotime_biomass <- biotime_sum_abund %>%
  filter(MEASURE=="Biomass") %>%
  group_by(STUDY_ID) %>%
  filter(length(unique(YEAR))>3) %>%
  group_by(STUDY_ID, YEAR) %>%
  summarize(Z = sum(SUM_VALUE==0),
            S = n_distinct(GENUS_SPECIES)-Z,
            B = sum(SUM_VALUE, na.rm=T),
            N_YEARS=N_YEARS[1]) %>%
  ungroup() 
 
biomass_timeseries <- ggplot(biotime_biomass,
         aes(x = YEAR, y = log(B+1), color = factor(STUDY_ID))) +
  geom_line(size = 0.7) +
  guides(color = "none") +
  xlab("") + ylab("Log Total Biomass + 1")
biomass_timeseries

# same as above, but add an overall trend line (locally linear in this case)
(biomass_timeseries_overall <- ggplot(biotime_biomass, aes(x = YEAR, y = log(B+1))) +
    geom_smooth(color = "black", alpha = 0.2, size = 1.5) + 
    geom_line(aes(color = factor(STUDY_ID)), size = 0.7) +
    guides(color = "none") +
    xlab("") + ylab("Log Total Biomass + 1")
)

#ggsave("../figures/biomass_change_over_time.jpg")

# saveRDS(list(biomass_timeseries = biomass_timeseries),
#         file = "../derived_data/6_biomass_timeseries.Rds")  
