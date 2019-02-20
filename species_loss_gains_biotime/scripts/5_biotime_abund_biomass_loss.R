#########################################
#' Process BIOTIME abundance data to get
#' species loss and change in relative abundance
#' 
#' 
#' @author Jarrett Byrnes
#' @details 
#' Changelog:
#########################################

#packages####
library(tidyverse)
library(data.table)
library(ggplot2)

theme_set(theme_classic(base_size = 16))

#biotime####
#read in full biotime dataset
#load(file="../derived_data/biotime_abund.Rdata")
load(file="../derived_data/biotime_sum_abund.Rdata")

biotime_duo <- biotime_sum_abund %>%
  group_by(STUDY_ID) %>%
  mutate(num_measures = length(unique(MEASURE))) %>%
  filter(num_measures==2) 


biotime_duo_reshape <- biotime_duo %>%
  dplyr::select(STUDY_ID, YEAR, GENUS_SPECIES, MEASURE, SUM_VALUE) %>%
  spread(MEASURE, SUM_VALUE) %>%
  group_by(STUDY_ID, YEAR) %>%
  mutate(abund_rank = rank(Abundance, ties = "random"),
         biomass_rank = rank(Biomass, ties = "random"),
         percapita_biomas_rank = rank(Biomass/Abundance)
  ) %>%
  group_by(STUDY_ID) %>%
  mutate(rel_abund = Abundance / max(Abundance, na.rm=TRUE),
         rel_biomass = Biomass / max(Biomass, na.rm=TRUE),
         rel_abund_rank = abund_rank / max(abund_rank, na.rm=TRUE),
         rel_biomass_rank = biomass_rank / max(biomass_rank, na.rm=TRUE),
         rel_percapita_biomas_rank = percapita_biomas_rank / max(percapita_biomas_rank, na.rm=TRUE),
  )

ggplot(biotime_duo_reshape, 
       mapping = aes(x = log(Abundance+1), y = log(Biomass+1), color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
  stat_smooth(method = "lm")
ggsave("../figures/biotime_log_abund_biomass.jpg")


#plot abundance rank by biomass rank 
ggplot(biotime_duo_reshape, 
       mapping = aes(x = rel_abund_rank, y = rel_biomass_rank, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
  stat_smooth(method = "lm", mapping = aes(group = 1))
  #stat_density_2d(aes(group = 1))
ggsave("../figures/rel_abund_rel_biomass.jpg")



#plot abundance rank by per capita biomass rank
ggplot(biotime_duo_reshape, 
       mapping = aes(x = rel_abund_rank, y = rel_percapita_biomas_rank, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
  stat_smooth(method="lm", mapping = aes(group = 1))
ggsave("../figures/rel_abund_rel_percapita_biomass.jpg")


#summed across dataset

#figure out who is plant, animal, alga, etc



#relationship between lost species biomass rank and abundance rank

#relationship between lost species per capita biomass rank and abundance rank

#slopes of species trajectories relative to 


#biomass timeseries
biotime_duo_reshape %>%
  group_by(STUDY_ID, YEAR) %>%
  summarize(biomass = sum(Biomass)) %>%
  ggplot(aes(x = YEAR, y = log(biomass), color = factor(STUDY_ID))) +
  geom_line() +
  guides(color = "none") 
ggsave("../figures/biomass_change_over_time.jpg")