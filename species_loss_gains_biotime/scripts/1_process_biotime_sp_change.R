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

#biotime####
#read in full biotime dataset
load(file="../derived_data/biotime_abund.Rdata")

#head(biotime_abund)

#### Make summed and relative abundances of each species in each year ####
biotime_sum_abund <- biotime_abund %>%
  #get values summed across SAMPLE_DESCs within a year
  group_by(STUDY_ID, YEAR, GENUS_SPECIES, MEASURE) %>%
  summarize(SUM_VALUE = sum(VALUE)) %>%
  ungroup() %>%
  
  #make sure we have 0s for missing species
  group_by(STUDY_ID, MEASURE) %>%
  nest() %>%
  mutate(reshaped_df = purrr::map(data, ~spread(., GENUS_SPECIES, SUM_VALUE, fill=0) %>%
                                       gather(GENUS_SPECIES, SUM_VALUE, -YEAR))) %>%
  unnest(reshaped_df) %>%
  group_by(STUDY_ID, MEASURE) %>%
  mutate(NUM_SPECIES = n_distinct(GENUS_SPECIES)) %>%
  ungroup() %>%
  
  #now create summed abundances, relative abundances,
  #absolute rank, and relative rank
  #both including and excluding absent species
  #in each year relative to community
  group_by(STUDY_ID, YEAR, MEASURE) %>%
  mutate(TOTAL_VALUE = sum(SUM_VALUE),
         REL_VALUE_COMM = SUM_VALUE/TOTAL_VALUE,
         ABS_RANK = rank(REL_VALUE_COMM, ties.method = "min"), 
         REL_RANK_ALL = scales::rescale(ABS_RANK, c(0,1)),
         REL_RANK_TIMEPOINT = scales::rescale(ifelse(ABS_RANK==1,NA, ABS_RANK), c(0,1))) %>%
  ungroup() %>%
  
  #now create total series abund for later calcultion
  group_by(STUDY_ID, MEASURE) %>%
  mutate( WHOLE_SERIES_TOTAL = sum(SUM_VALUE)) %>%
  ungroup() %>%
  
  #now create summed abundances and relative abundances 
  #in each year relative to self in first appearance
  #and max and whole series relative abundance
  group_by(STUDY_ID, MEASURE, GENUS_SPECIES) %>%
  arrange(YEAR) %>%
  mutate(LOG_SUM_VALUE = log(SUM_VALUE),
         REL_VALUE_TO_MAX = SUM_VALUE/max(SUM_VALUE),
         REL_VALUE_TO_FIRST = SUM_VALUE/SUM_VALUE[SUM_VALUE!=0][1],
         REL_VALUE_TO_T0 = SUM_VALUE/SUM_VALUE[1],
         REL_VALUE_TO_T0_T3 = SUM_VALUE/mean(SUM_VALUE[1:3]),
         REL_RANK_T0 = REL_RANK_TIMEPOINT[1],
         REL_RANK_T1 = REL_RANK_TIMEPOINT[2],
         REL_RANK_T3 = REL_RANK_TIMEPOINT[2],
         REL_RANK_MEAN_T0_TO_T3 = mean(REL_RANK_TIMEPOINT[1:3]),
         Z_TRANS_VALUE = (SUM_VALUE - mean(SUM_VALUE))/sd(SUM_VALUE),
         Z_TRANS_LOG_VALUE = (LOG_SUM_VALUE - mean(LOG_SUM_VALUE))/sd(LOG_SUM_VALUE),
         CV_LOG_VALUE = sd(LOG_SUM_VALUE)/mean(LOG_SUM_VALUE),
         REL_VALUE_COMM_T0 = REL_VALUE_COMM[1],
         REL_VALUE_COMM_T0_T3 = mean(REL_VALUE_COMM[1:3]),
         N_YEARS = length(YEAR),
         WHOLE_SERIES_REL_VALUE = sum(SUM_VALUE)/WHOLE_SERIES_TOTAL) %>%
  arrange(desc(YEAR)) %>%
  mutate(MISSING_FINAL_3_YRS = sum(SUM_VALUE[1:3])>0) %>%
  ungroup()

  #consec years absent?
  
#### Determine which species get lost ####
biotime_loss <- biotime_sum_abund %>%
  group_by(STUDY_ID, GENUS_SPECIES, MEASURE) %>%

  #cut out first year
  arrange(YEAR) %>%
  slice(2:n()) %>%
  
  #get loss information
  summarize(n_year = 1+n_distinct(YEAR),
            n_zeros = sum(SUM_VALUE==0),
            lost_at_end = as.numeric(SUM_VALUE[n()]==0)) %>%
  ungroup() %>%
  mutate(frac_missing_after_t0 = n_zeros/n_year)


#### Merge different derived datasets ####
biotime_loss_summary <- biotime_sum_abund %>%
  #get year 0 only
  group_by(STUDY_ID, GENUS_SPECIES, MEASURE) %>%
  arrange(YEAR) %>%
  slice(1L) %>%
  ungroup() %>%
  
  #join with loss data
  left_join(biotime_loss) %>%
  
  #get rid of rel abundance of 0 in t0
  filter(SUM_VALUE != 0) %>%
  
  #are you always present in the timeseries?
  mutate(absent_at_some_point = ifelse(frac_missing_after_t0 !=0, 1, 0)) %>%

  #don't need year, as it's first year only
  dplyr::select(-YEAR)

#### Save derive data ####
save(biotime_loss_summary, file="../derived_data/biotime_loss_summary.Rdata")

save(biotime_sum_abund, file="../derived_data/biotime_sum_abund.Rdata")
