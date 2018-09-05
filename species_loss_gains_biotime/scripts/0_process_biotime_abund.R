#########################################
#' Process BIOTIME data into abundances
#' per species in areas with fixed sampling
#' 
#' Code completely gakked from sChange working group
#'https://github.com/sChange-workshop/BiodivChange/blob/master/code/process%20biotime%20for%20scale%20analysis.R
#'
#' Generates ../derived_data/biotime_abundances_processed.RData
#' 
#' @author Jarrett Byrnes
#' @details 
#' Changelog:
#########################################

#packages####
library(tidyverse)
library(data.table)
library(iNEXT)

#biotime####
#read in full biotime dataset
#biotime_all<-fread("../raw_data/BioTIMEMarch.csv")
#saveRDS(biotime_all, file="../raw_data/BioTIMEMarch.rd")
biotime_all<-readRDS("../raw_data/BioTIMEMarch.rd")
biotime_meta<-read.csv("../raw_data/bioTIMEmetadataFeb.csv")

#calculate lat and long means and range
Biotime_latlong<-biotime_all %>% 
  group_by(YEAR,STUDY_ID) %>% 
  summarise(Lat_mean=mean(LATITUDE),Long_mean=mean(LONGITUDE),Lat_max=max(LATITUDE),Long_max=max(LONGITUDE),Lat_min=min(LATITUDE),Long_min=min(LONGITUDE)) %>% 
  arrange(STUDY_ID,YEAR)

#calculate sd of lat long means and range
Latlong_sd<-Biotime_latlong %>% 
  group_by(STUDY_ID) %>% 
  summarise_all(.funs = sd) %>% 
  arrange(desc(Lat_mean))

sd_cutoff<-0.3 #set this to whatever threshold of variation in lat long standard deviation you find acceptable

Fixed_locations<-Latlong_sd %>% 
  filter(Long_mean<sd_cutoff & 
           Lat_mean<sd_cutoff &
           Lat_max<sd_cutoff &
           Lat_min<sd_cutoff &
           Long_max<sd_cutoff &
           Long_min<sd_cutoff) %>% 
  arrange(desc(Lat_mean))

#subset datasets with fixed locations
biotime_fixed<-biotime_all %>% 
  filter(STUDY_ID %in% Fixed_locations$STUDY_ID)

#filter years with poor coverage####
names(biotime_fixed)[11:12]<-c("Abundance","Biomass")

#calculate abundance and biomass for all species in each plot in each year
biotime_summary <- biotime_fixed %>%
  group_by(STUDY_ID, YEAR,SAMPLE_DESC,GENUS_SPECIES) %>%
  summarise(
    Abundance = sum(Abundance, na.rm=TRUE),
    Biomass = sum(Biomass, na.rm=TRUE)) %>% 
  ungroup()

biotime_summary<-merge.data.frame(biotime_meta %>% 
                                    select(STUDY_ID,ABUNDANCE_TYPE),
                                  biotime_summary, by='STUDY_ID')


abund_coverage <- biotime_summary %>% 
  # get only rows representing count data (abundance)
  filter(ABUNDANCE_TYPE!='Presence/Absence' & ABUNDANCE_TYPE!='<NA>') %>%
  # remove zeroes and NAs(some studies are designated count, but record Biomass estimates)
  filter(Abundance > 0 & !is.na(Abundance)) %>%
  group_by(STUDY_ID, YEAR,SAMPLE_DESC) %>%
  summarise(
    # how many singletons
    singletons = sum(Abundance==1),
    # how many doubletons
    doubletons = sum(Abundance==2),
    # how many individuals in total sample
    N = sum(Abundance),
    # eqn 4a in Chao & Jost 2012 Ecology (==eqn 12 in Chao et al 2014 Ecol Monogr)
    Chat = 1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)),
    # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
    # correction for communities with no doubletons (this prevents NaN results for either singletons==0 or doubletons==0)
    Chat_corrected = ifelse(doubletons>0, 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)), 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2))),
    # the iNEXT coverage calculation has some extras in it that I don't understand
    # e.g., it corrects using f0.hat, the number of unobserved species (Chao1 estimator)? Why? Exactly how? Ref?
    f1_iNEXT = DataInfo(Abundance)$f1,
    f2_iNEXT = DataInfo(Abundance)$f2,
    #n_iNEXT = DataInfo(Abundance)$n,
    coverage = DataInfo(Abundance)$SC) %>%
  ungroup()

# are my estimates of singletons and doubletons equal to the iNEXT estimates?
sum(abund_coverage$singletons!=abund_coverage$f1_iNEXT)		# yes, good
sum(abund_coverage$doubletons!=abund_coverage$f2_iNEXT)		# yes, good

mn=mean(abund_coverage$Chat_corrected, na.rm=TRUE)
sd=sd(abund_coverage$Chat_corrected, na.rm=TRUE)

#histogram of corrected coverage
hist(abund_coverage$Chat_corrected)

# BIOMASS DATA coverage (as INCIDENCE)... by rarefyID (year as samples)
# need to check if I am doing this correctly - PLT - Feb 24, 2017
biotime_summary_biomass<-biotime_summary %>% 
  filter(is.na(ABUNDANCE_TYPE)) %>% 
  mutate(incidence = ifelse(Biomass==0, 0, 1)) %>% 
  group_by(STUDY_ID,YEAR,GENUS_SPECIES) %>% 
  summarize(incidence=sum(incidence))

head(biotime_summary_biomass)
biotime_summary_biomass %>% 
  group_by(STUDY_ID,YEAR)

biomass_coverage <- biotime_summary_biomass %>% 
  # get only rows representing count data (abundance)
  group_by(STUDY_ID, YEAR) %>%
  summarise(
    # how many singletons
    singletons = sum(incidence==1),
    # how many doubletons
    doubletons = sum(incidence==2),
    # how many individuals in total sample
    N = sum(incidence),
    # eqn 4a in Chao & Jost 2012 Ecology (==eqn 12 in Chao et al 2014 Ecol Monogr)
    Chat = 1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)),
    # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
    # correction for communities with no doubletons (this prevents NaN results for either singletons==0 or doubletons==0)
    Chat_corrected = ifelse(doubletons>0, 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)), 
                            1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2))),
    # the iNEXT coverage calculation has some extras in it that I don't understand
    # e.g., it corrects using f0.hat, the number of unobserved species (Chao1 estimator)? Why? Exactly how? Ref?
    f1_iNEXT = DataInfo(incidence)$f1,
    f2_iNEXT = DataInfo(incidence)$f2,
    #n_iNEXT = DataInfo(Abundance)$n,
    coverage = DataInfo(incidence)$SC) %>%
  ungroup()

mean(biomass_coverage$Chat_corrected,na.rm=T)
sd(biomass_coverage$Chat_corrected,na.rm=T)
hist(biomass_coverage$Chat_corrected)

## Make a list of the rarefyIDs to keep for analysis based on coverage threshold (>0.8)
##  NOTE: Count data will drop individual samples that don't meet the criteria, and Presence and Biomass data will drop entire years 

#Datasets to keep
biotime_abundance<-biotime_summary %>% 
  # get only rows representing count data (abundance)
  filter(ABUNDANCE_TYPE!='Presence/Absence' & ABUNDANCE_TYPE!='<NA>') %>% 
  merge(abund_coverage) %>% 
  filter(Chat_corrected>=0.8)

biotime_biomass<-biotime_summary %>% 
  filter(is.na(ABUNDANCE_TYPE)) %>% 
  merge(biomass_coverage) %>% 
  filter(Chat_corrected>=0.8)

#bind datasets back together - this is the dataset to use for rarefaction
biotime_filtered<-bind_rows(biotime_abundance,biotime_biomass)

#this is the proportion of the data that was kept after the coverage filter
dim(biotime_filtered)[1]/dim(biotime_summary)[1]
#this is the proportion of studies that were kept
length(unique(biotime_filtered$STUDY_ID))/length(unique(biotime_summary$STUDY_ID))

#rarefy####
#start with abundance
resample<-100
unique(biotime_filtered$STUDY_ID)

min_samples<-biotime_filtered %>% 
  group_by(STUDY_ID,YEAR,SAMPLE_DESC) %>% 
  summarise(abundance=sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(STUDY_ID,YEAR) %>% 
  summarise(samples=n()) %>% 
  summarise(min_samples=min(samples),
            max_samples=max(samples),
            prop=min_samples/max_samples)

biotime_filtered<-merge(biotime_filtered,min_samples)

biotime_filtered<-biotime_filtered %>% 
  filter(YEAR>1880)

save(biotime_filtered, file="../derived_data/biotime_filtered.Rdata")

#### New reshaping for nceas lter biodiv scaling here ####
biotime_abund <- biotime_filtered %>%
  gather(MEASURE, VALUE, Abundance:Biomass) %>% 
  group_by(STUDY_ID, MEASURE) %>%
  mutate(remove_me = !sum(VALUE)==0) %>%
  ungroup() %>%
  filter(remove_me) %>%
  select(-remove_me)

#### Save resulting file ####
save(biotime_abund, file="../derived_data/biotime_abund.Rdata")

