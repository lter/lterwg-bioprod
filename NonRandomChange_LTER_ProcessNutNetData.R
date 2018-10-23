##############################################################################################
### Process Nut Net Cover Data - Dominance Index and Other Abundance Attributes ############
## Laura Dee  Sept 14, 2018  #########################################################
## updated Oct 11, 2018 #########################################################
# to add the computation of the metric "site-level variable of ave relative abundance in plots per spp."
###################################################################################################

#Close graphics and clear local memory
graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(foreign)
library(rmarkdown)

setwd("~/Dropbox/IV in ecology/NutNet")
cover <- fread('full-cover-09-April-2018.csv',na.strings= 'NA')

## need to make max_cover NOT a character
cover$max_cover <- as.numeric(cover$max_cover)

###### If want to filter cover dataset to only the live cover ################# 
# cover <- cover[live == 1,]
# table(cover$live)

#########################################################################################
## Compare Taxon in Live (live==1) or non-live (live == 0) ################################
##########################################################################################
## determine how many Taxon are listed as live==0 ##
table(cover$live==1)
a <- table(cover$Taxon, cover$live==0)
head(a)

# It looks like these categories are primarily the "live == 0 " ones;
# seems like I should NOT include them inthe SR counts?
# OTHER ANIMAL DIGGING                            0   47
# OTHER ANIMAL DIGGINGS                           0   68
# OTHER ANIMAL DROPPINGS                          0   75
# OTHER ANT NEST                                  0    1
# OTHER CRUST                                     0   13
# OTHER DISTURBED SOIL                            0   16
# OTHER LIGNOTUBER                                0   11
# OTHER LITTER                                    0 1542
# OTHER ROCK                                      0  175
# OTHER SOIL BIOCRUST                             0    5
# OTHER STANDING DEAD                             0    2
# OTHER TRIODIA BASEDOWII (DEAD)                  0    5
# OTHER UNKNOWN SOIL_CRUST                        0     17
# OTHER WOOD                                      0     27
#       FUNGI                                       	0   	2
#       FUNGI SP.                                   	0	    2
#       GROUND                                      	0	   1039
#       OTHER WOODY OVERSTORY                       	0	    4


##### TOTAL & TOTAL LIVE COVER MEASURES ########################################
# make a total cover in a plot, site, year
cover[,totplotcover.yr := sum(max_cover, na.rm=T), by=.(plot, site_code, year)]

# same for LIVE species only
cover[,totplotcover.yr.live := sum(max_cover[live=="1"], na.rm=T), by=.(plot, site_code, year)]

#Make a relative cover for each species in each plot, site, year
# based on TOTAL cover (including dead cover)
cover[,relative_sp_cover.yr := max_cover/totplotcover.yr]

# same for LIVE species only
cover[,relative_sp_cover.yr.live := (max_cover*(live=="1"))/totplotcover.yr.live]

# in some cases, no live species in plot and year, so getting NA since totplotcover.yr.live is zero. Set these to zero.
cover[is.na(relative_sp_cover.yr.live),relative_sp_cover.yr.live:=0]

#Average relative abundance over time per species, plot, site considering: LIVE cover in a plot over time
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

######################################################################################################
### Simplify some of the names in the cover data #####################################################
#######################################################################################################

# covert Naturalised to INT (This use of Naturalised as a category is from one site)
cover[local_provenance=="Naturalised", local_provenance:="INT"]
# convert blank entries to NULL
cover[local_provenance=="", local_provenance:="NULL"]

## Combine Graminoid and Grass into a single category 
# Graminoid + Grass = Graminoids
cover[functional_group == "GRASS" , functional_group:= "GRAMINOID"]

####################################################################################################################
### Dominance Variables ############################################################################################
####################################################################################################################

### Do relative abundance of live cover only (discussed w Kim in June)

### Dominance indicator calculation ###
# Use the dominance indicator metric from Avoilio et al (seperates dominance indication from impact)
# DI = (average relative abundance + relative frequency)/2
#Relative abundance = abundance of a species a in sampling unit / total abundance of all species in a sampling unit
#Relative frequency = number of sampling units a species occurred / total number of sampling units
# note: ranges from 0-1; Relative abundance can be any measure of abundance. Does not incorporate a measure of impact.
# There is not a cutoff for "which range from 0-1 = dominant species versus subordinate or rare, in the Avolio et al paper # 
# so if we want to group species in these groups, we will need to make one (then also should test robustness to that decision)


# relative abundance per species and plot in pre-treatment year (year_trt == 0) 
# [NOTE its not the average per the Avolio et al metric bc only one plot and year]
cover[, rel_abundance_year0 := relative_sp_cover.yr.live[year_trt == 0], by = .(Taxon, plot, site_code)]

# average this, to get a site-level "relative abundance in a plot" for each species "Taxon"
cover[, Site.rel_abundance_year0 := ave(rel_abundance_year0), by = .(Taxon,  site_code)]


## Compute Relative Frequency per spp AT THE SITE (defining dominance in space, not time/across years)
# "Relative frequency = number of sampling units a species occurred / total number of sampling units" 
# this should be # of plots within a site that the species occured in / total # of plots within a site, for pre-treatment year 

#total # of plots within a site, for pre-treatment year 
# & to filter to do just on the live species:
cover[, tot.num.plots := length(unique(plot[year_trt == 0])), by =.(site_code)]

#number of plots within a site, in the pre-treatment year, a species occurred in:
# & to filter to do just on the live species:
cover[, tot.num.plots.with.spp := length(unique(plot[year_trt == 0 & live==1])), by =.(site_code, Taxon)]

# test to see if it works
abisko.test = cover[site_code == "abisko.se" , ]

#Compute Relative Frequench
# " Relative frequency = number of sampling units a species occurred / total number of sampling units" 
cover[, rel_freq.space :=  tot.num.plots.with.spp/tot.num.plots, by = .(plot, site_code)]
#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq.space)
summary(cover$rel_freq.space)

## Compute the DI per species 
# DI = (average relative abundance + relative frequency)/2

#FOR LIVE COVER -- the dead cover will be 0.
cover[, DI := (rel_abundance_year0 + rel_freq.space)/2 , by =.(Taxon, plot, site_code)]

## filter to only the live (the dead cover will be 0, which inflates the 0, bc of how we computed stuff above)
hist(cover[live == 1,DI])
summary(cover[live == 1,DI])




####### Make Categorical Variables of Dominance too ?########

# then given them a ranking into different categories? 
# to look at changes in those types of species and the impact on productivity?
## sub-ordinate and rare -- how can I compute this? ###########
cover[, DIgroup:=cut(DI, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group 
#*MAKE SURE ONLY TO DO FOR LIVE 
cover[, sr_domspp := length(unique(Taxon[DIgroup == "Dominant"])), by = .(plot, site_code, year)]
cover[, sr_rarespp := length(unique(Taxon[DIgroup == "Rare"])), by = .(plot, site_code, year)]
cover[, sr_subordspp := length(unique(Taxon[DIgroup == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$sr_domspp) # max is 2 species per plot that are dominant (makes sense!)
summary(cover$sr_subordspp) 
summary(cover$sr_rare) 

# compute change in richness in each group 
cover[order(year), change_sr_domspp := sr_domspp -shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), change_sr_rarespp := sr_rarespp -shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), change_sr_subordspp := sr_subordspp -shift(sr_subordspp), by =.(plot, site_code)]

summary(cover$change_sr_domspp) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -1.0000  0.0000  0.0000  0.0012  0.0000  1.0000    3005 
summary(cover$change_sr_subordspp) 
# > summary(cover$change_sr_subordspp) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0327   0.0000  10.0000     3005 

summary(cover$change_sr_rare) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# -13.0000   0.0000   0.0000  -0.0199   0.0000   7.0000     3005 





###############################################################################################
##### This computes DI over time, which is not what we want. ##############################
#### However, I am leaving this, in case we want to use any of the relative frequency over time measures
## Particularly to address questions like : ## For analysis-
# for rare species -- see if their DI score is more driven by the relative abundance or the low frequency ?#
# those might be different types of rare spp so it seems like we should not obscure them into one ## ] 
########################################################################################################################

#I did the following based on plot-years as sampling unit, but we decided it should be defined spatially.

## compute average relative abundance (averaged across years.)
# Relative abundance = abundance of a species in a sampling unit / total abundance of all species in a sampling unit
# here we need to average across the years.

# for TOTAL cover
cover[, ave_rel_abundance_over_time := ave(relative_sp_cover.yr), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time)
summary(cover$ave_rel_abundance_over_time) 

#for LIVE cover
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

## Compute Relative Frequency per spp per plot.
# "Relative frequency = number of sampling units a species occurred / total number of sampling units" 
# first calculate the total number of years for a given plot in a site "total # of sampling units"
cover[, tot.yr := length(unique(year)), by =.(plot, site_code)]
# to determine which is the most recent year: 
# cover[, max.yr := max(year), by =.(plot, site_code)]

# spot check
cover[site_code == "abisko.se", .(year, plot, tot.yr),]
cover[site_code == "barta.us", .(year, plot, tot.yr),]
cover[site_code == "smith.us", .(year, plot, tot.yr),]

#then calculate number of sampling units (years) a species occurred
# with a list of taxa per plot and year, count the number of records per (species, plot) pair. 
# so for ex., if the datatable is called dt  
# dt[,.N, by=.(plot, species)] 
#will count the records per plot x species and if they only record a species in a plot and year if it's there
# cover[, .N , by =.(plot, site_code, Taxon)]

# to store this count as part of the cover datatable:
cover[, Num_year_spp_occur := length(unique(year)), by = .(plot, site_code, Taxon)]
cover[, Num_recs_spp_plot_site := .N, by = .(plot, site_code, Taxon)]
#***send this to Ashley to find duplicates
# cover[Num_year_spp_occur!=Num_recs_spp_plot_site, ] #812 records with duplicates? 
#585 from site_name=="Antisana" & year==2013:
# cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013,]
# table(cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013, .(site_code, plot, Taxon)])
# cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013 & Taxon=="WERNERIA PUMILA", ]
# looks like this is repeated bc multiple subplots
#how many sites have multiple subplots and taxa recorded by it???? anything that would be recorded at the subplot level?
# the max cover for cover[Num_year_spp_occur!=Num_recs_spp_plot_site & site_name=="Antisana" & year==2013 & Taxon=="WERNERIA PUMILA", ] 
# is the same
cover[,Num_subplots:=length(unique(subplot)),.(site_code,plot, Taxon)]
cover[Num_subplots>1,length(unique(site_code))]
unique(cover[Num_subplots>1,.(site_code, year)])
#  site_code year
#    anti.ec 2013
#    elkh.us 2007
#    sgs.us  2007
#    ucsc.us 2007

#compute relative frequency
# " Relative frequency = number of sampling units a species occurred / total number of sampling units" 
cover[, rel_freq.time :=  Num_year_spp_occur/tot.yr, by = .(plot, site_code)]
#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq)
summary(cover$rel_freq)

## Compute the DI per species <-- at what scale does this apply? the plot over all years? 
# DI = (average relative abundance + relative frequency)/2
#FOR LIVE COVER:
cover[, DI.time := (ave_rel_abundance_over_time.live + rel_freq)/2 , by =.(Taxon, plot, site_code)]
hist(cover$DI.time)
summary(cover$DI.time)

## For analysis-
# for rare species -- see if their DI score is more driven by the relative abundance or the low frequency ?#
# those might be different types of rare spp so it seems like we should not obscure them into one ## ] 


