##################################################################
### Process Nut Net Cover Data - for Full Dataset ############
## Laura Dee for the LTER non-random loss paper############
# April 7 2020  ###########################

# some updates on April 23, 2019 to visualize different abundance metrics #
# updates on Aug 28 2019 to include other ways of processing the rare, dom, etc variables (based on frequency and relative cover)
# also other groups: rare & non-rare, and rare-non-native, rare-native, non-native/non-rare
#Updates on Dec 28 2019 to include other cut offs for rare groups
##########################################

#Close graphics and clear local memory
graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(foreign)
library(rmarkdown)

# data.table cheat sheets
# cover[, N_fixer_cover.yr := sum(relative_sp_cover.yr, by )]
# comb[,sum(trt=="Control"), by=.(newblockid, year)][V1>1,]

setwd("~/Dropbox/IV in ecology/NutNet")
cover <- fread('full-cover-09-April-2018.csv',na.strings= 'NA')

##  make max_cover NOT a character
cover$max_cover <- as.numeric(cover$max_cover)

#########################################################################################
## Compare Taxon in Live (live==1) or non-live (live == 0) ################################
##########################################################################################
## determine how many Taxon are listed as live==0 ##
a <- table(cover$Taxon, cover$live==0)
head(a)
# write.csv(a, "live_v_dead_spplist_NutNet.csv")

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

# # make a total cover in a plot, site, year
# cover[,totplotcover.yr := sum(max_cover, na.rm=T), by=.(plot, site_code, year)]
# 
# #Make a relative cover for each species in each plot, site, year
# # based on TOTAL cover (including dead cover)
# cover[,relative_sp_cover.yr := max_cover/totplotcover.yr]


###### Filter cover dataset to only the live cover before computing metrics ################# 
# as discussed with Kim Komatsu 

cover <- cover[live == 1,]

#Confirm this worked for max_cover:
hist(cover$max_cover)
hist(cover[live == 1,max_cover])
summary(cover$max_cover)
summary(cover[live == 1,max_cover])


##### TOTAL *LIVE* COVER MEASURES ########################################
#  make a total live cover in a plot, site, year
cover[,totplotcover.yr.live := sum(max_cover[live=="1"], na.rm=T), by=.(plot, site_code, year)]

#  make a *relative* live cover in a plot, site, year
cover[,relative_sp_cover.yr.live := (max_cover*(live=="1"))/totplotcover.yr.live]

# in some cases, no live species in plot and year, so getting NA since totplotcover.yr.live is zero. Set these to zero.
cover[is.na(relative_sp_cover.yr.live),relative_sp_cover.yr.live:=0]

## Compute Native species cover 
cover[, Native_cover.yr := sum(relative_sp_cover.yr.live[local_provenance=="NAT"]), by = .(plot, site_code, year)]

######################################################################################################### #################
## Make Dominant, Subordinate and Rarity Variables based on relative abundance - CUT OFF 1 breaks =  c(0.0,0.2,0.8,1.0) ###
##########################################################################################################################

#for LIVE cover - get relative cover per year per plot and species - which we can use later to see how this variable for the groups is changing
#through time.
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

# we will want to make our categorical groupings based on the pre-treament year, thus:
# Compute Relative abundance per species and plot in pre-treatment year (year_trt == 0)
cover[, rel_abundance_year0 := relative_sp_cover.yr.live[year_trt == 0], by = .(Taxon, plot, site_code)]
hist(cover$rel_abundance_year0, xlab = "Average relative abundance at a site", main ="")
summary(cover$rel_abundance_year0)

######################################################################################################### #############
## Compute for Average Relative Abundance: sub-ordinate, dominant and rare - Cut off1: breaks =  c(0.0,0.2,0.8,1.0)
######################################################################################################### #############

cover[, RelAbund_group:=cut(rel_abundance_year0, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, relabund_sr_domspp := length(unique(Taxon[RelAbund_group == "Dominant"])), by = .(plot, site_code, year)]
cover[, relabund_sr_rarespp := length(unique(Taxon[RelAbund_group == "Rare"])), by = .(plot, site_code, year)]
cover[, relabund_sr_subordspp := length(unique(Taxon[RelAbund_group == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$relabund_sr_domspp) # max 2 
summary(cover$relabund_sr_rarespp) # max 31 (mean ~9)
summary(cover$relabund_sr_subordspp) #max 3, min 0 ?? mean 2

# create a non-rare variable for relative abundance
cover[, non_rare_spp.RelA := RelAbund_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.RelA := length(unique(Taxon[non_rare_spp.RelA == "TRUE"])), by = .(plot, site_code, year)]

#rare species measures
cover[, sr_rare_spp.RelA := length(unique(Taxon[non_rare_spp.RelA == "FALSE"])), by = .(plot, site_code, year)]
# compute the rare native vs non group!
cover[, sr_rare_non.nat.RelA:= length(unique(Taxon[non_rare_spp.RelA == "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.RelA := length(unique(Taxon[non_rare_spp.RelA == "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.RelA:= length(unique(Taxon[non_rare_spp.RelA == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.RelA := length(unique(Taxon[non_rare_spp.RelA == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.RelA := length(unique(Taxon[RelAbund_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.RelA := length(unique(Taxon[RelAbund_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.RelA := length(unique(Taxon[RelAbund_group == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.RelA := length(unique(Taxon[RelAbund_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]


######################################################################################################################
## Compute for Average Relative Abundance: or Rare, Dom Subord Group Cut Off 2-- breaks=c(0.0,0.4,0.8,1.0)
# this expands the definition of rare to include more species (bottom 40% in terms of relative abundance)
######################################################################################################################

cover[, RelAbund_group2 :=cut(rel_abundance_year0, breaks=c(0.0,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, relabund_sr_domspp2 := length(unique(Taxon[RelAbund_group2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, relabund_sr_rarespp2 := length(unique(Taxon[RelAbund_group2 == "Rare"])), by = .(plot, site_code, year)]
cover[, relabund_sr_subordspp2 := length(unique(Taxon[RelAbund_group2 == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$relabund_sr_domspp2) # max 2 
summary(cover$relabund_sr_rarespp2) # max 43
summary(cover$relabund_sr_subordspp2) #max 3, min 0 ?? mean 2 

# create a non-rare variable for relative abundance
cover[, non_rare_spp.RelA2 := RelAbund_group2 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "TRUE"])), by = .(plot, site_code, year)]

#rare species measures
cover[, sr_rare_spp.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "FALSE"])), by = .(plot, site_code, year)]
# compute the rare native vs non group!
cover[, sr_rare_non.nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.RelA2 := length(unique(Taxon[non_rare_spp.RelA2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.RelA2 := length(unique(Taxon[RelAbund_group2 == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.RelA2 := length(unique(Taxon[RelAbund_group2 == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.RelA2 := length(unique(Taxon[RelAbund_group2 == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.RelA2 := length(unique(Taxon[RelAbund_group2 == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]


####################################################################################################################
### Dominance Variables using DI ############################################################################################
####################################################################################################################
### do this for live cover only (discussed w Kim in June)

### Do relative abundance of live cover only; and SR of live cover only.

### Dominance indicator calculation ###
# Use the dominance indicator metric from Avoilio et al (seperates dominance indication from impact)
# DI = (average relative abundance + relative frequency)/2
#Relative abundance = abundance of a species a in sampling unit / total abundance of all species in a sampling unit
#Relative frequency = number of sampling units a species occurred / total number of sampling units
# note: ranges from 0-1; Relative abundance can be any measure of abundance. Does not incorporate a measure of impact.
# There is not a cutoff for "which range from 0-1 = dominant species versus subordinate or rare, in the Avolio et al paper # 
# so if we want to group species in these groups, we will need to make one (then also should test robustness to that decision)

#for LIVE cover
cover[, ave_rel_abundance_over_time.live := ave(relative_sp_cover.yr.live), by = .(Taxon, plot, site_code)]
hist(cover$ave_rel_abundance_over_time.live)
summary(cover$ave_rel_abundance_over_time.live) 

# relative abundance per species and plot in pre-treatment year (year_trt == 0)
cover[, rel_abundance_year0 := relative_sp_cover.yr.live[year_trt == 0], by = .(Taxon, plot, site_code)]
hist(cover$rel_abundance_year0, xlab = "Average relative abundance at a site", main ="")
summary(cover$rel_abundance_year0)

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


#Compute Relative Frequency
# " Relative frequency = number of sampling units a species occurred / total number of sampling units" 
cover[, rel_freq.space :=  tot.num.plots.with.spp/tot.num.plots, by = .(plot, site_code)]
#check to make sure we took out duplicates, max should be 1
hist(cover$rel_freq.space, xlab = "Frequency at the site in pre-treatment year", main = "Frequency of occurrence")
summary(cover$rel_freq.space)


## Compute the DI per species defined as:  DI = (average relative abundance + relative frequency)/2

#FOR LIVE COVER -- the dead cover will be 0.
cover[, DI := (rel_abundance_year0 + rel_freq.space)/2 , by =.(Taxon, plot, site_code)]
# summary(cover$DI)
hist(cover$DI)

## filter to only the live (the dead cover will be 0, which inflates the 0, bc of how we computed stuff above)
hist(cover[live == 1,DI], xlab = "Dominance indicator (DI)", main = "Dominance indicator metric")
# summary(cover[live == 1,DI])

####### ####### ####### ####### ####### ####### ####### ####### # ####### ####### ####### # ####### ####### ####### 
####### Make Categorical Variables (Groups) of Rarity, Subordinate, Dominant Based on DI ########
####### ####### ####### ####### ####### ####### ####### ####### # ####### ####### ####### # ####### ####### ####### 

# then given them a ranking into different categories? 
# to look at changes in those types of species and the impact on productivity?
## sub-ordinate and rare -- how can I compute this? ###########
cover[, DIgroup:=cut(DI, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group 
#*MAKE SURE ONLY TO DO FOR LIVE 
cover[, sr_domspp := length(unique(Taxon[DIgroup == "Dominant"])), by = .(plot, site_code, year)]
cover[, sr_rarespp := length(unique(Taxon[DIgroup == "Rare"])), by = .(plot, site_code, year)]
cover[, sr_subordspp := length(unique(Taxon[DIgroup == "Subordinate"])), by = .(plot, site_code, year)]

### create a category that groups together all the non-rare
cover[, sr_non_rare_spp := length(unique(Taxon[DIgroup %in% c("Subordinate", "Dominant")])), by = .(plot, site_code, year)]

summary(cover$sr_domspp) # max is 2 species that are dominant 
summary(cover$sr_subordspp) 
summary(cover$sr_rare) 
summary(cover$sr_non_rare_spp)

# compute change in richness in each group 
cover[order(year), change_sr_domspp := sr_domspp -shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), change_sr_rarespp := sr_rarespp -shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), change_sr_subordspp := sr_subordspp -shift(sr_subordspp), by =.(plot, site_code)]
cover[order(year), change_sr_non_rare_spp := sr_non_rare_spp -shift(sr_non_rare_spp), by =.(plot, site_code)]

# lagged effects of sr in these groups
cover[order(year), lagged_sr_rarespp := shift(sr_rarespp), by =.(plot, site_code)]
cover[order(year), lagged_sr_domspp := shift(sr_domspp), by =.(plot, site_code)]
cover[order(year), lagged_sr_subordspp := shift(sr_subordspp), by =.(plot, site_code)]

plot(ave_rel_abundance_over_time.live ~ DI, data= cover)
abline(lm(ave_rel_abundance_over_time.live ~ DI, data= cover), col = "yellow")
plot(relative_sp_cover.yr.live ~ DI, data = cover)

## Make a variable that captures changes in *cover* of species classified as dominance, rare, subordinate each year 
# (so we can look at variations from year to year)
cover[, Dom_cover.yr := sum(relative_sp_cover.yr[DIgroup=="Dominant"]), by = .(plot, site_code, year)]
cover[, Subord_cover.yr := sum(relative_sp_cover.yr[DIgroup=="Subordinate"]), by = .(plot, site_code, year)]
cover[, Rare_cover.yr := sum(relative_sp_cover.yr[DIgroup=="Rare"]), by = .(plot, site_code, year)]

hist(cover$Dom_cover.yr)
hist(cover$Rare_cover.yr)
hist(cover$Subord_cover.yr)
plot(cover$Rare_cover.yr, cover$sr_rarespp)
plot(cover$Dom_cover.yr, cover$sr_rarespp)

# make a change in the relative_sp_cover.yr.live and then look at relationship with DI?
cover[order(year), change_rel_abundance_time.live := relative_sp_cover.yr.live - shift(relative_sp_cover.yr.live), by =.(plot, site_code)]

#########################################################################################################################
### Create variables that are combined groups of: #######################################################################
# non-native richness (dominant and rare) + native rare + native non-rare.  ################################################
##############################################################################################################################
# create combinations of all - as factor in a column
cover[, status.NN.RareDom := paste(DIgroup,local_provenance, sep = "_")]
# table(cover$status.NN.RareDom)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### SR variables by combined groupings ####### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#do SR for non-native, rare:
cover[, sr_non.nat_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover[, sr_nat_rare := length(unique(Taxon[DIgroup == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

## look a the data 
hist(cover$sr_non.nat_rare)
hist(cover$sr_nat_rare)

# create a non-rare variable
cover[, non_rare_spp := DIgroup %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat := length(unique(Taxon[non_rare_spp == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom := length(unique(Taxon[DIgroup == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom := length(unique(Taxon[DIgroup == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub := length(unique(Taxon[DIgroup == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub := length(unique(Taxon[DIgroup == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do for native, non-rare:
#cover[, sr_nat_non.rare := length(unique(Taxon[DIgroup %in% c("Subordinate", "Dominant") & local_provenance == "NAT"])), by = .(plot, site_code, year)]
# do for non-native, non-rare:
# cover[, sr_non.nat_non.rare := length(unique(Taxon[DIgroup %in%  c("Subordinate", "Dominant") & local_provenance == "INT"])), by = .(plot, site_code, year)]

### cover variables by each group and the combined groupings ###
#using relative_sp_cover.yr.live based on year 0 groups

#Compute overall sum of relative cover by groups over time variables
cover[, cover_tot_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" ]), by = .(plot, site_code, year)]
cover[, cover_tot_rare := sum(relative_sp_cover.yr.live[non_rare_spp == "FALSE" ]), by = .(plot, site_code, year)]
cover[, cover_tot_dom := sum(relative_sp_cover.yr.live[DIgroup == "Dominant" ]), by = .(plot, site_code, year)]
cover[, cover_tot_sub := sum(relative_sp_cover.yr.live[DIgroup == "Subordinate"]), by = .(plot, site_code, year)]
cover[, cover_tot_NAT := sum(relative_sp_cover.yr.live[local_provenance == "NAT"]), by = .(plot, site_code, year)]
cover[, cover_tot_INT := sum(relative_sp_cover.yr.live[local_provenance == "INT"]), by = .(plot, site_code, year)]

#for NON rare (grouping dominant and subordinate)
cover[, cover_non.nat_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" & local_provenance == "INT"]), by = .(plot, site_code, year)]
cover[, cover_nat_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" & local_provenance == "NAT"]), by = .(plot, site_code, year)]

#for rare
cover[, cover_nat_rare := sum(relative_sp_cover.yr.live[non_rare_spp == "FALSE" & local_provenance == "NAT"]), by = .(plot, site_code, year)]
cover[, cover_non.nat_rare := sum(relative_sp_cover.yr.live[non_rare_spp == "FALSE" & local_provenance == "INT"]), by = .(plot, site_code, year)]

# for dominant
cover[, cover_non.nat_dom := sum(relative_sp_cover.yr.live[DIgroup == "Dominant" & local_provenance == "INT"]), by = .(plot, site_code, year)]
cover[, cover_nat_dom := sum(relative_sp_cover.yr.live[DIgroup == "Dominant" & local_provenance == "NAT"]), by = .(plot, site_code, year)]

#native dominant cover is declining through time on average
trend1 <- lm(cover_nat_dom  ~ site_code + as.numeric(year_trt) , data = cover)
summary(trend1)

# for subordinate
cover[, cover_non.nat_sub := sum(relative_sp_cover.yr.live[DIgroup == "Subordinate" & local_provenance == "INT"]), by = .(plot, site_code, year)]
cover[, cover_nat_sub := sum(relative_sp_cover.yr.live[DIgroup == "Subordinate"& local_provenance == "NAT"]), by = .(plot, site_code, year)]

#cover[, cover_nat_non.rare := sum(relative_sp_cover.yr.live[non_rare_spp == "TRUE" & local_provenance == "NAT"]), by = .(plot, site_code, year)]

# plot cover
par(mfrow = c(2,2), pty = "s")
hist(cover$cover_non.nat_non.rare, main = "cover non-native non-rare", xlab = "cover non-native non-rare")
hist(cover$cover_nat_non.rare, main = "cover native non-rare", xlab = "cover native non-rare")
hist(cover$cover_nat_rare, main = "cover native rare")
hist(cover$cover_non.nat_rare, main = "cover non-native rare")

#plot changes 
plot(cover$cover_non.nat_non.rare ~ cover$year)
abline(lm(cover$cover_non.nat_non.rare ~ length(unique(cover$year))), col = "red")

trend1 <- lm(cover_nat_dom  ~ site_code +  as.numeric(year_trt), data = cover)
summary(trend1)
trend2 <- lm(cover_nat_dom  ~ site_code:as.numeric(year_trt) + as.numeric(year_trt) , data = cover)
summary(trend2)

# non-native rare trends - increaisng
trend1 <- lm(cover_non.nat_dom  ~ site_code +  as.numeric(year_trt), data = cover)

#non native dom cover also decreasing
trend1 <- lm(cover_non.nat_dom  ~ site_code +  as.numeric(year_trt), data = cover)

#EXAMPLE: explo[, pos_dev_thresh := dev_from_thresh[dev_from_thresh>1], by = .(site, year)]
#cover[, Fenced := trt[trt == "Fenced"], by = .(plot, site_code, year)]

#non native dom cover also decreasing
trend1 <- lm(cover_non.nat_rare  ~ site_code +  as.numeric(year_trt), data = cover)
summary(trend1)

# Dec 29 2019

## Do a different cut off for dominance, subordinate and rare:
cover[, DIgroup2:=cut(DI, breaks=c(0.0,0.2,0.7,1.0), labels=c("Rare","Subordinate","Dominant"))]
#richness in each group with different cutoff
cover[, sr_domspp2 := length(unique(Taxon[DIgroup2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, sr_rarespp2 := length(unique(Taxon[DIgroup2 == "Rare"])), by = .(plot, site_code, year)]
cover[, sr_subordspp2 := length(unique(Taxon[DIgroup2 == "Subordinate"])), by = .(plot, site_code, year)]
summary(cover$sr_domspp2)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### SR variables by combined groupings ####### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#do SR for non-native, rare:
cover[, sr_non.nat_rare2 := length(unique(Taxon[DIgroup2 == "Rare" & local_provenance == "INT"])), by = .(plot, site_code, year)]
#do SR for native, rare:
cover[, sr_nat_rare2 := length(unique(Taxon[DIgroup2 == "Rare" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# create a non-rare variable
cover[, non_rare_spp2 := DIgroup2 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat2 := length(unique(Taxon[non_rare_spp2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat2 := length(unique(Taxon[non_rare_spp2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom2 := length(unique(Taxon[DIgroup2 == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom2 := length(unique(Taxon[DIgroup2 == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub2 := length(unique(Taxon[DIgroup2 == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub2 := length(unique(Taxon[DIgroup2 == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

# # compute change in richness in each group for this diff cutoff (though only the dominance one should change)
# cover[order(year), change_sr_domspp2 := sr_domspp2 -shift(sr_domspp2), by =.(plot, site_code)]
# cover[order(year), change_sr_rarespp2 := sr_rarespp2 -shift(sr_rarespp2), by =.(plot, site_code)]
# cover[order(year), change_sr_subordspp2 := sr_subordspp2 -shift(sr_subordspp2), by =.(plot, site_code)]


###########################################################################################################################
######################################################################################################### #############
## Dominance and Rarity Variables based on relative abundance and frequency groups seperately (versus DI above) #############
######################################################################################################################

######################################################################################################### #############
## Compute for Frequency: sub-ordinate and rare - Cut Off 1 based on breaks=c(0.0,0.2,0.8,1.0)
######################################################################################################### #############

cover[, Freq_group:=cut(rel_freq.space, breaks=c(0.0,0.2,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group  #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, freq_sr_domspp := length(unique(Taxon[Freq_group == "Dominant"])), by = .(plot, site_code, year)]
cover[, freq_sr_rarespp := length(unique(Taxon[Freq_group == "Rare"])), by = .(plot, site_code, year)]
cover[, freq_sr_subordspp := length(unique(Taxon[Freq_group == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$freq_sr_domspp)  #max 23
summary(cover$freq_sr_rarespp) #max 13
summary(cover$freq_sr_subordspp) #max 24

# create a rare variable for frequency
cover[, sr_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "FALSE"])), by = .(plot, site_code, year)]
cover[, sr_rare_non.nat.Freq:= length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.Freq := length(unique(Taxon[non_rare_spp.Freq== "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# create a non-rare variable for frequency
cover[, non_rare_spp.Freq := Freq_group %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.Freq := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.Freq  := length(unique(Taxon[non_rare_spp.Freq == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.Freq := length(unique(Taxon[Freq_group == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.Freq := length(unique(Taxon[Freq_group == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]

######################################################################################################################
## Compute for Frequency: sub-ordinate and rare - Cut Off 1 based on breaks=c(0.0,0.4,0.8,1.0).
######################################################################################################################
cover[, Freq_group2:=cut(rel_freq.space, breaks=c(0.0,0.4,0.8,1.0), labels=c("Rare","Subordinate","Dominant"))]

#richness in each group  #*MAKE SURE ONLY TO DO FOR LIVE 
cover[, freq2_sr_domspp := length(unique(Taxon[Freq_group2 == "Dominant"])), by = .(plot, site_code, year)]
cover[, freq2_sr_rarespp := length(unique(Taxon[Freq_group2 == "Rare"])), by = .(plot, site_code, year)]
cover[, freq2_sr_subordspp := length(unique(Taxon[Freq_group2 == "Subordinate"])), by = .(plot, site_code, year)]

summary(cover$freq2_sr_domspp)  #previous was max 23; now: 23
summary(cover$freq2_sr_rarespp) #max 13; 19
summary(cover$freq2_sr_subordspp) #max 24; 18

# create a non-rare variable for frequency
cover[, non_rare_spp.Freq2 := Freq_group2 %in% c("Subordinate", "Dominant"), by = .(plot, site_code, year)]
cover[, sr_non_rare_spp.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "TRUE"])), by = .(plot, site_code, year)]

# non-rare native and non-native
cover[, sr_non.rare_non.nat.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "TRUE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_non.rare_nat.Freq2  := length(unique(Taxon[non_rare_spp.Freq2 == "TRUE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

# create a rare variable for frequency
cover[, sr_rare_spp.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "FALSE"])), by = .(plot, site_code, year)]
cover[, sr_rare_non.nat.Freq2:= length(unique(Taxon[non_rare_spp.Freq2 == "FALSE" & local_provenance == "INT"])), by = .(plot, site_code, year)]
cover[, sr_rare_nat.Freq2 := length(unique(Taxon[non_rare_spp.Freq2 == "FALSE" & local_provenance == "NAT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for dom
cover[, sr_nat_dom.Freq2 := length(unique(Taxon[Freq_group2 == "Dominant" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_dom.Freq2 := length(unique(Taxon[Freq_group2 == "Dominant" & local_provenance == "INT"])), by = .(plot, site_code, year)]

#do SR for native and non-native for subordinate
cover[, sr_nat_sub.Freq2 := length(unique(Taxon[Freq_group2 == "Subordinate" & local_provenance == "NAT"])), by = .(plot, site_code, year)]
cover[, sr_non.nat_sub.Freq2 := length(unique(Taxon[Freq_group2 == "Subordinate" & local_provenance == "INT"])), by = .(plot, site_code, year)]


##########################################################################################################
### Annual vs Perennial Variables ########################################################################
###########################################################################################################
### Now do for an annual and perennial cover
# Compute richness in annuals and perennials
cover[, sr_annual := length(unique(Taxon[local_lifespan=="ANNUAL"])), by = .(plot, site_code, year)]
cover[, sr_peren := length(unique(Taxon[local_lifespan=="PERENNIAL"])), by = .(plot, site_code, year)]
cover[, sr_null.lspan := length(unique(Taxon[local_lifespan=="NULL"])), by = .(plot, site_code, year)]  # this includes a lot non-live. 
cover[, sr_biennial := length(unique(Taxon[local_lifespan=="BIENNIAL"])), by = .(plot, site_code, year)]
cover[, sr_indeter := length(unique(Taxon[local_lifespan=="INDETERMINATE"])), by = .(plot, site_code, year)]

# Compute first differences 
cover[order(year), change_sr_annual := sr_annual-shift(sr_annual), by =.(plot, site_code)]
cover[order(year), change_sr_peren := sr_peren-shift(sr_peren), by =.(plot, site_code)]

# Plot first differences
hist(cover$change_sr_annual)
summary(cover$change_sr_annual)
hist(cover$change_sr_peren)
summary(cover$change_sr_peren)

# compute cover by annuals per plot and year
cover[, AnnualPercentcover.yr := sum(relative_sp_cover.yr[local_lifespan=="ANNUAL"]), by = .(plot, site_code, year)]
# compute cover by perennial per plot and year
cover[, PerenPercentcover.yr := sum(relative_sp_cover.yr[local_lifespan=="PERENNIAL"]), by = .(plot, site_code, year)]

hist(cover$PerenPercentcover.yr)
hist(cover$AnnualPercentcover.yr)

########################################################################################################
### Check for duplicates and write out file #############################################################
##########################################################################################################
# remove mistake/duplicate records from comp.pt as above 
#****** & likely need to check for others since there could be more in the full dataset!!!!
cover = cover[!(site_code == "comp.pt" & plot %in% c(5,19,34) & year %in% c(2013,2014,2015,2016) & year_trt==0),]

# write as csv datafile to use for R
write.csv(cover, "NutNet_FullCoverData_ProcessedApril2020.csv")


