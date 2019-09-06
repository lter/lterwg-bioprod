### NutNet and rare species #####
# June 13 2017

#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

setwd("~/Google Drive/LTER_Biodiversity_Productivity")

#kim's wd
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\nutrient network\\NutNet data')
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\nutrient network\\NutNet data')

library(grid)
library(lme4)
library(sjPlot)
library(merTools)
library(DHARMa)
library(car)
library(AICcmodavg)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  



#########
nutnetdf <-read.csv("full-cover-02-August-2019.csv")
# 
# nutnetpretreatdt <- data.table(nutnetpretreatdf)
# nutnetpretreatdt[, meanPTAbundance:=mean(rel_cover, na.rm=T), .(year, site_code, Taxon)]
# nutnetpretreatdt[, maxPTAbundance:=max(rel_cover, na.rm=T), .(year, site_code, Taxon)]
# nutnetpretreatdt[, PTfreq:=sum(live==1, na.rm=T), .(year, site_code, Taxon)]

#relative cover
nutnetSum <- nutnetdf%>%
  group_by(site_code, plot, year_trt)%>%
  summarise(total_cover=sum(max_cover))
nutnetRel <- nutnetdf%>%
  left_join(nutnetSum)%>%
  mutate(rel_cover=(max_cover/total_cover)*100)%>%
  filter(live==1)

## compute mean abundance, max abundance and frequency of each species in the pretreatment data 
#filter to pretreatment data 
nutnetpretreatdf <- nutnetRel[nutnetRel$year_trt == 0,]
meanAb_byspecies <- aggregate(nutnetpretreatdf$rel_cover, by = list(nutnetpretreatdf$year, nutnetpretreatdf$site_code, nutnetpretreatdf$Taxon), FUN = mean, na.rm = TRUE)
names(meanAb_byspecies) = c("year", "site_code", "Taxon", "meanPTAbundance")
meanAb_byspecies <- meanAb_byspecies%>%
  select(-year)
max_abund <-  aggregate(nutnetpretreatdf$rel_cover, by = list(nutnetpretreatdf$year, nutnetpretreatdf$site_code, nutnetpretreatdf$Taxon), FUN = max, na.rm = TRUE)
names(max_abund) = c("year", "site_code", "Taxon", "maxPTAbundance")
max_abund <- max_abund%>%
  select(-year)
#filter to species live in a plot and pretreatment year, then sum to get the # of plots the species appeared in
nutnetpretreatdf_live = nutnetpretreatdf[nutnetpretreatdf$live == 1,]
freq <- aggregate(nutnetpretreatdf_live$live, by = list(nutnetpretreatdf_live$year, nutnetpretreatdf_live$site_code, nutnetpretreatdf_live$Taxon), FUN = sum, na.rm = TRUE)
names(freq) = c("year", "site_code", "Taxon", "PTfreq")
tot_plots <- nutnetpretreatdf_live%>%
  select(site_code, plot)%>%
  unique()%>%
  group_by(site_code)%>%
  summarise(num_plots=max(plot))
freq <- freq%>%
  select(-year)%>%
  merge(tot_plots, by='site_code')%>%
  mutate(freq=PTfreq/num_plots)

#### Process data to create a max cover or 0 for each species, plot & year #### 
nutnetdf_allspp <- nutnetRel%>%
  select(-Family, -live:-total_cover)%>%
  group_by(site_name, site_code)%>%
  nest()%>%
  mutate(spread_df = purrr::map(data, ~spread(., key=Taxon, value=rel_cover, fill=0)%>%
                                  gather(key=Taxon, value=rel_cover,
                                         -year:-trt)))%>%
  unnest(spread_df)
nutnetdf_allspp <- as.data.frame(nutnetdf_allspp)


###make a column for presence absence of each species in a plot-year
nutnetdf_allspp$PA = ifelse(nutnetdf_allspp$rel_cover > 0, 1, 0)

nutnetdf_length <- as.data.frame(nutnetdf_allspp)%>%
  #make a column for max trt year
  group_by(site_code)%>%
  summarise(length=max(year_trt))

nutnetdf_allspp2 <- nutnetdf_allspp%>%
  merge(nutnetdf_length, by='site_code')%>%
  select(-year, -rel_cover)%>%
  filter(year_trt>0)%>%
  mutate(year_trt2=paste("yr", year_trt, sep=''))%>%
  select(-year_trt, -trt)%>%
  group_by(site_code, Taxon, site_name, block, plot, subplot, year_trt2, length)%>%
  summarise(PA2=mean(PA))%>%
  ungroup()%>%
  group_by(site_code, Taxon, plot, length)%>%
  summarise(PA3=sum(PA2))%>%
  ungroup()

nutnetpretrt <- nutnetpretreatdf_live%>%
  mutate(pretrt_cover=rel_cover)%>%
  select(site_code, plot, Taxon, pretrt_cover)

nutnetdf_allspp3 <- nutnetdf_allspp2%>%
  mutate(yrs_absent=(length-PA3)/length)%>%
  left_join(nutnetpretrt)%>%
  filter(!is.na(pretrt_cover))%>%
  merge(meanAb_byspecies, by=c('site_code', 'Taxon'))%>%
  merge(max_abund, by=c('site_code', 'Taxon'))%>%
  merge(freq, by=c('site_code', 'Taxon'))%>%
  mutate(abund_metric=((meanPTAbundance/100)+freq)/2)%>%
  select(-PTfreq)%>%
  filter(length>0)


#get back trt info
trt <- nutnetRel%>%
  select(year_trt, site_code, plot, trt)%>%
  filter(year_trt>0)%>%
  select(-year_trt)%>%
  unique()

nutnetdf_allspp3Trt <- nutnetdf_allspp3%>%
  merge(trt, by=c('site_code', 'plot'))


###count concecutive 0's: how long is a sp lost?
cumul_zeros <- function(x)  {
  x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the 
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}

nutnetOrder <- nutnetdf_allspp[with(nutnetdf_allspp, order(site_code, plot, Taxon, year)), ]

nutnetConsAbs <- nutnetOrder%>%
  filter(year_trt>0)%>%
  select(site_code, year_trt, plot, Taxon, trt, rel_cover, PA)%>%
  group_by(site_code, year_trt, plot, Taxon, trt, rel_cover, PA)%>%
  unique()%>%
  ungroup()%>%
  group_by(site_code, plot, Taxon)%>%
  nest()%>%
  mutate(mutate_df = purrr::map(data, ~mutate(., cons_abs=cumul_zeros(PA))))%>%
  unnest(mutate_df)

nutnetConsAbs2 <- nutnetConsAbs%>%
  group_by(site_code, plot, Taxon)%>%
  summarise(cons_abs_max=max(cons_abs))%>%
  ungroup()%>%
  merge(trt, by=c('site_code', 'plot'))%>%
  merge(meanAb_byspecies, by=c('site_code', 'Taxon'))%>%
  merge(max_abund, by=c('site_code', 'Taxon'))%>%
  merge(freq, by=c('site_code', 'Taxon'))%>%
  mutate(abund_metric=((meanPTAbundance/100)+PTfreq)/2)%>%
  merge(nutnetdf_length, by=c('site_code'))

# ggplot(nutnetConsAbs2, aes(x=abund_metric, y=cons_abs_max)) +
#   geom_point() +
#   xlab('Pre-Treatment Modified Importance Index') +
#   ylab('Consecutive Years Absent (max)') +
#   geom_smooth(method='loess') +
#   facet_wrap(~trt)

###pres/absence
nutnetPresAbs <- nutnetdf_allspp%>%
  merge(meanAb_byspecies, by=c('site_code', 'Taxon'))%>%
  merge(max_abund, by=c('site_code', 'Taxon'))%>%
  merge(freq, by=c('site_code', 'Taxon'))%>%
  mutate(abund_metric=((meanPTAbundance/100)+PTfreq)/2)%>%
  filter(year_trt>0)


###proportion of years absent
nutnetPropAbs <- nutnetPresAbs%>%
  group_by(site_code, Taxon, plot, trt, abund_metric, meanPTAbundance, PTfreq)%>%
  summarise(years_present=sum(PA))%>%
  ungroup()%>%
  left_join(nutnetdf_length)%>%
  mutate(prop_years_present=years_present/length)%>%
  mutate(prop_years_absent=1-prop_years_present)%>%
  filter(prop_years_present>0)

nutnetPropAbsCtl <- nutnetPropAbs%>%
  filter(trt=='Control')%>%
  rename(prop_years_absent_ctl=prop_years_absent)%>%
  select(site_code, Taxon, prop_years_absent_ctl)

nutnetPropAbsDiff <- nutnetPropAbs%>%
  filter(trt!='Control')%>%
  left_join(nutnetPropAbsCtl)%>%
  mutate(prop_years_absent_diff=prop_years_absent-prop_years_absent_ctl)


ggplot(nutnetPropAbsDiff, aes(x=abund_metric, y=prop_years_absent_diff)) +
  geom_point() +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Proportion Years Absent') +
  geom_smooth(method='loess') +
  facet_wrap(~trt)

n <- 100
newdata <- crossing(meanPTAbundance=seq(0,100,length.out=n), 
                    site_code = unique(nutnetPropAbs$site_code))
newdata_fixed <- data.frame(meanPTAbundance=seq(0,100,length.out=n), 
                            site_code = unique(nutnetPropAbs$site_code)[1])



#NPK plot model
modPropAbsentNPK <- glmer(prop_years_absent_diff ~ meanPTAbundance + (1 + as.factor(plot)|site_code),
                          data = subset(nutnetPropAbsDiff, trt=='NPK'))
pred_final_loss_fixed <- cbind(newdata_fixed, predictInterval(modPropAbsentNPK,
                                                              newdata=newdata_fixed,
                                                              which="fixed", type="probability",
                                                              include.resid.var = FALSE))
pred_final_loss_ranef <- cbind(newdata, predictInterval(modPropAbsentNPK,
                                                        newdata=newdata,
                                                        which="full", type="probability",
                                                        include.resid.var = F))

ggplot(subset(nutnetPropAbs, trt=='NPK'), aes(x=meanPTAbundance, y=prop_years_absent, group=as.character(site_code), color=PTfreq)) +
  geom_point() +
  geom_line(data = pred_final_loss_ranef,
            mapping=aes(y=fit, group=as.character(site_code)),
            lwd=0.4, alpha=0.4, color="lightgrey") +
  geom_line(data = pred_final_loss_fixed,
            mapping=aes(y=fit),
            color="black", lwd=1.3) +
  ylab("Proportion Years Absent") +
  xlab("Pre-Treatment Abundance")

#linear model
summary(modPropAbsentNPK <- glm(prop_years_absent_diff ~ meanPTAbundance,
                          data = subset(nutnetPropAbsDiff, trt=='NPK')))



#add in final year abundances to see if species do persist, how does their abundance change
nutnet_finalabund <- nutnetdf_allspp%>%
  filter(year_trt==9)%>%
  mutate(final_cover=rel_cover)%>%
  select(site_code, plot, trt, final_cover, Taxon)

nutnet_finalabund2 <- nutnetdf_allspp3Trt%>%
  group_by(site_code, plot, Taxon, length, PA3, yrs_absent, pretrt_cover, meanPTAbundance, maxPTAbundance, num_plots, freq, abund_metric, trt)%>%
  unique()%>%
  ungroup()%>%
  left_join(nutnet_finalabund)%>%
  filter(!is.na(final_cover))

ggplot(nutnet_finalabund2, aes(x=abund_metric, y=final_cover)) +
  geom_point(size=3) +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Final Year (9) Relative Abundance') + 
  facet_wrap(~trt)
#export at 1200x1200



###bring in biomass (function) component
biomass <- read.csv('full-biomass-22-February-2019.csv')%>%
  filter(live==1)%>%
  group_by(site_code, plot, subplot, year_trt)%>%
  summarise(anpp=sum(mass))%>%
  ungroup()

###BEF figure (traditional)
#notes, do we want to just include controls? and only 30 pretrt plots?

richness <- nutnetRel%>%
  filter(live==1, Family!='NULL')%>%
  group_by(year_trt, site_code, plot, trt)%>%
  summarise(richness=length(rel_cover))%>%
  ungroup()

biomassRichness <- biomass%>%
  left_join(richness)%>%
  filter(year_trt!=0)%>%
  filter(trt=='Control'|trt=='NPK'|trt=='Fence'|trt=='NPK+Fence')%>%
  mutate(logRich=log10(richness+1), logBio=log10(anpp+1))


bef_mod <- lmer(log10(anpp+1) ~ log10(richness+1)*year_trt + (log10(richness+1) | site_code), data=subset(biomassRichness, trt=='Control'))

#evaluate autocor
biomassRichnessModOut <- biomassRichness%>%
  ungroup()%>%
  filter(trt=='Control')%>%
  mutate(res_bef_mod = residuals(bef_mod))%>%
  group_by(site_code)%>%
  mutate(lag_res_bef_mod = lag(res_bef_mod))%>%
  ungroup()

qplot(lag_res_bef_mod, res_bef_mod, data = biomassRichnessModOut)

pred_bef_fix <- cbind(biomassRichnessModOut,
                      predictInterval(bef_mod, 
                                      newdata=biomassRichnessModOut%>%mutate(year_trt=1), 
                                      which="fixed"))
pred_bef_ranef <- cbind(biomassRichnessModOut,
                        predictInterval(bef_mod, newdata=biomassRichnessModOut, 
                                        which="full"))

BEFstrawmanFig <- ggplot(biomassRichness, aes(x=log10(richness+1), y=log10(anpp+1), color=as.character(site_code))) + 
  guides(color = "none") +
  geom_ribbon(pred_bef_fix, mapping=aes(ymin=lwr, ymax=upr), fill="lightgrey", color=NA, alpha=0.6) +
  geom_smooth(pred_bef_ranef, method='lm', mapping=aes(y=fit), alpha=0.7, lwd=1.3, se=F) +
  geom_smooth(pred_bef_fix, method='lm', mapping=aes(y=fit), color="black", lwd=1.5) +
  xlab("Log Richness + 1") + ylab("Log Biomass + 1")




###figures comparable to predicts database

# #biomass response by trt
# biomassResp <- biomass%>%
#   left_join(trt)%>%
#   filter(year_trt!=0)%>%
#   group_by(site_code, year_trt, trt)%>%
#   summarise(anpp_mean=mean(anpp))%>%
#   ungroup()%>%
#   spread(key=trt, value=anpp_mean)%>%
#   mutate(NPK_diff=(NPK-Control)/Control, Fence_diff=(Fence-Control)/Control, NPKfence_diff=(NPK+Fence-Control)/Control)%>%
#   select(site_code, year_trt, NPK_diff, Fence_diff, NPKfence_diff)%>%
#   na.omit()%>%
#   gather(key=trt, value=diff, NPK_diff:NPKfence_diff)
# 
# ggplot(data=barGraphStats(data=biomassResp, variable="diff", byFactorNames=c("year_trt", "trt")), aes(x=year_trt, y=mean, color=trt)) +
#   geom_point(size=5) +
#   stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
#   xlab('Treatment Year') + ylab('Biomass Difference (%)') +
#   geom_hline(yintercept=0)
# #export at 1200x800


###split species into groups within sites by dominance

#make a new dataframe with just the label
site_code=nutnetdf_allspp3Trt%>%
  select(site_code)%>%
  unique()

#makes an empty dataframe
nutnetSppGroups=data.frame(row.names=1) 


#loop through sites to get their species' dominance ranks (3 groups)
for(i in 1:length(site_code$site_code)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=nutnetdf_allspp3Trt[nutnetdf_allspp3Trt$site_code==as.character(site_code$site_code[i]),]%>%
    select(site_code, Taxon, abund_metric, meanPTAbundance)%>%
    unique()
  
  #group by Dominance Indicator index (DI)
  subset$DI_group <- ntile(subset$abund_metric, 3)
  
  #group by mean abundance in pretrt data
  subset$abund_group <- ntile(subset$meanPTAbundance, 3)
  
  #pasting dispersions into the dataframe made for this analysis
  nutnetSppGroups=rbind(subset, nutnetSppGroups)  
}



####do everything at site level (lost from entire trt, rather than just a plot)
nutnetAbsSite <- nutnetPresAbs%>%
  #get site level presence/absence by trt
  group_by(site_code, Taxon, year_trt, trt)%>%
  summarise(presence=sum(PA))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=richness)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_rich=(NPK-Control)/Control, Fence_rich=(Fence-Control)/Control, NPKfence_rich=(NPK+Fence-Control)/Control)%>%
  select(site_code, year_trt, DI, NPK_rich, Fence_rich, NPKfence_rich)%>%
  gather(key=trt, value=richness_diff, NPK_rich:NPKfence_rich)%>%
  na.omit()%>%
  filter(year_trt!=999999)

ggplot(data=barGraphStats(data=nutnetRichnessDI, variable="richness_diff", byFactorNames=c("year_trt", "trt", "DI")), aes(x=year_trt, y=mean, color=DI)) +
  mutate(PA=ifelse(presence>0, 1, 0))%>%
  #proportion of years absent for trts
  group_by(site_code, Taxon, trt)%>%
  summarise(years_present=sum(PA))%>%
  ungroup()%>%
  left_join(nutnetdf_length)%>%
  #calculate proportion of years absent
  mutate(prop_years_absent=1-(years_present/length))%>%
  #merge the species categories
  left_join(nutnetSppGroups)

#calculate differences -- this shows the propotional differences in the number of years a species is absent in the trt compared to years species is absent in the ctl plots (doesn't communicate loss from ctl to trt)
nutnetAbsSiteCtl <- nutnetAbsSite%>%
  filter(trt=='Control')%>%
  rename(prop_years_absent_ctl=prop_years_absent)%>%
  select(site_code, Taxon, prop_years_absent_ctl)
nutnetAbsSiteDiff <- nutnetAbsSite%>%
  filter(trt!='Control')%>%
  left_join(nutnetAbsSiteCtl)%>%
  mutate(prop_years_absent_diff=(prop_years_absent-prop_years_absent_ctl)/(prop_years_absent_ctl))%>%
  mutate(prop_years_absent_diff_corr=ifelse(is.nan(prop_years_absent_diff), 0, ifelse(is.infinite(prop_years_absent_diff), 1, prop_years_absent_diff)))%>%
  #drop spp that are always absent from controls, because those are gains not losses
  filter(prop_years_absent_ctl<1)%>%
  mutate(dom_group=ifelse(DI_group==1, 'rare', ifelse(DI_group==2, 'intermediate', 'common')))


ggplot(data=barGraphStats(data=subset(nutnetAbsSiteDiff, trt=='NPK'|trt=='N'|trt=='NP'), variable="prop_years_absent_diff_corr", byFactorNames=c("trt", "dom_group")), aes(x=trt, y=mean, color=as.factor(dom_group))) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment') + ylab('Difference in Proportion of Years Absent (%)') +
  geom_hline(yintercept=0)
#export at 800x800



#figure out when a spp is absent from trt but present in ctl
nutnetPASite <- nutnetPresAbs%>%
  #get site level presence/absence by trt
  group_by(site_code, Taxon, year_trt, trt)%>%
  summarise(presence=sum(PA))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=richness)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_rich=(NPK-Control)/Control, Fence_rich=(Fence-Control)/Control, NPKfence_rich=(NPK+Fence-Control)/Control)%>%
  select(site_code, year_trt, abund, NPK_rich, Fence_rich, NPKfence_rich)%>%
  gather(key=trt, value=richness_diff, NPK_rich:NPKfence_rich)%>%
  na.omit()%>%
  filter(year_trt!=999999)
  mutate(PA=ifelse(presence>0, 1, 0))

nutnetPACtl <- nutnetPASite%>%
  filter(trt=='Control')%>%
  rename(PA_ctl=PA)%>%
  select(site_code, Taxon, year_trt, PA_ctl)
nutnetPASiteTrt <- nutnetPASite%>%
  filter(trt!='Control')%>%
  left_join(nutnetPACtl)%>%
  #drop species that are never present in control in a year (because those are gains in trts)
  filter(PA_ctl>0)%>%
  #merge dominance categories
  left_join(nutnetSppGroups)

#proportion of years absent
nutnetAbsentSiteTrt <- nutnetPASiteTrt%>%
  group_by(site_code, Taxon, trt, abund_metric, DI_group)%>%
  summarise(present=sum(PA))%>%
  ungroup()%>%
  #merge in length of exp
  left_join(nutnetdf_length)%>%
  mutate(prop_years_absent=(length-present)/length)

ggplot(subset(nutnetAbsentSiteTrt, trt=='N'), aes(x=abund_metric, y=prop_years_absent, group=as.character(site_code))) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_smooth(method='loess', color='black', se=F, group=(1)) +
  ylab("Proportion Years Absent") +
  xlab("Dominance Indicator Index")

ggplot(subset(nutnetAbsentSiteTrt, trt=='NPK'), aes(x=abund_metric, y=prop_years_absent, group=as.character(site_code))) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_smooth(method='loess', color='black', se=F, group=(1)) +
  ylab("Proportion Years Absent") +
  xlab("Dominance Indicator Index")


###changes in mean abundances
#merge groupings with data
nutnetAbundDI <- nutnetPresAbs%>%
  select(site_code, Taxon, plot, year_trt, trt, rel_cover)%>%
  left_join(nutnetSppGroups)%>%
  mutate(DI=ifelse(DI_group==1, 'rare', ifelse(DI_group==2, 'int', 'common')))%>%
  filter(rel_cover>0)%>%
  group_by(site_code, year_trt, trt, DI)%>%
  summarise(mean_abund=mean(rel_cover))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=mean_abund)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_abund=(NPK-Control)/Control, Fence_abund=(Fence-Control)/Control, NPKfence_abund=(NPK+Fence-Control)/Control)%>%
  select(site_code, year_trt, DI, NPK_abund, Fence_abund, NPKfence_abund)%>%
  gather(key=trt, value=abund_diff, NPK_abund:NPKfence_abund)%>%
  na.omit()%>%
  filter(year_trt!=999999)


nutnetAbundAbund <- nutnetPresAbs%>%
  select(site_code, Taxon, plot, year_trt, trt, rel_cover)%>%
  left_join(nutnetSppGroups)%>%
  mutate(abund=ifelse(abund_group==1, 'rare', ifelse(abund_group==2, 'int', 'common')))%>%
  filter(rel_cover>0)%>%
  group_by(site_code, year_trt, trt, abund)%>%
  summarise(mean_abund=mean(rel_cover))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=mean_abund)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_abund=(NPK-Control)/Control, Fence_abund=(Fence-Control)/Control, NPKfence_abund=(NPK+Fence-Control)/Control)%>%
  select(site_code, year_trt, abund, NPK_abund, Fence_abund, NPKfence_abund)%>%
  gather(key=trt, value=abund_diff, NPK_abund:NPKfence_abund)%>%
  na.omit()%>%
  filter(year_trt!=999999)

#for each site, how many spp in each category lost?
nutnetLossSite <- nutnetPASiteTrt%>%
  filter(PA==0)%>%
  group_by(site_code, year_trt, trt, DI_group)%>%
  summarise(num_loss=length(PA))%>%
  ungroup()
nutnetNotLossSite <- nutnetPASiteTrt%>%
  filter(PA==1)%>%
  group_by(site_code, year_trt, trt, DI_group)%>%
  summarise(num_notloss=length(PA))%>%
  ungroup()
nutnetLossorNotSite <- nutnetLossSite%>%
  full_join(nutnetNotLossSite)%>%
  mutate(num_loss=replace_na(num_loss, 0))%>%
  mutate(num_notloss=replace_na(num_notloss, 0))%>%
  mutate(total_spp=num_loss+num_notloss)%>%
  mutate(prop_loss=num_loss/total_spp)%>%
  #remove sites with no DI groups (no pretrt data)
  filter(!is.na(DI_group))%>%
  #remove sites where trt year is not recorded
  filter(year_trt<11)

#across all years
ggplot(barGraphStats(data=subset(nutnetLossorNotSite, trt=='N'|trt=='NP'|trt=='NPK'), variable="prop_loss", byFactorNames=c("trt", "DI_group")), aes(x=as.factor(DI_group), y=mean, color=trt)) +
  geom_point(stat='identity', position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(width=0.7)) +
  scale_x_discrete(breaks=c("1", "2", "3"), labels=c("Rare", "Intermediate", "Dominant")) +
  xlab('') + ylab('Proportion of Species Lost')

#####FIX THIS -- put in model results instead of raw data, drop weird years beyond which most sites have data
#time series
badLossFig <- ggplot(barGraphStats(data=subset(nutnetLossorNotSite, trt=='N'|trt=='NP'|trt=='NPK'), variable="prop_loss", byFactorNames=c("trt", "DI_group", "year_trt")), aes(x=year_trt, y=mean, color=as.factor(DI_group))) +
  geom_point(stat='identity', position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(width=0.7)) +
  xlab('Year') + ylab('Proportion of Species Lost') +
  scale_color_discrete(breaks=c("1", "2", "3"), labels=c("Rare", "Intermediate", "Dominant")) +
  facet_wrap(~trt)

#time series with only 10 year datasets
nutnetLossorNotSite10yr <- nutnetLossorNotSite%>%
  left_join(nutnetdf_length)%>%
  filter(length==10)
ggplot(barGraphStats(data=subset(nutnetLossorNotSite10yr, trt=='N'|trt=='NP'|trt=='NPK'), variable="prop_loss", byFactorNames=c("trt", "DI_group", "year_trt")), aes(x=year_trt, y=mean, color=as.factor(DI_group))) +
  geom_point(stat='identity', position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(width=0.7)) +
  xlab('Year') + ylab('Proportion of Species Lost') +
  scale_color_discrete(breaks=c("1", "2", "3"), labels=c("Rare", "Intermediate", "Dominant")) +
  facet_wrap(~trt)
  





nutnetAbsSiteCtl <- nutnetAbsSite%>%
  filter(trt=='Control')%>%
  rename(prop_years_absent_ctl=prop_years_absent)%>%
  select(site_code, Taxon, prop_years_absent_ctl)
nutnetAbsSiteDiff <- nutnetAbsSite%>%
  filter(trt!='Control')%>%
  left_join(nutnetAbsSiteCtl)%>%
  mutate(prop_years_absent_diff=(prop_years_absent-prop_years_absent_ctl)/(prop_years_absent_ctl))%>%
  mutate(prop_years_absent_diff_corr=ifelse(is.nan(prop_years_absent_diff), 0, ifelse(is.infinite(prop_years_absent_diff), 1, prop_years_absent_diff)))%>%
  #drop spp that are always absent from controls, because those are gains not losses
  filter(prop_years_absent_ctl<1)%>%
  mutate(dom_group=ifelse(DI_group==1, 'rare', ifelse(DI_group==2, 'intermediate', 'common')))


ggplot(data=barGraphStats(data=subset(nutnetAbsSiteDiff, trt=='NPK'|trt=='N'|trt=='NP'), variable="prop_years_absent_diff_corr", byFactorNames=c("trt", "dom_group")), aes(x=trt, y=mean, color=as.factor(dom_group))) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment') + ylab('Difference in Proportion of Years Absent (%)') +
  geom_hline(yintercept=0)
#export at 800x800




















#########
#everything below is really messy and needs shortening
#########



###proportion of years absent
nutnetAbsentDI <- nutnetdf_allspp3Trt%>%
  select(site_code, Taxon, plot, trt, yrs_absent, length)%>%
  left_join(nutnetSppGroups)%>%
  mutate(DI=ifelse(DI_group==1, 'rare', ifelse(DI_group==2, 'int', 'common')))%>%
  mutate(prop_absent=yrs_absent/length)%>%
  filter(length>2)%>% #sets dataset length to be at least 3, so we don't have species 100% present or absent based on one year
  group_by(site_code, trt, DI)%>%
  summarise(prop_years_absent=mean(prop_absent))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=prop_years_absent)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_absent=(NPK-Control)/Control, Fence_absent=(Fence-Control)/Control, NPKfence_absent=(NPK+Fence-Control)/Control)%>%
  select(site_code, DI, NPK_absent, Fence_absent, NPKfence_absent)%>%
  gather(key=trt, value=absent_diff, NPK_absent:NPKfence_absent)%>%
  na.omit()%>%
  filter(absent_diff<1000)

ggplot(data=barGraphStats(data=nutnetAbsentDI, variable="absent_diff", byFactorNames=c("trt", "DI")), aes(x=trt, y=mean, color=DI)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment') + ylab('Difference in Proportion of Years Absent (%)') +
  geom_hline(yintercept=0)
#export at 800x800


nutnetAbsentAbund <- nutnetdf_allspp3Trt%>%
  select(site_code, Taxon, plot, trt, yrs_absent, length)%>%
  left_join(nutnetSppGroups)%>%
  mutate(abund=ifelse(abund_group==1, 'rare', ifelse(abund_group==2, 'int', 'common')))%>%
  mutate(prop_absent=yrs_absent/length)%>%
  filter(length>2)%>% #sets dataset length to be at least 3, so we don't have species 100% present or absent based on one year
  group_by(site_code, trt, abund)%>%
  summarise(prop_years_absent=mean(prop_absent))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=prop_years_absent)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_absent=(NPK-Control)/Control, Fence_absent=(Fence-Control)/Control, NPKfence_absent=(NPK+Fence-Control)/Control)%>%
  select(site_code, abund, NPK_absent, Fence_absent, NPKfence_absent)%>%
  gather(key=trt, value=absent_diff, NPK_absent:NPKfence_absent)%>%
  na.omit()%>%
  filter(absent_diff<1000)

ggplot(data=barGraphStats(data=nutnetAbsentAbund, variable="absent_diff", byFactorNames=c("trt", "abund")), aes(x=trt, y=mean, color=abund)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment') + ylab('Difference in Proportion of Years Absent (%)') +
  geom_hline(yintercept=0)
#export at 800x800




###cumulative number of years absent
nutnetConsAbsentDI <- nutnetConsAbs2%>%
  select(site_code, Taxon, plot, trt, cons_abs_max, length)%>%
  left_join(nutnetSppGroups)%>%
  mutate(DI=ifelse(DI_group==1, 'rare', ifelse(DI_group==2, 'int', 'common')))%>%
  mutate(cons_abs=cons_abs_max)%>%
  filter(length>2)%>% #sets dataset length to be at least 3, so we don't have species 100% present or absent based on one year
  group_by(site_code, trt, DI)%>%
  summarise(cons_abs_max=mean(cons_abs))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=cons_abs_max)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_absent=(NPK-Control)/Control, Fence_absent=(Fence-Control)/Control, NPKfence_absent=(NPK+Fence-Control)/Control)%>%
  select(site_code, DI, NPK_absent, Fence_absent, NPKfence_absent)%>%
  gather(key=trt, value=absent_diff, NPK_absent:NPKfence_absent)%>%
  na.omit()%>%
  filter(absent_diff<1000)

ggplot(data=barGraphStats(data=nutnetConsAbsentDI, variable="absent_diff", byFactorNames=c("trt", "DI")), aes(x=trt, y=mean, color=DI)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment') + ylab('Difference in Proportion of Years Absent (%)') +
  geom_hline(yintercept=0)
#export at 800x800


nutnetConsAbsentAbund <- nutnetConsAbs2%>%
  select(site_code, Taxon, plot, trt, cons_abs_max, length)%>%
  left_join(nutnetSppGroups)%>%
  mutate(abund=ifelse(abund_group==1, 'rare', ifelse(abund_group==2, 'int', 'common')))%>%
  mutate(cons_abs=cons_abs_max)%>%
  filter(length>2)%>% #sets dataset length to be at least 3, so we don't have species 100% present or absent based on one year
  group_by(site_code, trt, abund)%>%
  summarise(cons_abs_max=mean(cons_abs))%>%
  ungroup()%>%
  filter(trt!='NA')%>%
  spread(key=trt, value=cons_abs_max)%>%
  select(-N, -P, -K, -NP, -NK, -PK)%>%
  mutate(NPK_absent=(NPK-Control)/Control, Fence_absent=(Fence-Control)/Control, NPKfence_absent=(NPK+Fence-Control)/Control)%>%
  select(site_code, abund, NPK_absent, Fence_absent, NPKfence_absent)%>%
  gather(key=trt, value=absent_diff, NPK_absent:NPKfence_absent)%>%
  na.omit()%>%
  filter(absent_diff<1000)

ggplot(data=barGraphStats(data=nutnetConsAbsentAbund, variable="absent_diff", byFactorNames=c("trt", "abund")), aes(x=trt, y=mean, color=abund)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment') + ylab('Difference in Proportion of Years Absent (%)') +
  geom_hline(yintercept=0)
#export at 800x800











#spread(key=year_trt2, value=PA2)%>% ###I don't think we want to do this, because it is then unclear which are missing data and which are true 0s

#mutate(yr_present=(length-yr1-yr2-yr3-yr4-yr5-yr6-yr7-yr8-yr9)/length) ###problem: can't sum over cells with NA, but most experiments have short datasets or missing data

# library(data.table)
# nutnetdt <- data.table(nutnetRel)
# 
# # create a data.table of unique species x year combos by site
# site_spp_yr_dt <- nutnetdt[,expand.grid(unique(plot), unique(Taxon), unique(year)), by=site_code]
# names(site_spp_yr_dt) <- c("site_code", "plot", "Taxon", "year")
# filled_nutnetdt <- merge(nutnetdt, site_spp_yr_dt, by=c("site_code", "plot", "Taxon", "year"), all=T)
# filled_nutnetdt[,PA:=!is.na(rel_cover)]
# filled_nutnetdt[,tot_site_years:=length(unique(year)), by=site_code]
# filled_nutnetdt[,tot_spp_plot_years_present:=sum(PA), by=.(plot, site_code, Taxon)]
# filled_nutnetdt[,frac_absent:=(tot_site_years-tot_spp_plot_years_present)/tot_site_years]
# # make a variable that is just tot_site_years-tot_spp_plot_years_present
# head(filled_nutnetdt)
# summary(filled_nutnetdt$frac_absent)

## merge back into the full nutnet data (not only pretreatment year)
# nutnetdf_allspp  <- merge(nutnetdf_allspp ,meanAb_byspecies,  by = c("site_code", "Taxon"), all.x = T)
# nutnetdf_allspp  <- merge(max_abund, nutnetdf_allspp , by = c("site_code", "Taxon"), all = T)
# nutnetdf_allspp  <- merge(freq, nutnetdf_allspp ,  by = c("site_code", "Taxon"), all = T)


#nutnetdf_allspp2 <- reshape(nutnetdf_allspp, idvar=c('site_code', 'Taxon', 'PTfreq', 'maxPTAbundance', 'meanPTAbundance', 'site_name', 'trt'), timevar='year_trt', direction='wide')

# ### Code from Kim using Codyn 
# 
# library(codyn)
# 
# # codyn function modification ---------------------------------------------
# #modifying codyn functions to output integer numbers of species appearing and disappearing, plus total spp number over two year periods, rather than ratios
# turnover_allyears <- function(df, 
#                               time.var, 
#                               species.var, 
#                               abundance.var, 
#                               metric=c("total", "disappearance","appearance")) {
#   
#   # allows partial argument matching
#   metric = match.arg(metric) 
#   
#   # sort and remove 0s
#   df <- df[order(df[[time.var]]),]
#   df <- df[which(df[[abundance.var]]>0),]
#   
#   ## split data by year
#   templist <- split(df, df[[time.var]])
#   
#   ## create two time points (first year and each other year)
#   t1 <- templist[1]
#   t2 <- templist[-1]
#   
#   ## calculate turnover for across all time points
#   out <- Map(turnover_twoyears, t1, t2, species.var, metric)
#   output <- as.data.frame(unlist(out))
#   names(output)[1] = metric
#   
#   ## add time variable column
#   alltemp <- unique(df[[time.var]])
#   output[time.var] =  alltemp[2:length(alltemp)]
#   
#   # results
#   return(output)
# }
# 
# turnover_twoyears <- function(d1, d2, 
#                               species.var, 
#                               metric=c("total", "disappearance","appearance")){
#   
#   # allows partial argument matching
#   metric = match.arg(metric)
#   
#   # create character vectors of unique species from each df
#   d1spp <- as.character(unique(d1[[species.var]]))
#   d2spp <- as.character(unique(d2[[species.var]]))
#   
#   # ID shared species
#   commspp <- intersect(d1spp, d2spp)
#   
#   # count number not present in d2
#   disappear <- length(d1spp)-length(commspp)
#   
#   # count number that appear in d2
#   appear <- length(d2spp)-length(commspp)
#   
#   # calculate total richness
#   totrich <- sum(disappear, appear, length(commspp))
#   
#   # output based on metric 
#   if(metric == "total"){
#     output <- totrich
#   } else {
#     if(metric == "appearance"){
#       output <- appear
#     } else {
#       if(metric == "disappearance"){
#         output <- disappear
#       }
#     }
#   }
#   
#   # results
#   return(output)
# }
# 
# 
# # generating appearances and disappearances for each experiment ---------------------------------------------
# #make a new dataframe with just the label;
# site_code=nutnetdf_allspp%>%
#   select(site_code)%>%
#   unique()
# 
# #makes an empty dataframe
# for.analysis=data.frame(row.names=1)
# 
# for(i in 1:length(site_code$site_code)) {
#   
#   #creates a dataset for each unique site-year
#   subset=nutnetdf_allspp[nutnetdf_allspp$site_code==as.character(site_code$site_code[i]),]%>%
#     select(site_code, year, Taxon, rel_cover, plot)%>%
#     group_by(site_code, year, plot, Taxon)%>%
#     summarise(rel_cover=max(rel_cover))%>%
#     ungroup()%>%
#     filter(rel_cover>0)
#   
#   #need this to keep track of sites
#   labels=subset%>%
#     select(year, site_code)%>%
#     unique()
#   
#   #calculating appearances and disappearances (from previous year): for each year
#   appear<-turnover_allyears(df=subset, time.var='year', species.var='Taxon', abundance.var='rel_cover', metric='appearance')
#   disappear<-turnover_allyears(df=subset, time.var='year', species.var='Taxon', abundance.var='rel_cover', metric='disappearance')
#   total<-turnover_allyears(df=subset, time.var='year', species.var='Taxon', abundance.var='rel_cover', metric='total')
#   
#   #merging back with labels to get back experiment labels
#   turnover<-merge(appear, disappear, by=c('year'))
#   turnoverAll<-merge(turnover, total, by=c('year'))
#   turnoverLabel<-merge(turnoverAll, labels, by=c('year'), all=T)
#   
#   #pasting into the dataframe made for this analysis
#   for.analysis=rbind(turnoverLabel, for.analysis)  
# }


#biomass response by trt
biomassResp <- biomass%>%
  left_join(trt)%>%
  filter(year_trt!=0)%>%
  group_by(site_code, year_trt, trt)%>%
  summarise(anpp_mean=mean(anpp))%>%
  ungroup()%>%
  spread(key=trt, value=anpp_mean)%>%
  mutate(NPK_diff=(NPK-Control)/Control, N_diff=(N-Control)/Control, NP_diff=(NP+Fence-Control)/Control)%>%
  select(site_code, year_trt, NPK_diff, N_diff, NP_diff)%>%
  na.omit()%>%
  gather(key=trt, value=diff, NPK_diff:NP_diff)

badBiomassFig <- ggplot(data=barGraphStats(data=biomassResp, variable="diff", byFactorNames=c("year_trt", "trt")), aes(x=year_trt, y=mean, color=trt)) +
  geom_point(size=5) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  xlab('Treatment Year') + ylab('Biomass Difference (%)') +
  geom_hline(yintercept=0)
#export at 1200x800


#figure placeholder - needs updates
pushViewport(viewport(layout=grid.layout(1,3)))
print(BEFstrawmanFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(badLossFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(badBiomassFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
