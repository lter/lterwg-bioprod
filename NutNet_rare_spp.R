### NutNet and rare species #####
# June 13 2017

#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

setwd("~/Google Drive/LTER_Biodiversity_Productivity")

#kim's wd
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\nutrient network\\NutNet data')

library(tidyverse)
library(grid)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

nutnetdf <-read.csv("full-cover-31-August-2018.csv")
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
  mutate(rel_cover=(max_cover/total_cover)*100)

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


# make a column for presence absence of each species in a plot-year
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
  group_by(site_code, Taxon, site_name, block, plot, subplot, year_trt2, length)%>% ###there's duplicate entries for some spp; need to fix this some day
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
  mutate(abund_metric=(2*meanPTAbundance+freq)/3)%>%
  select(-PTfreq)%>%
  filter(length>0)


###figures of proportion of years absent for the various abundance metrics
metricPlot <- ggplot(nutnetdf_allspp3, aes(x=abund_metric, y=yrs_absent, color=length)) +
  geom_point(alpha=0.1, size=3) +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Proportion of Years Absent')

freqPlot <- ggplot(nutnetdf_allspp3, aes(x=freq, y=yrs_absent, color=length)) +
  geom_point(alpha=0.1, size=3) +
  xlab('Pre-Treatment Frequency') +
  ylab('Proportion of Years Absent')

avgAbundPlot <- ggplot(nutnetdf_allspp3, aes(x=meanPTAbundance, y=yrs_absent, color=length)) +
  geom_point(alpha=0.1, size=3) +
  xlab('Pre-Treatment Mean Abundance') +
  ylab('Proportion of Years Absent')

maxAbundPlot <- ggplot(nutnetdf_allspp3, aes(x=maxPTAbundance, y=yrs_absent, color=length)) +
  geom_point(alpha=0.1, size=3) +
  xlab('Pre-Treatment Max Abundance') +
  ylab('Proportion of Years Absent')

pushViewport(viewport(layout=grid.layout(2,2)))
print(metricPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(freqPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(avgAbundPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(maxAbundPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600


#get back trt info
trt <- nutnetRel%>%
  select(year_trt, site_code, plot, trt)%>%
  filter(year_trt>0)%>%
  select(-year_trt)%>%
  unique()

nutnetdf_allspp3Trt <- nutnetdf_allspp3%>%
  merge(trt, by=c('site_code', 'plot'))

ggplot(nutnetdf_allspp3Trt, aes(x=abund_metric, y=yrs_absent, color=length)) +
  geom_point(size=3) +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Proportion of Years Absent') +
  facet_wrap(~trt)
#export at 1200x1200


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
  mutate(abund_metric=(2*meanPTAbundance+PTfreq)/3)%>%
  merge(nutnetdf_length, by=c('site_code'))

ggplot(subset(nutnetConsAbs2, length==9), aes(x=abund_metric, y=cons_abs_max)) +
  geom_point(alpha=0.1, size=3) +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Consecutive Years Absent') +
  facet_wrap(~trt)
#export at 1200x1200

# ggplot(subset(nutnetConsAbs2, length==9), aes(x=meanPTAbundance, y=cons_abs_max)) +
#   geom_point(alpha=0.1) +
#   xlab('Pre-Treatment Mean Abundance') +
#   ylab('Consecutive Years Absent') +
#   facet_wrap(~trt)





###figures for year and pres/absence
nutnetPresAbs <- nutnetdf_allspp%>%
  merge(meanAb_byspecies, by=c('site_code', 'Taxon'))%>%
  merge(max_abund, by=c('site_code', 'Taxon'))%>%
  merge(freq, by=c('site_code', 'Taxon'))%>%
  mutate(abund_metric=(2*meanPTAbundance+PTfreq)/3)%>%
  filter(year_trt>0)

ggplot(nutnetPresAbs, aes(x=abund_metric, y=PA, color=trt)) +
  geom_point() +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Presence (1)/Absence (0)') +
  facet_wrap(~year_trt)



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
biomass <- read.csv('full-biomass-09-June-2017.csv')%>%
  filter(live==1)%>%
  group_by(site_code, plot, subplot, year_trt)%>%
  summarise(mass2=sum(mass))%>%
  ungroup()%>%
  group_by(site_code, plot, year_trt)%>%
  summarise(anpp=mean(mass2))

biomassSpp <- nutnetdf_allspp%>%
  filter(rel_cover>0)%>%
  merge(biomass, by=c('site_code', 'plot', 'year_trt'))%>%
  merge(meanAb_byspecies, by=c('site_code', 'Taxon'))%>%
  merge(max_abund, by=c('site_code', 'Taxon'))%>%
  merge(freq, by=c('site_code', 'Taxon'))%>%
  mutate(abund_metric=(2*meanPTAbundance+freq)/3)%>%
  select(-PTfreq)%>%
  group_by(site_code, plot, year_trt)%>%
  summarise(biomass=mean(anpp), importance=mean(abund_metric), min_importance=min(abund_metric))

ggplot(biomassSpp, aes(x=importance, y=biomass)) +
  geom_point(size=3, alpha=0.1) +
  xlab('Pre-Treatment Modified Importance Index') +
  ylab('Aboveground Biomass')
#export at 900x900



















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