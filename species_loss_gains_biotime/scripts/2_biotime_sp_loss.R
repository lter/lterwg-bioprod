#########################################
#' Analyze BIOTIME to see if relative abundance 
#' influences loss over time
#' 
#' 
#' @author Jarrett Byrnes
#' @details 
#' Changelog:
#########################################

library(tidyverse)
library(ggplot2)
library(lme4)
library(sjPlot)
library(merTools)
library(DHARMa)
library(car)
library(piecewiseSEM)
library(AICcmodavg)

###

#### Set some graphical parameters

theme_set(theme_bw(base_size=17) +
            theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black")))

#for a labeller
add_year <- function(x) paste(x,"years", sep=" ")

#### Load data ####
load(file="../derived_data/biotime_loss_summary_filtered.Rdata")

load(file="../derived_data/biotime_sum_abund.Rdata")

biotime_loss_summary_filtered <- biotime_loss_summary %>%
  filter(N_YEARS>=3) %>%
  mutate(STUDY_ID = paste(STUDY_ID, MEASURE, sep="_"))

#------------------------
# How important was relative rank for being missing at the end? ####
#------------------------

#Plot the raw data, curve for each site, and a pooled curve
ggplot(biotime_loss_summary_filtered,
       aes(x=REL_RANK_T0, y=lost_at_end)) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_smooth(method="glm", mapping=aes(group=as.character(STUDY_ID)),
              method.args = list(family=binomial),
               fill=NA, lwd=0.4, alpha=0.1, color="lightgrey") +
  geom_smooth(method="glm", mapping= aes(group=1),
              method.args = list(family=binomial),
              color="black", lwd=1.3) +
  ylab("Missing at End?") +
  xlab("Relative Rank at T1")


#Build Models
mod_final_loss <- glmer(lost_at_end ~ REL_RANK_T0 + (1 + REL_RANK_T0|STUDY_ID),
                       data = biotime_loss_summary_filtered,
                       family = binomial())

mod_final_loss_nyear <- glmer(lost_at_end ~ REL_RANK_T0+N_YEARS + (1 + REL_RANK_T0|STUDY_ID),
                        data = biotime_loss_summary_filtered,
                        family = binomial())

mod_final_loss_nyear_int <- glmer(lost_at_end ~ REL_RANK_T0*N_YEARS + (1 + REL_RANK_T0|STUDY_ID),
                              data = biotime_loss_summary_filtered,
                              family = binomial())

aictab(list(mod_final_loss,mod_final_loss_nyear,mod_final_loss_nyear_int))

#Check model assumptions with DHarma
simulationOutput <- simulateResiduals(fittedModel = mod_final_loss, n = 250)
plot(simulationOutput)


#Assuming we've passed, get LRT for fixed effect
Anova(mod_final_loss)
summary(mod_final_loss)
rsquared(mod_final_loss)

#sjPlot the results
sjp.glmer(mod_final_loss)

#get fixef and CI prediction
n <- 100
newdata <- crossing(REL_RANK_T0=seq(0,1,length.out=n), 
                    STUDY_ID = unique(biotime_loss_summary_filtered$STUDY_ID))

newdata_fixed <- data.frame(REL_RANK_T0=seq(0,1,length.out=n), 
                    STUDY_ID = unique(biotime_loss_summary_filtered$STUDY_ID)[1])

pred_final_loss_fixed <- cbind(newdata_fixed, predictInterval(mod_final_loss,
                          newdata=newdata_fixed,
                          which="fixed", type="probability",
                          include.resid.var = FALSE))


#get ranef predictions
pred_final_loss_ranef <- cbind(newdata, predictInterval(mod_final_loss,
                                                              newdata=newdata,
                                                              which="full", type="probability",
                                                              include.resid.var = F))



#replot with raw data and model predictions
jpeg("../figures/effect_of_rarity_on_final_loss.jpg")
ggplot(biotime_loss_summary_filtered,
       aes(x=REL_RANK_T0, y=lost_at_end, group=as.character(STUDY_ID))) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_line(data = pred_final_loss_ranef,
            mapping=aes(y=fit, group=as.character(STUDY_ID)),
              lwd=0.4, alpha=0.4, color="lightgrey") +
  geom_line(data = pred_final_loss_fixed,
              mapping=aes(y=fit),
              color="black", lwd=1.3) +
  ylab("Probability of Being\n Missing at End") +
  xlab("Relative Rank at 1st Time Step\n
       (Rare to Common)")
dev.off()

#---------------------------------------------------------------------------------------
# Does initial rarity affect fraction of the timeseries you are missing ? ####
#---------------------------------------------------------------------------------------


#Plot - it's a hurdle model, so, first, are you missing at all?
ggplot(biotime_loss_summary_filtered,
       aes(x=REL_RANK_T0, y=absent_at_some_point,
           group=as.character(STUDY_ID))) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_smooth(method="glm", mapping=aes(group=as.character(STUDY_ID)),
              method.args = list(family=binomial),
              fill=NA, lwd=0.4, alpha=0.1, color="lightgrey") +
  geom_smooth(method="glm", mapping= aes(group=1),
              method.args = list(family=binomial),
              color="black", lwd=1.3) +
  scale_fill_viridis()+
  ylab("Lost during the time series?") +
  xlab("Relative Rank at T1") 

#For those that are present some of the time, fit a curve

#Plot - it's a hurdle model, so, first, are you missing at all?
ggplot(biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0),
       aes(x=REL_RANK_T0, y=frac_missing_after_t0,
           group=as.character(STUDY_ID))) +
  geom_point(alpha=0.3, color="grey") +
  guides(colour=FALSE) +
  geom_smooth(method="glm", mapping=aes(weight=n_year),
              method.args = list(family=binomial),
              fill = NA,
              alpha=0.4,
              color="darkgrey",
              lwd=0.4)+
  geom_smooth(method="glm", mapping=aes(weight=n_year, group=1),
              method.args = list(family=binomial),
              fill = NA,
              alpha=0.9,
              color="black",
              lwd=1.5,
              fill="lightblue")+
  ylab("Fraction of Time Series\n Not Found") +
  xlab("Relative Rank at T0") 


#Build Models

hurdle_mod <- glmer(absent_at_some_point ~ REL_RANK_T0 + (1 + REL_RANK_T0 | STUDY_ID),
                    data = biotime_loss_summary_filtered,
                    family=binomial(), weights=biotime_loss_summary_filtered$N_YEARS)

hurdle_mod_time <- glmer(absent_at_some_point ~ REL_RANK_T0+N_YEARS + (1 + REL_RANK_T0 | STUDY_ID),
                    data = biotime_loss_summary_filtered,
                    family=binomial(), weights=biotime_loss_summary_filtered$N_YEARS)

hurdle_mod_time_int <- glmer(absent_at_some_point ~ REL_RANK_T0*N_YEARS + (1 + REL_RANK_T0 | STUDY_ID),
                         data = biotime_loss_summary_filtered,
                         family=binomial(), weights=biotime_loss_summary_filtered$N_YEARS)

aictab(list(hurdle_mod, hurdle_mod_time, hurdle_mod_time_int))
#additive time model seems best

rsquared(hurdle_mod)
rsquared(hurdle_mod_time)


#### Frac missing
frac_missing_mod <- glmer(frac_missing_after_t0 ~ REL_RANK_T0 + (1 + REL_RANK_T0 | STUDY_ID),
                    data = biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0),
                    family=binomial(),
                    weights = (biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0))$n_year)


frac_missing_mod_time <- glmer(frac_missing_after_t0 ~ REL_RANK_T0+I(scale(N_YEARS, scale=F)) + 
                                 (1 + REL_RANK_T0 | STUDY_ID),
                          data = biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0),
                          family=binomial(),
                          weights = (biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0))$n_year)



frac_missing_mod_time_int <- glmer(frac_missing_after_t0 ~ REL_RANK_T0*I(scale(N_YEARS, scale=F)) + 
                                     (1 + REL_RANK_T0 | STUDY_ID),
                          data = biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0),
                          family=binomial(),
                          weights = (biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0))$n_year)

aictab(list(frac_missing_mod, frac_missing_mod_time, frac_missing_mod_time_int))
rsquared(frac_missing_mod)
rsquared(frac_missing_mod_time_int)

#Check model assumptions
simulationOutputHurdle <- simulateResiduals(fittedModel = hurdle_mod, n = 250)
plot(simulationOutputHurdle)

simulationOutputFrac<- simulateResiduals(fittedModel = frac_missing_mod, n = 250)
plot(simulationOutputFrac)

#Assuming we've passed, get LRT for fixed effect
Anova(hurdle_mod)
Anova(frac_missing_mod)
#get fixef and CI prediction
pred_hurdle_fixed <- cbind(newdata_fixed, predictInterval(hurdle_mod,
                                                              newdata=newdata_fixed,
                                                              which="fixed", type="probability",
                                                              include.resid.var = FALSE))


#get ranef predictions
pred_hurdle_ranef <- cbind(newdata, predictInterval(hurdle_mod,
                                                        newdata=newdata,
                                                        which="full", type="probability",
                                                        include.resid.var = FALSE))

pred_frac_missing_fixed <- cbind(newdata_fixed, predictInterval(frac_missing_mod,
                                                              newdata=newdata_fixed,
                                                              which="fixed", type="probability",
                                                              include.resid.var = FALSE))


#get ranef predictions
pred_frac_missing_ranef <- cbind(newdata, predictInterval(frac_missing_mod,
                                                        newdata=newdata,
                                                        which="full", type="probability",
                                                        include.resid.var = FALSE))


##replot with raw data and model predictions
jpeg("../figures/effect_of_rarity_on_missing_hurdle_notime.jpg")
ggplot(biotime_loss_summary_filtered,
       aes(x=REL_RANK_T0, y=absent_at_some_point, group=as.character(STUDY_ID))) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_line(data = pred_hurdle_ranef,
            mapping=aes(y=fit, group=as.character(STUDY_ID)),
            lwd=0.4, alpha=0.4, color="lightgrey") +
  geom_line(data = pred_hurdle_fixed,
            mapping=aes(y=fit),
            color="black", lwd=1.3) +
  ylab("Missing at some point\nduring the time series?") +
  xlab("Relative Rank at 1st Time Step\n (Rare to Common)")
dev.off()

jpeg("../figures/effect_of_rarity_on_missing_fraction_notime.jpg")
ggplot(biotime_loss_summary_filtered %>% filter(frac_missing_after_t0 !=0),
       aes(x=REL_RANK_T0, y=frac_missing_after_t0, group=as.character(STUDY_ID))) +
  geom_point(position=position_jitter(width=0.05, height=0.05),
             alpha=0.2, color="grey") +
  geom_line(data = pred_frac_missing_ranef,
            mapping=aes(y=fit, group=as.character(STUDY_ID)),
            lwd=0.5, alpha=0.4, color="darkgrey") +
  geom_line(data = pred_frac_missing_fixed,
            mapping=aes(y=fit),
            color="black", lwd=1.3) +
  ylab("Fraction of Time Series\n Not Found") +
  xlab("Relative Rank at 1st Time Step\n (Rare to Common)")
dev.off()

## Let's think about visualizing the hurdle and fraction analysis with time

cut_interval(biotime_loss_summary_filtered$N_YEARS, n=4)
times <- c(5, 15, 35, 50)

newdata_time <- crossing(REL_RANK_T0=seq(0,1,length.out=n/4),
                         N_YEARS=times,
                    STUDY_ID = unique(biotime_loss_summary_filtered$STUDY_ID))

newdata_fixed_time <- crossing(data.frame(REL_RANK_T0=seq(0,1,length.out=n), 
                            STUDY_ID = unique(biotime_loss_summary_filtered$STUDY_ID)[1]),
                            N_YEARS=times)

hurdle_time_fixed <- cbind(newdata_fixed_time, predictInterval(hurdle_mod_time,
                                                              newdata=newdata_fixed_time,
                                                              which="fixed", type="probability",
                                                              include.resid.var = TRUE))


#get ranef predictions
hurdle_time_ranef <- cbind(newdata_time, predictInterval(hurdle_mod_time,
                                                        newdata=newdata_time,
                                                        which="full", type="probability",
                                                        include.resid.var = FALSE))


jpeg("../figures/effect_of_rarity_on_missing_hurdle_time.jpg")
ggplot(hurdle_time_fixed) + 
  geom_line(hurdle_time_ranef, alpha=0.1,
            mapping=aes(x=REL_RANK_T0, y=fit, 
                                   color=as.factor(N_YEARS),
                                   group=paste(as.character(N_YEARS), STUDY_ID))) +
  geom_line(hurdle_time_fixed, mapping=aes(x=REL_RANK_T0, y=fit, 
                                   
                                   color=as.factor(N_YEARS))) +
  geom_ribbon(hurdle_time_fixed, alpha=0.3, color=NA, 
              mapping=aes(x=REL_RANK_T0, ymin=lwr, ymax=upr, fill=as.factor(N_YEARS))) +
  facet_wrap(~N_YEARS, labeller=labeller(N_YEARS = add_year))  +
  ylab("Missing at some point\nduring the time series?") +
  xlab("Relative Rank at 1st Time Step\n (Rare to Common)") +
  scale_color_discrete(guide = guide_legend(title = "Duration of\nTimeseries"))+
  scale_fill_discrete(guide = guide_legend(title = "Duration of\nTimeseries"))
dev.off()

## Fraction of loss
frac_time_fixed <- cbind(newdata_fixed_time, predictInterval(frac_missing_mod_time_int,
                                                               newdata=newdata_fixed_time,
                                                               which="fixed", type="probability",
                                                               include.resid.var = TRUE))

ggplot(frac_time_fixed, aes(x=REL_RANK_T0, y=fit, 
                              ymin=lwr, ymax=upr,
                              color=as.factor(N_YEARS),
                              fill=as.factor(N_YEARS))) + 
  geom_line() +
  geom_ribbon(alpha=0.3, color=NA)

#get ranef predictions
frac_time_ranef <- cbind(newdata_time, predictInterval(frac_missing_mod_time_int,
                                                         newdata=newdata_time,
                                                         which="full", type="probability",
                                                         include.resid.var = FALSE))


jpeg("../figures/effect_of_rarity_on_missing_frac_time.jpg", width=800, height=600)
ggplot(frac_time_fixed) + 
  geom_line(frac_time_ranef, alpha=0.1,
            mapping=aes(x=REL_RANK_T0, y=fit, 
                        color=as.factor(N_YEARS),
                        group=paste(as.character(N_YEARS), STUDY_ID))) +
  geom_line(frac_time_fixed, mapping=aes(x=REL_RANK_T0, y=fit, 
                                           color=as.factor(N_YEARS))) +
  geom_ribbon(frac_time_fixed, alpha=0.3, color=NA, 
              mapping=aes(x=REL_RANK_T0, ymin=lwr, ymax=upr, fill=as.factor(N_YEARS))) +
  facet_wrap(~N_YEARS, labeller=labeller(N_YEARS = add_year))  +
  ylab("Fraction of Time Series\n Not Found") +
  xlab("Relative Rank at 1st Time Step\n (Rare to Common)") +
  scale_color_discrete(guide = guide_legend(title = "Duration of\nTimeseries"))+
  scale_fill_discrete(guide = guide_legend(title = "Duration of\nTimeseries"))
dev.off()