#########################################
#' Analyze BIOTIME to look for BEF
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
#load(file="../derived_data/biotime_loss_summary_filtered.Rdata")
load(file="../derived_data/biotime_loss_summary.Rdata")
load(file="../derived_data/biotime_sum_abund.Rdata")


biotime_biomass <- biotime_sum_abund %>%
  filter(MEASURE=="Biomass") %>%
  group_by(STUDY_ID, MEASURE, GENUS_SPECIES) %>%
  mutate(N_YEARS = n_distinct(YEAR)) %>%
  ungroup() %>%
  filter(N_YEARS>=3) %>%
  group_by(STUDY_ID, YEAR) %>%
  summarize(Z = sum(SUM_VALUE==0),
            S = n_distinct(GENUS_SPECIES)-Z,
            B = sum(SUM_VALUE, na.rm=T),
            N_YEARS=N_YEARS[1]) %>%
  ungroup() %>%
  group_by(STUDY_ID) %>%
  mutate(STD_YEAR = YEAR - min(YEAR)) %>%
  ungroup() %>%
  group_by(STUDY_ID) %>%
  arrange(STD_YEAR) %>%
  mutate(DS = S - lag(S),
         DB = B - lag(B),
         LRDS = log(S+1) - log(lag(S)+1),
         LRDB = log(B+1) - log(lag(B)+1)) %>%
  ungroup()

#A plot of the S-biomass relationship
ggplot(biotime_biomass, aes(x=log(S), y=log(B), color=as.character(STUDY_ID))) + 
  geom_point() + 
  guides(color = "none") + 
  stat_smooth(method="lm", fill=NA)

ggplot(biotime_biomass, aes(x=DS, y=DB, color=as.character(STUDY_ID))) + 
  geom_point() + 
  guides(color = "none") 

ggplot(biotime_biomass, aes(x=LRDS, y=LRDB, color=as.character(STUDY_ID))) + 
  geom_point() + 
  guides(color = "none") +
  xlab("Log Ratio of Richness Change\nLog(S(t)/S(t-1) +1)")+
  ylab("Log Ratio of Biomass Change\nLog(B(t)/B(t-1) +1)")

####---------
#Let's do some modeling
####-----------
bef_mod <- lmer(log(B+1) ~ log(S+1)*STD_YEAR + (log(S+1) | STUDY_ID), data=biotime_biomass)

#evaluate autocor
biotime_biomass <- biotime_biomass %>%
  mutate(res_bef_mod = residuals(bef_mod)) %>%
  group_by(STUDY_ID) %>%
  mutate(lag_res_bef_mod = lag(res_bef_mod)) %>%
  ungroup()

qplot(lag_res_bef_mod, res_bef_mod, data = biotime_biomass) #looks good!

pred_bef_fix <- cbind(biotime_biomass,
                      predictInterval(bef_mod, 
                                      newdata=biotime_biomass %>% mutate(STD_YEAR=1), 
                                      which="fixed"))
pred_bef_ranef <- cbind(biotime_biomass,
                      predictInterval(bef_mod, newdata=biotime_biomass, 
                                      which="full"))

ggplot(biotime_biomass, aes(x=log(S+1), y=log(B+1), color=as.character(STUDY_ID))) + 
  geom_point() + 
  guides(color = "none") +
  geom_ribbon(pred_bef_fix, mapping=aes(ymin=lwr, ymax=upr), fill="lightgrey", color=NA, alpha=0.8) +
  geom_line(pred_bef_ranef, mapping=aes(y=fit), alpha=0.7, lwd=1.3) +
  geom_line(pred_bef_fix, mapping=aes(y=fit), color="black", lwd=1.5) +
  xlab("Log Richness + 1") + ylab("Log Biomass + 1")


####---------
#Let's do some modeling  of deltas
####-----------
change_mod <- lmer(LRDB ~ LRDS + (LRDS | STUDY_ID) , data=biotime_biomass)

#assumptions
_nyear,mod_final_loss_nyear_int))

#Check model assumptions with DHarma
simulationOutput <- simulateResiduals(fittedModel = change_mod, n = 250, 
                                      plot=T)



pred_lr_fix <- cbind(biotime_biomass,
                      predictInterval(change_mod, 
                                      newdata=biotime_biomass, 
                                      which="fixed"))
pred_lr_ranef <- cbind(biotime_biomass,
                        predictInterval(change_mod, newdata=biotime_biomass, 
                                        which="full"))

ggplot(biotime_biomass, aes(x=LRDS, y=LRDB, color=as.character(STUDY_ID))) + 
  geom_point() + 
  guides(color = "none") +
  geom_ribbon(pred_lr_fix, mapping=aes(ymin=lwr, ymax=upr), fill="lightgrey", color=NA, alpha=0.8) +
  geom_line(pred_lr_ranef, mapping=aes(y=fit), alpha=0.5, lwd=0.9) +
  geom_line(pred_lr_fix, mapping=aes(y=fit), color="black", lwd=1.5)+
  xlab("Log Ratio of Richness Change\nLog(S(t)/S(t-1) +1)")+
  ylab("Log Ratio of Biomass Change\nLog(B(t)/B(t-1) +1)")


#-------------
# Slope versus CV? From bef_mod
#-------------
biotime_biomass_cv <- biotime_biomass %>%
  group_by(STUDY_ID) %>%
  summarize(CV_B = sd(B, na.rm=T)/mean(B, na.rm=T)) %>%
  ungroup() %>%
  left_join(tibble(STUDY_ID = as.numeric(rownames(ranef(bef_mod)$STUDY_ID)),
                       slope = ranef(bef_mod)$STUDY_ID[,2]))

ggplot(biotime_biomass_cv, aes(x=slope, y=CV_B)) +
  geom_point() +
  stat_smooth(method="lm") +
  xlab("Individual Study BLUP of\nSlope Random Effect") +
  ylab("Temporal Coefficient of Variation\nin Biomass")
