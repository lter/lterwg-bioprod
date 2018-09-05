#########################################
#' Analyze BIOTIME to see if relative abundance 
#' influences change in relative abundance over time
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
library(glmmTMB)


i#### Set some graphical parameters

theme_set(theme_bw(base_size=17) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black")))

#### Load data ####
load(file="../derived_data/biotime_loss_summary.Rdata")

load(file="../derived_data/biotime_sum_abund.Rdata")


biotime_sum_abund_filtered <- biotime_sum_abund %>%
  filter(REL_VALUE_COMM_T0_T3>0,
         N_YEARS>=9) %>%
  group_by(STUDY_ID)%>%
  mutate(YEAR_STD = YEAR - min(YEAR)+1) %>%
  ungroup()

#### Plots of Raw Data ####

ggplot(biotime_sum_abund_filtered , 
       aes(x=YEAR_STD,
           y = REL_VALUE_COMM,
           color=as.character(STUDY_ID),
           group = paste(STUDY_ID, GENUS_SPECIES))) +
  geom_line(alpha=0.3)+
  guides(colour=FALSE) +
  facet_wrap( ~  cut(REL_VALUE_COMM_T0_T3, c(0, 0.01, 0.2, 0.5, 1),
                     include.lowest=TRUE),
              scales = "free_y")  +
  theme_bw(base_size=17) +
  ylab("Relative Abundance in Community") +
  xlab("Years Since First Sample") 

ggplot(biotime_sum_abund_filtered , 
       aes(x=YEAR_STD,
           y = REL_VALUE_COMM,
           color=as.character(STUDY_ID),
           group = paste(STUDY_ID, GENUS_SPECIES))) +
  geom_line(alpha=0.3)+
  guides(colour=FALSE) +
  facet_wrap( ~ cut(REL_RANK_T0, c(0, 0.01, 0.25, 0.5, 0.75, 1),
                    include.lowest=TRUE))  +
  theme_bw(base_size=17) +
  ylab("Relative Abundance in Community") +
  xlab("Years Since First Sample") 

ggplot(biotime_sum_abund_filtered , 
       aes(x=YEAR_STD,
           y = REL_RANK_TIMEPOINT,
           color=as.character(STUDY_ID),
           group = paste(STUDY_ID, GENUS_SPECIES))) +
  geom_line(alpha=0.3)+
  guides(colour=FALSE) +
  facet_wrap( ~  cut(REL_RANK_T0, c(0, 0.001, 0.1, 0.25, 0.5, 1),
                     include.lowest=TRUE))  +
  theme_bw(base_size=17) +
  ylab("Relative Rank in Community") +
  xlab("Years Since First Sample") 


#### Model with lme4 ####

mod_com_abund_mixed <- lmer(REL_RANK_TIMEPOINT ~ 
                              YEAR_STD*REL_RANK_T0 +
                              (1|STUDY_ID) +
                              (1|GENUS_SPECIES),
                            data = biotime_sum_abund_filtered)

mod_com_abund <- glmmTMB(REL_VALUE_COMM ~ 
                              YEAR_STD*REL_RANK_T0 +
                              (1|STUDY_ID) +
                              (1|GENUS_SPECIES),
                         family=betar(link = "probit"),
                         data = biotime_sum_abund_filtered)

glmmTMB(REL_VALUE_COMM ~ 
          YEAR_STD*REL_RANK_T0, 
        family=list(),
        data = biotime_sum_abund_filtered)


#checks of model
simulationOutputAbund <- simulateResiduals(mod_com_abund_mixed)
plotSimulatedResiduals(simulationOutputAbund, quantreg = FALSE)

broom::tidy(mod_com_abund_mixed)

#Betareg?

#Make a prediction data frame
#we're only lookign at fixed
#so, study and species don't matter
pred_df <- crossing(REL_VALUE_COMM_T0_T3 = seq(0.025, 0.975, length.out=100),
                    YEAR_STD=c(5,25,50), 
                    STUDY_ID=biotime_sum_abund_filtered$STUDY_ID[1],
                    GENUS_SPECIES=biotime_sum_abund_filtered$GENUS_SPECIES[100]) 

pred_df <- cbind(pred_df,
                 predictInterval(mod_com_abund_mixed, 
                                 pred_df, 
                                 which="fixed",
                                 include.resid.var = TRUE)
) %>%
  mutate(fit = arm::invlogit(fit),
         upr = arm::invlogit(upr),
         lwr = arm::invlogit(lwr))
  
ggplot(pred_df, 
       aes(x=REL_VALUE_COMM_T0_T3,
           y = fit,
           ymin=lwr,
           ymax=upr,
           color = factor(YEAR_STD),
           fill = factor(YEAR_STD))) +
  geom_ribbon(color=NA, alpha=0.7) +
  geom_line() +
  theme_bw(base_size=17) +
  ylab("Relative Abundance in Community") +
  xlab("Mean Relative Abundance in Years 1-3") +
  guides(color = guide_legend(title = "Years Since\nFirst Sample\n"),
         fill = guide_legend(title = "Years Since\nFirst Sample\n"))
