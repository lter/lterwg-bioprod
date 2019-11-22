#########################################
#' Process BIOTIME abundance data to get
#' species loss and change in relative abundance
#' 
#' 
#' @author Jarrett Byrnes
#' @details 
#' Changelog:
#########################################


#figure out who is plant, animal, alga, etc
#rel rank is wonky

#packages####
library(tidyverse)
library(data.table)
library(ggplot2)
library(lme4)
library(DHARMa)
library(merTools)
library(modelr)

theme_set(theme_bw(base_size = 16))

#biotime####
#read in full biotime dataset
#load(file="../derived_data/biotime_abund.Rdata")
load(file="../derived_data/biotime_sum_abund.Rdata")

biotime_duo <- biotime_sum_abund %>%
  group_by(STUDY_ID) %>%
  mutate(num_measures = length(unique(MEASURE))) %>%
  filter(num_measures==2)  %>%
  filter(length(unique(YEAR))>3)


biotime_duo_reshape <- biotime_duo %>%
  dplyr::select(STUDY_ID, YEAR, GENUS_SPECIES, MEASURE, SUM_VALUE) %>%
  spread(MEASURE, SUM_VALUE) %>%
  group_by(STUDY_ID, YEAR) %>%
  #make ranks, etc
  mutate(abund_rank = rank(Abundance,  ties.method = "min"),
         biomass_rank = rank(Biomass, ties.method = "min"),
         percapita_biomass_rank = rank(Biomass/Abundance, ties.method = "min")
  ) %>%
  #deal with zeroes - make them rank 0 and recenter ranks
  mutate(abund_rank = ifelse(Abundance==0,0,abund_rank),
         biomass_rank = ifelse(Abundance==0,0,biomass_rank),
         percapita_biomass_rank = ifelse(Abundance==0,0,percapita_biomass_rank)
  ) %>%
  #rescaling
  mutate(abund_rank = abund_rank - min(abund_rank[abund_rank!=0]) + 1,
         biomass_rank = biomass_rank - min(biomass_rank[biomass_rank!=0]) + 1,
         percapita_biomass_rank = percapita_biomass_rank - min(percapita_biomass_rank[percapita_biomass_rank!=0]) + 1) %>%
  mutate(abund_rank = ifelse(abund_rank<0,0,abund_rank),
         biomass_rank = ifelse(biomass_rank<0,0,biomass_rank),
         percapita_biomass_rank = ifelse(percapita_biomass_rank<0,0,percapita_biomass_rank)) %>%
  
  #make relative ranks, etc
  #group_by(STUDY_ID) %>%
  mutate(rel_abund = Abundance / max(Abundance, na.rm=TRUE),
         rel_biomass = Biomass / max(Biomass, na.rm=TRUE),
         rel_abund_rank = abund_rank / max(abund_rank, na.rm=TRUE),
         rel_biomass_rank = biomass_rank / max(biomass_rank, na.rm=TRUE),
         rel_percapita_biomass_rank = percapita_biomass_rank / max(percapita_biomass_rank, na.rm=TRUE),
  ) %>%
  # put together losses and initial relative ranks 
  group_by(STUDY_ID, GENUS_SPECIES) %>%
  arrange(YEAR) %>%
  mutate(rel_abund_rank = rel_abund_rank[1],
         rel_biomass_rank = rel_biomass_rank[1],
         rel_percapita_biomass_rank = rel_percapita_biomass_rank[1],
         rel_abund_rank_t1_t3 = mean(rel_abund_rank[1:3]),
         rel_biomass_rank_t1_t3 = mean(rel_biomass_rank[1:3]),
         rel_percapita_biomass_rank_t1_t3 = mean(rel_percapita_biomass_rank[1:3]),
         lag_rel_abund_rank = lag(rel_abund_rank),
         lag_rel_biomass_rank = lag(rel_biomass_rank),
         lag_rel_percapita_biomass_rank = lag(rel_percapita_biomass_rank),
         lost = Abundance==0 & lag_rel_abund_rank !=0,
         present_t1_t3 = sum(Abundance[1:3]!=0) == 3) %>%
  ungroup()
    
#BEF ####
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
       mapping = aes(x = rel_abund_rank, y = rel_percapita_biomass_rank, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
  stat_smooth(method="lm", mapping = aes(group = 1))
ggsave("../figures/rel_abund_rel_percapita_biomass.jpg")




#relationship between lost species biomass rank and abundance rank

ggplot(biotime_duo_reshape %>% filter(lost), 
       mapping = aes(x = lag_rel_abund_rank, y = lag_rel_biomass_rank, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
  stat_smooth(method="lm", mapping = aes(group = 1))


ggplot(biotime_duo_reshape %>% filter(lost) %>% filter(present_t1_t3), 
       mapping = aes(x = rel_abund_rank_t1_t3, y = rel_biomass_rank_t1_t3, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
#  stat_smooth(method="lm", mapping = aes(group = 1)) +
  stat_density_2d(aes(group = 1))

#relationship between lost species per capita biomass rank and abundance rank

ggplot(biotime_duo_reshape %>% filter(lost), 
       mapping = aes(x = lag_rel_abund_rank, y = lag_rel_percapita_biomass_rank, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.1) +
  guides(color = "none") +
  stat_smooth(method="lm", mapping = aes(group = 1)) +
#  stat_density_2d(aes(group = 1, fill = ..density..), geom = "polygon") +
  scale_fill_viridis_c()

ggplot(biotime_duo_reshape %>% filter(lost) %>% filter(present_t1_t3), 
       mapping = aes(x = rel_abund_rank_t1_t3, y = rel_percapita_biomass_rank_t1_t3, color = factor(STUDY_ID))) +
  geom_point(alpha = 0.5) +
  guides(color = "none")

#models ####
loss_mod_lag <- glmer(lost ~ lag_rel_abund_rank*lag_rel_percapita_biomass_rank +
                    (1|GENUS_SPECIES) +
                    (1 |STUDY_ID),
                  data = biotime_duo_reshape,
                  family = "binomial")

piecewiseSEM::rsquared(loss_mod_lag)


car::Anova(loss_mod_lag)
summary(loss_mod_lag)

simulationOutput <- simulateResiduals(fittedModel = loss_mod_lag, n = 250)
plot(simulationOutput)

# prediction plots ####

n <- 100
f <- factor(levels(cut_interval(seq(0,1,.01),4)), levels = levels(cut_interval(seq(0,1,.01),4)))

newdata <- crossing(lag_rel_abund_rank = seq(0,1,length.out=n),
                     tibble(lag_rel_percapita_biomass_rank = c(0.125, 0.375, 0.625, 0.875),
                            lagsize_split = f),
                    biotime_duo %>%
                      group_by(STUDY_ID, GENUS_SPECIES) %>%
                      slice(1L) %>%
                      ungroup() %>%
                      dplyr::select(STUDY_ID, GENUS_SPECIES) %>% slice(1L))
  
loss_lag_pred_fix <- predictInterval(loss_mod_lag,
                                 newdata=newdata,
                                 which="fixed", type="probability",
                                 include.resid.var = FALSE)

loss_lag_pred_fix <- cbind(newdata, loss_lag_pred_fix)


ggplot(data = biotime_duo_reshape %>% 
         mutate(lagsize_split = cut_interval(lag_rel_percapita_biomass_rank, 4)) %>%
         filter(!is.na(lagsize_split))
) +
  geom_jitter(mapping = aes(x = lag_rel_abund_rank, y = as.numeric(lost)),
              position = position_jitter(width = 0.1, height = 0.1),
              alpha = 0.1) +
  facet_wrap(~lagsize_split) +
   geom_line(data = loss_lag_pred_fix,
            mapping = aes(x = lag_rel_abund_rank, y = fit, group = lagsize_split,
                          color = factor(lag_rel_percapita_biomass_rank)),
            size = 2) +
  scale_color_viridis_d(guide = guide_legend(title = "Size Rank"),
                        option = "D") +
  scale_fill_viridis_d(guide = guide_legend(title = "Size Rank"),
                        option = "D") +
  geom_ribbon(data = loss_lag_pred_fix,
              mapping = aes(x = lag_rel_abund_rank,  group = lagsize_split,
                            fill = factor(lag_rel_percapita_biomass_rank),
                            ymin = lwr, ymax = upr), alpha = 0.5) 
  


#Show no data
plot_fit_loss_mod <- ggplot() +
  geom_line(data = loss_lag_pred_fix,
            mapping = aes(x = lag_rel_abund_rank, y = fit, group = lagsize_split,
                          color = factor(lag_rel_percapita_biomass_rank)),
            size = 2) +
  scale_color_viridis_d(guide = guide_legend(title = "Relative\nSize Rank\nAt Time T-1"),
                        option = "D") +
  scale_fill_viridis_d(guide = guide_legend(title = "Relative\nSize Rank\nAt Time T-1"),
                       option = "D") +
  geom_ribbon(data = loss_lag_pred_fix,
              mapping = aes(x = lag_rel_abund_rank,  group = lagsize_split,
                            fill = factor(lag_rel_percapita_biomass_rank),
                            ymin = lwr, ymax = upr), alpha = 0.3) +
  xlab("Relative Abundance Rank At Time T-1") +
  ylab("Probability of Loss at Time T")
plot_fit_loss_mod
ggsave("../figures/loss_size_abund_rank_model.jpg")


#The things we save...

saveRDS(list(loss_mod_lag = loss_mod_lag,
             loss_lag_pred_fit = loss_lag_pred_fix, 
             plot_fit_loss_mod = plot_fit_loss_mod),
        file = "../derived_data/5_biotime_abund_biomass_loss.Rds")
