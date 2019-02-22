#'----------------------------
#' Make 4 panel plot
#' 
#' 
#' 
#'----------------------------

library(ggplot2)
library(patchwork)


bef <- readRDS("../derived_data/4_bef_biotime.Rds")
rarity <- readRDS("../derived_data/2_biotime_sp_loss.Rds")
size <- readRDS("../derived_data/5_biotime_abund_biomass_loss.Rds")
timeseries <- readRDS("../derived_data/6_biomass_timeseries.Rds")

#make the 4 panel plot
jpeg("../figures/four_panel_plot.jpg")
bef$log_biomass_plot +
  rarity$missing_at_end_plot +
  size$plot_fit_loss_mod +
  timeseries$biomass_timeseries +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a")
dev.off()