#########################################
#' What LTERs are involved in BIOTIME
#' 
#' 
#' @author Jarrett Byrnes
#' @details 
#' Changelog:
#########################################

#packages####
library(tidyverse)
library(data.table)

#biotime####
#read in full biotime dataset
biotime <- read_csv("../raw_data/bioTIMEmetadataFeb.csv")
names(biotime_abund)


bt <- unique(biotime$TITLE)
bt[str_detect(bt, "LTER")]
