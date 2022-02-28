################################################################################
# CALCULATE DATE DISTRIBUTIONS FOR SITES
#
################################################################################

library(tidyverse)
library(lubridate)
library(viridis)
library(sf)

# load data filtered to just within site polygons
bandedbirds_resights <- readRDS("./data/bandedbirds_resights_subset.rds")

# check if any BirdIDs are NAs
any(is.na(bandedbirds_resights$BirdID))

# get just columns needed
bb_data <- bandedbirds_resights %>% 
  select(BirdID, site, ResightDate) %>% 
  st_drop_geometry()

# drop year from dates
bb_dates <- bb_data %>% 
  mutate(daymonth = ResightDate) %>% 
  select(-ResightDate)

bb_dates$daymonth <- ymd(bb_dates$daymonth)
bb_dates$daymonth <- format(bb_dates$daymonth, format = "%d-%b")

# get min/max dates for each site
minmax_dates <- bb_dates %>% 
  mutate(julianday = yday(daymonth)) %>% 
  group_by(site) %>% 
  summarize(min = min(julianday),
            max = max(julianday)) %>% 
  ungroup()
