################################################################################
# PREPARE MULTISTATE ENCOUNTER HISTORIES
#
################################################################################

library(tidyverse)
library(readxl)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(fuzzyjoin)

# map data
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
lakes <- ne_load(type = "lakes", scale = "medium", category = "physical",
                 returnclass = "sf",
                 destdir = "./map-data/lakes")

# load data
all_data <- read_excel("./data/all_projects_resights.xlsx", sheet = "REKNresights")
jb_data <- read_excel("./data/jb_project_resights.xlsx", sheet = "Sheet1")
cc_data <- read_excel("./data/cc_project_resights.xlsx", sheet = "CapeCod_REKNresights_JEL")
br_data <- read_excel("./data/br_project_resights.xlsx", sheet = "proj87REKNresights")

all_resights <- bind_rows(all_data, jb_data)
all_resights <- bind_rows(all_resights, br_data)

# summarize number of resights per project, year
resights_year <- all_resights %>% 
  mutate(year = year(ResightDate)) %>% 
  group_by(year) %>% 
  summarize(n = n())

resights_project <- all_resights %>% 
  mutate(year = year(ResightDate)) %>% 
  filter(year %in% 2009:2018) %>%
  group_by(ProjectID) %>% 
  summarize(n = n())

# start filtering process
all_resights <- all_resights %>% 
  mutate(year = year(ResightDate)) %>% 
  filter(year %in% 2009:2018) %>% 
  filter(ProjectID %in% c(3, 6, 8, 15, 21, 22, 27, 31, 41, 69, 87)) %>% 
  mutate(Longitude = ifelse(Longitude > 0, Longitude*(-1), Longitude))

# remove uncertain resights
all_resights_clean <- all_resights %>% 
  mutate(ResightCertainty = str_replace(ResightCertainty, "^75% \\(COULD BE TV9\\)\\?$", "75")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^100%$", "100")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.95$", "95")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.75$", "75")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^NOT 100% SURE$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^POSSIBLY H$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^90%$", "90")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^1$", "100")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.8$", "80")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^75%$", "75")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^NOT 100% - COULD HAVE BEEN AHU\\?$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^N VERY FADED$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.5$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^80%$", "80")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^0.9$", "90")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^COULD BE N9V$", "50")) %>% 
  mutate(ResightCertainty = str_replace(ResightCertainty, "^Y$", "100")) %>% 
  mutate(ResightCertainty = str_replace(ResightCertainty, "^C$", "100")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^\\?$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^N$", "50")) %>%
  mutate(ResightCertainty = str_replace(ResightCertainty, "^L$", "50"))

all_resights_clean$ResightCertainty <- as.numeric(all_resights_clean$ResightCertainty)  
  
all_resights_clean <- all_resights_clean %>% 
  filter(!str_detect(FlagCode, "Q")) %>% 
  filter(is.na(ResightCertainty) | ResightCertainty > 94 | ResightCertainty == 100 )

# remove/correct manually identified suspect detections
unique(all_resights_clean$FlagID)

suspect_resights <- all_resights_clean %>%
  filter(FlagID %in% c("FEY", "FEBK", "EY", "FEDP", "NA", "CB", "FEGY"))

all_resights_clean <- all_resights_clean %>% 
  filter(!FlagID %in% c("FEY", "FEBK", "EY", "FEDP", "CB", "FEGY"))

# remove resights with no banding records (except Argentina and Brazil)

# load banding records
banding_records <- read_excel("./data/banding_records.xlsx", sheet = "REKNcaptures")

band_numbers <- banding_records %>% 
  select(MetalID) %>% 
  distinct() %>% pull()

# separate out brazil and argentina - don't have the banding records for them so keep all
brar_resights <- all_resights_clean %>% 
  filter(FlagID %in% c("FDB", "FEDB", "FO", "FEO"))

other_resights <- all_resights_clean %>% 
  filter(!FlagID %in% c("FDB", "FEDB", "FO", "FEO"))

# remove resights that don't have a banding record
other_resights <- other_resights %>% 
  filter(MetalID %in% band_numbers)

# add brazil and argentina resights back in
all_resights_clean <- bind_rows(other_resights, brar_resights)

# format CC data
cc_data$ResightDate <- mdy(cc_data$ResightDate)

cc_data <- cc_data %>% 
  mutate(year = year(ResightDate))

# set map limits
xmin <- min(all_resights_clean$Longitude) - 5
xmax <- max(all_resights_clean$Longitude) + 5
ymin <- min(all_resights_clean$Latitude) - 5
ymax <- max(all_resights_clean$Latitude) + 5

# plot map
all_resights_clean$ProjectID <- as.factor(all_resights_clean$ProjectID)

resights_map <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = all_resights_clean,
             aes(x = Longitude, y = Latitude, colour = ProjectID),
             alpha = 0.5) +
  scale_color_viridis_d() +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

png(filename = paste0("figures/included-projects-map.png"),
    width=5, height=8, units="in", res=600)

print(resights_map)

dev.off()

# find date ranges
date_ranges <- all_resights_clean %>% 
  mutate(JulianDay = yday(ResightDate)) %>% 
  group_by(ProjectID) %>% 
  summarize(minday = min(JulianDay),
            maxday = max(JulianDay)) %>% 
  ungroup()

# check for Jim's proofed CC data - same as bandedbirds data
cc_dates <- cc_data %>% 
  mutate(JulianDay = yday(ResightDate)) %>% 
  group_by(ProjectID) %>% 
  summarize(minday = min(JulianDay),
            maxday = max(JulianDay)) %>% 
  ungroup()

# check date ranges for each SE project
se_dates <- all_resights_clean %>% 
  filter(ProjectID %in% c(41, 69, 15, 27)) %>% 
  mutate(JulianDay = yday(ResightDate)) %>% 
  select(ProjectID, JulianDay) %>% 
  distinct()

date_counts <- all_resights_clean %>% 
  filter(!ProjectID == 13) %>% 
  mutate(JulianDay = yday(ResightDate)) %>% 
  group_by(ProjectID, JulianDay) %>% 
  mutate(n = n()) %>% 
  select(ProjectID, JulianDay, n) %>% 
  distinct() %>% arrange(ProjectID, JulianDay)

# test encounter histories
all_resights_clean$ProjectID <- as.character(all_resights_clean$ProjectID)

all_resights_clean <- all_resights_clean %>% 
  filter(!ProjectID == 13) %>% 
  mutate(Region = case_when(ProjectID %in% c(3, 6, 8) ~ "DB",
                            ProjectID == 22 ~ "JB",
                            ProjectID == 21 ~ "MI",
                            ProjectID == 31 ~ "CC",
                            ProjectID %in% c(15, 27, 41, 69) ~ "SE",
                            ProjectID == 87 ~ "BR"))

set_periods <- tibble(Period = c(1:39),
                      Start = c("2009-04-01", "2009-07-01", "2009-10-01", "2009-11-16",
                                "2010-04-01", "2010-07-01", "2010-10-01", "2010-11-16",
                                "2011-04-01", "2011-07-01", "2011-10-01", "2011-11-16",
                                "2012-04-01", "2012-07-01", "2012-10-01", "2012-11-16",
                                "2013-04-01", "2013-07-01", "2013-10-01", "2013-11-16",
                                "2014-04-01", "2014-07-01", "2014-10-01", "2014-11-16",
                                "2015-04-01", "2015-07-01", "2015-10-01", "2015-11-16",
                                "2016-04-01", "2016-07-01", "2016-10-01", "2016-11-16",
                                "2017-04-01", "2017-07-01", "2017-10-01", "2017-11-16",
                                "2018-04-01", "2018-07-01", "2018-10-01"),
                      End = c("2009-05-30", "2009-09-30", "2009-11-15", "2010-03-31",
                              "2010-05-30", "2010-09-30", "2010-11-15", "2011-03-31",
                              "2011-05-30", "2011-09-30", "2011-11-15", "2012-03-31",
                              "2012-05-30", "2012-09-30", "2012-11-15", "2013-03-31",
                              "2013-05-30", "2013-09-30", "2013-11-15", "2014-03-31",
                              "2014-05-30", "2014-09-30", "2014-11-15", "2015-03-31",
                              "2015-05-30", "2015-09-30", "2015-11-15", "2016-03-31",
                              "2016-05-30", "2016-09-30", "2016-11-15", "2017-03-31",
                              "2017-05-30", "2017-09-30", "2017-11-15", "2018-03-31",
                              "2018-05-30", "2018-09-30", "2018-11-15"))

set_periods$Start <- as_date(set_periods$Start)
set_periods$End <- as_date(set_periods$End)

# remove bandedbirds CC data and add Jim's proofed CC data
cc_resights <- cc_data %>% 
  mutate(Region = "CC") %>% 
  select(ResightDate, Region, BirdID)

all_resights_clean <- all_resights_clean %>% 
  filter(!Region == "CC") %>% 
  select(ResightDate, Region, BirdID)

all_resights_clean <- bind_rows(all_resights_clean, cc_resights)

# assign each resight to a time period
resights_enchist_full <- all_resights_clean %>% 
  distinct()

# get number of resights per state
resights_state <- resights_enchist_full %>% 
  group_by(Region) %>% 
  summarize(n = n())

resights_enchist_full <- fuzzy_left_join(resights_enchist_full, set_periods,
                                    by = c("ResightDate" = "Start",
                                           "ResightDate" = "End"),
                                    match_fun = list(`>=`, `<=`))

# back to formatting for encounter histories
resights_enchist <- resights_enchist_full %>% 
  arrange(ResightDate) %>% 
  filter(!is.na(Period)) %>% 
  select(Period, Region, BirdID) %>% 
  distinct() %>% 
  arrange(Period)

# dataframe to look up resight dates of duplicate states
resights_check <- resights_enchist_full %>% 
  filter(!is.na(Period)) %>% 
  distinct() %>% 
  arrange(BirdID, ResightDate)

# pull out birds detected in >1 state in a period
check <- resights_enchist %>%
  group_by(Period, BirdID) %>%
  count() %>% filter(n > 1) %>% 
  ungroup()

duplicates <- check %>% 
  select(BirdID) %>% 
  distinct() %>% pull()

# create dummy df to ensure all periods included
dummy <- tibble(Period = 1:38,
                Region = "XX",
                BirdID = 99999)

resights_enchist_dummy <- bind_rows(resights_enchist, dummy)

# create encounter history
enchist <- resights_enchist_dummy %>% 
  filter(!BirdID %in% duplicates) %>% 
  distinct() %>% 
  mutate(Region = str_replace(Region, "DB", "1")) %>%
  mutate(Region = str_replace(Region, "JB", "2")) %>%
  mutate(Region = str_replace(Region, "MI", "3")) %>%
  mutate(Region = str_replace(Region, "CC", "4")) %>%
  mutate(Region = str_replace(Region, "SE", "5")) %>%
  mutate(Region = str_replace(Region, "BR", "6")) %>%
  arrange(Period) %>% 
  tidyr::pivot_wider(id_cols = BirdID, names_from = Period, values_from = Region) %>% 
  filter(!BirdID == 99999)
  #column_to_rownames(var = "BirdID")

# get number of rekn resighted per period
apply(enchist, 2, function(x) length(which(!is.na(x))))

sum(apply(enchist, 2, function(x) length(which(!is.na(x)))))

# convert to numeric matrix
enchist <- sapply(enchist, as.numeric)
enchist[is.na(enchist)] <- 0

# get regions for each period
enchist_df <- as.data.frame(enchist)
regions <- lapply(enchist_df, unique)

# get number resights per state
resights_state <- resights_enchist %>% 
  filter(!BirdID %in% duplicates) %>%
  group_by(Region) %>% 
  summarize(n = n())

# export encounter histories
write.csv(enchist, file = "./processed-data/rekn-multistate-enchist.csv", row.names = FALSE)
saveRDS(enchist, file = "./processed-data/rekn-multistate-enchist.rds")

# summarize number observations per state per year
state_year_sum <- resights_enchist_dummy %>% 
  group_by(Period, Region) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = Period, values_from = n) %>% 
  filter(!Region == "XX") %>% 
  ungroup() %>% 
  replace(is.na(.), 0) %>% 
  arrange(factor(Region, levels = c("DB", "JB", "MI", "CC", "SE", "BR")))

write.csv(state_year_sum, file = "./processed-data/resights-per-period.csv", row.names = FALSE)

# convert encounter histories to m-array
# function to create multistate m-array
marray <- function(CH, unobs = 0){ # unobs = number of unobservable states
  n.states <- max(CH) + unobs
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  out <- matrix(0, ncol = n.states*(n.occasions-1)+1, nrow = n.states*(n.occasions-1))
  
  # remove capture histories of individuals marked in last occasion
  get.first <- function(x) min(which(x!=0))
  first <- apply(CH, 1, get.first)
  last.only <- which(first==n.occasions)
  if (length(last.only) > 0) CH <- CH[-last.only,]
  
  # create m-array
  for (i in 1:dim(CH)[1]){
    cap.occ <- which(CH[i,]!=0)
    state <- CH[i,cap.occ]
    if (length(state) == 1) {
      out[state[1] + n.states*(cap.occ[1]-1), n.states*(n.occasions-1)+1] <- 
        out[state[1] + n.states*(cap.occ[1]-1), n.states*(n.occasions-1)+1] + 1
    }
    
    if (length(state) > 1) {
      for (t in 2:length(cap.occ)){
        out[(cap.occ[t-1]-1)*n.states + state[t-1], (cap.occ[t]-2)*n.states + state[t]] <-
          out[(cap.occ[t-1]-1)*n.states + state[t-1], (cap.occ[t]-2)*n.states + state[t]] + 1
      } # t
      
      if (max(cap.occ) < n.occasions){
        out[(cap.occ[t]-1)*n.states + state[t], n.states*(n.occasions-1)+1] <-
          out[(cap.occ[t]-1)*n.states + state[t], n.states*(n.occasions-1)+1] + 1
      } # i
    }
  }
  
  return(out)
  
}

enchist <- enchist[,2:39]

marr <- marray(enchist)

write.csv(marr, file = "./processed-data/rekn-multistate-marray.csv", row.names = FALSE)
saveRDS(marr, file = "./processed-data/rekn-multistate-marray.rds")

# create single-state encounter histories for CJS seasonal survival
resights_enchist_cjs <- resights_enchist_full %>% 
  arrange(ResightDate) %>% 
  filter(!is.na(Period)) %>% 
  select(Period, BirdID) %>% 
  distinct() %>% 
  arrange(Period)

# create dummy df to ensure all periods included
dummy <- tibble(Period = 1:38,
                BirdID = 99999)

resights_enchist_cjs_dummy <- bind_rows(resights_enchist_cjs, dummy) %>% 
  mutate(Seen = 1)

# create encounter history
enchist_cjs <- resights_enchist_cjs_dummy %>% 
  distinct() %>% 
  arrange(Period) %>% 
  tidyr::pivot_wider(id_cols = BirdID, names_from = Period, values_from = Seen) %>% 
  replace(is.na(.), 0) %>%
  filter(!BirdID == 99999)
#column_to_rownames(var = "BirdID")

# get number of rekn resighted per period
apply(enchist_cjs, 2, function(x) length(which(!is.na(x))))

sum(apply(enchist_cjs, 2, function(x) length(which(!is.na(x)))))

# export encounter histories
write.csv(enchist_cjs, file = "./processed-data/rekn-cjs-enchist.csv", row.names = FALSE)
saveRDS(enchist_cjs, file = "./processed-data/rekn-cjs-enchist.rds")

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  nind <- dim(CH)[1]
  
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

enchist_cjs <- enchist_cjs[,2:39]

# remove capture histories of individuals marked in last occasion
get.first <- function(x) min(which(x!=0))
first <- apply(enchist_cjs, 1, get.first)
last.only <- which(first==n.occasions)
if (length(last.only) > 0) enchist_cjs <- enchist_cjs[-last.only,]

marr <- marray(enchist_cjs)

write.csv(marr, file = "./processed-data/rekn-cjs-marray.csv", row.names = FALSE)
saveRDS(marr, file = "./processed-data/rekn-cjs-marray.rds")
