################################################################################
# PREPARE CJS ENCOUNTER HISTORIES
#
# Season 1 = Apr-Jun
# Season 2 = Jul-Sep
# Season 3 = Oct-Mar
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

# load data
all_data <- read_excel("./data/all_projects_resights.xlsx", sheet = "REKNresights")
jb_data <- read_excel("./data/jb_project_resights.xlsx", sheet = "Sheet1")
cc_data <- read_excel("./data/cc_project_resights.xlsx", sheet = "CapeCod_REKNresights_JEL")
br_data <- read_excel("./data/br_project_resights.xlsx", sheet = "proj87REKNresights")

all_resights <- bind_rows(all_data, jb_data)
all_resights <- bind_rows(all_resights, br_data)

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
  distinct() %>% 
  filter(ResightDate > "2009-03-31") %>% 
  mutate(Year = year(ResightDate)) %>% 
  mutate(Month = month(ResightDate))

set_periods <- tibble(Period = c(1:29),
                      Start = c("2009-04-01", "2009-07-01", "2009-10-01",
                                "2010-04-01", "2010-07-01", "2010-10-01",
                                "2011-04-01", "2011-07-01", "2011-10-01",
                                "2012-04-01", "2012-07-01", "2012-10-01",
                                "2013-04-01", "2013-07-01", "2013-10-01",
                                "2014-04-01", "2014-07-01", "2014-10-01",
                                "2015-04-01", "2015-07-01", "2015-10-01",
                                "2016-04-01", "2016-07-01", "2016-10-01",
                                "2017-04-01", "2017-07-01", "2017-10-01",
                                "2018-04-01", "2018-07-01"),
                      End = c("2009-06-30", "2009-09-30", "2010-03-31",
                              "2010-06-30", "2010-09-30", "2011-03-31",
                              "2011-06-30", "2011-09-30", "2012-03-31",
                              "2012-06-30", "2012-09-30", "2013-03-31",
                              "2013-06-30", "2013-09-30", "2014-03-31",
                              "2014-06-30", "2014-09-30", "2015-03-31",
                              "2015-06-30", "2015-09-30", "2016-03-31",
                              "2016-06-30", "2016-09-30", "2017-03-31",
                              "2017-06-30", "2017-09-30", "2018-03-31",
                              "2018-06-30", "2018-09-30"))

set_periods$Start <- as_date(set_periods$Start)
set_periods$End <- as_date(set_periods$End)

resights_enchist_full <- fuzzy_left_join(resights_enchist_full, set_periods,
                                         by = c("ResightDate" = "Start",
                                                "ResightDate" = "End"),
                                         match_fun = list(`>=`, `<=`))

# create single-state encounter histories for CJS seasonal survival
resights_enchist_cjs <- resights_enchist_full %>% 
  arrange(ResightDate) %>% 
  filter(!is.na(Period)) %>% 
  select(Period, BirdID) %>% 
  distinct() %>% 
  arrange(Period) %>% 
  mutate(Seen = 1)

# create encounter history
enchist_cjs <- resights_enchist_cjs %>% 
  distinct() %>% 
  arrange(Period) %>% 
  tidyr::pivot_wider(id_cols = BirdID, names_from = Period, values_from = Seen) %>% 
  replace(is.na(.), 0)
#column_to_rownames(var = "BirdID")

# get number of rekn resighted per period
apply(enchist_cjs, 2, function(x) length(which(!is.na(x))))

sum(apply(enchist_cjs, 2, function(x) length(which(!is.na(x)))))

# export encounter histories
write.csv(enchist_cjs, file = "./processed-data/rekn-cjs-enchist.csv", row.names = FALSE)
saveRDS(enchist_cjs, file = "./processed-data/rekn-cjs-enchist.rds")

# create INP file to check assumptions
# now do formatting of encounter history file for .inp file for MARK
enchist <- enchist_cjs

enchist$eh <- apply(enchist[2:ncol(enchist)],1,paste,collapse="") # concatenates encounter columns into eh
enchist[2:(ncol(enchist)-1)] <- NULL # drops individual encounter columns

# create commented tag
enchist$BirdID <- paste("/*", enchist$BirdID, "*/", sep=" ")

# sort by descending encounter histories
enchist <- enchist[order(enchist$eh,decreasing=TRUE),]

# tack on the frequency for the individual
enchist$end <- "1;"

write.table(enchist, file = "./processed-data/rekn-cjs-enchist.inp", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
