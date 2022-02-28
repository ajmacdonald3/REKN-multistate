################################################################################
# MAP PUBLIC RESIGHTS
#
################################################################################

library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(nngeo)

# load data
raw_data <- read_excel("./data/bandedbirds_publicresights.xlsx")
saveRDS(raw_data, file = "./data/bandedbirds_publicresights.rds")

public_resights <- readRDS("./data/bandedbirds_publicresights.rds")

public_resights<- public_resights %>% 
  filter(!is.na(Longitude)) %>% 
  filter(!is.na(Latitude)) %>% 
  filter(!Longitude < (-1000)) %>% 
  mutate(Longitude = ifelse(Longitude > 0 & Latitude > 0, (Longitude*(-1)), Longitude))

# map data
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
lakes <- ne_load(type = "lakes", scale = "medium", category = "physical",
                 returnclass = "sf",
                 destdir = "./map-data/lakes")
states <- ne_states(country = c("United States of America", "Argentina", "Brazil",
                                   "Uruguay", "Chile"), returnclass = "sf")

# set map limits
xmin <- min(public_resights$Longitude)
xmax <- max(public_resights$Longitude)
ymin <- min(public_resights$Latitude) - 1
ymax <- max(public_resights$Latitude) + 1

# plot map
resights_map <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = public_resights,
             aes(x = Longitude, y = Latitude),
             colour = "#440154FF", alpha = 0.5)

# convert to spatial data
public_resights_sf <- public_resights %>% 
  filter(Longitude > (-104)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"),
           crs = 4326)

world_sf <- world %>% 
  st_transform(crs = 4326)

# find points within polygons
resights_country <- st_join(public_resights_sf, world_sf, join = st_within)

resights_country_count <- count(as_tibble(resights_country), sovereignt)

resights_country <- resights_country %>% 
  filter(!(sovereignt %in% c("Australia", "Canada", "Costa Rica", "France",
                             "Mexico", "Panama", "Peru", "Venezuela")))

resights_map2 <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = resights_country,
          aes(colour = is.na(sovereignt)), alpha = 0.5) +
  scale_colour_manual(values = c("#440154FF", "#35B779FF"))

# find points within and near polygons
resights_country_buf <- st_join(public_resights_sf, world_sf, join = st_nn, k = 1, maxdist = 1000)

# region level
states_sf <- states %>% 
  st_transform(crs = 4326)

resights_state <- st_join(public_resights_sf, states_sf, join = st_within)

resights_state_count <- count(as_tibble(resights_state), name)

states_list <- c("Ceará", "Delaware", "Florida", "Georgia", "Massachusetts",
                 "New Jersey", "Pará", "Rio Grande do Sul", "Río Negro",
                 "Rocha", "Santa Catarina", "Santa Cruz", "São Paulo", "South Carolina")

states_sf <- states_sf %>% 
  filter(name %in% states_list)

resights_map3 <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = resights_state,
          aes(colour = is.na(name)), alpha = 0.5) +
  scale_colour_manual(values = c("#440154FF", "#35B779FF"))

resights_states_buf <- st_join(public_resights_sf, states_sf, join = st_nn, k = 1, maxdist = 5000)

resights_state_count2 <- count(as_tibble(resights_states_buf), name)

resights_map4 <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = resights_states_buf,
          aes(colour = is.na(name)), alpha = 0.5) +
  scale_colour_manual(values = c("#440154FF", "#35B779FF"))

resights_map4
