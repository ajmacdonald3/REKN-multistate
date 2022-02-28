################################################################################
# PREPARE BBL RESIGHTS
#
################################################################################

library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(nngeo)
library(viridis)
library(leaflet)

# load data
#raw_data <- read_csv("./data/BBL_REKNresights.csv")
#saveRDS(raw_data, file = "./data/BBL_REKNresights.rds")

bbl_raw <- readRDS("./data/BBL_REKNresights.rds")

# filter out years and select needed columns
bbl_resights <- bbl_raw %>% 
  filter(ENCOUNTER_YEAR %in% (2009:2018)) %>% 
  select(B_MARKER_LONG_DESC, B_SPECIES_NAME, BAND_NUM, E_COMMENTS, E_LAT_DECIMAL_DEGREES, E_LON_DECIMAL_DEGREES,
         E_MARKER_LONG_DESC, E_SPECIES_NAME, ENCOUNTER_DATE, ENCOUNTER_DAY, ENCOUNTER_MONTH, ENCOUNTER_YEAR,
         ORIGINAL_BAND, OTHER_BANDS, REC_SOURCE, ENC_ERROR)

# load bandedbirds public resights
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

# set map limits
xmin <- min(bbl_resights$E_LON_DECIMAL_DEGREES) - 1
xmax <- max(bbl_resights$E_LON_DECIMAL_DEGREES) + 1
ymin <- min(bbl_resights$E_LAT_DECIMAL_DEGREES) - 1
ymax <- max(bbl_resights$E_LAT_DECIMAL_DEGREES) + 1

# plot map
resights_map <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = bbl_resights,
             aes(x = E_LON_DECIMAL_DEGREES, y = E_LAT_DECIMAL_DEGREES),
             colour = "#440154FF", alpha = 0.5) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)
  
# get cape cod shoreline shapefile
cape_cod <- st_read("map-data/shorelines/CapeCod_shorelines/CapeCod_shorelines.shp")

cape_cod <- cape_cod %>% 
  filter(Year_ == 2000)

cape_cod <- st_zm(cape_cod)

# transform to CRS that uses metres so can add buffers around polygons
cape_cod_buf <- st_transform(cape_cod, crs = 3395)

# add a 50 km buffer to each polygon
cape_cod_buf <- st_buffer(cape_cod_buf, dist = 50000)

# transform back to WGS 84 CRS
cape_cod_buf <- st_transform(cape_cod_buf, crs = 4326)

cape_cod_buf <- st_union(cape_cod_buf)

# add site label
cape_cod_df <- data.frame(site = "CC")

cape_cod_buf <- st_sf(cape_cod_df, geometry = cape_cod_buf)

xmin <- st_bbox(cape_cod_buf)[1] - 0.5
xmax <- st_bbox(cape_cod_buf)[3] + 1.5
ymin <- st_bbox(cape_cod_buf)[2] - 0.5
ymax <- st_bbox(cape_cod_buf)[4] + 0.5

# plot buffered coastline
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = bbl_resights,
             aes(x = E_LON_DECIMAL_DEGREES, y = E_LAT_DECIMAL_DEGREES),
             colour = "#440154FF", alpha = 0.5) +
  geom_point(data = public_resights,
             aes(x = Longitude, y = Latitude),
             colour = "#35B779FF", alpha = 0.5) +
  geom_sf(data = cape_cod_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# get south carolina shoreline shapefile
s_carolina <- st_read("map-data/shorelines/SouthCarolina_shorelines/sc2000.shp")

# transform to CRS that uses metres so can add buffers around polygons
s_carolina_buf <- st_transform(s_carolina, crs = 3395)

# add a 50 km buffer to each polygon
s_carolina_buf <- st_buffer(s_carolina_buf, dist = 50000)

# transform back to WGS 84 CRS
s_carolina_buf <- st_transform(s_carolina_buf, crs = 4326)

s_carolina_buf <- st_union(s_carolina_buf)

xmin <- st_bbox(s_carolina_buf)[1] - 0.5
xmax <- st_bbox(s_carolina_buf)[3] + 1.5
ymin <- st_bbox(s_carolina_buf)[2] - 0.5
ymax <- st_bbox(s_carolina_buf)[4] + 0.5

# plot buffered coastline
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = s_carolina_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# get georgia shoreline shapefile
georgia <- st_read("map-data/shorelines/Georgia_shorelines/ga1999.shp")

# transform to CRS that uses metres so can add buffers around polygons
georgia_buf <- st_transform(georgia, crs = 3395)

# add a 50 km buffer to each polygon
georgia_buf <- st_buffer(georgia_buf, dist = 50000)

# transform back to WGS 84 CRS
georgia_buf <- st_transform(georgia_buf, crs = 4326)

georgia_buf <- st_union(georgia_buf)

xmin <- st_bbox(georgia_buf)[1] - 0.5
xmax <- st_bbox(georgia_buf)[3] + 1.5
ymin <- st_bbox(georgia_buf)[2] - 0.5
ymax <- st_bbox(georgia_buf)[4] + 0.5

# plot buffered coastline
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = georgia_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# get florida shoreline shapefile
florida <- st_read("map-data/shorelines/Florida_shorelines/fl1999.shp")

# transform to CRS that uses metres so can add buffers around polygons
florida_buf <- st_transform(florida, crs = 3395)

# add a 50 km buffer to each polygon
florida_buf <- st_buffer(florida_buf, dist = 50000)

# transform back to WGS 84 CRS
florida_buf <- st_transform(florida_buf, crs = 4326)

florida_buf <- st_union(florida_buf)

xmin <- st_bbox(florida_buf)[1] - 0.5
xmax <- st_bbox(florida_buf)[3] + 1.5
ymin <- st_bbox(florida_buf)[2] - 0.5
ymax <- st_bbox(florida_buf)[4] + 0.5

# plot buffered coastline
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = florida_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# plot SC, GA, FL together
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = florida_buf, alpha = 0.3) +
  geom_sf(data = georgia_buf, alpha = 0.3) +
  geom_sf(data = s_carolina_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# combine polygons
se_states_buf <- st_union(s_carolina_buf, georgia_buf)
se_states_buf <- st_union(se_states_buf, florida_buf)

# add site label
se_states_df <- data.frame(site = "SE")

se_states_buf <- st_sf(se_states_df, geometry = se_states_buf)

# plot combined polygon
xmin <- st_bbox(se_states_buf)[1] - 0.5
xmax <- st_bbox(se_states_buf)[3] + 1.5
ymin <- st_bbox(se_states_buf)[2] - 0.5
ymax <- st_bbox(se_states_buf)[4] + 0.5

ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = bbl_resights,
             aes(x = E_LON_DECIMAL_DEGREES, y = E_LAT_DECIMAL_DEGREES),
             colour = "#440154FF", alpha = 0.5) +
  geom_point(data = public_resights,
             aes(x = Longitude, y = Latitude),
             colour = "#35B779FF", alpha = 0.5) +
  geom_sf(data = se_states_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# get brazil shoreline shapefile
brazil <- st_read("map-data/shorelines/Brazil_shorelines/cd01.shp")

#brazil <- st_as_sf(brazil)

# crop to just northern brazil coast
brazil <- st_crop(brazil, xmin = -51, xmax = -33.3, ymin = -6.5, ymax = 4.6)

# transform to CRS that uses metres so can add buffers around polygons
brazil_buf <- st_transform(brazil, crs = 3395)

# add a 50 km buffer to each polygon
brazil_buf <- st_buffer(brazil_buf, dist = 50000)

# transform back to WGS 84 CRS
brazil_buf <- st_transform(brazil_buf, crs = 4326)

brazil_buf <- st_union(brazil_buf)

brazil_df <- data.frame(site = "BR")

brazil_buf <- st_sf(brazil_df, geometry = brazil_buf)

xmin <- st_bbox(brazil_buf)[1] - 0.5
xmax <- st_bbox(brazil_buf)[3] + 1.5
ymin <- st_bbox(brazil_buf)[2] - 0.5
ymax <- st_bbox(brazil_buf)[4] + 0.5

# plot buffered coastline
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = bbl_resights,
             aes(x = E_LON_DECIMAL_DEGREES, y = E_LAT_DECIMAL_DEGREES),
             colour = "#440154FF", alpha = 0.5) +
  geom_point(data = public_resights,
             aes(x = Longitude, y = Latitude),
             colour = "#35B779FF", alpha = 0.5) +
  geom_sf(data = brazil_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# create sf dataframe of other sites as point locations
jb <- st_point(c(-80.571, 51.6557))
mi <- st_point(c(-63.818, 50.2117))
db <- st_point(c(-75.145, 39.0583))
ar <- st_point(c(-64.945, -40.7816))

add_sites <- st_sfc(jb, mi, db, ar, crs = 4326)

# transform to CRS that uses metres so can add buffers around polygons
add_sites_buf <- st_transform(add_sites, crs = 3395)

# add a 50 km buffer to each polygon
add_sites_buf <- st_buffer(add_sites_buf, dist = 50000)

# transform back to WGS 84 CRS
add_sites_buf <- st_transform(add_sites_buf, crs = 4326)

# add site labels
add_sites_df <- data.frame(site = c("JB", "MI", "DB", "AR"))

add_sites_buf <- st_sf(add_sites_df, geometry = add_sites_buf)

xmin <- st_bbox(add_sites_buf)[1] - 0.5
xmax <- st_bbox(add_sites_buf)[3] + 1.5
ymin <- st_bbox(add_sites_buf)[2] - 0.5
ymax <- st_bbox(add_sites_buf)[4] + 0.5

# plot buffered coastline
ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = bbl_resights,
             aes(x = E_LON_DECIMAL_DEGREES, y = E_LAT_DECIMAL_DEGREES),
             colour = "#440154FF", alpha = 0.5) +
  geom_point(data = public_resights,
             aes(x = Longitude, y = Latitude),
             colour = "#35B779FF", alpha = 0.5) +
  geom_sf(data = add_sites_buf, alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE)

# combine polygons
all_sites_buf <- rbind(cape_cod_buf, se_states_buf, brazil_buf, add_sites_buf)

# save as shapefile
st_write(all_sites_buf, "./map-data/regions/regions-polygons.shp")

xmin <- st_bbox(all_sites_buf)[1] - 0.5
xmax <- st_bbox(all_sites_buf)[3] + 1.5
ymin <- st_bbox(all_sites_buf)[2] - 0.5
ymax <- st_bbox(all_sites_buf)[4] + 0.5

sites <- c(AR = "Argentina", BR = "Brazil", CC = "Cape Cod", DB = "Delaware Bay",
           JB = "James Bay", MI = "Mingan", SE = "Southeast")

# plot buffered coastline
regions_map <- ggplot(data = world) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_point(data = bbl_resights,
             aes(x = E_LON_DECIMAL_DEGREES, y = E_LAT_DECIMAL_DEGREES),
             colour = "#440154FF", alpha = 0.5) +
  geom_point(data = public_resights,
             aes(x = Longitude, y = Latitude),
             colour = "#35B779FF", alpha = 0.5) +
  geom_sf(data = all_sites_buf, aes(colour = site, fill = site),
          alpha = 0.3) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  scale_color_viridis(discrete = TRUE, labels = sites) +
  scale_fill_viridis(discrete = TRUE, labels =  sites)

# leaflet interactive map
pal_fun <- colorFactor(viridis(7), all_sites_buf$site)

int_resight_map <- leaflet(all_sites_buf) %>% 
  addCircleMarkers(data = bbl_resights,
                   lng = ~E_LON_DECIMAL_DEGREES, lat = ~E_LAT_DECIMAL_DEGREES,
    radius = 3,
    color = "blue",
    stroke = FALSE, fillOpacity = 0.5
  ) %>% 
  addCircleMarkers(data = public_resights,
                   #lng = ~E_LON_DECIMAL_DEGREES, lat = ~E_LAT_DECIMAL_DEGREES,
                   radius = 3,
                   color = "red",
                   stroke = FALSE, fillOpacity = 0.5
  ) %>% 
  addPolygons(stroke = FALSE,
              fillColor = ~pal_fun(site), fillOpacity = 0.3) %>% 
  addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}') %>% 
  addLegend("bottomright",
            colors = viridis(7),
            labels = sites) %>% 
  addLegend("bottomleft",
            colors = c("blue", "red"),
            labels = c("BBL", "BandedBirds"))

htmlwidgets::saveWidget(int_resight_map, file = "./figures/int_resight_map.html")

# filter out resights within study regions
bbl_resights_sf <- st_as_sf(bbl_resights, coords = c("E_LON_DECIMAL_DEGREES", "E_LAT_DECIMAL_DEGREES"),
                            crs = 4326)

bbl_subset <- st_intersection(bbl_resights_sf, all_sites_buf)

bandedbirds_resights_sf <- st_as_sf(public_resights, coords = c("Longitude", "Latitude"),
                            crs = 4326)

bandedbirds_subset <- st_intersection(bandedbirds_resights_sf, all_sites_buf)

leaflet(all_sites_buf) %>% 
  addCircleMarkers(data = bbl_subset,
                   #lng = ~E_LON_DECIMAL_DEGREES, lat = ~E_LAT_DECIMAL_DEGREES,
                   radius = 3,
                   color = "blue",
                   stroke = FALSE, fillOpacity = 0.5
  ) %>% 
  addCircleMarkers(data = bandedbirds_subset,
                   #lng = ~E_LON_DECIMAL_DEGREES, lat = ~E_LAT_DECIMAL_DEGREES,
                   radius = 3,
                   color = "red",
                   stroke = FALSE, fillOpacity = 0.5
  ) %>% 
  addPolygons(stroke = FALSE,
              fillColor = ~pal_fun(site), fillOpacity = 0.3) %>% 
  addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}') %>% 
  addLegend("bottomright",
            colors = viridis(7),
            labels = sites) %>% 
  addLegend("bottomleft",
            colors = c("blue", "red"),
            labels = c("BBL", "BandedBirds"))

# save subsetted data
saveRDS(bbl_subset, "./data/BBL_resights_subset.rds")
writexl::write_xlsx(bbl_subset, "./data/BBL_resights_subset.xlsx")

saveRDS(bandedbirds_subset, "./data/bandedbirds_resights_subset.rds")
writexl::write_xlsx(bandedbirds_subset, "./data/bandedbirds_resights_subset.xlsx")


