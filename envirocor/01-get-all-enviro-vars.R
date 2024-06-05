# get roms covariates from pacea-data-main and using their spatial_average function
# updated here with a default polygon for whole outer coast because I couldn't install the package
# remotes::install_github("pbs-assess/pacea") ## not working so retrieved files from github
library(tidyverse)
theme_set(ggsidekick::theme_sleek())


# Trim grid to depths of interest ----
# depth ranges for Love 2011
species_min <- 20
species_max <- 450
spawn_min <- 300
spawn_max <- 450

# depth range from our survey data encompassing 95% of this species’ biomass
summer_min <- 50
summer_max <- 200

# this uses a new_grid created with 00-pacea-grid.R
load("data-generated/grid.rda")

sp_grid <- new_grid |> filter(depth_min >= species_min & depth_max <= species_max & max_depth < 0)
# plot(sp_grid)
spawn_grid <- new_grid |> filter(depth_min >= spawn_min & depth_max <= spawn_max & max_depth < 0)
# plot(spawn_grid)
summer_grid <- new_grid |> filter(depth_min >= summer_min & depth_max <= summer_max & max_depth < 0)
# plot(summer_grid)


# Calculate spatial average from pacea objects and defined by these grids ----
# default area is roughly whole coast
spatial_average <- function(pacea_st_obj,
                            area = list(matrix(c(-134, -131, -128, -124, -125.5, -134, -134,
                                                 54.4, 54.4, 50.5, 48, 47.5, 52, 54.4),
                                               ncol = 2))
){
  stopifnot("pacea_st_obj must be of class pacea_st" =
              ("pacea_st" %in% class(pacea_st_obj)))

  if(is(area, "sfc_POLYGON")) {
    area_sf <- area |> sf::st_transform(crs = 3005)
  } else {
  # convert area to a simple features object, with the correct projection
  area_sf <- sf::st_sfc(sf::st_polygon(area),
                        crs = 4326) %>%
    sf::st_as_sf() %>%
    sf::st_transform(crs = 3005)
  }
  # this filters to just the required area
  obj_area <- pacea_st_obj[area_sf, ]

  obj_area_drop <- sf::st_drop_geometry(obj_area) %>%
    as_tibble()

  avg <- colMeans(obj_area_drop)

  obj_area_tib <- tibble::tibble(value = avg)

  obj_area_tib$year <- as.numeric(substr(names(avg),
                                         1,
                                         4))
  obj_area_tib$month <- as.numeric(substr(names(avg),
                                          6,
                                          7))

  obj_area_tib <- dplyr::relocate(obj_area_tib,
                                  year,
                                  month)
  obj_area_tib
}


load("../../pacea-data-main/data/bccm_primaryproduction_01.rds")
# pp <- spatial_average(bccm_primaryproduction, area = area)
pp <- spatial_average(bccm_primaryproduction, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_primary_production.rds")
pp <- spatial_average(bccm_primaryproduction, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_primary_production_spawn.rds")
pp <- spatial_average(bccm_primaryproduction, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_primary_production_summer.rds")

load("../../pacea-data-main/data/bccm_phytoplankton_01.rds")
# pp <- spatial_average(bccm_phytoplankton, area = area)
pp <- spatial_average(bccm_phytoplankton, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_phytoplankton.rds")
pp <- spatial_average(bccm_phytoplankton, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_phytoplankton_spawn.rds")
pp <- spatial_average(bccm_phytoplankton, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_phytoplankton_summer.rds")

load("../../pacea-data-main/data/bccm_bottom_oxygen_01.rds")
# pp <- spatial_average(bccm_bottom_oxygen, area = area)
pp <- spatial_average(bccm_bottom_oxygen, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_oxygen.rds")
pp <- spatial_average(bccm_bottom_oxygen, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_oxygen_spawn.rds")
pp <- spatial_average(bccm_bottom_oxygen, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_oxygen_summer.rds")

load("../../pacea-data-main/data/bccm_surface_oxygen_01.rds")
# pp <- spatial_average(bccm_surface_oxygen, area = area)
pp <- spatial_average(bccm_surface_oxygen, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_oxygen.rds")
pp <- spatial_average(bccm_surface_oxygen, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_oxygen_spawn.rds")
pp <- spatial_average(bccm_surface_oxygen, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_oxygen_summer.rds")

load("../../pacea-data-main/data/bccm_bottom_temperature_01.rds")
# pp <- spatial_average(bccm_bottom_temperature, area = area)
pp <- spatial_average(bccm_bottom_temperature, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_temperature.rds")
pp <- spatial_average(bccm_bottom_temperature, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_temperature_spawn.rds")
pp <- spatial_average(bccm_bottom_temperature, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_temperature_summer.rds")

load("../../pacea-data-main/data/bccm_surface_temperature_01.rds")
# pp <- spatial_average(bccm_surface_temperature, area = area)
pp <- spatial_average(bccm_surface_temperature, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_temperature.rds")
pp <- spatial_average(bccm_surface_temperature, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_temperature_spawn.rds")
pp <- spatial_average(bccm_surface_temperature, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_temperature_summer.rds")


load("../../pacea-data-main/data/bccm_surface_salinity_01.rds")
# pp <- spatial_average(bccm_surface_salinity, area = area)
pp <- spatial_average(bccm_surface_salinity, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_salinity.rds")
pp <- spatial_average(bccm_surface_salinity, area = spawn_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_salinity_spawn.rds")
pp <- spatial_average(bccm_surface_salinity, area = summer_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_salinity_summer.rds")

load("../../pacea-data-main/data/bccm_surface_ph_01.rds")
# pp <- spatial_average(bccm_surface_salinity, area = area)
pp <- spatial_average(bccm_surface_pH, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_surface_pH.rds")

load("../../pacea-data-main/data/bccm_bottom_ph_01.rds")
# pp <- spatial_average(bccm_surface_salinity, area = area)
pp <- spatial_average(bccm_bottom_pH, area = sp_grid$geometry)
saveRDS(pp, "data-generated/cw_bottom_pH.rds")



# Get mean for a specific set of months and standardize ----
# contains filter function
source("envirocor/utils.R")

load("data-generated/npgo.rda")
load("data-generated/oni.rda")
load("data-generated/pdo.rda")
load("data-generated/soi.rda")
load("data-generated/npi_monthly.rda")

# these are now all potentially occupied depths year-round
pp_monthly <- readRDS("data-generated/cw_primary_production.rds")
pt_monthly <- readRDS("data-generated/cw_phytoplankton.rds")
sst_monthly <- readRDS("data-generated/cw_surface_temperature.rds")
tob_monthly <- readRDS("data-generated/cw_bottom_temperature.rds")
do_monthly <- readRDS("data-generated/cw_bottom_oxygen.rds")
so2_monthly <- readRDS("data-generated/cw_surface_oxygen.rds")
ssa_monthly <- readRDS("data-generated/cw_surface_salinity.rds")
ph_monthly <- readRDS("data-generated/cw_bottom_pH.rds")

npgoA <- npgo |> rename(value = anomaly) |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "NPGO")
soiA <- soi |> rename(value = anomaly) |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "SOI")
oniA <- oni |> mutate(value = anomaly) |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "ENSO")
npiA <- npi_monthly |> filter_months(c(1,2,3), "NPI")
pdoA <- pdo |> rename(value = anomaly) |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "PDO")
ppA <- pp_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Primary production")
ptA <- pt_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Phytoplankton")
sstA <- sst_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Sea surface temperature")
tobA <- tob_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Sea floor temperature")
o2A <- do_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), type = "Sea floor O2")
so2A <- so2_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Sea surface O2")
ssaA <- ssa_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Sea surface salinity")
phA <- ph_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9,10,11,12), "Sea bottom pH")


# Explore covariates ----

(ev1 <- bind_rows(pdoA, oniA, npgoA) |>
   filter(year >= 1990) |>
   mutate(type = factor(type, levels = c("ENSO", "PDO", "NPGO"))) %>%
   ggplot() +
   geom_line(aes(year, value, colour = type), alpha = 0.7, linewidth = 1) +
   scale_colour_manual(values = RColorBrewer::brewer.pal(n = 8, name = "Paired")[c(1,2,8)]) +
   scale_x_continuous(limits = c(1990,2025), breaks = seq(1990, 2025, 5) ) +
   theme(
     axis.title = element_blank(),
     legend.justification=c(0, 1)) +
   labs(x = "Year", y = "Standardized index", colour = "Climate Indices"))

# ggsave("figs/climate-indices.png", width = 4, height = 2)


(ev2 <- bind_rows(sstA, tobA, ppA, ptA, phA, o2A, ssaA) |>
    ggplot() +
    geom_line(aes(year, value, colour = type), alpha = 0.7, linewidth = 1) +
    scale_colour_manual(values = RColorBrewer::brewer.pal(n = 10, name = "Paired")[c(3,4,8,9,5,10,6)])  +
    scale_x_continuous(limits = c(1990,2025), breaks = seq(1990, 2025, 5) ) +
    theme(
      axis.title.y = element_blank(),
      legend.justification=c(0, 1)) +
    labs(x = "Year", y = "Standardized value", colour = "BCCM Variables"))
ev2


y_lab_big <- ggplot() +
  annotate(geom = "text", x = 1, y = 1, size = 4,
           colour = "grey30",
           label = "Standardized annual values", angle = 90) +
  coord_cartesian(clip = "off")+
  theme_void()

y_lab_big + (ev1/ev2) + patchwork::plot_layout(width = c(0.1,1))

ggsave("envirocor/figs/ev-indices.png", width =7.5, height = 4.5)


# Plot correlations -----
library(GGally)

dw <- bind_rows(pdoA, npgoA, sstA, tobA, ppA, ptA, o2A) |> select(year, type, value) |>
  pivot_wider(names_from = type, values_from = value)

ggpairs(dw, columns = c(2:8),
        upper = list(continuous = wrap(cor_func, method = 'spearman', symbol = expression('\u03C1 ='))),
        progress = FALSE)


ggsave("envirocor/figs/clim-variables-correlations.png", width = 11, height = 11)


