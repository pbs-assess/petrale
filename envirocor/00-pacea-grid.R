# add lat lons to new depth grid in PACea package
library(pacea)
library(sf)
library(tidyverse)
# grid26_depth comes with pacea but at the moment lacks lat lons
# plot(grid26_depth)
nd <- st_as_sf(grid26_depth)

nd_utm <- st_transform(nd, crs = 3156)
nd_utm$geometry <- st_centroid(nd_utm) %>%
  # since you want the centroids in a second geometry col:
  st_geometry()

nd_ll <- st_transform(nd, crs = 4326)
nd_ll$geometry <- st_centroid(nd_ll) %>%
  # since you want the centroids in a second geometry col:
  st_geometry()
nd_ll$x <- NULL

nd_ll_sf <- st_as_sf(nd_ll)

# plot(st_geometry(nd_ll))
# plot(st_geometry(nd_utm))
# plot(st_geometry(nd_utm), cex = 0.2, add = T)

ret <- sf::st_coordinates(nd)
ret <- tibble::as_tibble(ret)

sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") #&& inherits(sf::st_geometry(x),"sfc_POINT")
  )
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

nd_ll$x <- NULL
nd_ll_sf <- st_as_sf(nd_ll)
nd <- sfc_as_cols(nd_ll_sf, c("X", "Y")) |> rename(longitude = X, latitude = Y)
nd$geometry <- grid26_depth$x
nd$area <- st_area(nd)
new_grid <- nd

new_grid <- new_grid |>
  mutate(
    area_m2 = as.numeric(area),
    depth_m = -mean_depth,
    depth_min = -max_depth, depth_max = -min_depth
  )

save(new_grid, file = "data-generated/grid.rda")
load("data-generated/grid.rda")

offshore <- new_grid |> filter(area_m2 > 4000000, max_depth < 0)
plot(offshore)

inshore <- new_grid |> filter(area_m2 <= 4000000, max_depth < 0)
plot(inshore)

save(inshore, file = "data-generated/grid-inshore.rda")

