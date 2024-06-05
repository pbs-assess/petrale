# 1. Create a total biomass distribution model (to use for prediction grid)
# 2. Create prediction grid with density and depth
options(scipen = 100, digits = 4)

# gfenvirocor is currently a private repository that contains plotting and model refining functions
library(gfenvirocor)
library(tidyverse)
library(sdmTMB)
library(gfplot)
library(ggsidekick)
library(patchwork)

species <- "Petrale Sole"


options(scipen = 100, digits = 4)
theme_set(theme_sleek())
fig_height <- 4 * 2
fig_width <- 5 * 2

dens_model_name_long <- "depth, DOY, and survey type"

## these need to be in the global environment to allow up of update()
knot_distance <<- 20
set_priors <<- sdmTMBpriors(
  matern_s = pc_matern(
    range_gt = knot_distance,
    sigma_lt = 2
  ),
  matern_st = pc_matern(
    range_gt = knot_distance,
    sigma_lt = 2
  )
)

set_family <- delta_lognormal()
set_family2 <- delta_gamma()

dens_model_name0 <- "dln-1998-new"
dens_model_name <- "dln-1998-good-split"


# subset of surveys used for density models
surveys_included <- c(
  "HS MSA",
  "SYN HS", "SYN QCS",
  "SYN WCHG", "SYN WCVI",
  "MSSM QCS", "MSSM WCVI"
)

set_spatial <- "on"
set_spatiotemporal <- list("rw", "rw")

custom_maturity_code <- NULL
custom_length_threshold <- NULL


spp <- gsub(" ", "-", gsub("\\/", "-", tolower(species)))

# # Load data ----
dset <- readRDS("data-generated/all-sets-used.rds") %>%
  filter(species_common_name == tolower(species)) %>%
  filter(survey_abbrev %in% surveys_included) %>%
  distinct()

dsamp <- readRDS("data-generated/all-samples-used.rds") %>%
  filter(species_common_name == tolower(species)) %>%
  distinct()

check_for_duplicates <- dset[duplicated(dset$fishing_event_id), ]
if (nrow(check_for_duplicates) > 0) {
  stop(paste(species, "has duplicate event ids."))
}

check_for_duplicates2 <- dsamp[duplicated(dsamp$specimen_id), ]
if (nrow(check_for_duplicates2) > 0) {
  stop(paste(species, "has duplicate specimen ids."))
}

# I'm not using sample is so shouldn't be any unless something else can cause this
# test_event <- dset[ dset$fishing_event_id == 329270, ]

maturity_possible <- TRUE

sets_mat_m <- filter(dsamp, !is.na(maturity_code), !is.na(length), maturity_code != 0, sex == 1)
sets_mat_f <- filter(dsamp, !is.na(maturity_code), !is.na(length), maturity_code != 0, sex == 2)

#
if (min(length(unique(sets_mat_m$fishing_event_id)), length(unique(sets_mat_f$fishing_event_id))) >= 20) {
  set_sample_id_re <- TRUE
} else {
  set_sample_id_re <- FALSE
}

# Split by maturity ----
# use SD to choose between aggregating by year or survey first, then survey, then all
set_min_sample_number <- 6
mat_threshold <- 0.5

dss <- gfplot::split_catch_by_sex(dset, dsamp,
  split_by_weight = TRUE, # automatically switches to TRUE for default weight-based catch variable
  sample_id_re = set_sample_id_re, # used for maturity ogives
  immatures_pooled = TRUE,
  min_sample_number = set_min_sample_number,
  split_by_maturity = maturity_possible,
  custom_maturity_at = custom_maturity_code,
  custom_length_thresholds = custom_length_threshold,
  p_threshold = mat_threshold,
  plot = maturity_possible
)

dir.create(paste0("data-generated/split-catch-data/"), showWarnings = FALSE)


meaneffort1 <- dss$data %>%
  filter(group_name %in% c("Mature", "Females", "Mature females") &
    usability_code == 1) %>%
  mutate(
    doorspread_m = ifelse(doorspread_m == 0 | doorspread_m > 150, NA_real_, doorspread_m),
    doorspread_m = ifelse(survey_series_id == 68 & is.na(doorspread_m), mouth_width_m * 0.74 + 40, doorspread_m),
    speed_mpm = ifelse(speed_mpm == 0, NA_real_, speed_mpm)
  ) %>%
  group_by(year, survey_id, survey_series_id) %>%
  summarise(
    mean_doorspread = mean(doorspread_m, na.rm = TRUE),
    mean_speed = mean(speed_mpm, na.rm = TRUE)
  )

meaneffort2 <- dss$data %>%
  filter(group_name %in% c("Mature females") &
    usability_code == 1) %>%
  mutate(
    doorspread_m = ifelse(doorspread_m == 0, NA_real_, doorspread_m),
    speed_mpm = ifelse(speed_mpm == 0, NA_real_, speed_mpm)
  ) %>%
  group_by(survey_series_id) %>%
  summarise(
    n = n(),
    mean_series_doorspread = mean(doorspread_m, na.rm = TRUE),
    mean_series_speed = mean(speed_mpm, na.rm = TRUE)
  )

ds <- dss$data %>%
  left_join(meaneffort1) %>%
  left_join(meaneffort2) %>%
  mutate(
    DOY = as.numeric(strftime(time_deployed, format = "%j")),
    days_to_solstice = DOY - 172,
    fishing_event_id = as.factor(fishing_event_id),
    usability_name = paste(usability_code, "-", usability_desc),
    usability_name = fct_reorder(usability_name, usability_code),
    # correct outlier duration (probably only need for Hake survey), looks like date problem
    duration_min = ifelse(duration_min > 1000, 30, duration_min),
    # get rid of false 0 and huge outliers for doorspread
    doorspread_m = ifelse(doorspread_m == 0 | doorspread_m > 150, NA_real_, doorspread_m),
    mouth_width_m = ifelse(mouth_width_m == 0, NA_real_, mouth_width_m),
    # use survey or survey series means with missing doorspread and speed
    mean_doorspread = ifelse(is.na(mean_doorspread), mean_series_doorspread, mean_doorspread),
    doorspread_m = ifelse(is.na(doorspread_m), mean_doorspread, doorspread_m),
    mean_speed = ifelse(is.na(mean_speed), mean_series_speed, mean_speed),
    speed_mpm = ifelse(speed_mpm == 0 | is.na(speed_mpm), mean_speed, speed_mpm),
    log_depth = log(depth_m),
    log_depth_c = log_depth - 5, # mean and median for whole data set
    area_swept = ifelse(is.na(tow_length_m),
      doorspread_m * duration_min * speed_mpm,
      doorspread_m * tow_length_m
    )
  )

ds <- ds[, which(unlist(lapply(ds, function(x) !all(is.na(x))))), with = FALSE]

# Select what data to include ----
ds <- ds %>%
  filter(usability_code %in% c(0, 1, 22, 16, 6)) %>%
  # if speed recorded, it isn't too slow
  filter(!is.na(speed_mpm) & speed_mpm != 0 & speed_mpm >= 50) %>%
  # if time recorded, it was at least 10 min
  filter(!is.na(duration_min) & duration_min != 0 & duration_min >= 10) %>%
  # if tow length available it's at least 200m
  filter(is.na(tow_length_m) | tow_length_m > 500) %>%
  filter(area_swept > 0) %>%
  # early MSSM not reliable for small fish
  filter(survey_type != "MSSM <03" & year >= 1998)

ds <- ds %>% filter(!is.na(catch_weight))
ds <- ds %>% filter(!is.na(depth_m))
ds <- ds %>% filter(!is.na(area_swept))
ds <- ds %>% filter(!is.na(latitude))
ds <- ds %>% filter(!is.na(longitude))
ds$area_swept_km2 <- ds$area_swept / 1000000 # converts m2 to km2
ds$offset <- log(ds$area_swept_km2 * 100) # converts km2 to ha
ds$survey_type <- relevel(ds$survey_type, "SYN")


# Plot usabilities ----
ds %>%
  mutate(
    group_name = ifelse(is.na(group_catch_est), "Unsampled", group_name),
    group_catch_est = ifelse(group_name == "Unsampled", catch_weight, group_catch_est)
  ) %>%
  filter(group_name %in% c("Mature", "Females", "Mature females", "Unsampled")) %>%
  distinct() %>%
  ggplot() +
  geom_point(aes(area_swept, usability_name, size = log(catch_weight + 1)), alpha = 0.3) +
  scale_y_discrete(limits = rev) +
  facet_grid(~survey_type, scales = "free") +
  theme(axis.title.y = element_blank())
ggsave(paste0("envirocor/figs/all-usabilities-", spp, ".png"),
  width = 17, height = 3.5
)



unique(ds$survey_abbrev)
unique(ds$survey_type)
sort(unique(ds$year))

# check that offset doesn't contain NAs of Inf
range(ds$offset)
range(ds$area_swept)
mean(ds$offset)

ggplot(ds) +
  geom_histogram(aes(offset)) +
  facet_wrap(~year)

# Make grid ----
ds$X <- NULL
ds$Y <- NULL

d <- sdmTMB::add_utm_columns(ds,
  ll_names = c("longitude", "latitude"),
  utm_crs = 32609
)

dir.create(paste0("data-generated/density-data/"), showWarnings = FALSE)
saveRDS(d, paste0("data-generated/density-data/", spp, ".rds"))


# Set naming conventions ----
dir.create(paste0("data-generated/density-models/", dens_model_name, "/"))
dir.create(paste0("data-generated/density-models/", dens_model_name0, "/"))

dir0 <- paste0("data-generated/density-models/", dens_model_name0, "/total/")
dir1 <- paste0("data-generated/density-models/", dens_model_name, "/mat-fem/")
dir2 <- paste0("data-generated/density-models/", dens_model_name, "/mat-m/")
dir3 <- paste0("data-generated/density-models/", dens_model_name, "/imm/")

### folder names for all density models
dir.create(dir0, showWarnings = FALSE)
dir.create(dir1, showWarnings = FALSE)
dir.create(dir2, showWarnings = FALSE)
dir.create(dir3, showWarnings = FALSE)

### directories for all generated indices
dir.create(paste0("data-generated/density-index/"), showWarnings = FALSE)
dir.create(paste0("data-generated/density-index/", dens_model_name0), showWarnings = FALSE)
dir.create(paste0("data-generated/density-index/", dens_model_name), showWarnings = FALSE)
dir.create(paste0(
  "data-generated/density-index/", dens_model_name0,
  "/total/"
), showWarnings = FALSE)
dir.create(paste0(
  "data-generated/density-index/", dens_model_name,
  "/mat-fem/"
), showWarnings = FALSE)
dir.create(paste0(
  "data-generated/density-index/", dens_model_name,
  "/mat-m/"
), showWarnings = FALSE)
dir.create(paste0(
  "data-generated/density-index/", dens_model_name,
  "/imm/"
), showWarnings = FALSE)
dir.create(paste0("data-generated/density-split-ind/"), showWarnings = FALSE)

m0 <- paste0(spp, "-total-", dens_model_name0, "-", knot_distance, "-km")
m1 <- paste0(spp, "-mat-fem-", dens_model_name, "-", knot_distance, "-km")
m2 <- paste0(spp, "-mat-m-", dens_model_name, "-", knot_distance, "-km")
m3 <- paste0(spp, "-imm-", dens_model_name, "-", knot_distance, "-km")

fm <- paste0(dir0, m0, ".rds")
fmf <- paste0(dir1, m1, ".rds")
fmm <- paste0(dir2, m2, ".rds")
fmi <- paste0(dir3, m3, ".rds")

pfn <- paste0("data-generated/density-predictions/p-", m0, ".rds")
pmfn <- paste0("data-generated/density-predictions/p-", m1, ".rds")
pmfn2 <- paste0("data-generated/density-predictions/p-", m2, ".rds")
pifn <- paste0("data-generated/density-predictions/p-", m3, ".rds")

i0 <- paste0("data-generated/density-index/", dens_model_name0, "/total/i-", m0, ".rds")
i1 <- paste0("data-generated/density-index/", dens_model_name, "/mat-fem/i-", m1, ".rds")
i2 <- paste0("data-generated/density-index/", dens_model_name, "/mat-m/i-", m2, ".rds")
i3 <- paste0("data-generated/density-index/", dens_model_name, "/imm/i-", m3, ".rds")

# Generate synoptic grid for all modelled years ----
grid <- replicate_df(gfplot::synoptic_grid,
  time_name = "year",
  time_values = min(d$year):max(d$year)
) %>%
  mutate(
    fishing_event_id = as.factor(paste(d$fishing_event_id[1])),
    days_to_solstice = 0,
    log_depth = log(depth),
    log_depth_c = log_depth - 5, # mean and median for whole data set
    survey_type = as.factor("SYN")
  )

# Total density dataframe ----

# simplify dataframe to be included with models
d <- d %>% select(
  species_common_name,
  fishing_event_id,
  catch_weight,
  group_name, group_catch_est,
  n_events_sampled, n_fish_sampled,
  proportion, median_prop_ann,
  X, Y, latitude, longitude,
  year, days_to_solstice, DOY,
  survey_abbrev, survey_type,
  depth_m, log_depth, log_depth_c,
  area_swept_km2, offset,
  duration_min, speed_mpm,
  doorspread_m, tow_length_m,
  usability_name
)

d1 <- d %>%
  filter(group_name %in% c("Mature", "Females", "Mature females"))

# if(nrow(d1a)!= nrow(d1)) { print("")}

## check what years we have data for ----
## used to filter grid
survey_years <- d1 %>%
  select(survey_abbrev, year) %>%
  distinct() %>%
  mutate(survey = ifelse(survey_abbrev == "HS MSA", "SYN HS",
    ifelse(survey_abbrev == "MSSM QCS", "SYN QCS",
      ifelse(survey_abbrev %in% c("OTHER", "MSSM WCVI"), "SYN WCVI",
        survey_abbrev
      )
    )
  ))

all_years <- unique(d1$year)
extra_years <- sdmTMB:::find_missing_time(d1$year)

# check that my data and grid are on the same XY scale
range(grid$year)
range(d$year)
range(d$X)
range(grid$X)
range(d$Y)
range(grid$Y)

dp <- d1 %>% filter(catch_weight > 0)

hist(log(dp$catch_weight))
range(d$catch_weight)

rm(dsamp, dset, dss, dp, ds)

# Make mesh for total density ----

mesh <- make_mesh(d1, c("X", "Y"), cutoff = knot_distance)

# browser()
## plot mesh ----

plot_mesh <- function(
    mesh_obj = mesh,
    data_obj = d1,
    catch_var = "catch_weight",
    group = "Total") {
  g <- ggplot() +
    inlabru::gg(mesh_obj$mesh) +
    coord_fixed() +
    geom_point(aes(X, Y),
      shape = "x",
      size = 0.75,
      data = data_obj
    ) +
    geom_point(
      aes(X, Y,
        fill = .data[[catch_var]],
        colour = .data[[catch_var]],
        size = .data[[catch_var]],
        shape = survey_type
      ),
      data = filter(data_obj, .data[[catch_var]] != 0)
    ) +
    facet_wrap(~year, ncol = 5) +
    scale_fill_viridis_c(trans = "fourth_root_power") +
    scale_color_viridis_c(trans = "fourth_root_power") +
    scale_size_continuous(guide = NULL) +
    ggtitle(paste0(group, " biomass")) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
    )
  if (catch_var == "density_kgha") {
    g + labs(
      fill = "Density (kg/ha)",
      colour = "Density (kg/ha)",
      shape = "Survey"
    )
  } else {
    g
  }
}

d1$density_kgha <- d1$catch_weight / (d1$area_swept_km2 * 100) # converts to ha

plot_mesh(
  mesh_obj = mesh,
  data_obj = d1,
  catch_var = "density_kgha"
)

ggsave(paste0(
  "envirocor/figs/density-mesh-",
  spp,
  "-total-",
  dens_model_name,
  ".png"
), width = 9, height = 11)


d1 %>% ggplot() +
  geom_histogram(aes(depth_m)) +
  geom_histogram(aes(depth_m), fill = "red", data = filter(d1, catch_weight > 0)) +
  facet_wrap(~survey_type, scales = "free_y")

d1 %>% ggplot() +
  geom_histogram(aes(days_to_solstice)) +
  geom_histogram(aes(days_to_solstice), fill = "red", data = filter(d1, catch_weight > 0)) +
  facet_wrap(~survey_type, scales = "free_y")


# Start modeling ----

## total abundance model ----
if (!file.exists(fm)) {
  m <- sdmTMB(
    catch_weight ~ 1 + survey_type +
      poly(log_depth_c, 2) +
      poly(days_to_solstice, 2),
    offset = "offset",
    mesh = mesh,
    data = d1,
    spatial = set_spatial,
    spatiotemporal = set_spatiotemporal,
    share_range = FALSE,
    silent = FALSE,
    time = "year",
    extra_time = seq(min(d1$year), max(d1$year)),
    family = set_family,
    priors = set_priors
  )

  saveRDS(m, fm)

  if (!all(sanity(m, gradient_thresh = 0.005))) {
    m <- refine_model(m, alternate_family = set_family2, use_priors = set_priors)
  }
  saveRDS(m, fm)
} else {
  m <- readRDS(fm)
  m <- sdmTMB:::update_version(m)
  # browser()
  if (!all(sanity(m, gradient_thresh = 0.005))) {
    m <- refine_model(m, alternate_family = set_family2, use_priors = set_priors)
  }
  saveRDS(m, fm)
}

m
m$sd_report
tidy(m, "ran_pars", conf.int = TRUE, model = 1)
try(tidy(m, "ran_pars", conf.int = TRUE, model = 2))

if (is.null(extra_years)) {
  grid <- filter(grid, year %in% c(sort(unique(m$data$year))))
} else {
  grid <- filter(grid, year %in% sort(union(m$data$year, m$extra_time)))
}


if (file.exists(pfn) & file.exists(i0)) {
  p <- readRDS(pfn)
} else {
  p <- predict(m, re_form_iid = NA, newdata = grid, return_tmb_object = TRUE)

  map_density(p, pfn, variable = "density_trimmed") +
    labs(
      title = paste0("Total biomass (", paste(unique(m$data$survey_type), collapse = ", "), ")"),
      # subtitle = paste0("Variance explained:", TODO: r2_total$R2[1])
      fill = "Density (kg/ha)",
      colour = "Density (kg/ha)"
    ) +
    theme(
      legend.position = c(0.9, 0.1),
      legend.direction = "vertical", legend.box = "horizontal"
    )

  ggsave(paste0("envirocor/figs/density-map-", m0, ".png"),
    height = 10, width = fig_width
  )

  plot_index(p, species, "Total", dens_model_name0, i0) +
    ggtitle(paste0(species, ": total biomass (", dens_model_name, ")"))
}

if (!file.exists(i0)) {
  plot_index(p, species, "Total", dens_model_name0, i0) +
    ggtitle(paste0(species, ": total biomass (", dens_model_name, ")"))
}

# certain model formulation functions (like formula() and terms() ) grab the ENTIRE global environment
# need to purge previous models before building new ones to save HD space when saving models
set_spatial <- as.list(m[["spatial"]])
set_spatiotemporal <- as.list(m[["spatiotemporal"]])
set_family <- m$family
rm(m, p)

## mature female model ----
d2 <- d1 %>%
  filter(!is.na(group_catch_est))

mesh2 <- make_mesh(d2, c("X", "Y"), cutoff = knot_distance)

d2$density_kgha <- d2$group_catch_est / (d2$area_swept_km2 * 100)

plot_mesh(
  mesh_obj = mesh2,
  data_obj = d2,
  catch_var = "density_kgha",
  group = "Mature females"
)

ggsave(paste0(
  "envirocor/figs/density-mesh-",
  spp,
  "-mat-fem-",
  dens_model_name,
  ".png"
), width = 9, height = 11)




if (!file.exists(fmf)) {
  mf <- sdmTMB(
    group_catch_est ~ 1 + survey_type +
      poly(log_depth_c, 2) +
      poly(days_to_solstice, 2),
    offset = "offset",
    spatial = set_spatial,
    spatiotemporal = set_spatiotemporal,
    share_range = FALSE,
    time = "year",
    family = set_family,
    extra_time = seq(min(d2$year), max(d2$year)),
    priors = set_priors,
    silent = FALSE,
    mesh = mesh2, data = d2
  )

  saveRDS(mf, fmf)

  if (!all(sanity(mf, gradient_thresh = 0.005))) {
    mf <- refine_model(mf, alternate_family = set_family2, use_priors = set_priors)
  }
  saveRDS(mf, fmf)
} else {
  mf <- readRDS(fmf)
  mf <- sdmTMB:::update_version(mf)
  if (!all(sanity(mf, gradient_thresh = 0.005))) {
    mf <- refine_model(mf, alternate_family = set_family2, use_priors = set_priors)
  }
  saveRDS(mf, fmf)
}


if (file.exists(pmfn) & file.exists(i1)) {
  pf <- readRDS(pmfn)
} else {
  pf <- predict(mf,
    re_form_iid = NA, # only needed if random intercepts
    newdata = filter(grid, year %in% sort(union(mf$data$year, mf$extra_time))),
    return_tmb_object = TRUE
  )

  map_density(pf, pmfn) +
    labs(
      title = paste0("Mature female biomass (", paste(unique(mf$data$survey_type), collapse = ", "), ")"),
      # subtitle = paste0("Variance explained:", TODO: r2_total$R2[1])
      fill = "Density (kg/ha)",
      colour = "Density (kg/ha)"
    ) +
    theme(
      legend.position = c(0.9, 0.1),
      legend.direction = "vertical", legend.box = "horizontal"
    )

  ggsave(paste0("envirocor/figs/density-map-", m1, ".png"),
    height = 10, width = fig_width
  )
}

if (!file.exists(i1)) {
  plot_index(pf, species, "Mature female", dens_model_name, i1) +
    ggtitle(paste0(species, ": mature female biomass (", dens_model_name, ")"))
}

rm(mf, pf, d1, d2)


## mature male model ----

d2b <- d %>%
  filter(group_name %in% c("Mature males")) %>%
  filter(!is.na(group_catch_est))

mesh2b <- make_mesh(d2b, c("X", "Y"), cutoff = knot_distance)

d2b$density_kgha <- d2b$group_catch_est / (d2b$area_swept_km2 * 100)

plot_mesh(
  mesh_obj = mesh2b,
  data_obj = d2b,
  catch_var = "density_kgha",
  group = "Mature males"
)

ggsave(paste0(
  "envirocor/figs/density-mesh-",
  spp,
  "-mat-m-",
  dens_model_name,
  ".png"
), width = 9, height = 11)


if (!file.exists(fmm)) {
  mm <- sdmTMB(
    group_catch_est ~ 1 + survey_type +
      poly(log_depth_c, 2) +
      poly(days_to_solstice, 2),
    offset = "offset",
    spatial = set_spatial,
    spatiotemporal = set_spatiotemporal,
    share_range = FALSE,
    time = "year",
    family = set_family,
    priors = set_priors,
    extra_time = sdmTMB:::find_missing_time(d2b$year),
    mesh = mesh2b,
    data = d2b
  )

  saveRDS(mm, fmm)

  if (!all(sanity(mm, gradient_thresh = 0.005))) {
    mm <- refine_model(mm, alternate_family = set_family2, use_priors = set_priors)
  }
  saveRDS(mm, fmm)
} else {
  mm <- readRDS(fmm)
  mm <- sdmTMB:::update_version(mm)
  if (!all(sanity(mm, gradient_thresh = 0.005))) {
    mm <- refine_model(mm, alternate_family = set_family2, use_priors = set_priors)
  }
  saveRDS(mm, fmm)
}

if (file.exists(pmfn2) & file.exists(i2)) {
  pm <- readRDS(pmfn2)
} else {
  pm <- predict(mm,
    re_form_iid = NA,
    newdata = filter(grid, year %in% sort(union(mm$data$year, mm$extra_time))),
    return_tmb_object = TRUE
  )

  map_density(pm, pmfn2) +
    labs(
      title = paste0("Mature male biomass (", paste(unique(mm$data$survey_type), collapse = ", "), ")"),
      fill = "Density (kg/ha)",
      colour = "Density (kg/ha)"
    ) +
    theme(
      legend.position = c(0.9, 0.1),
      legend.direction = "vertical", legend.box = "horizontal"
    )

  ggsave(paste0("envirocor/figs/density-map-", m2, ".png"),
    height = 10, width = fig_width
  )
}

if (!file.exists(i2)) {
  plot_index(pm, species, "Mature male", dens_model_name, i2) +
    ggtitle(paste0(species, ": mature male biomass (", dens_model_name, ")"))
}

rm(mm, pm, d2b)

## immature model ----
d3 <- d %>%
  filter(group_name %in% c("Immature")) %>%
  filter(!is.na(group_catch_est))

mesh3 <- make_mesh(d3, c("X", "Y"), cutoff = knot_distance)

d3$density_kgha <- d3$group_catch_est / (d3$area_swept_km2 * 100)

plot_mesh(
  mesh_obj = mesh3,
  data_obj = d3,
  catch_var = "density_kgha",
  group = "Immatures"
)
ggsave(paste0(
  "envirocor/figs/density-mesh-",
  spp,
  "-imm-",
  dens_model_name,
  ".png"
), width = 9, height = 11)


if (!file.exists(fmi)) {
  mi <- sdmTMB(
    group_catch_est ~ 1 + survey_type +
      poly(log_depth_c, 2) +
      poly(days_to_solstice, 2),
    offset = "offset",
    spatial = set_spatial,
    spatiotemporal = set_spatiotemporal,
    share_range = FALSE,
    time = "year",
    family = set_family,
    extra_time = sdmTMB:::find_missing_time(d3$year),
    priors = set_priors,
    # priors = set_priors,
    silent = FALSE,
    mesh = mesh3, data = d3
  )

  saveRDS(mi, fmi)
  if (!all(sanity(mi, gradient_thresh = 0.005))) {
    mi <- refine_model(mi,
      alternate_family = set_family2,
      use_priors = set_priors
    )
  }
  saveRDS(mi, fmi)
} else {
  mi <- readRDS(fmi)

  if (!all(sanity(mi, gradient_thresh = 0.005))) {
    mi <- refine_model(mi,
      alternate_family = set_family2,
      use_priors = set_priors
    )
  }
  saveRDS(mi, fmi)
}

if (exists("mi")) {
  s <- sanity(mi, gradient_thresh = 0.005)
  if (all(s)) {
    if (file.exists(pifn) & file.exists(i3)) {
      pi <- readRDS(pifn)
    } else {
      pi <- predict(mi,
        re_form_iid = NA,
        newdata = filter(grid, year %in% sort(union(mi$data$year, mi$extra_time))),
        return_tmb_object = TRUE
      )

      map_density(pi, pifn) +
        labs(
          title = paste0("Immature biomass (", paste(unique(mi$data$survey_type), collapse = ", "), ")"),
          fill = "Density (kg/ha)",
          colour = "Density (kg/ha)"
        ) +
        theme(
          legend.position = c(0.9, 0.1),
          legend.direction = "vertical", legend.box = "horizontal"
        )

      ggsave(paste0("envirocor/figs/density-map-", m3, ".png"),
        height = 10, width = fig_width
      )
    }
  }
}

if (!file.exists(i3)) {
  try(plot_index(pi, species, "Immature", dens_model_name, i3) +
    ggtitle(paste0(species, ": immature biomass (", dens_model_name, ")")))
}

ind0 <- readRDS(i0) %>% mutate(index = "Total")
ind1 <- readRDS(i1) %>% mutate(index = "Mature female")
ind2 <- readRDS(i2) %>% mutate(index = "Mature male")
try(ind3 <- readRDS(i3) %>% mutate(index = "Immature"))

bc_inds <- bind_rows(ind0, ind1, ind2)
try(bc_inds <- bind_rows(ind0, ind1, ind2, ind3))

## Plot coastwide indices ----
m <- readRDS(fm)

(p1 <- bc_inds %>%
  mutate(index = fct_relevel(index, rev)) %>%
  ggplot(aes(year, est / 1000)) +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000, fill = index), alpha = 0.3) +
  geom_line(aes(colour = index), linewidth = 0.7) +
  scale_colour_viridis_d(direction = 1, end = 0.8, option = "A") +
  scale_fill_viridis_d(direction = 1, end = 0.8, option = "A") +
  labs(colour = "Biomass Index", fill = "Biomass Index") +
  xlab("Year") +
  ylab("Relative biomass estimate (tonnes)") +
  ggtitle(paste0(""),
    subtitle = paste0(
      "Model: ",
      ifelse(isTRUE(m$family$delta), m$family$clean_name, paste0(m$family[1], "(link = 'log')")),
      ", spatial (", m[["spatial"]][1], ", ", m[["spatial"]][2],
      ") with st RW and ", dens_model_name_long,
      " (", paste(unique(m$data$survey_type), collapse = ", "), ")"
    )
  ))

p1 + theme(legend.position = c(0.1, 0.75))
ggsave(
  paste0(
    "envirocor/figs/density-index-", spp, "-all-",
    dens_model_name, "-", knot_distance, "-km.png"
  ),
  height = fig_height / 2, width = 9.5
)

# Generate split indices ----
fsi <- paste0(
  "data-generated/density-split-ind/", spp, "-split-",
  dens_model_name, "-", knot_distance, "-km.rds"
)

if (!file.exists(fsi)) {
  mf <- readRDS(fmf)
  mm <- readRDS(fmm)
  inds0 <- split_index_by_survey(m, grid, species, "Total")
  inds1 <- split_index_by_survey(mf, grid, species, "Mature female")
  inds2 <- split_index_by_survey(mm, grid, species, "Mature male")

  all_split_inds <- bind_rows(inds0, inds1, inds2)
  saveRDS(all_split_inds, fsi)

  if (file.exists(i3)) {
    inds3 <- split_index_by_survey(mi, grid, species, "Immature")
    all_split_inds <- bind_rows(inds0, inds1, inds2, inds3)
  }
  saveRDS(all_split_inds, fsi)
} else {
  all_split_inds <- readRDS(fsi)
}

## Plot split indices ----
p2 <- all_split_inds %>%
  left_join(survey_years, ., relationship = "many-to-many") %>%
  select(-survey_abbrev) %>%
  filter(!is.na(est)) %>%
  distinct() %>%
  filter(group == "Total") %>%
  mutate(index = fct_relevel(index, rev)) %>%
  ggplot(aes(year, est / 1000)) +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000, fill = survey), alpha = 0.3) +
  geom_line(aes(colour = survey)) +
  scale_colour_viridis_d(direction = -1, end = 0.9) +
  scale_fill_viridis_d(direction = -1, end = 0.9) +
  labs(colour = "Area", fill = "Area") +
  coord_cartesian(ylim = c(0, max(all_split_inds$est / 1000) * 1.5)) +
  ggtitle("Total") +
  xlab("Year") +
  ylab("Relative biomass estimate (tonnes)")

p2

p3dat <- filter(all_split_inds, group != "Total") %>%
  left_join(survey_years, ., relationship = "many-to-many") %>%
  select(-survey_abbrev) %>%
  filter(!is.na(est)) %>%
  distinct()

p3 <- p3dat %>%
  mutate(index = fct_relevel(index, rev)) %>%
  ggplot(aes(year, est / 1000)) +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000, fill = survey), alpha = 0.3) +
  geom_line(aes(colour = survey)) +
  facet_wrap(~index, scales = "free_x") +
  scale_colour_viridis_d(direction = -1, end = 0.9) +
  scale_fill_viridis_d(direction = -1, end = 0.9) +
  labs(colour = "Area", fill = "Area") +
  coord_cartesian(ylim = c(0, max(p3dat$est / 1000) * 1.25)) +
  xlab("Year") +
  ylab("Relative biomass estimate (tonnes)")

if (length(unique(all_split_inds$model)) > 1) {
  p3 <- p3 + geom_text(aes(label = model),
    x = 2006, y = max(p3dat$est / 1000) * 1.1, size = 2.5, hjust = 0
  )
}

p1a <- p1 + theme(axis.text.x = element_blank(), axis.title = element_blank())
p2a <- p2 + theme(axis.title.x = element_blank())
p3a <- p3 + theme(axis.title.y = element_blank(), legend.position = "none")

p1a + p2a + p3a + plot_layout(ncol = 1)

ggsave(
  paste0(
    "envirocor/figs/density-all-indices-", spp, "-",
    dens_model_name, "-", knot_distance, "-km.png"
  ),
  height = fig_height, width = fig_width
)


# Effect plots ----

variable <- "days_to_solstice"
# variable <- "log_depth_c"

m <- readRDS(fm)
mf <- readRDS(fmf)
mm <- readRDS(fmm)
mi <- readRDS(fmi)

plot_effects_by_class <- function(m, group_label = "Total biomass") {
  if (variable == "log_depth_c") {
    nd <- data.frame(
      days_to_solstice = 0,
      survey_type = as.factor("SYN"),
      log_depth_c = seq(min(m$data$log_depth_c),
        max(m$data$log_depth_c),
        length.out = 50
      ),
      year = max(m$data$year) # a chosen year
    )
    xlabel <- "Log depth (centered on mean)"
  }

  if (variable == "days_to_solstice") {
    nd <- data.frame(
      days_to_solstice = seq(min(m$data$days_to_solstice),
        max(m$data$days_to_solstice),
        length.out = 50
      ),
      survey_type = as.factor("SYN"),
      log_depth_c = 0,
      year = max(m$data$year) # a chosen year
    )
    xlabel <- "Days from solstice"
  }

  pd <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA)
  pd$depth <- exp(pd$log_depth_c + 5)
  m$data$depth <- exp(m$data$log_depth_c + 5)

  set_alpha <- 0.8

  if (variable == "log_depth_c") {
    variable <- "depth"
    xlabel <- "Depth (m)"
  }

  (g <- ggplot(pd, aes(.data[[variable]], exp(est),
    ymin = exp(est - 1.96 * est_se),
    ymax = exp(est + 1.96 * est_se)
  )) +
    geom_rug(
      data = m$data, aes(.data[[variable]], y = 0.9),
      sides = "b", alpha = 0.1, inherit.aes = FALSE
    ) +
    geom_line(alpha = set_alpha) +
    geom_ribbon(alpha = set_alpha / 2) +
    scale_x_continuous() +
    coord_cartesian(expand = F, ylim = c(NA, max(exp(pd$est)) + max(exp(pd$est) * 0.2))) +
    labs(
      x = xlabel,
      y = "Conditional effect"
    ) +
    ggtitle(paste0(group_label)))
}

if (variable == "log_depth_c") {
  xlabel <- "Depth (m)"
}

if (variable == "days_to_solstice") {
  xlabel <- "Days from solstice"
}
y_lab_big <- ggplot() +
  annotate(
    geom = "text", x = 1, y = 1, size = 5,
    label = paste0("Conditional effect"), angle = 90
  ) +
  coord_cartesian(clip = "off") +
  theme_void()

x_lab_big <- ggplot() +
  annotate(
    geom = "text", x = 1, y = 1, size = 5,
    label = paste("", xlabel)
  ) +
  coord_cartesian(clip = "off") +
  theme_void()

p1 <- plot_effects_by_class(m)
p2 <- plot_effects_by_class(mi, "Immatures")
p3 <- plot_effects_by_class(mm, "Mature males")
p4 <- plot_effects_by_class(mf, "Mature females")

p1a <- p1 + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title = element_blank()
)
p2a <- p2 + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title = element_blank()
)
p3a <- p3 + theme(axis.title = element_blank())
p4a <- p3 + theme(axis.title = element_blank())

p <- list(p1a, p2a, p3a, p4a)

design <- "
    AAAAAAAA
    #BBBBBBB
    "

(g <- ((y_lab_big |
  wrap_plots(gglist = p, ncol = 2)) +
  plot_layout(widths = c(0.05, 1)))
/ x_lab_big + plot_layout(heights = c(1, 0.05), design = design)
)

ggsave(paste0("envirocor/figs/dens-effects-", variable, "-", dens_model_name, "-total-filtered.png"),
  height = 5, width = 7
)
