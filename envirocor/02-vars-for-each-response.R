# prepare environmental variables for recruitment
library(tidyverse)
library(pacea)
source("envirocor/utils.R")

### Climate indices -----
npgoP <- npgo |> rename(value = anomaly) |> filter_months(c(1,2,3,4,5,6,7,8,9), "NPGO (Q1-3: all pelagic)")

# load("data-generated/pdo.rda")
pdoW <- pdo |> rename(value = anomaly) |> filter_months(c(1,2,3), "PDO (Q1: spawning)")
pdoSp <- pdo |> rename(value = anomaly) |> filter_months(c(4,5,6),"PDO (Q2: spring pelagic)")

### production variable ----
pt_monthly <- readRDS("data-generated/cw_phytoplankton.rds")
ptP <- pt_monthly |>
  filter_months(c(1,2,3,4,5,6,7,8,9), "Phytoplankton (Q1-3: all pelagic)")

## SST ----
sst_monthly <- readRDS("data-generated/cw_surface_temperature.rds")
sstW <- readRDS("data-generated/cw_surface_temperature_spawn.rds") |>
  filter_months(c(1,2,3), "SST (Q1: spawning)")
sstSp <- sst_monthly |>
  filter_months(c(4,5,6), "SST (Q2: spring pelagic)")
sstSu <- sst_monthly |>
  filter_months(c(7,8,9), "SST (Q3: summer pelagic)")

## TOB -----
tob_monthly <- readRDS("data-generated/cw_bottom_temperature.rds")
tobW <- readRDS("data-generated/cw_bottom_temperature_spawn.rds") |>
  filter_months(c(1,2,3), "Sea floor temperature (Q1: spawning)")

## O2 ----
do_monthly <- readRDS("data-generated/cw_bottom_oxygen.rds")
o2W <- readRDS("data-generated/cw_bottom_oxygen_spawn.rds") |>
  filter_months(c(1,2,3), "Sea floor O2 (Q1: spawning)")

load("Euphausiids/fit_ph_all_index_bioregion.Rdata")

ggplot(filter(index_bioregion,
       bioregion != "Offshore bioregion")) +
  geom_line(aes(year, log_est, colour = forcats::fct_rev(bioregion))) +
  facet_wrap(~model)

euphausiids <- index_bioregion |>
  filter(bioregion != "Offshore bioregion") |>
  rename(type = model) |>
  group_by(year, type) |> summarise(value = sum(est)) |> ungroup() |>
  group_by(type) |>
  mutate(
    value = value/1000000,
    time = seq_along(year),
         value_raw = value,
         mean = mean(value),
         sd = sd(value),
         value = (value - mean(value))/ sd(value)
         )

euph_regional <- index_bioregion |>
  rename(type = model) |>
  select(year, type, bioregion, est) |>
  group_by(type, bioregion) |>
  mutate(
    value = est/1000000,
    time = seq_along(year),
    value_raw = value,
    mean = mean(value),
    sd = sd(value),
    value = (value - mean(value))/ sd(value)
  ) |> ungroup()

ggplot(euphausiids) + geom_line(aes(year, value, colour = type))

epac <- filter(euphausiids, type == "E. pacifica")
tspin <- filter(euphausiids, type == "T. spinifera")
euph <- filter(euphausiids, type == "Total euphausiid")
euphN <- filter(euph_regional, type == "Total euphausiid"&bioregion == "Northern bioregion")|> select(-bioregion, -est)
euphS <- filter(euph_regional, type == "Total euphausiid"&bioregion == "Southern bioregion")|> select(-bioregion, -est)
euphO <- filter(euph_regional, type == "Total euphausiid"&bioregion == "Offshore bioregion")|> select(-bioregion, -est)

epac$type <- "Euphausiid index (E. pacifica)"
tspin$type <- "Euphausiid index (T. spinifera)"
euph$type <- "Euphausiid index (Total)"

euphN$type <- "Euphausiid index (Total QCS)"
euphS$type <- "Euphausiid index (Total WCVI)"
euphO$type <- "Euphausiid index (Total offshore)"


dvs <- bind_rows(
  npgoP,
  euphO,
  euphN,
  euphS,
  pdoW,
  pdoSp,
  o2W,
  tobW,
  ptP,
  sstW,
  sstSp,
  sstSu
)
saveRDS(dvs, "data-generated/envrio-vars-for-rdevs.rds")

dvs2 <- bind_rows(
  euphN,
  pdoW,
  npgoP,
  sstW,
  sstSp,
  sstSu
)
saveRDS(dvs2, "data-generated/envrio-vars-for-rdevs-shortlist.rds")

# for condition
pdo0 <- pdo |>rename(value = anomaly) |>
  filter_months(c(4,5,6), "PDO (Q2: Apr-Jun)")
npgo0 <- npgo |>rename(value = anomaly) |>
  filter_months(c(4,5,6), "NPGO (Q2: Apr-Jun)")

### production variables ----
pp0 <- readRDS("data-generated/cw_primary_production_summer.rds") |>
  filter_months(c(4,5,6), "Primary production (Q2: Apr-Jun)")

## TOB -----
tob0 <- readRDS("data-generated/cw_bottom_temperature_summer.rds") |>
  filter_months(c(4,5,6), "Sea floor temperature (Q2: Apr-Jun)")

## O2 ----
o20 <- readRDS("data-generated/cw_bottom_oxygen_summer.rds") |>
  filter_months(c(4,5,6),  "Sea floor O2 (Q2: Apr-Jun)")


dvc <- bind_rows(
  euph,
  epac,
  pdo0,
  npgo0,
  o20,
  tob0,
  pp0
)

saveRDS(dvc, "data-generated/envrio-vars-for-condition.rds")

# more correlation plots

library(GGally)

dw <- bind_rows(pdoW, sstW, pdoSp, sstSp,  pdoSu, sstSu) |> select(year, type, value) |>
  pivot_wider(names_from = type, values_from = value)

ggpairs(dw, columns = c(2:7),
        upper = list(continuous = wrap(cor_func, method = 'spearman', symbol = expression('\u03C1 ='))),
        progress = FALSE)

ggsave("envirocor/figs/pdo-sst-variables-correlations.png", width = 11, height = 11)


dw <- bind_rows(
  pdo0,
  ssa0,
  pt0,
  pp0,
  epac,
  tspin,
  euph) |> select(year, type, value) |> pivot_wider(names_from = type, values_from = value)

ggpairs(dw, columns = c(2:8),
        upper = list(continuous = wrap(cor_func, method = 'spearman', symbol = expression('\u03C1 ='))),
        progress = FALSE)

ggsave("envirocor/figs/pdo-biotic-variables-correlations.png", width = 14, height = 14)

