# 1. Calculate Le Cren’s relative condition factor for each fish sampled
# 2. Remove the most extreme outliers that are likely errors
# 3. Split samples into immature and mature and into length bins

library(tidyverse)
library(gfplot)
theme_set(ggsidekick::theme_sleek())

species <- "Petrale Sole"

spp <- gsub(" ", "-", gsub("\\/", "-", tolower(species)))

# # # Load data ----
dat <- readRDS("data-generated/all-samples-used.rds") %>%
  filter(
    !is.na(longitude), !is.na(latitude),
    species_common_name == tolower(species)
  )

check_for_duplicates <- dat[duplicated(dat$specimen_id), ]

if (nrow(check_for_duplicates) > 0) {
  stop(paste(species, "has duplicate specimen ids."))
}

sort(unique(dat$survey_abbrev))

dset <- readRDS("data-generated/all-sets-used.rds") %>%
  filter(species_common_name == tolower(species))

# TODO: temporary because false zeros were in an older data pull
dat$catch_count <- ifelse(dat$catch_weight > 0 & dat$catch_count == 0, NA, dat$catch_count)
dat$catch_weight <- ifelse(dat$catch_count > 0 & dat$catch_weight == 0, NA, dat$catch_weight)


datf <- filter(dat, !is.na(length))
datf <- filter(datf, !is.na(survey_abbrev))
unique(datf$survey_abbrev)
unique(dat$survey_abbrev)

fish <- dat

if (length(unique(fish$length_type)) > 1) {
  stop("Stop. Two different length types.")
}

ggplot(
  filter(fish, !is.na(latitude) & !is.na(longitude) & !is.na(length) & !is.na(weight)),
  aes(longitude, latitude, colour = survey_abbrev)
) +
  geom_point() +
  facet_wrap(~year)

# 2. Split samples into immature and mature and into length bins ----
# Starting with code lifted from split by maturity function
split_by_maturity <- TRUE
split_by_sex <- TRUE
immatures_pooled <- TRUE

## get the same maturity model from the density split
m <- readRDS(paste0("data-generated/split-catch-data/", spp, ".rds"))$m

p_threshold <- mat_threshold <- 0.5
custom_maturity <- NULL
custom_length_threshold <- NULL

## could use separate estimates for each year
# year_re <- TRUE
## discovered that petrale length at maturity was unusually high in WCVI 2004 and 2006
year_re <- FALSE
sample_id_re <- TRUE

# -----
# does maturity data exist at all for this species?
if (split_by_maturity) {
  maturity_codes <- unique(fish$maturity_code)
  if (length(maturity_codes) < 3 & is.null(custom_length_threshold)) {
    return("Fewer than 3 maturity codes. Please provide custom length thresholds from literature.")
  }

  # Check if only some years without maturity data, and set year_re = FALSE in that case
  years_w_maturity <- fish %>%
    group_by(year) %>%
    summarise(maturity_levels = length(unique(maturity_code)))

  levels_per_year <- unique(years_w_maturity$maturity_levels)

  if (min(levels_per_year) < 3) { # some years lack maturity data
    years_w_maturity <- fish %>%
      group_by(year) %>%
      summarise(maturity_levels = length(unique(maturity_code)))

    levels_per_year <- unique(years_w_maturity$maturity_levels)
  }

  f_fish <- fish %>%
    filter(sex == 2) %>%
    mutate(year_f = as.character(year))

  m_fish <- fish %>%
    filter(sex == 1) %>%
    mutate(year_f = as.character(year))

  if (min(levels_per_year) < 3) { # some years lack maturity data
    if (length(levels_per_year) < 1) {
      warning("Maturity data not recorded, so catch not split.", call. = FALSE)
    } else {
      warning("Some years lack maturity data, but catch still split.", call. = FALSE)

      if (!is.null(custom_length_threshold)) {
        f_fish$threshold <- custom_length_threshold[2]
        m_fish$threshold <- custom_length_threshold[1]
      } else {
        if (p_threshold == 0.5) {
          f_fish$threshold <- m$mat_perc$f.p0.5
          m_fish$threshold <- m$mat_perc$m.p0.5
        }
        if (p_threshold == 0.05) {
          f_fish$threshold <- m$mat_perc$f.p0.05
          m_fish$threshold <- m$mat_perc$m.p0.05
        }
        if (p_threshold == 0.95) {
          f_fish$threshold <- m$mat_perc$f.p0.95
          m_fish$threshold <- m$mat_perc$m.p0.95
        }
      }
    }
  } else {
    if (year_re) {
      if (!is.null(custom_length_threshold)) {
        f_fish$threshold <- custom_length_threshold[2]
        m_fish$threshold <- custom_length_threshold[1]
      } else {
        if (p_threshold == 0.5) {
          f_fish$threshold <- lapply(f_fish$year_f, function(x) m$mat_perc[[x]]$f.p0.5)
          m_fish$threshold <- lapply(m_fish$year_f, function(x) m$mat_perc[[x]]$m.p0.5)
        }
        if (p_threshold == 0.05) {
          f_fish$threshold <- lapply(f_fish$year_f, function(x) m$mat_perc[[x]]$f.p0.05)
          m_fish$threshold <- lapply(m_fish$year_f, function(x) m$mat_perc[[x]]$m.p0.05)
        }
        if (p_threshold == 0.95) {
          f_fish$threshold <- lapply(f_fish$year_f, function(x) m$mat_perc[[x]]$f.p0.95)
          m_fish$threshold <- lapply(m_fish$year_f, function(x) m$mat_perc[[x]]$m.p0.95)
        }
      }
    } else {
      if (!is.null(custom_length_threshold)) {
        f_fish$threshold <- custom_length_threshold[2]
        m_fish$threshold <- custom_length_threshold[1]
      } else {
        # apply global estimates to all catches
        if (p_threshold == 0.5) {
          f_fish$threshold <- m$mat_perc$f.p0.5
          m_fish$threshold <- m$mat_perc$m.p0.5
        }
        if (p_threshold == 0.05) {
          f_fish$threshold <- m$mat_perc$f.p0.05
          m_fish$threshold <- m$mat_perc$m.p0.05
        }
        if (p_threshold == 0.95) {
          f_fish$threshold <- m$mat_perc$f.p0.95
          m_fish$threshold <- m$mat_perc$m.p0.95
        }
      }
    }
  }

  if (spp %in% c(" ")) {
    f_fish <- mutate(f_fish, mature = 1)
    m_fish <- mutate(m_fish, mature = 1)
  } else {
    # classify each fish as immature or mature based on above thresholds
    f_fish <- mutate(f_fish, mature = if_else(length >= threshold, 1, 0, missing = NULL))
    m_fish <- mutate(m_fish, mature = if_else(length >= threshold, 1, 0, missing = NULL))
  }

  # get unsexed immature fish
  imm_fish <- fish %>%
    filter(!(sex %in% c(1, 2)) &
      length < min(c(
        min(unique(unlist(f_fish$threshold)), na.rm = TRUE),
        min(unique(unlist(m_fish$threshold)), na.rm = TRUE)
      ), na.rm = TRUE)) %>%
    mutate(
      mature = 0,
      year_f = as.character(year)
    )

  # create groups
  if (split_by_sex) {
    if (immatures_pooled) {
      # since not splitting by sex for immatures, the unsexed imm can be added on
      fish_groups <- bind_rows(f_fish, m_fish, imm_fish) %>%
        mutate(group_name = ifelse(mature == 1,
          paste("Mature", ifelse(sex == 1, "males", "females")),
          "Immature"
        ))
    } else {
      fish_groups <- rbind(f_fish, m_fish) %>%
        mutate(group_name = ifelse(mature == 1,
          paste("Mature", ifelse(sex == 1, "males", "females")),
          paste("Immature", ifelse(sex == 1, "males", "females"))
        ))
    }
  } else {
    fish_groups <- rbind(f_fish, m_fish, imm_fish) %>%
      mutate(group_name = ifelse(mature == 1, "Mature", "Immature"))
  }
} else {
  # just split by sex
  fish_groups <- rbind(f_fish, m_fish) %>%
    mutate(group_name = ifelse(sex == 1, "Males", "Females"))
}

# 2. Le Cren’s relative condition factor ----
# defaults to a TMB model with a robust likelihood using t distribution with df = 3
mf <- gfplot::fit_length_weight(fish_groups, sex = "female", usability_codes = NULL)
mm <- gfplot::fit_length_weight(fish_groups, sex = "male", usability_codes = NULL)

## Length-weight plot ----
gfplot::plot_length_weight(object_female = mf, object_male = mm)

## Remove black swan outliers ----
sd_threshold <- 2
is_heavy_tail <- TRUE

filter_lw_outliers <- function(model,
                               numsd = sd_threshold,
                               heavy_tailed = is_heavy_tail) {
  if (heavy_tailed) {
    l <- qt(0.025, 3) * exp(model$pars$log_sigma) * numsd
    u <- qt(0.975, 3) * exp(model$pars$log_sigma) * numsd
  } else {
    l <- qnorm(0.025, 0, sd = exp(model$pars$log_sigma) * numsd)
    u <- qnorm(0.975, 0, sd = exp(model$pars$log_sigma) * numsd)
  }
  model$data$resids <- log(model$data$weight) - (model$pars$log_a + model$pars$b * log(model$data$length))
  out <- model$data |> filter(!(resids > u) & !(resids < l))
  out
}

df <- filter_lw_outliers(mf)
dm <- filter_lw_outliers(mm)

df$wbar <- exp(mf$pars$log_a) * df$length^mf$pars$b
dm$wbar <- exp(mm$pars$log_a) * dm$length^mm$pars$b

# include unknown sex individuals for now, because immature individuals can be difficult to sex and differences in growth rate may be slim
du <- dplyr::filter(fish_groups, sex %in% c(0, 3), !is.na(weight), !is.na(length))
# Apply an intermediate slope and intercept to these individuals
# weight is in grams, so convert to kg
du$weight <- du$weight / 1000

du$wbar <- exp((mm$pars$log_a + mf$pars$log_a) / 2) * du$length^((mm$pars$b + mf$pars$b) / 2)

dd <- bind_rows(df, dm)
dd$cond_fac <- dd$weight / dd$wbar
du$cond_fac <- du$weight / du$wbar

dd2 <- filter(du, cond_fac >= min(dd$cond_fac) & cond_fac <= max(dd$cond_fac)) |> bind_rows(dd)
# plot(cond_fac~length, data = dd2)

# 4. Add some useful variables ----
ds <- dd2 %>%
  group_by(fishing_event_id, group_name) %>%
  mutate(
    log_depth = log(depth_m),
    group_sampled_weight = sum(weight, na.rm = T),
    group_num_sampled = n()
  ) %>%
  ungroup() %>%
  group_by(fishing_event_id) %>%
  mutate(
    sampled_weight = sum(weight, na.rm = T),
    num_sampled = n(),
    mean_weight = mean(weight, na.rm = T),
    total_weight = ifelse(
      is.na(catch_weight) & catch_count > 0, catch_count * mean_weight, ifelse(
        catch_weight > sampled_weight, catch_weight, sampled_weight
      )
    ),
    est_count = round(total_weight / mean_weight),
    est_count = ifelse(num_sampled > est_count, num_sampled, est_count),
    prop_n_in_group = (group_num_sampled / num_sampled), # this is proportion by count
    prop_w_in_group = (group_sampled_weight / sampled_weight), # this is proportion by weight
    group_catch_weight = total_weight * prop_w_in_group, # this is est group catch weight
    est_num_unsampled_group_members = (est_count - num_sampled) * prop_n_in_group
  ) %>%
  unique()

# 5. Plotted Le Cren's----

ggplot(
  dat |> filter(
    sex %in% c(1, 2)
  ) |>
    mutate(
      sex_label = ifelse(sex == 1, "Male", "Female"),
      weight = weight / 1000,
      vline = ifelse(sex == 1, m_fish$threshold[1], f_fish$threshold[1])
    ),
  aes(length, weight, shape = as.factor(sex))
) +
  geom_point(colour = "red") +
  geom_point(colour = "white", data = ds) +
  geom_point(aes(colour = cond_fac),
    data = filter(ds, sex %in% c(1, 2)) |> mutate(sex_label = ifelse(sex == 1, "Male", "Female")),
    alpha = 0.4
  ) +
  facet_wrap(~sex_label) +
  geom_vline(aes(xintercept = vline), linetype = "dashed") +
  labs(
    colour = "Le Cren's", shape = "Sex",
    x = "Length (cm)", y = "Weight (kg)"
  ) +
  scale_colour_viridis_c() +
  scale_shape_discrete(guide = NULL) +
  ggsidekick::theme_sleek() +
  theme(legend.position = c(0.1, 0.8))

ggsave(paste0("envirocor/figs/cond-black-swan-",
  ifelse(is_heavy_tail, "t", "norm"), "-", sd_threshold, "sd-", spp, ".png"
), width = 6, height = 4)


# save data
dir.create(paste0("data-generated/condition-data-black-swan/"), showWarnings = FALSE)
saveRDS(ds, paste0("data-generated/condition-data-black-swan/", spp, "-mat-", mat_threshold, "-condition.rds"))


# sampling date stats and figure?
d2 <- ds |> mutate(DOY = as.numeric(strftime(time_deployed, format = "%j")))

# hist(d2$DOY)
# median(d2$DOY) # june 6
# mean(d2$DOY) # june 14
# range(d2$DOY) # April 25 - Sept 22nd

ggplot(d2) +
  geom_vline(xintercept = 172, colour = "darkgrey", linetype = "dashed") +
  geom_histogram(aes(DOY, fill = survey_abbrev)) +
  coord_cartesian(expand = FALSE) +
  scale_fill_viridis_d(option = "G", direction = -1, end = 0.95, begin = 0.15) +
  labs(y = "Specimens sampled", fill = "Survey")

ggsave("envirocor/figs/sample-date-hist.png", width = 5, height = 3)
