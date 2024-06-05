# Generate condition models and indices
library(tidyverse)
library(sdmTMB)
library(ggsidekick)
# remotes::install_github("seananderson/ggeffects", ref = "sdmTMB")
library(gfenvirocor)

index_list <- expand.grid(
  species = "Petrale Sole",
  maturity = c(
    "mat",
    "imm"
  ),
  males = c(TRUE, FALSE)
) %>%
  mutate(
    females = ifelse(males == FALSE & maturity == "mat", TRUE, FALSE),
    males = ifelse(maturity == "imm", FALSE, males)
  ) %>%
  distinct()


# Function for generating condition indices ----

calc_condition_indices <- function(species, maturity, males, females) {
  spp <- gsub(" ", "-", gsub("\\/", "-", tolower(species)))

  # define models density model used
  dens_model_total <- "dln-1998-new"
  dens_model_name <- "dln-1998-good-split"
  mat_threshold <- 0.5
  knot_distance <- 20

  # some plotting options
  fig_height <- 8
  fig_width <- 10
  theme_set(ggsidekick:::theme_sleek())

  # set other labels and switches

  # stop_early <- TRUE
  stop_early <- FALSE

  plot_covariates <- FALSE
  # plot_covariates <- TRUE

  add_density <- FALSE
  # add_density <- TRUE

  get_mvn_sims <- TRUE

  mat_class <- maturity
  just_males <- males
  just_females <- females

  print(paste(species, maturity, "(males:", just_males, ")"))


  if (stop_early) {
    return(itest)
  }

  # Load condition data and attach lagged density estimates ----
  dir.create(paste0("data-generated/condition-data-w-dens/"),
    showWarnings = FALSE
  )

  f <- paste0(
    "data-generated/condition-data-w-dens/",
    spp, "-", mat_class, "-condition-1998.rds"
  )

  if (!file.exists(f)) {
    ds <- readRDS(paste0(
      "data-generated/condition-data-black-swan/",
      spp, "-mat-", mat_threshold, "-condition.rds"
    )) %>%
      ungroup() %>%
      mutate(
        log_depth_c = log_depth - 5,
        DOY = as.numeric(strftime(time_deployed, format = "%j")),
        days_to_solstice = DOY - 172
      )

    # remove NAs in necessary variables
    ds <- ds %>% filter(!is.na(depth_m))
    ds <- ds %>% filter(!is.na(days_to_solstice))
    ds <- ds %>% filter(!is.na(latitude))
    ds <- ds %>% filter(!is.na(longitude))

    ds$X <- NULL
    ds$Y <- NULL
    ds2 <- sdmTMB::add_utm_columns(ds,
      ll_names = c("longitude", "latitude"),
      utm_crs = 32609
    )

    nd0 <- ds2 %>%
      select(year, survey_abbrev, fishing_event_id, X, Y, log_depth) %>%
      mutate(
        survey_group = "TRAWL",
        survey_type = "SYN",
        days_to_solstice = 0,
        log_depth_c = log_depth - 5
      )

    if (mat_class == "mat") {
      mdf <- readRDS(paste0(
        "data-generated/density-models/",
        dens_model_name, "/mat-fem/",
        spp, "-mat-fem-", dens_model_name, "-",
        knot_distance,
        "-km.rds"
      ))
      mdm <- readRDS(paste0(
        "data-generated/density-models/",
        dens_model_name, "/mat-m/",
        spp, "-mat-m-", dens_model_name, "-",
        knot_distance,
        "-km.rds"
      ))

      mdf <- sdmTMB:::update_version(mdf)
      mdm <- sdmTMB:::update_version(mdm)

      nd0 <- nd0 %>% filter(year >= max(min(mdf$data$year), min(mdm$data$year)))
      # can't predict prior to first year of density model

      pd0f <- predict(mdf, newdata = nd0)
      pd0m <- predict(mdm, newdata = nd0)

      if (length(mdf$family$family) > 1) {
        pd0f <- pd0f %>% mutate(
          density = mdf$family[[1]]$linkinv(est1) * mdf$family[[2]]$linkinv(est2)
        )
      } else {
        pd0f <- pd0f %>% mutate(density = exp(est))
      }

      if (length(mdm$family$family) > 1) {
        pd0m <- pd0m %>% mutate(
          density = mdm$family[[1]]$linkinv(est1) * mdm$family[[2]]$linkinv(est2)
        )
      } else {
        pd0m <- pd0m %>% mutate(density = exp(est))
      }

      pd0 <- pd0f
      pd0$density <- pd0f$density + pd0m$density
    } else {
      dens_model_total <- dens_model_name

      fmi <- paste0(
        "data-generated/density-models/",
        dens_model_name, "/imm/",
        spp, "-imm-", dens_model_name, "-",
        knot_distance,
        "-km.rds"
      )

      if (file.exists(fmi)) {
        md <- readRDS(fmi)
      } else {
        return(NA)
      }

      nd0 <- nd0 %>% filter(year >= min(md$data$year))
      # can't predict prior to first year of density model

      pd0 <- predict(md, newdata = nd0)

      if (length(md$family$family) > 1) {
        pd0 <- pd0 %>% mutate(
          density = md$family[[1]]$linkinv(est1) * md$family[[2]]$linkinv(est2)
        )
      } else {
        pd0 <- pd0 %>% mutate(density = exp(est))
      }
    }

    pd0 <- pd0 %>%
      mutate(
        log_density = log(density),
        model = dens_model_name
      ) %>%
      select(
        year, survey_abbrev, fishing_event_id, X, Y, log_depth,
        density, log_density, model
      ) %>%
      distinct()

    d2 <- left_join(ds2, pd0)
    saveRDS(d2, f)
  } else {
    d2 <- readRDS(f)
  }


  # Select relevant data and grid ----
  if (mat_class == "mat") {
    if (just_males) {
      pf <- paste0(
        "data-generated/density-predictions/p-", spp,
        "-mat-m-", dens_model_name, "-", knot_distance, "-km.rds"
      )
      if (file.exists(pf)) {
        d <- d2 %>% filter(group_name == "Mature males")
        group_tag <- "mat-m"
        group_label <- "mature males"

        # get current year density to scale condition index with
        gridA <- readRDS(pf) %>%
          select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
          group_by(year) %>%
          mutate(
            sum_density = sum(density),
            prop_density = density / sum_density
          ) %>%
          ungroup() %>%
          group_by(year, survey) %>%
          mutate(
            survey_density = sum(density),
            prop_density_by_survey = density / survey_density
          )
      } else {
        print(paste("No density predictions for mature male ", species, "densities."))
        return(NA)
      }
    } else {
      if (just_females) {
        pf <- paste0(
          "data-generated/density-predictions/p-", spp,
          "-mat-fem-", dens_model_name, "-", knot_distance, "-km.rds"
        )
        if (file.exists(pf)) {
          d <- d2 %>% filter(group_name == "Mature females")
          group_tag <- "mat-fem"
          group_label <- "mature females"

          # get current year density to scale condition index with
          gridA <- readRDS(pf) %>%
            select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
            group_by(year) %>%
            mutate(
              sum_density = sum(density),
              prop_density = density / sum_density
            ) %>%
            ungroup() %>%
            group_by(year, survey) %>%
            mutate(
              survey_density = sum(density),
              prop_density_by_survey = density / survey_density
            )
        } else {
          print(paste("No density predictions for mature female ", species, "densities."))
          return(NA)
        }
      } else {
        pf <- paste0(
          "data-generated/density-predictions/p-", spp,
          "-all-mat-", dens_model_name, "-", knot_distance, "-km.rds"
        )
        if (file.exists(pf)) {
          d <- d2 %>%
            filter(group_name %in% c("Mature females", "Mature males")) %>%
            rename(group_catch_weight_split = group_catch_weight)
          d3 <- d %>%
            group_by(fishing_event_id, group_name) %>%
            select(fishing_event_id, group_name, group_catch_weight_split) %>%
            distinct() %>%
            ungroup() %>%
            group_by(fishing_event_id) %>%
            summarise(group_catch_weight = sum(group_catch_weight_split))
          d <- left_join(d, d3)
          group_tag <- "mat"
          group_label <- "mature (females and males)"

          # get current year density to scale condition index with
          gridA <- readRDS(pf) %>%
            select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
            group_by(year) %>%
            mutate(
              sum_density = sum(density),
              prop_density = density / sum_density
            ) %>%
            ungroup() %>%
            group_by(year, survey) %>%
            mutate(survey_density = sum(density), prop_density_by_survey = density / survey_density)
        } else {
          print(paste("No density predictions for all mature ", species, "densities."))
          return(NA)
        }
      }
    }
  } else {
    if (mat_class == "imm") {
      pf <- paste0(
        "data-generated/density-predictions/p-", spp,
        "-imm-", dens_model_name, "-", knot_distance, "-km.rds"
      )
      if (file.exists(pf)) {
        d <- d2 %>% filter(group_name %in% c("Immature"))
        group_tag <- "imm"
        group_label <- "immatures"

        # get current year density to scale condition index with
        gridA <- readRDS(pf) %>%
          select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
          group_by(year) %>%
          mutate(
            sum_density = sum(density),
            prop_density = density / sum_density
          ) %>%
          ungroup() %>%
          group_by(year, survey) %>%
          mutate(
            survey_density = sum(density),
            prop_density_by_survey = density / survey_density
          )

        # browser()
      } else {
        print(paste("No density predictions for immature ", species, "densities."))
        return(NA)
      }
    } else {
      pf <- paste0(
        "data-generated/density-predictions/p-", spp,
        "-total-", dens_model_total, "-", knot_distance, "-km.rds"
      )
      if (file.exists(pf)) {
        # model everything together
        # haven't tested this option recently
        d <- d2

        # get current year density to scale condition index with
        gridA <- readRDS(pf) %>%
          select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
          group_by(year) %>%
          mutate(
            sum_density = sum(density),
            prop_density = density / sum_density
          ) %>%
          ungroup() %>%
          group_by(year, survey) %>%
          mutate(
            survey_density = sum(density),
            prop_density_by_survey = density / survey_density
          )
        ungroup()
      } else {
        print(paste("No density predictions for total ", species, "densities."))
        return(NA)
      }
    }
  }

  # Add density covariates to grid ----
  if (add_density) {
    if (mat_class == "mat") {
      # load lagged density predictions for full survey grid if going to be used as covariates condition
      gridBf <- readRDS(paste0(
        "data-generated/density-predictions/p-", spp,
        "-mat-fem-", dens_model_name, "-", knot_distance, "-km.rds"
      )) %>%
        select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
        rename(
          density_f = density
        )

      gridBm <- readRDS(paste0(
        "data-generated/density-predictions/p-", spp,
        "-mat-m-", dens_model_name, "-", knot_distance, "-km.rds"
      )) %>%
        select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) %>%
        rename(
          density_m = density
        )

      gridB <- full_join(gridBf, gridBm) %>% mutate(density = density_f + density_m)

      gridB <- gridB |>
        select(year, X, Y, survey, depth, days_to_solstice, log_depth, density, density_f, density_m) |>
        group_by(year) |>
        mutate(
          sum_density_f = sum(density_f),
          prop_density_f = density_f / sum_density_f,
          sum_density_m = sum(density_m),
          prop_density_m = density_m / sum_density_m
        ) |>
        ungroup() |>
        group_by(X, Y) |>
        mutate(
          log_density = log(density),
          cell_mean_density = mean(density),
          log_mean_density = log(cell_mean_density)
        ) |>
        ungroup()
    } else {
      gridB <- readRDS(paste0(
        "data-generated/density-predictions/p-", spp,
        "-imm-", dens_model_name, "-", knot_distance, "-km.rds"
      )) |>
        select(year, X, Y, survey, depth, days_to_solstice, log_depth, density) |>
        group_by(year) |>
        mutate(
          sum_density = sum(density),
          prop_density = density / sum_density
        ) |>
        ungroup() |>
        group_by(X, Y) |>
        mutate(
          log_density = log(density),
          cell_mean_density = mean(density),
          log_mean_density = log(cell_mean_density)
        ) |>
        ungroup()
    }
  }

  # Make mesh ----
  hist(d$cond_fac)
  hist(log(d$cond_fac))

  d <- d %>% mutate(
    survey_group = as.factor(
      case_when(
        survey_abbrev %in% c("HBLL OUT N", "HBLL OUT S") ~ "HBLL",
        survey_abbrev %in% c("MSSM QCS", "MSSM WCVI") ~ "MSSM",
        survey_abbrev %in% c("OTHER", "HS MSA", "SYN HS", "SYN QCS", "SYN WCHG", "SYN WCVI") ~ "TRAWL",
        survey_series_id == 68 ~ "HAKE",
        TRUE ~ survey_abbrev
      )
    )
  )
  d <- d %>% filter(year >= min(gridA$year))

  mesh <- make_mesh(d, c("X", "Y"), cutoff = knot_distance)

  ggplot() +
    inlabru::gg(mesh$mesh) +
    coord_fixed() +
    geom_point(aes(X, Y, size = total_weight), data = d2) +
    geom_point(aes(X, Y, colour = group_catch_weight), data = d) +
    facet_wrap(~year) +
    scale_color_viridis_c(trans = ggsidekick::fourth_root_power_trans())

  ggplot() +
    inlabru::gg(mesh$mesh) +
    coord_fixed() +
    geom_point(aes(X, Y, size = group_catch_weight), data = d2) +
    geom_point(aes(X, Y, colour = log(cond_fac)), size = 0.5, data = d) +
    facet_wrap(~year) +
    scale_color_gradient2()

  # load refine function
  refine_cond_model <- function(m, set_formula = cond_formula, dist = knot_distance) {
    s <- sanity(m)
    t <- tidy(m, "ran_pars", conf.int = TRUE)
    if (!all(s) & !m$call$share_range) {
      t <- tidy(m, "ran_pars", conf.int = TRUE)
      if (abs(diff(t$estimate[t$term == "range"])) < dist | !s$range_ok) {
        m <- update(m,
          formula = set_formula,
          share_range = TRUE,
          weights = m$data$sample_multiplier,
          data = m$data, mesh = m$spde
        )
        s <- sanity(m)
        t <- tidy(m, "ran_pars", conf.int = TRUE)
      }
    }
    if (!all(s)) {
      if (t$estimate[t$term == "sigma_O"] < 0.001) {
        m <- update(m,
          formula = set_formula,
          weights = m$data$sample_multiplier,
          spatial = "off",
          data = m$data, mesh = m$spde
        )
      }
    }
    sanity(m)
    return(m)
  }

  # Estimate condition model ----
  d %>%
    group_by(survey_group) %>%
    summarise(n = n())

  model_name <- "1998"

  dir.create(
    paste0(
      "data-generated/condition-models-",
      group_tag, "/", model_name, "/"
    ),
    showWarnings = FALSE
  )
  dir.create(paste0("data-generated/cond-index/", model_name, "/"),
    showWarnings = FALSE
  )


  mf <- paste0(
    "data-generated/condition-models-", group_tag, "/", model_name, "/", spp, "-c-",
    group_tag, "-", model_name, "-", knot_distance, "-km.rds"
  )

  if (diff(range(d$days_to_solstice)) < 60) {
    if (length(unique(d$survey_group)) == 1) {
      cond_formula <- cond_fac ~ 1 + days_to_solstice
    } else {
      cond_formula <- cond_fac ~ survey_group + days_to_solstice
    }
  } else {
    if (length(unique(d$survey_group)) == 1) {
      cond_formula <- cond_fac ~ 1 + poly(days_to_solstice, 2)
    } else {
      cond_formula <- cond_fac ~ survey_group +
        poly(days_to_solstice, 2)
    }
  }

  # if commented out, will be forced to rerun model
  if (file.exists(mf)) {
    try(m <- readRDS(mf))
  }

  if (!exists("m")) {
    sort(unique(d$year))

    m <- sdmTMB(cond_formula,
      mesh = mesh,
      data = d,
      spatial = "on",
      spatiotemporal = "rw",
      extra_time = sdmTMB:::find_missing_time(d$year),
      share_range = FALSE,
      silent = FALSE,
      time = "year",
      offset = NULL,
      family = lognormal(link = "log"),
      priors = sdmTMBpriors(
        matern_s = pc_matern(
          range_gt = knot_distance,
          sigma_lt = 2
        ),
        matern_st = pc_matern(
          range_gt = knot_distance,
          sigma_lt = 2
        )
      )
    )

    print(paste(species, maturity, "(males:", just_males, ")"))

    m <- refine_cond_model(m,
      set_formula = cond_formula,
      dist = knot_distance
    )
    saveRDS(m, mf)
  } else {
    ## make sure it converged last time
    m <- refine_cond_model(m, set_formula = cond_formula, dist = knot_distance)
    saveRDS(m, mf)
  }

  m
  m$sd_report
  tidy(m, "ran_pars", conf.int = TRUE)


  # Add density dependence ----
  ## don't do this for now, but can be used to explore utility of covariates
  if (add_density) {
    rm(m)

    d$log_density_c <- d$log_density - mean(d$log_density, na.rm = TRUE)
    model_name <- "1998-density"

    dir.create(paste0(
      "data-generated/condition-models-", group_tag,
      "/", model_name, "/"
    ), showWarnings = FALSE)

    dir.create(paste0("data-generated/cond-index/", model_name, "/"),
      showWarnings = FALSE
    )

    if (diff(range(d$days_to_solstice)) < 60) {
      if (length(unique(d$survey_group)) == 1) {
        cond_formula2 <- cond_fac ~ 1 +
          days_to_solstice + log_density_c
      } else {
        cond_formula2 <- cond_fac ~ survey_group +
          days_to_solstice + log_density_c
      }
    } else {
      if (length(unique(d$survey_group)) == 1) {
        cond_formula2 <- cond_fac ~ 1 +
          poly(days_to_solstice, 2) +
          log_density_c
      } else {
        cond_formula2 <- cond_fac ~ survey_group +
          poly(days_to_solstice, 2) +
          log_density_c
      }
    }

    mesh2 <- make_mesh(d, c("X", "Y"), cutoff = 15)

    mf2 <- paste0(
      "data-generated/condition-models-", group_tag, "/", model_name, "/",
      spp, "-c-", group_tag, "-", model_name, "-", knot_distance, "-km.rds"
    )

    if (file.exists(mf2)) {
      try(m2 <- readRDS(mf2))
    }

    if (!exists("m2")) {
      m2 <- sdmTMB(cond_formula2,
        mesh = mesh2,
        data = d,
        spatial = "on",
        spatiotemporal = "rw",
        extra_time = seq(min(d$year), max(d$year)),
        share_range = FALSE,
        silent = FALSE,
        time = "year",
        family = lognormal(link = "log"),
        priors = sdmTMBpriors(
          matern_s = pc_matern(
            range_gt = knot_distance,
            sigma_lt = 2
          ),
          matern_st = pc_matern(
            range_gt = knot_distance,
            sigma_lt = 2
          )
        )
      )

      saveRDS(m2, mf2)
      m2 <- refine_cond_model(m2,
        set_formula = cond_formula2,
        dist = knot_distance
      )
      saveRDS(m2, mf2)
    } else {
      m2 <- refine_cond_model(m2,
        set_formula = cond_formula2,
        dist = knot_distance
      )
      saveRDS(m2, mf2)
    }

    sanity(m2)
    m2
    m2$sd_report
    tidy(m2, "ran_pars", conf.int = TRUE)

    ## check if effect of density is negative by more than the SE on the estimate
    ## if not, revert to m1 for index generation
    t <- tidy(m2, conf.int = TRUE)

    if (t$estimate[t$term == "log_density_c"] < -t$std.error[t$term == "log_density_c"]) {
      m <- m2
    } else {
      m <- readRDS(mf)
    }
  }

  # Filter grid ----
  if (add_density) {
    grid <- gridB |> mutate(
      log_density_c = log_mean_density - mean(d$log_density, na.rm = TRUE)
    )
    if (mat_class == "mat") {
      if (just_males) {
        grid <- grid |> mutate(prop_density = prop_density_m)
      } else {
        grid <- grid |> mutate(prop_density = prop_density_f)
      }
    }
  } else {
    grid <- gridA
  }

  grid$survey_group <- "TRAWL"

  sort(unique(m$data$year))
  sort(unique(grid$year))

  # might be redundant
  grid <- filter(grid, year %in% sort(union(unique(m$data$year), m$extra_time)))

  # Has this model been run before? ----
  i1 <- paste0(
    "data-generated/cond-index/", model_name,
    "/cond-index-", group_tag, "-", spp, "-",
    model_name, "-", knot_distance, "-km.rds"
  )

  # remove if wishing to rebuild maps
  if (!file.exists(i1)) {
    pc <- predict(m, newdata = grid, return_tmb_object = TRUE)

    p2 <- pc$data %>%
      mutate(cond = exp(est))

    # filter to plot only cells representing 99% of mean predicted biomass
    # cells must be defined by "X", "Y", time by "year", and biomass/abundance stored as "density"
    p2 <- trim_predictions_by_year(p2, 0.001)


    # Map model predictions ----

    dset <- readRDS(paste0("data-generated/density-data/", spp, ".rds")) %>%
      filter(year >= min(m$data$year), year <= max(m$data$year))

    unique(dset$survey_abbrev)

    # all sets
    set_list <- dset %>%
      select(fishing_event_id, longitude, latitude, X, Y, year, catch_weight, catch_count) %>%
      distinct() %>%
      mutate(
        fishing_event_id = as.factor(fishing_event_id),
        lon = longitude, lat = latitude
      )

    # # just sampled sets
    # set_list <- d2 %>% select(fishing_event_id, longitude, latitude, X, Y) %>% distinct() %>%
    #   mutate(lon = longitude, lat = latitude)

    model_dat <- d %>%
      group_by(fishing_event_id) %>%
      mutate(
        fishing_event_id = as.factor(fishing_event_id),
        count = n()
      )

    model_dat <- left_join(set_list, model_dat, multiple = "all") %>% mutate(
      density = group_catch_weight,
      caught = ifelse(catch_count > 0 | catch_weight > 0, 1, 0),
      count = ifelse(is.na(count), 0, count),
      present = ifelse(count > 0, 1, ifelse(caught == 1, 0, NA_integer_))
    )

    p2$log_cond <- log(p2$cond)
    p2 <- p2 %>% mutate(cond_trim = ifelse(cond > quantile(p2$cond, 0.99),
      quantile(p2$cond, 0.99), cond
    ))

    # remove years with no data for petrale
    p3 <- filter(p2, !(year %in% c(1999, 2001)))
    model_dat2 <- filter(model_dat, !(year %in% c(1999, 2001)))

    # just plot predictions
    g2 <- gfenvirocor::plot_predictions(p3, model_dat2,
      fill_column = "cond_trim",
      fill_label = "Condition \nfactor",
      pt_column = "count",
      pt_label = "Fish \nsampled",
      pt_size_range = c(0.5, 1),
      pos_pt_fill = NA,
      bin_pt_col = NA,
      pos_pt_col = NA,
      fill_scale =
        ggplot2::scale_fill_viridis_c(),
      bounds = grid,
      rotation_angle = 30, show_raw_data = FALSE
    )


    g2 <- g2 + facet_wrap(~year, ncol = 8) +
      ggtitle(paste0(toupper(group_label)))

    ggsave(
      paste0(
        "envirocor/figs/cond-", model_name, "/condition-map-wide-", spp, "-",
        group_tag, "-", model_name, "-", knot_distance, "-km-p.png"
      ),
      height = fig_height, width = fig_width
    )


    # with sampling locations
    g3 <- gfenvirocor::plot_predictions(p3, model_dat2,
      fill_column = "cond_trim",
      fill_label = "Condition \nfactor",
      pt_column = "count",
      pt_label = "Fish \nsampled",
      pt_size_range = c(0.5, 2),
      pos_pt_fill = NA,
      bin_pt_col = "black",
      pos_pt_col = "red",
      fill_scale =
        ggplot2::scale_fill_viridis_c(),
      bounds = grid,
      rotation_angle = 30, show_raw_data = TRUE
    )

    g3 <- g3 + facet_wrap(~year, ncol = 8) +
      ggtitle(paste0(toupper(group_label)))

    ggsave(
      paste0(
        "envirocor/figs/cond-", model_name, "/condition-map-wide-", spp, "-",
        group_tag, "-", model_name, "-", knot_distance, "-km-n.png"
      ),
      height = fig_height, width = fig_width
    )
  }

  # Get coastwide index ----
  if (!file.exists(i1)) {
    ind2 <- get_index(pc, area = grid$prop_density, bias_correct = TRUE)
    ind2$species <- species
    ind2$group <- group_label
    ind2$model_string <- model_name
    ind2$dens_model_name <- dens_model_name

    saveRDS(ind2, i1)
  } else {
    ind2 <- readRDS(i1)
  }

  ggplot(ind2, aes(year, est)) +
    geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
    xlab("Year") +
    ylab("Predicted average condition factor") +
    labs(title = paste0(species, ": ", group_label, " ", model_name))

  ggsave(
    paste0(
      "figs/cond-", model_name, "/condition-index-", spp, "-",
      group_tag, "-", model_name, "-", knot_distance, "-km.png"
    ),
    height = fig_height / 2, width = fig_width / 2
  )

  if (plot_covariates) {
    nd <- data.frame(
      days_to_solstice = seq(min(d$days_to_solstice),
        max(d$days_to_solstice),
        length.out = 50
      ),
      survey_group = "TRAWL",
      log_density_c = 0,
      year = 2021L # a chosen year
    )
    pd <- predict(m, newdata = nd, se_fit = TRUE, re_form = NA)

    ggplot(pd, aes(days_to_solstice, exp(est),
      ymin = exp(est - 1.96 * est_se),
      ymax = exp(est + 1.96 * est_se)
    )) +
      geom_line() +
      geom_ribbon(alpha = 0.4) +
      scale_x_continuous() +
      coord_cartesian(expand = F) +
      labs(x = "Days to solstice", y = "condition") +
      ggtitle(paste0(species, ": ", group_label, " ", model_name))

    ggsave(
      paste0(
        "figs/cond-", model_name, "/effect-plot-doy-", spp, "-",
        group_tag, "-", model_name, "-", knot_distance, "-km.png"
      ),
      height = fig_height / 2, width = fig_width / 2
    )


    # Effect plots ----
    if (add_density) {
      nd2 <- data.frame(
        days_to_solstice = 0,
        survey_group = "TRAWL",
        log_density_c = seq(min(d$log_density_c),
          max(d$log_density_c),
          length.out = 50
        ),
        year = 2021L # a chosen year
      )
      pd <- predict(m, newdata = nd2, se_fit = TRUE, re_form = NA)


      ggplot(pd, aes(log_density_c, exp(est),
        ymin = exp(est - 1.96 * est_se),
        ymax = exp(est + 1.96 * est_se)
      )) +
        geom_line() +
        geom_ribbon(alpha = 0.4) +
        scale_x_continuous() +
        coord_cartesian(expand = F) +
        labs(x = "log_density_c", y = "condition") +
        ggtitle(paste0(species, ": ", group_label, " ", model_name))


      ggsave(
        paste0(
          "figs/cond-", model_name, "/effect-plot-density-", spp, "-",
          group_tag, "-", model_name, "-", knot_distance, "-km.png"
        ),
        height = fig_height / 2, width = fig_width / 2
      )
    }
  }

  # Get MVN simulation samples ----
  if (get_mvn_sims) {
    psims <- predict(m, newdata = grid, nsim = 100)
    isims <- get_index_sims(psims, area = grid$prop_density, return_sims = TRUE)

    saveRDS(isims, paste0(
      "data-generated/cond-index/", model_name,
      "/cond-index-sims-", group_tag, "-", spp, "-",
      model_name, "-", knot_distance, "-km.rds"
    ))
  }
}

# Run with pmap -----
pmap(index_list, calc_condition_indices)
