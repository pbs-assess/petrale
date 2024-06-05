# time series regressions for body condition
library(sdmTMB)
library(glmmTMB)
library(ggeffects)
library(tidyverse)
library(patchwork)
library(brms)
source("envirocor/utils.R")
theme_set(ggsidekick::theme_sleek())

which_cond_model <- "1998"

scenario <- 300
start_year <- 1998
end_year <- 2018

ts <- readRDS(paste0("data-generated/df-", scenario, ".RData"))
ts <- ts |> filter(year >= start_year & year <= end_year)

model_name1 <- "1998"

f1 <- list.files(paste0(
  "data-generated/cond-index/",
  model_name1
), pattern = ".rds", full.names = TRUE)

d1 <- purrr::map_dfr(f1, readRDS) %>%
  mutate(group = factor(group,
    levels = c("immatures", "mature males", "mature females"),
    labels = c("Immatures", "Mature males", "Mature females")
  ))
d1$model <- "Density-agnostic"

min_est <- min(d1$est)
max_est <- max(d1$est)

yrs <- sort(unique(d1$year))

dd1 <- d1 |>
  group_by(group) |>
  mutate(
    value = (est - 1) / sd(est),
    type = ifelse(group == "Mature females", "Female condition",
      ifelse(group == "Mature males", "Male condition",
        "Immature condition"
      )
    ),
    group = factor(group,
      levels = c("immatures", "mature males", "mature females"),
      labels = c("Immatures", "Mature males", "Mature females")
    )
  ) |>
  ungroup() |>
  select(year, value, type)
dd1 <- na.omit(dd1)

dd2 <- dd1 |> filter(type != "Immature condition")

ts2 <- left_join(ts, dd2)
ts2 <- na.omit(ts2)



# correlations between condition and recruitment ----

poly <- TRUE
lag <- TRUE

p <- list()
m <- list()
coefs <- list()

n_draws <- 100
colours <- c(
  "#fde725",
  "#21918c"
)


if (lag) {
  data <- ts2 |>
    group_by(type) |>
    mutate(year = lead(year)) |>
    ungroup() |> ## lag condition
    filter(year >= start_year & year <= end_year)
} else {
  data <- ts2 |> filter(year >= start_year & year <= end_year)
}

data$response <- data[["rdev"]]
data$var_names <- data[["type"]]

for (i in seq_along(sort(unique(data$var_names)))) {
  dat <- filter(data, var_names == sort(unique(data$var_names))[[i]])
  dat <- dat |> mutate(time = as.integer(seq_along(year)))

  if (poly) {
    m[[i]] <- brm(
      bf(response ~ poly(value, 2) + ar(time = time)),
      data = dat,
      iter = 1000,
      chains = 3,
      prior =
        c(
          set_prior("normal(0, 1)", class = "ar"),
          set_prior("normal(0, 10)", class = "b"),
          set_prior("student_t(3, 0, 2)", class = "sigma"),
          set_prior("normal(0, 10)", class = "Intercept")
        ),
      backend = "cmdstan"
    )
  } else {
    m[[i]] <- brm(
      bf(response ~ value + ar(time = time)),
      data = dat,
      iter = 1000,
      chains = 3,
      prior =
        c(
          set_prior("normal(0, 1)", class = "ar"),
          set_prior("normal(0, 10)", class = "b"),
          set_prior("student_t(3, 0, 2)", class = "sigma"),
          set_prior("normal(0, 10)", class = "Intercept")
        ),
      backend = "cmdstan"
    )
  }

  if (unique(dat$var_names) == "Immature condition") group <- "imm"
  if (unique(dat$var_names) == "Female condition") group <- "mat-fem"
  if (unique(dat$var_names) == "Male condition") group <- "mat-m"

  # retrieve MCMC/MVN sample data frames for 'response' and covariate:
  dd <- purrr::map_dfr(seq_len(n_draws), \(j) {
    .draw <- readRDS(paste0("data-generated/dfs/df", j, ".RData"))
    .cov <- readRDS(paste0(
      "data-generated/cond-index/",
      which_cond_model, "/cond-index-sims-", group,
      "-petrale-sole-", which_cond_model, "-20-km.rds"
    )) |>
      filter(.iteration == j) |>
      mutate(
        value = (.value - 1) / sd(.value),
        type = unique(dat$var_names)
      ) |>
      select(year, value, type)

    if (lag) {
      # only one iteration and type, so shouldn't be grouped
      .cov <- .cov |>
        mutate(year = lead(year)) |>
        ungroup() ## lag condition
    }

    .d <- left_join(.draw, .cov) |> filter(year >= start_year & year <= end_year)
    .d <- na.omit(.d)
    .d$response <- .d[["rdev"]]
    .d$var_names <- .d[["type"]]
    .d <- filter(.d, var_names == unique(dat$var_names))
    .d <- .d |> mutate(time = as.integer(seq_along(year)))
    .d$original_iter <- j
    .d
  })


  dd_sum <- select(dat, year, value, response) |>
    rename(value_med = value, response_med = response) |>
    right_join(dd) |>
    group_by(year, value_med, response_med) |>
    summarise(
      y_max = quantile(response, probs = 0.975),
      y_min = quantile(response, probs = 0.025),
      x_max = quantile(value, probs = 0.975),
      x_min = quantile(value, probs = 0.025),
      response_new_med = median(response),
      value_new_med = median(value)
    )

  if (poly) {
    fits <- dd |>
      split(dd$original_iter) |>
      lapply(do_fit)
  } else {
    fits <- dd |>
      split(dd$original_iter) |>
      lapply(do_fit, poly = FALSE)
  }

  nd <- data.frame(value = seq(min(dd$value), max(dd$value), length.out = 200), time = NA)

  # make predictions for each:
  preds <- fits |> lapply(\(x) {
    posterior_epred(
      x,
      incl_autocor = FALSE,
      re.formula = NA,
      newdata = nd
    )
  })
  pred <- do.call(rbind, preds)

  # plot
  nd$est <- apply(pred, 2, median)
  nd$lwr <- apply(pred, 2, quantile, probs = 0.025)
  nd$upr <- apply(pred, 2, quantile, probs = 0.975)

  # nd$value_raw <- nd$value*dd$sd[1]+dd$mean[1]
  preds_full_posterior <- fits |> lapply(\(x) {
    posterior_predict(
      x,
      incl_autocor = FALSE,
      re.formula = NA,
      newdata = nd
    )
  })
  pred2 <- do.call(rbind, preds_full_posterior)

  nd$lwr2 <- apply(pred2, 2, quantile, probs = 0.025)
  nd$upr2 <- apply(pred2, 2, quantile, probs = 0.975)

  (p[[i]] <- ggplot() +
    geom_ribbon(
      data = nd, aes(value, ymin = lwr, ymax = upr),
      alpha = 0.5, fill = colours[[i]]
    ) +
    geom_ribbon(
      data = nd, aes(value, ymin = lwr2, ymax = upr2),
      alpha = 0.25, fill = colours[[i]]
    ) +
    # this is median of simulation draws
    geom_point(
      data = dd_sum, aes(value_new_med, response_new_med),
      colour = colours[[i]]
      # colour = "white"
    ) +
    # a place holder to set the axes correctly
    geom_linerange(
      data = dd_sum, aes(x = value_new_med, ymin = y_min, ymax = y_max),
      colour = colours[[i]]
    ) +
    geom_linerange(
      data = dd_sum, aes(y = response_new_med, xmin = x_min, xmax = x_max),
      colour = colours[[i]]
    ) +
    geom_line(
      data = nd, aes(value, est),
      colour = colours[[i]]
    ) +
    labs(
      x = unique(dd$var_names), y = "",
      colour = "", fill = ""
    ) +
    ggtitle("") +
    ggsidekick::theme_sleek()
  )

  if (poly) {
    coefs[[i]] <- fits |> purrr::map_dfr(\(x) {
      a <- as.data.frame(x)
      data.frame(poly1 = a$b_polyvalue21, poly2 = a$b_polyvalue22, p = a$`ar[1]`, sigma = a$sigma)
    })
  } else {
    coefs[[i]] <- fits |> purrr::map_dfr(\(x) {
      a <- as.data.frame(x)
      data.frame(slope = a$b_value, p = a$`ar[1]`, sigma = a$sigma)
    })
  }
  coefs[[i]]$type <- unique(dd$var_names)[1]
}

if (lag) {
  if (poly) {
    saveRDS(coefs, paste0("data-generated/rdev-condition-corr-coefs-", n_draws, "-poly-lag.rds"))
    saveRDS(p, paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, "-poly-lag.rds"))
    saveRDS(m, paste0("data-generated/rdev-condition-corr-model-list-", n_draws, "-poly-lag.rds"))
  } else {
    saveRDS(coefs, paste0("data-generated/rdev-condition-corr-coefs-", n_draws, "-lag.rds"))
    saveRDS(p, paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, "-lag.rds"))
    saveRDS(m, paste0("data-generated/rdev-condition-corr-model-list-", n_draws, "-lag.rds"))
  }
} else {
  if (poly) {
    saveRDS(coefs, paste0("data-generated/rdev-condition-corr-coefs-", n_draws, "-poly.rds"))
    saveRDS(p, paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, "-poly.rds"))
    saveRDS(m, paste0("data-generated/rdev-condition-corr-model-list-", n_draws, "-poly.rds"))
  } else {
    saveRDS(coefs, paste0("data-generated/rdev-condition-corr-coefs-", n_draws, ".rds"))
    saveRDS(p, paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, ".rds"))
    saveRDS(m, paste0("data-generated/rdev-condition-corr-model-list-", n_draws, ".rds"))
  }
}

# # reload saved versions
# if(lag){
#   if(poly){
#     coefs <- readRDS(paste0("data-generated/rdev-condition-corr-coefs-", n_draws, "-poly-lag.rds"))
#     p <- readRDS(paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, "-poly-lag.rds"))
# m <- readRDS(paste0("data-generated/rdev-condition-corr-model-list-", n_draws, "-poly-lag.rds"))
#   } else{
#     coefs <- readRDS(paste0("data-generated/rdev-condition-corr-coefs-", n_draws, "-lag.rds"))
#     p <- readRDS(paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, "-lag.rds"))
#     m <- readRDS(paste0("data-generated/rdev-condition-corr-model-list-", n_draws, "-lag.rds"))
#   }
# } else{
#   if(poly){
#     coefs <- readRDS(paste0("data-generated/rdev-condition-corr-coefs-", n_draws, "-poly.rds"))
#     p <- readRDS(paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, "-poly.rds"))
#     m <- readRDS(paste0("data-generated/rdev-condition-corr-model-list-", n_draws, "-poly.rds"))
#   } else{
#     coefs <- readRDS(paste0("data-generated/rdev-condition-corr-coefs-", n_draws, ".rds"))
#     p <- readRDS(paste0("data-generated/rdev-condition-corr-plot-list-", n_draws, ".rds"))
#     m <- readRDS(paste0("data-generated/rdev-condition-corr-model-list-", n_draws, ".rds"))
#   }
# }

# check for convergence
lapply(m, get_ess)
lapply(m, max_rhat)


y_lab_big <- ggplot() +
  annotate(
    geom = "text", x = 1, y = 1, size = 4.5, colour = "grey30",
    label = paste0("Recruitment deviations"), angle = 90
  ) +
  coord_cartesian(clip = "off") +
  theme_void()

(pp <- ((y_lab_big |
  wrap_plots(gglist = p, ncol = 2) &
    theme(
      text = element_text(size = 12), legend.position = "none",
      plot.tag.position = c(.2, .77)
    )) +
  plot_annotation(tag_levels = list(c("", "A", "B"))) +
  plot_layout(widths = c(0.04, 1)))
)


if (lag) {
  if (poly) {
    ggsave(
      paste0(
        "envirocor/figs/rdev-condition-corr-timeseries-", scenario,
        "-start", start_year, "-", n_draws, "-draws-bmrs-poly-lag.png"
      ),
      width = 7, height = 3
    )
  } else {
    ggsave(
      paste0(
        "envirocor/figs/rdev-condition-corr-timeseries-", scenario,
        "-start-", start_year, "-", n_draws, "-draws-bmrs-lag.png"
      ),
      width = 7, height = 3
    )
  }
} else {
  if (poly) {
    ggsave(
      paste0(
        "envirocor/figs/rdev-condition-corr-timeseries-", scenario,
        "-start-", start_year, "-", n_draws, "-draws-bmrs-poly2.png"
      ),
      width = 7, height = 3
    )
  } else {
    ggsave(
      paste0(
        "envirocor/figs/rdev-condition-corr-timeseries-", scenario,
        "-start-", start_year, "-", n_draws, "-draws-bmrs2.png"
      ),
      width = 7, height = 3
    )
  }
}


coefs2 <- do.call(rbind, coefs)
head(coefs2)

coefs2 |>
  pivot_longer(1:(ncol(coefs2) - 1), values_to = "est", names_to = "coef") |>
  mutate(coef = factor(coef, levels = c("poly1", "poly2", "slope", "p", "sigma"))) |>
  ggplot() +
  geom_violin(aes(forcats::fct_rev(type), est, fill = type), colour = NA, alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  facet_grid(~coef, scales = "free") +
  labs(x = "", y = "Estimate", colour = "Variable", fill = "Variable") +
  theme(legend.position = "none")


if (lag) {
  if (poly) {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start", start_year, "-", n_draws, "-draws-bmrs-poly-lag.png"), width = 6, height = 1.5)
  } else {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start-", start_year, "-", n_draws, "-draws-bmrs-lag.png"), width = 6, height = 1.5)
  }
} else {
  if (poly) {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start-", start_year, "-", n_draws, "-draws-bmrs-poly.png"), width = 6, height = 1.5)
  } else {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start-", start_year, "-", n_draws, "-draws-bmrs.png"), width = 6, height = 1.5)
  }
}



coefs2 |>
  pivot_longer(1:2, values_to = "est", names_to = "coef") |>
  mutate(coef = factor(coef, levels = c("poly1", "poly2"))) |>
  ggplot() +
  geom_violin(aes(forcats::fct_rev(type), est, fill = type), colour = NA, alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  facet_grid(rows = vars(coef), scales = "free") +
  labs(x = "", y = "Estimate", colour = "Variable", fill = "Variable") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


if (lag) {
  if (poly) {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start", start_year, "-", n_draws, "-draws-bmrs-poly-lag-inset.png"), width = 1.5, height = 2)
  } else {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start-", start_year, "-", n_draws, "-draws-bmrs-lag-inset.png"), width = 1.5, height = 2)
  }
} else {
  if (poly) {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start-", start_year, "-", n_draws, "-draws-bmrs-poly-inset.png"), width = 1.5, height = 2)
  } else {
    ggsave(paste0("envirocor/figs/rdev-condition-corr-coef-violins-", scenario, "-start-", start_year, "-", n_draws, "-draws-bmrs-inset.png"), width = 1.5, height = 2)
  }
}


# correlations with environment ----
dvc <- readRDS("data-generated/envrio-vars-for-condition.rds")
data <- dd1 |>
  rename(group = type, est_st = value) |>
  left_join(dvc, relationship = "many-to-many")

data$response <- data[["est_st"]]
data$var_names <- data[["type"]]

poly <- TRUE

p <- list()
m <- list()
coefs <- list()

n_draws <- 100
colours <- c(7, 5, 3, 2, 4, 9, 6)
pal <- RColorBrewer::brewer.pal(n = 12, name = "Paired")

for (i in seq_along(sort(unique(data$var_names)))) {
  # for(i in 1) {
  dat0 <- filter(data, var_names == sort(unique(data$var_names))[[i]])

  for (g in seq_along(sort(unique(data$group)))) {
    idx <- length(unique(data$group)) * (i - 1) + g

    dat <- filter(dat0, group == unique(data$group)[[g]])

    if (unique(dat$group) == "Immature condition") {
      group <- "imm"
    }
    if (unique(dat$group) == "Female condition") {
      group <- "mat-fem"
    }
    if (unique(dat$group) == "Male condition") {
      group <- "mat-m"
    }

    # retrieve a bunch of `.d` data frames above as MCMC samples from 'response' posterior:
    dd <- purrr::map_dfr(seq_len(n_draws), \(j) {
      .draw <- readRDS(paste0(
        "data-generated/cond-index/", which_cond_model, "/cond-index-sims-", group,
        "-petrale-sole-", which_cond_model, "-20-km.rds"
      )) |>
        filter(
          .iteration == j
        ) |>
        mutate(
          est_st = (.value - 1) / sd(.value),
          group = unique(dat$var_names)
        ) |>
        select(year, est_st)
      .d <- left_join(.draw, dvc)
      .d <- na.omit(.d)

      .d$response <- .d[["est_st"]]
      .d$var_names <- .d[["type"]]
      .d <- filter(.d, var_names == sort(unique(data$var_names))[[i]])
      .d <- .d |> mutate(time = as.integer(seq_along(year)))
      .d$original_iter <- j
      .d
    })

    dd_sum <- dd |>
      group_by(year, value_raw) |>
      summarise(
        max = quantile(response, probs = 0.025),
        min = quantile(response, probs = 0.975),
        response_new_med = median(response)
      )

    m[[idx]] <- tryCatch(brm(
      bf(response ~ poly(value, 2) + ar(time = time)),
      data = dat,
      iter = 1000,
      chains = 3,
      prior =
        c(
          set_prior("normal(0, 1)", class = "ar"),
          set_prior("normal(0, 10)", class = "b"),
          set_prior("student_t(3, 0, 2)", class = "sigma"),
          set_prior("normal(0, 10)", class = "Intercept")
        ),
      backend = "cmdstan"
    ))

    summary(m[[idx]])


    if (max(rhat(m[[idx]])) > 1.01) {
      p[[idx]] <- NULL
    } else {
      nd <- data.frame(value = seq(min(dat$value), max(dat$value), length.out = 200), time = NA)
      using_posteriors <- TRUE
      if (using_posteriors) {
        fits <- dd |>
          split(dd$original_iter) |>
          lapply(do_fit)

        # make predictions for each:
        preds <- fits |> lapply(\(x) {
          posterior_epred(
            x,
            incl_autocor = FALSE,
            re.formula = NA,
            newdata = nd
          )
        })
        pred <- do.call(rbind, preds)

        # plot
        nd$est <- apply(pred, 2, median)
        nd$lwr <- apply(pred, 2, quantile, probs = 0.025)
        nd$upr <- apply(pred, 2, quantile, probs = 0.975)

        nd$value_raw <- nd$value * dd$sd[1] + dd$mean[1]

        # make predictions for each:
        preds_full_posterior <- fits |> lapply(\(x) {
          posterior_predict(
            x,
            incl_autocor = FALSE,
            re.formula = NA,
            newdata = nd
          )
        })
        pred2 <- do.call(rbind, preds_full_posterior)

        nd$lwr2 <- apply(pred2, 2, quantile, probs = 0.025)
        nd$upr2 <- apply(pred2, 2, quantile, probs = 0.975)


        (p[[idx]] <- ggplot() +
          # a place holder to set the axes correctly
          geom_linerange(
            data = dd_sum, aes(value_raw, ymin = min, ymax = max),
            colour = pal[colours[i]]
          ) +
          geom_point(
            data = dat, aes(value_raw, response),
            colour = pal[colours[i]]
          ))

        # combine coefs:
        coefs[[idx]] <- fits |> purrr::map_dfr(\(x) {
          a <- as.data.frame(x)
          data.frame(poly1 = a$b_polyvalue21, poly2 = a$b_polyvalue22)
        })
      } else {
        pred <- posterior_epred(
          m[[idx]],
          incl_autocor = FALSE,
          re.formula = NA,
          newdata = nd
        )

        # plot
        nd$est <- apply(pred, 2, median)
        nd$lwr <- apply(pred, 2, quantile, probs = 0.025)
        nd$upr <- apply(pred, 2, quantile, probs = 0.975)

        nd$value_raw <- nd$value * dd$sd[1] + dd$mean[1]

        # make predictions for each:
        pred2 <- posterior_predict(
          m[[idx]],
          incl_autocor = FALSE,
          re.formula = NA,
          newdata = nd
        )

        nd$lwr2 <- apply(pred2, 2, quantile, probs = 0.025)
        nd$upr2 <- apply(pred2, 2, quantile, probs = 0.975)

        (p[[idx]] <- ggplot() +
          geom_point(data = dat, aes(value_raw, response), pal[colours[i]]))

        # combine coefs:
        a <- as.data.frame(m[[i]])
        coefs[[idx]] <- data.frame(poly1 = a$b_polyvalue21, poly2 = a$b_polyvalue22, ar1 = a$`ar[1]`, sigma = a$sigma)
      }


      (p[[idx]] <- p[[idx]] +
        geom_line(
          data = nd, aes(value_raw, est),
          colour = RColorBrewer::brewer.pal(n = 12, name = "Paired")[colours[i]]
        ) +
        geom_ribbon(
          data = nd, aes(value_raw, ymin = lwr, ymax = upr),
          alpha = 0.5, fill = pal[colours[i]]
        ) +
        geom_ribbon(
          data = nd, aes(value_raw, ymin = lwr2, ymax = upr2),
          alpha = 0.25, fill = pal[colours[i]]
        ) +
        labs(
          x = unique(dat$var_names), y = "",
          colour = "", fill = ""
        ) +
        ggsidekick::theme_sleek()
      )

      if (idx %in% seq_along(unique(data$group))) {
        p[[idx]] <- p[[idx]] + ggtitle(unique(dat$group))
      }

      if (!(idx %in% c(seq(from = 2, to = length(unique(data$group)) * length(unique(data$var_names)), by = 3)))) {
        p[[idx]] <- p[[idx]] + theme(axis.title.x = element_blank())
      }
    }

    # combine coefs:
    coefs[[idx]] <- fits |> purrr::map_dfr(\(x) {
      a <- as.data.frame(x)
      data.frame(poly1 = a$b_polyvalue21, poly2 = a$b_polyvalue22, p = a$`ar[1]`, sigma = a$sigma)
    })

    coefs[[idx]]$group <- unique(dat$group)[1]
    coefs[[idx]]$var_names <- unique(dat$var_names)[1]
  }
}



saveRDS(p, paste0(
  "data-generated/cond-enviro-corr-plot-list-poly-", n_draws, "-draws-",
  length(unique(data$var_names)), "-all-years3.rds"
))

saveRDS(m, paste0(
  "data-generated/cond-enviro-corr-model-list-poly-", n_draws, "-draws-",
  length(unique(data$var_names)), "-all-years3.rds"
))

saveRDS(coefs, paste0("data-generated/cond-enviro-corr-coefs-", n_draws, "-draws-",
                      length(unique(data$var_names)), "-all-years3.rds"))


p <- readRDS(paste0(
  "data-generated/cond-enviro-corr-plot-list-poly-", n_draws, "-draws-",
  length(unique(data$var_names)), "-all-years3.rds"
))

# check convergence
lapply(m, max_rhat)
lapply(m, get_ess)

y_lab_big <- ggplot() +
  annotate(
    geom = "text", x = 1, y = 1, size = 4.5, colour = "grey30",
    label = paste0("Condition index (scaled)"), angle = 90
  ) +
  coord_cartesian(clip = "off") +
  theme_void()


(pp <- ((y_lab_big |
  wrap_plots(gglist = p, ncol = 3) &
    scale_color_viridis_d(option = "D", direction = 1) &
    theme(
      text = element_text(size = 12), plot.background = element_rect(fill = "transparent", color = NA),
      legend.position = "none"
    )) +
  plot_layout(widths = c(0.05, 2)))
)

ggsave(paste0(
  "envirocor/figs/cond-enviro-corr-timeseries-", scenario, "-", n_draws, "-draws-brms-",
  length(unique(data$var_names)), ".png"
), width = 6, height = 10)


coefs2 <- do.call(rbind, coefs)
head(coefs2)

coefs2 |>
  pivot_longer(1:4, values_to = "est", names_to = "coef") |>
  mutate(
    group = factor(group, levels = c("Immature condition", "Male condition", "Female condition"))
  ) |>
  mutate(coef = factor(coef, levels = c("poly1", "poly2", "slope", "p", "ar1", "sigma"))) |>
  ggplot() +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  geom_violin(aes(forcats::fct_rev(var_names), est, fill = var_names), colour = NA, alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = pal[colours]) +
  scale_colour_manual(values = pal[colours]) +
  facet_grid(group ~ coef, scales = "free") +
  labs(x = "", y = "Estimate", colour = "Variable", fill = "Variable") +
  theme(legend.position = "none")

ggsave(paste0(
  "envirocor/figs/cond-enviro-corr-coef-violins-", scenario, "-", n_draws, "-draws-brms-",
  length(unique(data$var_names)), ".png"
), width = 7, height = 4.5)

coefs2 |>
  pivot_longer(1:2, values_to = "est", names_to = "coef") |>
  mutate(
    group = factor(group, levels = c("Immature condition", "Male condition", "Female condition"))
  ) |>
  mutate(coef = factor(coef, levels = c("poly1", "poly2", "slope", "p", "ar1", "sigma"))) |>
  ggplot() +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  geom_violin(aes(forcats::fct_rev(group), est, fill = group), colour = NA, alpha = 0.7) +
  coord_flip() +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  facet_grid(var_names ~ coef, switch = "y", scales = "free") +
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y.left = element_text(angle = 0)
  ) +
  labs(x = NULL, y = "Estimate", colour = "", fill = "")

ggsave(paste0(
  "envirocor/figs/cond-enviro-corr-coef-violins-by-group-", scenario, "-", n_draws, "-draws-brms-",
  length(unique(data$var_names)), ".png"
), width = 7, height = 4.5)



library(GGally)

dvcw <- dvc |>
  select(year, type, value) |>
  pivot_wider(names_from = type, values_from = value)
unique(dvc$type)

ggplot(dvcw) +
  geom_point(aes(`Primary production (Q2: Apr-Jun)`, `Sea floor O2 (Q2: Apr-Jun)`))

ggpairs(dvcw, columns = 2:8)
ggsave("envirocor/figs/condition-variable-correlations.png", width = 15, height = 15)

ggpairs(dvcw, columns = 5:8)
ggsave("envirocor/figs/condition-ROMS-variable-correlations.png", width = 9, height = 9)
