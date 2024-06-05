# time series regressions for recruitment deviations
library(glmmTMB)
library(brms)
library(tidyverse)
library(patchwork)

# install cmdstanr to use instead of rstan as the backend:
if (FALSE) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  cmdstanr::install_cmdstan()
}

source("envirocor/utils.R")
theme_set(ggsidekick::theme_sleek())

scenario <- 300

# little environmental and age data pre-1995
start_year <- 1995

# last ages are from 2019, so assuming last informed recruitment deviations were 2018
end_year <- 2018
shortlist <- TRUE
# shortlist <- FALSE

p <- list()
m <- list()
coefs <- list()
n_draws <- 100

if (shortlist) {
  colours <- c(5, 3, 2, 7, 8, 6)
  pal <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
} else {
  colours <- c(11, 5, 12, 3, 2, 1, 4, 9, 10, 7, 8, 6)
  pal <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
  # plot(1:length(pal), pch = 20, cex = 4, col = pal)
  # change yellow to be more visible
  pal[11] <- "#E5E74C"
  # plot(1:length(pal), pch = 20, cex = 4, col = pal)
}

if (shortlist) {
  dvs <- readRDS("data-generated/envrio-vars-for-rdevs-shortlist.rds")
} else {
  dvs <- readRDS("data-generated/envrio-vars-for-rdevs.rds")
}

ts <- readRDS(paste0("data-generated/df-", scenario, ".RData"))
data <- left_join(ts, dvs) |> filter(year >= start_year & year <= end_year)
data <- na.omit(data) # not needed but kept as a precaution

data$response <- data[["rdev"]]
data$var_names <- data[["type"]]

for (i in seq_along(sort(unique(data$type)))) {
  # for(i in 1) {
  dat <- filter(data, var_names == sort(unique(data$type))[[i]])

  # retrieve a bunch of `.d` data frames above as MCMC samples from 'response' posterior:
  dd <- purrr::map_dfr(seq_len(n_draws), \(j) {
    .draw <- readRDS(paste0("data-generated/dfs/df", j, ".RData"))
    .d <- left_join(.draw, dvs) |> filter(year >= start_year & year <= end_year)
    .d <- na.omit(.d)
    .d <- .d |>
      group_by(type) |>
      mutate(time = as.integer(seq_along(year))) |>
      ungroup()

    .d$response <- .d[["rdev"]]
    .d$var_names <- .d[["type"]]
    .d <- filter(.d, var_names == sort(unique(data$type))[[i]])
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

  summary(m[[i]])

  if (max(rhat(m[[i]])) > 1.01) {
    p[[i]] <- NULL
  } else {
    fits <- dd |>
      split(dd$original_iter) |>
      lapply(do_fit)

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


    (p[[i]] <- ggplot() +
      # a place holder to set the axes correctly
      geom_linerange(
        data = dd_sum, aes(value_raw, ymin = min, ymax = max),
        colour = pal[colours[i]]
      ) +
      geom_point(
        data = dat, aes(value_raw, response),
        colour = pal[colours[i]]
      ) +
      geom_line(
        data = nd, aes(value_raw, est),
        colour = pal[colours[i]]
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
        x = unique(dd$var_names), y = "",
        colour = "", fill = ""
      ) +
      ggtitle("") +
      ggsidekick::theme_sleek()
    )

    # combine coefs:
    coefs[[i]] <- fits |> purrr::map_dfr(\(x) {
      a <- as.data.frame(x)
      data.frame(poly1 = a$b_polyvalue21, poly2 = a$b_polyvalue22, p = a$`ar[1]`, sigma = a$sigma)
    })
    coefs[[i]]$var_names <- sort(unique(data$var_names))[[i]]
  }
}

if (shortlist) {
  saveRDS(coefs, paste0(
    "data-generated/rdev-enviro-corr-coefs-", n_draws, "-draws-",
    length(unique(data$type)), "-short.rds"
  ))
  saveRDS(p, paste0(
    "data-generated/rdev-enviro-corr-plot-list-", n_draws, "-draws-",
    length(unique(data$type)), "-short.rds"
  ))
  saveRDS(m, paste0(
    "data-generated/rdev-enviro-corr-model-list-", n_draws, "-draws-",
    length(unique(data$type)), "-short.rds"
  ))
} else {
  saveRDS(coefs, paste0(
    "data-generated/rdev-enviro-corr-coefs-", n_draws, "-draws-",
    length(unique(data$type)), ".rds"
  ))
  saveRDS(p, paste0(
    "data-generated/rdev-enviro-corr-plot-list-", n_draws, "-draws-",
    length(unique(data$type)), ".rds"
  ))
  saveRDS(m, paste0(
    "data-generated/rdev-enviro-corr-model-list-", n_draws, "-draws-",
    length(unique(data$type)), ".rds"
  ))
}

# ## load saved
# if(shortlist){
#   coefs <- readRDS(paste0("data-generated/rdev-enviro-corr-coefs-", n_draws, "-draws-",
#                         length(unique(data$type)), "-short.rds"))
#   p <- readRDS(paste0("data-generated/rdev-enviro-corr-plot-list-", n_draws, "-draws-",
#                     length(unique(data$type)), "-short.rds"))
#   m <- readRDS(paste0("data-generated/rdev-enviro-corr-model-list-", n_draws, "-draws-",
#                     length(unique(data$type)), "-short.rds"))
# } else{
#   coefs <- readRDS(paste0("data-generated/rdev-enviro-corr-coefs-", n_draws, "-draws-",
#                         length(unique(data$type)), ".rds"))
#   p <- readRDS(paste0("data-generated/rdev-enviro-corr-plot-list-", n_draws, "-draws-",
#                     length(unique(data$type)), ".rds"))
#   m <- readRDS(paste0("data-generated/rdev-enviro-corr-model-list-", n_draws, "-draws-",
#                     length(unique(data$type)), ".rds"))
# }


lapply(m, get_ess)
lapply(m, max_rhat)


y_lab_big <- ggplot() +
  annotate(
    geom = "text", x = 1, y = 1, size = 5, colour = "grey30",
    label = paste0("Recruitment deviations"), angle = 90
  ) +
  coord_cartesian(clip = "off") +
  theme_void()


(pp <- ((y_lab_big |
  wrap_plots(gglist = p, ncol = 3) &
    scale_color_viridis_d(option = "D", direction = 1) &
    theme(
      text = element_text(size = 12),
      plot.title = element_blank(),
      legend.position = "none",
      plot.tag.position = c(.2, .85)
    )) +
  plot_annotation(tag_levels = list(c(
    "", "A", "B", "C", "D",
    "E", "F", "G", "H", "I",
    "J", "K", "L"
  ))) +
  plot_layout(widths = c(0.05, 2)))
)

if (shortlist) {
  ggsave(paste0(
    "envirocor/figs/rdev-enviro-corr-timeseries-", scenario, "-", n_draws, "-draws-",
    length(unique(data$type)), "-short2.png"
  ), width = 9, height = 6)
} else {
  ggsave(paste0(
    "envirocor/figs/rdev-enviro-corr-timeseries-", scenario, "-", n_draws, "-draws-",
    length(unique(data$type)), "2.png"
  ), width = 10, height = 10)
}

if (!shortlist) {
  coefs2 <- do.call(rbind, coefs)
  head(coefs2)

  g <- coefs2 |>
    pivot_longer(1:4, values_to = "est", names_to = "coef") |>
    mutate(coef = factor(coef, levels = c("poly1", "poly2", "slope", "p", "sigma"))) |>
    ggplot() +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    geom_violin(aes(forcats::fct_rev(var_names), est, fill = var_names), colour = NA, alpha = 0.7) +
    coord_flip() +
    scale_fill_manual(values = pal[colours]) +
    scale_colour_manual(values = pal[colours]) +
    facet_grid(~coef, scales = "free_x") +
    labs(x = "", y = "Estimate", colour = "Variable", fill = "Variable") +
    theme(legend.position = "none")
  g

  ggsave(paste0(
    "envirocor/figs/rdev-enviro-corr-coef-violins-", scenario, "-", n_draws, "-draws-brms-",
    length(unique(data$type)), "-2.png"
  ), width = 8, height = 3)

  coefs2 |>
    pivot_longer(1:2, values_to = "est", names_to = "coef") |>
    mutate(coef = factor(coef, levels = c("poly1", "poly2", "slope", "p", "sigma"))) |>
    ggplot() +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    geom_violin(aes(forcats::fct_rev(var_names), est, fill = var_names), colour = NA, alpha = 0.7) +
    coord_flip() +
    scale_fill_manual(values = pal[colours]) +
    scale_colour_manual(values = pal[colours]) +
    facet_grid(~coef, scales = "free_x") +
    labs(x = "", y = "Estimate", colour = "Variable", fill = "Variable") +
    theme(legend.position = "none")

  ggsave(paste0(
    "envirocor/figs/rdev-enviro-corr-coef-violins-", scenario, "-", n_draws, "-draws-brms-",
    length(unique(data$type)), "-just-poly2.png"
  ), width = 5, height = 3)
}

saveRDS(dd_sum, "data-generated/rdev-uncertainty-range.rds")


if (shortlist) {

  library(GGally)

  dvsw <- dvs |>
    select(year, type, value) |>
    mutate(type = ifelse(type == "Sea surface salinity (Q1-3: all pelagic)", "Surface salinity (all pelagic)", type)) |>
    pivot_wider(names_from = type, values_from = value) ## |> View()
  unique(dvs$type)

  ggpairs(dvsw,
    columns = c(2:7),
    upper = list(continuous = wrap(cor_func, method = "spearman", symbol = expression("\u03C1 ="))),
    progress = FALSE
  )

  ggsave("envirocor/figs/rdev-clim-variables-correlations.png", width = 11, height = 11)
}

