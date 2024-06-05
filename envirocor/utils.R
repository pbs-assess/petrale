# Get environmental means for a specific set of months and standardize
filter_months <- function(data, months, type_label){
  data |>
    filter(month %in% months) |>
    group_by(year) |> summarise(value = mean(value, na.rm = TRUE)) |>
    filter(!is.na(value)) |>
    mutate(time = seq_along(year),
           value_raw = value,
           mean = mean(value),
           sd = sd(value),
           value = (value - mean(value))/ sd(value),
           type = type_label)
}

get_model_table <- function(model){
  .m <- summary(model)
  .t <- as.data.frame(round(.m$coefficients$cond, 3))
  df <- tibble::rownames_to_column(.t, "Parameter")
  df
}

# brms helper functions

# fit once per response draw:
do_fit <- function(dat, poly = TRUE) {
  if(poly){
    tryCatch(brm(
    bf(response ~  poly(value, 2) + ar(time = time)),
    data = dat,
    iter = 500,
    chains = 1,
    prior =
      # c(set_prior("normal(0, 0.5)", class = "ar"),
      c(set_prior("normal(0, 1)", class = "ar"),
        set_prior("normal(0, 10)", class = "b"),
        set_prior("student_t(3, 0, 2)", class = "sigma"),
        set_prior("normal(0, 10)", class = "Intercept")
      ),
    backend = "cmdstan"
  ))
  } else{
   tryCatch(brm(
      bf(response ~  value + ar(time = time)),
      data = dat,
      iter = 500,
      chains = 1,
      prior =
        # c(set_prior("normal(0, 0.5)", class = "ar"),
        c(set_prior("normal(0, 1)", class = "ar"),
          set_prior("normal(0, 10)", class = "b"),
          set_prior("student_t(3, 0, 2)", class = "sigma"),
          set_prior("normal(0, 10)", class = "Intercept")
        ),
      backend = "cmdstan"
    ))
  }
}

max_rhat <- function(m){max(rhat(m))}

get_ess <- function(m){
  ms <- summary(m)
  min(min(ms$spec_pars$Bulk_ESS),min(ms$spec_pars$Tail_ESS))
}


# hack found on stackoverflow https://stackoverflow.com/questions/37889222/change-colors-in-ggpairs-now-that-params-is-deprecated
cor_func <- function(data, mapping, method, symbol, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)

  corr <- cor(x, y, method=method, use='complete.obs')

  colFn <- colorRampPalette(c("brown1", "white", "dodgerblue"),
                            interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length = 100))]

  ggally_text(
    label = paste(symbol, as.character(round(corr, 2))),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...) +
    theme_void() +
    theme(panel.background = element_rect(fill = fill))
}

