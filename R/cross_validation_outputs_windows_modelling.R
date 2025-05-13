#additional function for the simulation and block cross validation
cross_validation_outputs_windows_modelling = function(
  z,
  bio_data = bio_data,
  window_ranges_df = window_ranges_df,
  rolling.data = rolling.data,
  refday = 305,
  rollwin = 1,
  formula_model = formula('log.seed ~ mean.temperature'),
  model_type = 'lm'
) {
  covariates.of.interest = as.character(formula_model)[3]
  #en gros apres je veux predire avec les valeurs de parametres que jai eu
  window.open <- window_ranges_df$window.open[z]
  window.close <- window_ranges_df$window.close[z]
  window_method <- window_ranges_df$method[z]
  intercept = window_ranges_df$intercept.estimate[z]
  slope = window_ranges_df$slope.estimate[z]

  # Filter the rolling temperature data according to the current window range
  climate_windows_best <- rolling.data %>%
    filter(days.reversed <= window.open & days.reversed >= window.close) %>%
    group_by(LONGITUDE, LATITUDE, year) %>%
    dplyr::summarise(
      covariates.of.interest.aggregate = mean(
        !!sym(covariates.of.interest),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    mutate(window_method = window_method)

  colnames(climate_windows_best)[4] <- covariates.of.interest

  for (col in colnames(bio_data)) {
    if (col == "Year") {
      colnames(bio_data)[colnames(bio_data) == "Year"] <- "year"
    }
  }

  bio.with.window.climate <- bio_data %>%
    left_join(climate_windows_best, by = "year")

  #make predictions
  predictions <- intercept +
    slope * bio.with.window.climate[, covariates.of.interest]

  #get MAE
  absolute_errors = abs(bio.with.window.climate$log.seed - predictions[, 1])
  mae = mean(absolute_errors, na.rm = T) #avoid it is a dataframe
  mae.scaled = mae / sd(bio.with.window.climate$log.seed, na.rm = T)
  r2.p = cor(predictions, bio.with.window.climate$log.seed)^2
  pbiais = Metrics::percent_bias(
    bio.with.window.climate$log.seed,
    predictions[, 1]
  )
  rmse = Metrics::rmse(bio.with.window.climate$log.seed, predictions[, 1])
  rmse.scaled = rmse / sd(bio.with.window.climate$log.seed, na.rm = T)
  rmse.normalized <- rmse / mean(bio.with.window.climate$log.seed, na.rm = T)
  ccc = DescTools::CCC(bio.with.window.climate$log.seed, predictions[, 1])

  out = data.frame(
    sitenewname = unique(bio_data$sitenewname),
    plotname.lon.lat = unique(bio_data$plotname.lon.lat),
    reference.day = refday,
    windows.size = rollwin,
    window.open = window.open,
    window.close = window.close,
    mae = mae,
    scaled.mae = mae.scaled,
    r2.validation = r2.p,
    pbiais = pbiais,
    rmse = rmse,
    scaled.rmse = rmse.scaled,
    rmse.normalized = rmse.normalized,
    ccc.rhos.est = unname(ccc$rho.c[1]),
    ccc.rhos.lowCI = unname(ccc$rho.c[2]),
    ccc.rhos.highCI = unname(ccc$rho.c[3]),
    method.cues = window_method
  )
}
