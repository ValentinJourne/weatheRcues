#additional function for the simulation and block cross validation 
cross_validation_outputs_windows_modelling = function(z, 
                                                      tible.sitelevel = tible.sitelevel, 
                                                      window_ranges_df = window_ranges_df, 
                                                      rolling.data = rolling.data, 
                                                      refday = 305, 
                                                      rollwin = 1,
                                                      myform.fin = formula('log.seed ~ mean.temperature'),
                                                      model_type = 'lm') {
  covariates.of.interest = as.character(myform.fin)[3]
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
    dplyr::summarise(variableofinterest.aggregate = mean(!!sym(covariates.of.interest), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(window_method = window_method)
  
  colnames(climate_windows_best)[4]<-covariates.of.interest
  
  for (col in colnames(tible.sitelevel)) {
    if (col == "Year") {
      colnames(tible.sitelevel)[colnames(tible.sitelevel) == "Year"] <- "year"
    }
  }
  
  seed.with.window.climate <- tible.sitelevel %>%
    left_join(climate_windows_best, by = "year")
  
  
  #make predictions
  predictions <- intercept + slope * seed.with.window.climate[,covariates.of.interest] 
  
  #get MAE 
  absolute_errors = abs(seed.with.window.climate$log.seed - predictions)
  mae = mean(absolute_errors, na.rm = T)
  r2.p = cor(predictions, seed.with.window.climate$log.seed)^2
  pbiais = Metrics::percent_bias(seed.with.window.climate$log.seed, predictions)
  
  
  out = data.frame(sitenewname = unique(tible.sitelevel$sitenewname),
                   plotname.lon.lat = unique(tible.sitelevel$plotname.lon.lat),
                   reference.day = refday,
                   windows.size = rollwin,
                   window.open = window.open,
                   window.close = window.close,
                   mae = mae,
                   r2.validation = r2.p,
                   pbiais = pbiais,
                   method.cues = window_method)
}