#additional function for the simulation and block cross validation 
cross_validation_outputs_windows_modelling = function(z, 
                                                      tible.sitelevel = tible.sitelevel, 
                                                      window_ranges_df = window_ranges_df, 
                                                      rolling.temperature.data = rolling.temperature.data, 
                                                      refday = 305, 
                                                      rollwin = 1) {
  
  #en gros apres je veux predire avec les valeurs de parametres que jai eu 
  window.open <- window_ranges_df$window.open[z]
  window.close <- window_ranges_df$window.close[z]
  window_method <- window_ranges_df$method[z]
  intercept = window_ranges_df$intercept.estimate[z]
  slope = window_ranges_df$slope.estimate[z]
  
  # Filter the rolling temperature data according to the current window range
  climate_windows_best <- rolling.temperature.data %>%
    filter(days.reversed <= window.open & days.reversed >= window.close) %>%
    group_by(LONGITUDE, LATITUDE, year) %>%
    summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(window_method = window_method)
  
  seed.with.window.climate <- tible.sitelevel %>%
    mutate(year = Year) %>%
    left_join(climate_windows_best, by = "year")
  
  
  #make predictions
  predictions <- intercept + slope * seed.with.window.climate$mean.temperature 
  
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


#for me to change the name from climwin output 
rename_columns_if_needed <- function(df) {
  if ("WindowOpen" %in% colnames(df)) {
    df <- df %>%
      rename(window.open = WindowOpen)
  }
  if ("WindowClose" %in% colnames(df)) {
    df <- df %>%
      rename(window.close = WindowClose)
  }
  return(df)
}

#do the faction split analysis
#random selection of the data 
#Take time to run ! 
run_sampling_fraction_all_methods = function(fraction = .3, 
                                             range, 
                                             Fagus.seed , 
                                             beech.site.all, 
                                             climate.path){
  
  seed.data = Fagus.seed %>% 
    ungroup() %>% 
    group_by(plotname.lon.lat) %>%   
    sample_frac(fraction) %>%        
    ungroup() %>% 
    as.data.frame()
  
  #run climwin 
  library(furrr)#make it in parallele, https://cran.r-project.org/web/packages/future/vignettes/future-1-overview.html
  plan(multisession)#background R sessions (on current machine) , the computer is flying (alt cluster)
  #but because of that you might need to restart R before ... 
  statistics_absolute_climwin_frac <- future_map_dfr(
    beech.site.all,
    ~ {
      # format the climate, .x is the current site
      climate_data <- format_climate_data(site = .x ,
                                          path = climate.beech.path,
                                          scale.climate = T)
      # Run climwin per site
      climwin_site_days(
        climate_data = climate_data,
        data = seed.data %>% filter(plotname.lon.lat == .x),
        site.name = .x,
        range = range,
        cinterval = 'day',
        refday = c(01, 11),
        optionwindows = 'absolute',
        climate_var = 'TMEAN'
      )
    }
  )
  
  results.moving.site_frac = FULL.moving.climate.analysis(seed.data.all = seed.data,
                                                          climate.path = climate.path,
                                                          refday = 305,
                                                          lastdays = max(range),
                                                          rollwin = 1)
  
  
  Results_daily_frac = results.moving.site_frac %>%
    arrange(plotname.lon.lat) %>% 
    group_by(plotname.lon.lat) %>%
    group_split()
  
  name_daily_frac =   results.moving.site_frac %>% 
    arrange(plotname.lon.lat) %>% 
    group_by(plotname.lon.lat) %>%
    group_keys()
  
  # Set names of the list based on group keys - run csp
  names(Results_daily_frac) <- apply(name_daily_frac, 1, paste)
  
  statistics_csp_method_frac = map_dfr(
    1:length(Results_daily_frac), 
    ~CSP_function_site(
      Results_daily_frac[[.]], 
      unique(Results_daily_frac[[.]]$plotname.lon.lat),
      seed.data = seed.data, 
      climate.path = climate.path,
      refday = 305,
      lastdays = max(range),
      rollwin = 1
    ))
  
  
  #run psr 
  statistics_psr_method_frac = map_dfr(
    beech.site.all, 
    ~ PSR_function_site(
      site = .x,
      climate.path = climate.path,
      seed.data = seed.data,  # Use .x correctly here
      tot_days = max(range),
      refday = 305,
      matrice = c(3,1),
      knots = NULL
    )
  )
  
  #run basic 
  statistics_basic_method_frac = map_dfr(
    1:length(Results_daily_frac), 
    ~basiccues_function_site(
      Results_daily_frac[[.]], 
      unique(Results_daily_frac[[.]]$plotname.lon.lat),
      seed.data = seed.data, 
      climate.path = climate.path,
      refday = 305,
      lastdays = 600,
      rollwin = 1
    ))
  
  #add method character for each tible 
  list.all.mm = list(statistics_absolute_climwin_frac%>% dplyr::mutate(method = 'climwin'),
                     statistics_csp_method_frac %>% dplyr::mutate(method = 'csp'),
                     statistics_psr_method_frac %>% dplyr::mutate(method = 'psr'),
                     statistics_basic_method_frac %>% dplyr::mutate(method = 'signal'))
  
  #just add the percentage used 
  list.all = lapply(list.all.mm, function(df) {
    df <- df %>%
      dplyr::mutate(sample_fraction = fraction)  # Add sample_fraction column
    return(df)
  })
  
  return(list.all)
}