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

#for simulating fake data bio 
simulated_fake_bio_data = function(fakeclimatewindow,
                                   setup.param, 
                                   num_simulations, 
                                   save_fake_data = TRUE,
                                   overwrite=FALSE){
  results_list <- vector("list", num_simulations) #save in list
  for (i in seq_len(num_simulations)) {
    
    # Randomly sample alpha, beta, and sigma from the provided ranges to the uniform
    alpha.random.value <- runif(1, min = setup.param$alpha[1], max = setup.param$alpha[2])
    beta.random.value <- runif(1, min = setup.param$beta[1], max = setup.param$beta[2])
    sigma.random.value <- runif(1, min = setup.param$sigma[1], max = setup.param$sigma[2])
    
    # Simulate seed production, need same june week data, coming from the same climate data simulated
    seed_production_simulated <- fakeclimatewindow %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(TMEAN = mean(TMEAN)) %>%
      # Standardize the mean temperature for proper correlation, as in the main analysis, already scaled before..!
      # Generate log.seed with controlled correlation
      mutate(log.seed = alpha.random.value + beta.random.value * TMEAN + 
               rnorm(n(), mean = 0, sd = sigma.random.value)) %>% #add the random noise based on sigma 
      dplyr::filter(year > startyear) %>%
      dplyr::mutate(year = as.numeric(as.character(year))) %>%
      dplyr::mutate(Date = paste0("15/06/", year)) %>%
      dplyr::mutate(Date2 = strptime(as.character(Date), format = "%d/%m/%Y"))
    
    # Run regression to get R^2 before climwin
    reg.summary.simulation <- summary(lm(log.seed ~ TMEAN, data = seed_production_simulated))
    r2sim.before.climwin <- reg.summary.simulation$adj.r.squared
    
    
    # Bind results and store them
    result <- seed_production_simulated %>% 
      dplyr::bind_cols(r2.before.climwin = r2sim.before.climwin,
                alpha.random.value.setup = alpha.random.value,
                beta.random.value.setup = beta.random.value,
                sigma.random.value.setup = sigma.random.value) %>% 
      dplyr::mutate(sitenewname = as.character(i),
             plotname.lon.lat = as.character(i))
    
    #generate X simulated datasets
    results_list[[i]] <- result
  }
  
  
  
  if (file.exists(here(paste0('outputs/simulated_bio_data_nsim_', num_simulations,'.qs')))) {
    if (overwrite) {
      message("File exists. Overwriting. Hope it was not too bad")
      qs::qsave(results_list, here(paste0('outputs/simulated_bio_data_nsim_', num_simulations,'.qs')))
    } else {
      message("File exists. Skipping overwrite as 'overwrite = FALSE'.")
    }
  } else {
    message("File does not exist. Creating a new file.")
    qs::qsave(results_list, here(paste0('outputs/simulated_bio_data_nsim_', num_simulations,'.qs')))
  }
  
  return(results_list)
  
}
