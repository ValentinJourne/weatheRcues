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