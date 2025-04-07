#do the faction split analysis
#random selection of the data 
#Take time to run ! 
run_sampling_fraction_all_methods <- function(years_to_extract = c(5, 10, 15, 20),
                                              range, 
                                              Fagus.seed, 
                                              beech.site.all, 
                                              climate.path,
                                              matrix.psr = c(3,1)) { #usually 3,1
  
  # Initialize a list to store results for different year blocks
  results_all_years <- list()
  
  # Maximum years we need to extract
  max_year_block <- max(years_to_extract)
  
  library(dplyr)
  library(purrr)
  
  # Function to get a valid consecutive block
  get_valid_consecutive_block <- function(years_available, max_year_block) {
    # Sort and ensure unique years
    years_available <- sort(unique(years_available))
    
    # If the total years available ≤ 30, return all years
    if (length(years_available) <= 30) {
      return(years_available)  # Use full data range
    }
    
    # Find valid start positions for a `max_year_block`-long sequence
    possible_starts <- which((years_available + max_year_block - 1) %in% years_available)
    
    # If no valid block is found, return the full dataset (fallback)
    if (length(possible_starts) == 0) {
      return(years_available)  # Use full data range
    }
    
    # Select a random valid start year
    selected_start <- sample(possible_starts, 1)
    
    # Extract the valid block
    return(years_available[selected_start:(selected_start + max_year_block - 1)])
  }
  
  
  
  # Compute valid site year ranges ensuring consecutive observations
  site_year_ranges <- Fagus.seed %>%
    dplyr::group_by(plotname.lon.lat) %>%
    dplyr::summarise(
      min_year = min(Year),
      max_year = max(Year),
      years_available = list(unique(Year))
    ) %>%
    dplyr::mutate(valid_years = purrr::map(years_available, function(years) {
      get_valid_consecutive_block(years, max(years_to_extract))  # Ensure consistent block size
    })) %>%
    dplyr::filter(!is.null(valid_years)) %>%  # Remove sites without valid blocks
    dplyr::mutate(random_start_year = purrr::map_dbl(valid_years, min))  # Get start year
  
  
  
  
  #do extract for consisntent year
  # Iterate over each number of years to extract
  for (year_block in 2:length(years_to_extract)) {
    
    # Filter data to get consistent blocks
    seed.data <- Fagus.seed %>%
      dplyr::left_join(site_year_ranges, by = "plotname.lon.lat") %>%
      dplyr::group_by(plotname.lon.lat) %>%
      dplyr::arrange(Year) %>%  # Ensure correct order
      dplyr::mutate(start_row = which(Year == random_start_year)[1]) %>%  # Find starting row
      dplyr::filter(row_number() >= start_row & row_number() < start_row + years_to_extract[year_block]) %>%  # Select next `year_block` rows
      dplyr::select(-random_start_year, -years_available,-valid_years) %>%  # Drop extra columns
      dplyr::ungroup() %>%
      as.data.frame()
    
    
    # Debugging: Print number of rows per site
    #print(paste("Processing:", years_to_extract[year_block], "years"))
    #print(seed.data %>% dplyr::group_by(plotname.lon.lat) %>% dplyr::reframe(Years = range(Year)))
    
    # Run climwin in parallel
    library(furrr) #when using future_map_dfr
    plan(multisession)
    devtools::load_all(here::here())
    #for package in progress
    #devtools::install()
    library(weatheRcues)
    
    statistics_absolute_climwin_frac <- map_dfr(
      beech.site.all,
      ~ {
        # Format the climate data for the site
        climate_data <- format_climate_data(site = .x, path = climate.path, scale.climate = TRUE)
        
        seed_subset <- seed.data %>% dplyr::filter(plotname.lon.lat == .x) %>% as.data.frame()
        # Print debugging info
        print(paste("Processing site:", .x))
        print(paste("Number of rows in seed subset:", nrow(seed_subset)))
        
        # Run climwin per site
        climwin_site_days(
          climate_data = climate_data,
          data = seed_subset,
          site.name = .x,
          range = range,
          cinterval = 'day',
          refday = c(01, 11),
          optionwindows = 'absolute',
          climate_var = 'TMEAN'
        )
      }
    )
    
    # Run moving climate analysis
    results.moving.site_frac <- FULL.moving.climate.analysis(
      seed.data.all = seed.data,
      climate.path = climate.path,
      refday = 305,
      lastdays = max(range),
      rollwin = 1
    )
    
    # Process moving climate results
    Results_daily_frac <- results.moving.site_frac %>%
      dplyr::arrange(plotname.lon.lat) %>% 
      dplyr::group_by(plotname.lon.lat) %>%
      dplyr::group_split()
    
    name_daily_frac <- results.moving.site_frac %>%
      dplyr::arrange(plotname.lon.lat) %>%
      dplyr::group_by(plotname.lon.lat) %>%
      dplyr::group_keys()
    
    names(Results_daily_frac) <- apply(name_daily_frac, 1, paste)
    
    # Run CSP method
    statistics_csp_method_frac <- map_dfr(
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
    
    #Run PSR method - to few row is not working
    statistics_psr_method_frac <- map_dfr(
      beech.site.all,
      ~PSR_function_site(
        site = .x,
        climate.path = climate.path,
        seed.data = seed.data,
        lastdays = max(range),
        refday = 305,
        matrice = c(2,1),
        knots = NULL
      )
    )
    
    # Run Basic method
    statistics_basic_method_frac <- map_dfr(
      1:length(Results_daily_frac), 
      ~basiccues_function_site(
        Results_daily_frac[[.]], 
        unique(Results_daily_frac[[.]]$plotname.lon.lat),
        seed.data = seed.data, 
        climate.path = climate.path,
        refday = 305,
        lastdays = 600,
        rollwin = 1, lag = 100, threshold = 3
      ))
    
    # Add method type
    list.all.mm <- list(
      statistics_absolute_climwin_frac %>% dplyr::mutate(method = 'climwin'),
      statistics_csp_method_frac %>% dplyr::mutate(method = 'csp'),
      statistics_psr_method_frac %>% dplyr::mutate(method = 'psr'),
      statistics_basic_method_frac %>% dplyr::mutate(method = 'signal')
    )
    
    # Add the year block used in sampling
    list.all <- lapply(list.all.mm, function(df) {
      df <- df %>% dplyr::mutate(year_cum_sample = year_block)  # Track which year block was used
      return(df)
    })
    
    # Store results for this year block
    results_all_years[[paste0("years_", year_block)]] <- list.all
  }
  
  return(results_all_years)
}