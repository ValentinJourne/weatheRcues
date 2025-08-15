#function to format MASTREE data focusing on fagus syvlatica. It is a standard procedure for filtering
#if you want to know more about the reason, please check in Journe et al, New Phytologist, 2023 and Journe et al, Nature Communication, 2023
formatting_mastree_fagus = function(mastreedata) {
  mastreedata %>%
    filter(
      !Variable == "flower" &
        VarType == "C" &
        !Unit == "index" &
        !Variable == "pollen"
    ) %>%
    group_by(Species, VarType) %>%
    mutate(nMax = n()) %>%
    filter(Species == "Fagus sylvatica") %>%
    filter(Year > 1952 & Year < 2023) %>%
    mutate(
      sitenewname = paste0(Alpha_Number, "_", Site_number, "_", Species_code)
    ) %>%
    group_by(sitenewname) %>%
    mutate(
      log.seed = log(1 + Value),
      scale.seed = scale(Value),
      n = n(),
      scaling01 = (Value - min(Value, na.rm = T)) /
        (max(Value) - min(Value, na.rm = T)),
      scalingbeta = y_transformation_betareg(scaling01)
    ) %>%
    filter(n > 19) %>% #initially 14
    mutate() %>%
    mutate(Date = paste0("15/06/", Year)) %>%
    mutate(
      Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
      plotname.lon.lat = paste0(
        "longitude=",
        Longitude,
        "_",
        "latitude=",
        Latitude
      )
    ) %>%
    group_by(plotname.lon.lat) %>%
    ungroup() %>%
    as.data.frame()
}

#format to get collection method of seed production
obtained_cleaned_method_collection = function(data) {
  data %>%
    group_by(sitenewname, Country, Collection_method, Length) %>%
    summarise(
      average.log.seed = mean(log.seed, na.rm = T),
      sd.log.seed = sd(log.seed, na.rm = T)
    ) %>%
    dplyr::select(
      sitenewname,
      Country,
      Collection_method,
      Length,
      average.log.seed,
      sd.log.seed
    ) %>%
    mutate(
      Collection_method = factor(
        dplyr::recode(
          as_factor(Collection_method),
          "seed count" = "Seed count",
          "seed trap" = "Seed trap",
          "visual crop assessment" = "Visual crop",
        ),
        levels = c("Seed count", "Seed trap", "Visual crop")
      )
    ) %>%
    ungroup() %>%
    distinct()
}

#function to make figure 1 plot
Figure1.main = function(Fagus.seed) {
  seed.production.plot = Fagus.seed %>%
    group_by(plotname.lon.lat) %>%
    complete(Year = full_seq(Year, 1)) %>%
    drop_na(Collection_method) %>%
    ggplot(aes(
      x = Year,
      y = log.seed,
      group = plotname.lon.lat
    )) +
    geom_line(na.rm = FALSE, alpha = .2, col = "#008EA0FF") +
    ggpubr::theme_pubr() +
    ylab('Seed production (log)') +
    theme(legend.position = 'none') +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14)
    )

  plot_locations <- st_as_sf(
    Fagus.seed %>%
      dplyr::select(
        Longitude,
        Latitude,
        plotname.lon.lat,
        Collection_method
      ) %>%
      distinct(),
    coords = c("Longitude", "Latitude")
  ) %>%
    st_set_crs(4326)

  mapBase <- getMap(resolution = "high") %>%
    st_as_sf() %>%
    st_make_valid() %>%
    st_crop(mapBase, xmin = -12, xmax = 35, ymin = 35, ymax = 70)

  MapsMastree <- ggplot(data = mapBase$geometry) +
    geom_sf(fill = "grey90", colour = "black", alpha = .8, linewidth = .25) + #plot map of France
    xlab(" ") +
    ylab(" ") + #white or aliceblue
    geom_point(
      data = plot_locations %>%
        mutate(x = unlist(map(geometry, 1)), y = unlist(map(geometry, 2))) %>%
        as.data.frame(),
      aes(x = x, y = y, col = Collection_method, fill = Collection_method),
      stroke = 1,
      alpha = .8,
      shape = 21,
      size = 4.3,
      col = "#008EA0FF",
      fill = "#008EA0FF",
    ) + #size = rmse,
    scale_size_continuous(
      breaks = c(0.1, 0.3, 0.8),
      range = c(0, 7)
    ) +
    coord_sf(xlim = c(-10, 25), ylim = c(40, 60)) +
    scale_y_continuous(breaks = c(40, 50, 60, 70)) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20)) +
    ylab("Latitude") +
    xlab("Longitude") +
    ggpubr::theme_pubr() +
    #scale_color_futurama() +
    #scale_fill_futurama() +
    guides(
      col = "none",
      fill = guide_legend(title = NULL, override.aes = list(size = 8))
    ) +
    theme(
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14)
    )

  print(
    MapsMastree +
      seed.production.plot +
      plot_annotation(tag_levels = 'a') &
      theme(plot.tag = element_text(size = 12))
  )
  cowplot::save_plot(
    here(here("Application_MASTREE/figures/Figure1.png")),
    MapsMastree +
      seed.production.plot +
      plot_annotation(tag_levels = 'a') &
      theme(plot.tag = element_text(size = 12)),
    ncol = 2,
    nrow = 2,
    dpi = 300
  )
}

# Fagus.seed %>%
#   group_by(plotname.lon.lat) %>%
#   complete(Year = full_seq(Year, 1)) %>%
#   drop_na(Collection_method) %>%
#   ggplot(aes(
#     x = Year,
#     y = log.seed,
#     group = plotname.lon.lat,
#     col = Collection_method
#   )) +
#   geom_line(na.rm = FALSE, alpha = .1) +
#   ggpubr::theme_pubr() +
#   ylab('Seed production (log)') +
#   facet_grid(Collection_method ~ ., scales = 'free') +
#   scale_color_futurama() +
#   scale_fill_futurama() +
#   theme(legend.position = 'none')

#function to get block of 30 years for block cross validation
get_30_year_block <- function(years) {
  years <- sort(unique(years))
  for (i in 1:(length(years) - 29)) {
    block <- years[i:(i + 29)]
    if (diff(range(block)) == 29) {
      return(block)
    }
  }
  return(NULL)
}

#do the faction split analysis
#random selection of the data
#Take time to run !
run_sampling_fraction_all_methods <- function(
  years_to_extract = c(5, 10, 15, 20),
  range,
  bio.data,
  site.needed,
  climate.path,
  matrix.psr = c(3, 1)
) {
  #usually 3,1

  # Initialize a list to store results for different year blocks
  results_all_years <- list()

  # Maximum years we need to extract
  max_year_block <- max(years_to_extract)

  library(dplyr)
  library(purrr)

  # Function to get a valid consecutive block
  # get_valid_consecutive_block <- function(years_available, max_year_block) {
  #   # Sort and ensure unique years
  #   years_available <- sort(unique(years_available))
  #
  #   # If the total years available ≤ 30, return all years
  #   if (length(years_available) <= 30) {
  #     return(years_available) # Use full data range
  #   }
  #
  #   # Find valid start positions for a `max_year_block`-long sequence
  #   possible_starts <- which(
  #     (years_available + max_year_block - 1) %in% years_available
  #   )
  #
  #   # If no valid block is found, return the full dataset (fallback)
  #   if (length(possible_starts) == 0) {
  #     return(years_available) # Use full data range
  #   }
  #
  #   # Select a random valid start year
  #   selected_start <- sample(possible_starts, 1)
  #
  #   # Extract the valid block
  #   return(years_available[
  #     selected_start:(selected_start + max_year_block - 1)
  #   ])
  # }

  get_valid_consecutive_block <- function(years_available, block_size) {
    years_available <- sort(unique(years_available))

    if (length(years_available) < block_size) return(NULL)

    # Find all valid start years
    possible_starts <- which(
      (years_available + block_size - 1) %in% years_available
    )

    if (length(possible_starts) == 0) return(NULL)

    # Randomly pick one
    selected_start <- sample(possible_starts, 1)

    return(years_available[selected_start:(selected_start + block_size - 1)])
  }

  # get_valid_consecutive_block <- function(years_available, block_size) {
  #   years_available <- sort(unique(years_available))
  #   if (length(years_available) < block_size) return(NULL)
  #
  #   for (i in 1:(length(years_available) - block_size + 1)) {
  #     block <- years_available[i:(i + block_size - 1)]
  #     if (length(block) == block_size && diff(range(block)) <= block_size + 5) {
  #       return(block)
  #     }
  #   }
  #
  #   return(NULL)
  # }

  # Compute valid site year ranges ensuring consecutive observations
  site_year_ranges <- bio.data %>%
    dplyr::group_by(plotname.lon.lat) %>%
    dplyr::summarise(
      min_year = min(Year),
      max_year = max(Year),
      years_available = list(unique(Year))
    ) %>%
    dplyr::mutate(
      valid_years = purrr::map(years_available, function(years) {
        get_valid_consecutive_block(years, max(years_to_extract)) # Ensure consistent block size
      })
    ) %>%
    dplyr::filter(!is.null(valid_years)) %>% # Remove sites without valid blocks
    dplyr::mutate(random_start_year = purrr::map_dbl(valid_years, min)) # Get start year
  #do extract for consisntent year
  # Iterate over each number of years to extract
  for (year_block in 1:length(years_to_extract)) {
    # Filter data to get consistent blocks
    # seed.data <- bio.data %>%
    #   dplyr::left_join(site_year_ranges, by = "plotname.lon.lat") %>%
    #   dplyr::group_by(plotname.lon.lat) %>%
    #   dplyr::arrange(Year) %>% # Ensure correct order
    #   dplyr::mutate(start_row = which(Year == random_start_year)[1]) %>% # Find starting row
    #   dplyr::filter(
    #     row_number() >= start_row &
    #       row_number() < start_row + years_to_extract[year_block]
    #   ) %>% # Select next `year_block` rows
    #   dplyr::select(-random_start_year, -years_available, -valid_years) %>% # Drop extra columns
    #   dplyr::ungroup() %>%
    #   as.data.frame()
    # Loop across sites and extract site-specific data

    seed.data <- site_year_ranges %>%
      dplyr::select(plotname.lon.lat, valid_years) %>%
      dplyr::mutate(
        selected_years = purrr::map(
          valid_years,
          ~ .x[1:years_to_extract[year_block]]
        )
      ) %>%
      tidyr::unnest(cols = selected_years) %>%
      dplyr::rename(Year = selected_years) %>%
      dplyr::left_join(bio.data, by = c("plotname.lon.lat", "Year"))

    #do it with future map
    statistics_absolute_climwin_frac <- future_map_dfr(
      site.needed,
      ~ {
        # Format the climate data for the site
        climate_data <- format_climate_data(
          site = .x,
          path = climate.path,
          scale.climate = TRUE,
          date_column = "DATEB"
        )

        seed_subset <- seed.data %>%
          dplyr::filter(plotname.lon.lat == .x) %>%
          as.data.frame()
        # Print debugging info
        print(paste("Processing site:", .x))
        print(paste("Number of rows in seed subset:", nrow(seed_subset)))

        # Run climwin per site
        runing_climwin(
          climate_data = climate_data,
          bio_data = seed_subset,
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
    results.moving.site_frac <- ATS_moving_climate(
      bio_data_all = seed.data,
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
      ~ ATS_CSP(
        Results_daily_frac[[.]],
        unique(Results_daily_frac[[.]]$plotname.lon.lat),
        bio_data_all = seed.data,
        climate.path = climate.path,
        refday = 305,
        lastdays = max(range),
        rollwin = 1
      )
    )

    #Run PSR method - to few row is not working
    statistics_psr_method_frac <- map_dfr(
      site.needed,
      ~ ATS_PSR(
        site = .x,
        climate.path = climate.path,
        bio_data_all = seed.data,
        lastdays = max(range),
        refday = 305,
        matrice = matrix.psr,
        knots = NULL
      )
    )

    # Run Basic method
    statistics_basic_method_frac <- map_dfr(
      1:length(Results_daily_frac),
      ~ ATS_peak_detection(
        Results_daily_frac[[.]],
        unique(Results_daily_frac[[.]]$plotname.lon.lat),
        bio_data_all = seed.data,
        climate.path = climate.path,
        refday = 305,
        lastdays = max(range),
        rollwin = 1,
        lag = 100,
        threshold = 3
      )
    )

    # Add method type
    list.all.mm <- list(
      statistics_absolute_climwin_frac %>% dplyr::mutate(method = 'climwin'),
      statistics_csp_method_frac %>% dplyr::mutate(method = 'csp'),
      statistics_psr_method_frac %>% dplyr::mutate(method = 'psr'),
      statistics_basic_method_frac %>% dplyr::mutate(method = 'signal')
    )

    # Add the year block used in sampling
    list.all <- lapply(list.all.mm, function(df) {
      df <- df %>% dplyr::mutate(year_cum_sample = year_block) # Track which year block was used
      return(df)
    })

    # Store results for this year block
    results_all_years[[paste0(
      "years_",
      years_to_extract[year_block]
    )]] <- list.all
  }

  return(results_all_years)
}


#function to run simulation in parallel
simulation_runing_window_parallel <- function(
  climate_data_simulated = climate_data_simulated,
  results_list_fakedata = results_list_fakedata,
  save = TRUE
) {
  library(doParallel)
  library(foreach)
  print(paste0("Number of simulations: ", length(results_list_fakedata)))
  if (
    length(results_list_fakedata) > 10 & length(results_list_fakedata) < 100
  ) {
    print("It might take a few hours.")
  }
  if (length(results_list_fakedata) > 100) {
    print("It might take some days (considering this is parallelized).")
  }

  # Set up parallel backend
  num_cores <- parallel::detectCores() - 3 #just to avoid too much memory :(
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

  # Use foreach for parallel processing
  #needed to add all function and everything ...
  fin.sim <- foreach::foreach(
    k = 1:length(results_list_fakedata),
    .combine = rbind,
    .packages = c("weatheRcues", 'tidyverse', "mgcv"), #export all the new functions :') and the obect
    .export = c(
      "runing_daily_relationship",
      "runing_csp",
      "runing_peak_detection",
      "runing_psr",
      'reformat_climate_backtothepast',
      'reruning_windows_modelling',
      'correlation.spearman.se',
      "optimize_and_fit_gam",
      'get_predictions_windows',
      'extract_consecutive_sequences',
      'save_window_ranges',
      "find_concurrent_period",
      "Thresholding_algorithm",
      "replace_0",
      ls(globalenv()) #just use environment .. because some parameters included in env
    )
  ) %dopar%
    {
      library(tidyverse)
      library(weatheRcues)
      library(mgcv)
      # climwin method
      statistics_climwin_method.simulated <- runing_climwin(
        climate_data = climate_data_simulated,
        bio_data = results_list_fakedata[[k]] %>% dplyr::select(-TMEAN),
        site.name = as.character(k),
        range = c(365, 0),
        cinterval = 'day',
        refday = c(1, 11),
        optionwindows = 'absolute',
        climate_var = 'TMEAN'
      ) %>%
        dplyr::mutate(method = 'climwin')

      # csp and basic methods
      run.sim.day.res <- runing_daily_relationship(
        bio_data = results_list_fakedata[[k]] %>% dplyr::select(-TMEAN),
        climate_data = climate_data_simulated,
        lastdays = 365,
        formula_model = formula('log.seed ~ TMEAN'),
        refday = 305,
        yearneed = 1
      )

      statistics_csp_method.simulated <- runing_csp(
        Results_days = run.sim.day.res,
        bio_data = results_list_fakedata[[k]] %>% dplyr::select(-TMEAN),
        siteneame.forsub = as.character(k),
        option.check.name = TRUE,
        climate_data = climate_data_simulated,
        refday = 305,
        lastdays = 365,
        rollwin = 1,
        optim.k = F,
        formula_model = formula('log.seed~TMEAN'),
        yearneed = 1
      ) %>%
        dplyr::slice(which.max(r2)) %>%
        dplyr::mutate(method = 'csp')

      statistics_basic_method.simulated <- runing_peak_detection(
        lag = 100,
        threshold = 3,
        influence = 0,
        tolerancedays = 7,
        refday = 305,
        lastdays = 365,
        rollwin = 1,
        siteforsub = as.character(k),
        climate_data = climate_data_simulated,
        Results_days = run.sim.day.res,
        formula_model = formula('log.seed ~ TMEAN'),
        bio_data = results_list_fakedata[[k]] %>% dplyr::select(-TMEAN),
        yearneed = 1
      ) %>%
        dplyr::slice(which.max(r2)) %>%
        dplyr::mutate(method = 'signal')

      #psr method
      statistics_psr_method.simulated <- runing_psr(
        bio_data = results_list_fakedata[[k]] %>% dplyr::select(-TMEAN),
        site = as.character(k),
        climate_data = climate_data_simulated,
        lastdays = 365,
        refday = 305,
        rollwin = 1,
        formula_model = formula('log.seed ~ TMEAN'),
        matrice = c(2, 1),
        knots = NULL,
        tolerancedays = 7,
        yearneed = 1
      ) %>%
        dplyr::slice(which.max(r2)) %>%
        dplyr::mutate(method = 'psr')

      # Combine results into a temporary dataframe
      fin.sim.temp <- dplyr::bind_rows(
        statistics_climwin_method.simulated,
        statistics_csp_method.simulated,
        statistics_basic_method.simulated,
        statistics_psr_method.simulated
      ) %>%
        dplyr::select(
          sitenewname,
          reference.day,
          method,
          window.open,
          window.close,
          slope.estimate,
          intercept.estimate,
          AIC,
          r2
        ) %>%
        dplyr::mutate(
          r2.before.simulate = unique(
            results_list_fakedata[[k]]$r2.before.climwin
          ),
          alpha.random.value.setup = unique(
            results_list_fakedata[[k]]$alpha.random.value.setup
          ),
          beta.random.value.setup = unique(
            results_list_fakedata[[k]]$beta.random.value.setup
          ),
          sigma.random.value.setup = unique(
            results_list_fakedata[[k]]$sigma.random.value.setup
          ),
          nb.year.biosimulate = nrow(results_list_fakedata[[k]]),
          simulation = as.character(k)
        )

      fin.sim.temp
    }

  # Stop parallel processing
  parallel::stopCluster(cl)

  if (save == TRUE) {
    qs::qsave(
      fin.sim,
      here(paste0(
        'Application_MASTREE/outputs/result_simulations',
        lubridate::today(),
        '.qs'
      ))
    )
  }
  return(fin.sim)
}

#alternative figure for michal

# datatype20.windows = merged_df %>%
#   group_by(method, source, sitenewname) %>%
#   mutate(
#     r2_highest = ifelse(method %in% c('csp', 'psr', 'signal'), r2 == max(r2), T)
#   ) %>%
#   dplyr::filter(r2_highest) %>%
#   mutate(
#     method = recode(
#       method,
#       "climwin" = "Sliding time window",
#       "csp" = "Climate sensitivity profile",
#       "psr" = "P-spline regression",
#       'signal' = 'Peak signal identification'
#     ),
#     source = recode(
#       source,
#       "years_1" = 5,
#       "years_2" = 10,
#       "years_3" = 15,
#       'years_4' = 20
#     )
#   ) %>%
#   filter(source == 20) %>% #do it for the 20 years
#   left_join(methods.collection.mv2) %>%
#   dplyr::select(method, Collection_method, window.open, window.close) %>%
#   group_by(method, Collection_method) %>%
#   dplyr::select(method, source, window.open, window.close) %>%
#   summarise(
#     wind.open_median = median(window.open, na.rm = TRUE),
#     wind.close_median = median(window.close, na.rm = TRUE),
#     wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
#     wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
#     wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
#     wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
#     wind.open_mean = mean(window.open, na.rm = TRUE),
#     wind.close_mean = mean(window.close, na.rm = T),
#     wind.open_se = se(window.open),
#     wind.close_se = se(window.close)
#   )
#
# datatype20.windows.figure = datatype20.windows %>%
#   dplyr::select(method, Collection_method, wind.open_median:wind.open_q75) %>%
#   pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
#   separate_wider_delim(
#     typewind,
#     delim = '_',
#     names = c('windows.type', 'metric')
#   ) %>%
#   pivot_wider(names_from = 'metric', values_from = value) %>%
#   ggplot(aes(x = Collection_method, group = windows.type, col = windows.type)) +
#   geom_point(
#     aes(y = median),
#     size = 2,
#     position = position_dodge(width = 0.5)
#   ) +
#   geom_pointrange(
#     mapping = aes(y = median, ymin = q25, ymax = q75),
#     position = position_dodge(width = 0.5)
#   ) +
#   coord_flip() +
#   facet_grid(. ~ method, scales = 'free_y') +
#   xlab('') +
#   ylab('Days reversed') +
#   theme(legend.position = 'bottom', legend.title = element_blank()) +
#   geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
#   ggpubr::theme_pubr() +
#   theme(legend.position = 'bottom', legend.title = element_blank())
# datatype20.windows.figure
#
#
#
#   pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
#   ggplot(aes(x = Collection_method      , col = typewind)) +
#   geom_boxplot(aes(y = value)) +
#   coord_flip() +
#   facet_grid(.~method)+
#   xlab('') +
#   ylab('Days reversed') +
#   theme(legend.position = 'bottom',
#         legend.title = element_blank()) +
#   geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')+
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5,6)])+
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5,6)])+
#   ggpubr::theme_pubr()+
#   theme(legend.position = 'bottom',
#         legend.title = element_blank())
# cowplot::save_plot(
#   here("figures/Figure.20yrsdata.subset.windows.png"),
#   datatype20.windows.figure,
#   ncol = 2,
#   nrow = 1.4,
#   dpi = 300
# )

####################################################################################
##############################################################################
#######################

#now test the frac dataset
#take more than 12 hours if 10 iteration (climwin takes some times)
# Initialize a list to store all results
# result_list_all <- list()
#
# # Define the sample fractions you want to iterate over
# fractions <- c(0.3, 0.5, 0.7, 0.9)
#
# # Define the number of iterations for each fraction
# n_iter <- 10 #take > 12 and < 24 hours on mac os
#
# # Outer loop to iterate over the fractions
# for (fraction in fractions) {
#   # Initialize a list to store results for each fraction
#   result_list <- list()
#
#   # Inner loop to perform the operation multiple times for the given fraction
#   for (i in 1:n_iter) {
#     # Get the sample list for the given fraction
#     sample_result <- run_sampling_fraction_all_methods(
#       fraction = fraction,
#       range = range,
#       Fagus.seed = Fagus.seed,
#       beech.site.all = beech.site.all,
#       climate.path = climate.beech.path
#     )
#
#     # Append this sample's results to the overall list for this fraction
#     result_list[[i]] <- sample_result
#   }
#
#   # Store the results for this fraction in the overall list, using the fraction as the key
#   result_list_all[[paste0("fraction_", fraction * 100)]] <- result_list
# }
#
#
# sampling_data_effect <- map_dfr(names(result_list_all), function(sublist) {
#   # extract my sublist
#   new_list <- result_list_all[[sublist]]
#   new_list[[1]] <- lapply(new_list[[1]], rename_columns_if_needed) # Apply renaming to the first list
#
#   # rbind dplyr
#   datatible.sim <- bind_rows(new_list)
#   return(datatible.sim)
# })
#
# #qs::qsave(sampling_data_effect, here('outputs/sampling_data_effect.qs'))
# final_combined_df = qs::qread(here('final_combined_df.qs'))
#
# summary.test.sampling = final_combined_df %>%
#   mutate(
#     method = recode(
#       method,
#       "climwin" = "Sliding time window",
#       "csp" = "Climate sensitivity profile",
#       "psr" = "P-spline regression",
#       'signal' = 'Signal processing'
#     )
#   ) %>%
#   group_by(method, sample_fraction) %>%
#   summarise(
#     wind.open_median = median(window.open, na.rm = TRUE),
#     wind.close_median = median(window.close, na.rm = TRUE),
#     wind.close_q25 = quantile(window.close, 0.3, na.rm = TRUE),
#     wind.close_q75 = quantile(window.close, 0.7, na.rm = TRUE),
#     wind.open_q25 = quantile(window.open, 0.3, na.rm = TRUE),
#     wind.open_q75 = quantile(window.open, 0.7, na.rm = TRUE),
#     wind.open_mean = mean(window.open),
#     wind.close_mean = mean(window.close),
#     wind.open_se = se(window.open),
#     wind.close_se = se(window.close)
#   )
#
#
# sensi.data.plot = summary.test.sampling %>%
#   dplyr::select(method, sample_fraction, wind.open_median:wind.open_q75) %>%
#   pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
#   separate_wider_delim(
#     typewind,
#     delim = '_',
#     names = c('windows.type', 'metric')
#   ) %>%
#   pivot_wider(names_from = 'metric', values_from = value) %>%
#   mutate(
#     sample_fraction = as.numeric(sample_fraction),
#     method_sample = paste(method, sample_fraction, sep = " (n=") %>%
#       paste0(")"),
#     sample_fraction_factor = as_factor(paste0(sample_fraction * 100, '%'))
#   ) %>%
#   ggplot(aes(
#     x = sample_fraction_factor,
#     group = windows.type,
#     col = windows.type
#   )) +
#   geom_point(
#     aes(y = median),
#     size = 2,
#     position = position_dodge(width = 0.5)
#   ) +
#   geom_pointrange(
#     mapping = aes(y = median, ymin = q25, ymax = q75),
#     position = position_dodge(width = 0.5)
#   ) +
#   coord_flip() +
#   facet_wrap(. ~ method, scales = 'free_y') +
#   xlab('Method (Sample Size)') +
#   ylab('Days reversed') +
#   theme(legend.position = 'bottom', legend.title = element_blank()) +
#   geom_hline(yintercept = 498, color = 'black', linetype = 'dashed') +
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
#   scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Paired")[c(5, 6)]) +
#   ggpubr::theme_pubr() +
#   theme(legend.position = 'bottom', legend.title = element_blank())
#
# sensi.data.plot
#
# cowplot::save_plot(
#   here("figures/sensitivity.method.png"),
#   sensi.data.plot,
#   ncol = 1.4,
#   nrow = 1.4,
#   dpi = 300
# )

# window.block.cv = formatting.data.block.cross %>%
#   group_by(method.cues, Collection_method) %>%
#   summarise(
#     wind.open_median = median(window.open, na.rm = TRUE),
#     wind.close_median = median(window.close, na.rm = TRUE),
#     wind.close_q25 = quantile(window.close, 0.25, na.rm = TRUE),
#     wind.close_q75 = quantile(window.close, 0.75, na.rm = TRUE),
#     wind.open_q25 = quantile(window.open, 0.25, na.rm = TRUE),
#     wind.open_q75 = quantile(window.open, 0.75, na.rm = TRUE),
#     wind.open_mean = mean(window.open),
#     wind.close_mean = mean(window.close),
#     wind.open_se = se(window.open),
#     wind.close_se = se(window.close),
#     n.obs = n()
#   ) %>%
#   dplyr::select(
#     method.cues,
#     Collection_method,
#     wind.open_median:wind.open_q75
#   ) %>%
#   pivot_longer(starts_with('wind'), names_to = 'typewind') %>%
#   separate_wider_delim(
#     typewind,
#     delim = '_',
#     names = c('windows.type', 'metric')
#   ) %>%
#   pivot_wider(names_from = 'metric', values_from = value) %>%
#   mutate(
#     method.cues = factor(
#       method.cues,
#       levels = rev(c(
#         "Climate sensitivity profile",
#         "P-spline regression",
#         "Peak signal detection",
#         "Sliding time window"
#       ))
#     ),
#     windows.type = factor(
#       dplyr::recode(
#         as_factor(windows.type),
#         "wind.open" = "Open",
#         "wind.close" = "Close"
#       ),
#       levels = c("Close", "Open")
#     )
#   ) %>%
#   #filter(method.cues=="Climate sensitivity profile") %>%
#   ggplot(aes(
#     x = method.cues,
#     group = windows.type,
#     col = windows.type,
#     shape = windows.type
#   )) +
#   geom_rect(
#     aes(ymin = 519, ymax = 428, xmin = -Inf, xmax = Inf),
#     fill = 'grey',
#     col = 'white',
#     alpha = 0.05
#   ) +
#   geom_point(
#     aes(y = median),
#     size = .001,
#     position = position_dodge(width = 0.2)
#   ) +
#   geom_pointrange(
#     mapping = aes(y = median, ymin = q25, ymax = q75),
#     position = position_dodge(width = 0.2),
#     size = .5
#   ) +
#   facet_grid(. ~ Collection_method) +
#   coord_flip() +
#   xlab('') +
#   ylab('Days reversed') +
#   #scale_color_brewer(palette = "Paired") +
#   #scale_fill_brewer(palette = "Paired") +
#   scale_color_manual(values = c("Open" = "#56B4E9", "Close" = "#D55E00")) +
#   scale_shape_manual(values = c("Open" = 15, "Close" = 19)) +
#   ggpubr::theme_pubr() +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12),
#     axis.text = element_text(size = 13),
#     axis.title = element_text(size = 14)
#   ) +
#   geom_hline(yintercept = 498, color = 'black', linetype = 'dashed')
#
#
# window.block.cv
# cowplot::save_plot(
#   here("figures/methodblock.window.png"),
#   window.block.cv,
#   ncol = 1.6,
#   nrow = 1.4,
#   dpi = 300
# )
