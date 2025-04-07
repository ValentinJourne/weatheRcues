#' Identify Basic Weather Cues Using a Thresholding Algorithm
#'
#' This function identifies weather cues based on the `ThresholdingAlgo` for a specific site. It computes the relationship between slope and R-squared of a CSP method result, applies a thresholding algorithm to detect signal changes, and analyzes the climate data using a rolling window.
#'
#' @param lag Integer. The number of data points to consider for calculating thresholds in the thresholding algorithm (default is 100).
#' @param threshold Numeric. The threshold for triggering signals in the `ThresholdingAlgo` function (default is 3).
#' @param influence Numeric. The influence of new signals on the algorithm; set to 0 to avoid influence (default is 0).
#' @param tolerancedays Integer. The number of days for gap tolerance when identifying sequences (default is 7).
#' @param refday Integer. The reference day for climate data (default is 305).
#' @param lastdays Integer. The total number of days used for the analysis (default is 600).
#' @param rollwin Integer. The size of the rolling window for climate data smoothing (default is 1).
#' @param siteforsub String. The unique identifier for the site (e.g., longitude and latitude).
#' @param climate.beech.path String. The path to the directory containing the climate data.
#' @param Results_CSPsub A data frame containing the results of a CSP method for a specific site, including slope and R-squared values.
#' @param data A data frame containing biological data (e.g., seed production) with columns such as `plotname.lon.lat` and `Year`.
#'
#' @details
#' The function applies the thresholding algorithm to detect signals from the product of slope and R-squared values from a CSP analysis. Identified signals are used to extract sequences of days, and these sequences are analyzed with respect to climate data using a rolling window approach.
#'
#' @return
#' A data frame summarizing the fitted model for each identified weather cue window. If no signals are detected, the function will return a message indicating that no windows were identified.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' basic_cues <- runing_basic_cues(lag = 100, 
#'                                 threshold = 3, 
#'                                 siteforsub = "longitude=-0.15_latitude=50.85",
#'                                 climate.beech.path = "path/to/climate/data", 
#'                                 Results_CSPsub = Results_CSPsub,
#'                                 data = bio_data)
#' }
#'
#' @export
runing_basic_cues = function(lag  = 100,
                             threshold = 3,
                             influence = 0,
                             tolerancedays = 7,
                             refday = 305,
                             lastdays = 600,
                             rollwin = 1,
                             siteforsub = "longitude=-0.15_latitude=50.85",
                             climate_csv = climate_csv,
                             Results_CSPsub = Results_CSPsub,
                             data = data,
                             variablemoving = 'TMEAN',
                             yearneed = 2,
                             myform.fin=formula('log.seed ~ mean.temperature'),
                             model_type = 'lm'){
  
  #mostly similar strucutre to CSP methods 
  list_slope <- as.list(Results_CSPsub$estimate)
  list_rs <- as.list(Results_CSPsub$r.squared)
  day = seq(rollwin,(lastdays-1),1)
  slope = list_slope
  r_s = list_rs
  
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), 
                          day = day, 
                          r_s = unlist(r_s))
  
  data_beforesign = temporary %>% 
    mutate(combo_r2_slope = r_s*slope) 
  
  y = data_beforesign$combo_r2_slope
  
  #use algo found on stakoverflow 
  result <- ThresholdingAlgo(y,lag,threshold,influence)
  
  #now made them together
  formatzek = data.frame(day = seq(rollwin,(lastdays-1),1),
                         signal = replace_0((result$signals), threshold = tolerancedays)) %>% 
    dplyr::filter(signal !=0)
  
  window.basic = as.vector(formatzek$day)
  
  #tolerance of 20 days gaps between 0 
  sequences_days = extract_sequences_auto(window.basic, tolerance = tolerancedays) 
  #i guess now will do the same shit as the other methods 
  window_ranges_df <- save_window_ranges(sequences_days) %>% 
    mutate(windows.sequences.number = 1:nrow(.))
  
  yearneed <- yearneed#2
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  
  # Apply the function across all years in yearperiod and combine results
  rolling.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                          climate = climate_csv, 
                          yearneed, 
                          refday = refday, 
                          lastdays = lastdays, 
                          rollwin = 1, 
                          variablemoving = variablemoving)
  
  output_fit_summary.temp.basic <- map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = data, 
                                                                                                 window_ranges_df = window_ranges_df,
                                                                                                 rolling.data = rolling.data,
                                                                                                 myform.fin = myform.fin,
                                                                                                 yearneed,
                                                                                                 model_type = model_type))
  
  
}