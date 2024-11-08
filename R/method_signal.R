#' Thresholding Algorithm for Peak and Valley Detection
#' obtained from Stackoverflow : https://stackoverflow.com/questions/72784873/conditional-peak-valley-signal-detection-in-realtime-timeseries-data-r
#' 
#' This function detects peaks and valleys in a time series using a threshold-based approach. It identifies significant deviations from a moving average using a threshold based on standard deviations and an influence factor to control the response to new signals.
#'
#' @param y A numeric vector representing the time series data.
#' @param lag An integer specifying the number of observations to use for the moving window average and standard deviation. Defaults to 100.
#' @param threshold A numeric value defining the number of standard deviations that a data point must deviate from the moving average to be considered a signal (either peak or valley). Defaults to 3.
#' @param influence A numeric value between 0 and 1 that controls how much influence a new signal has on the filtered time series. If `influence = 1`, the algorithm reacts fully to new signals, while `influence = 0` ignores new signals. Defaults to 0.
#'
#' @details
#' The algorithm works by:
#' \itemize{
#'   \item Calculating a moving average and standard deviation over a specified `lag`.
#'   \item Comparing new data points with the previous moving average.
#'   \item Flagging points that deviate more than `threshold` times the standard deviation from the average as a signal (1 for peak, -1 for valley).
#'   \item Updating the moving average and standard deviation based on the influence of the new signals.
#' }
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item `signals`: A vector where 1 represents a peak, -1 represents a valley, and 0 indicates no signal.
#'   \item `avgFilter`: The filtered moving average at each point.
#'   \item `stdFilter`: The filtered moving standard deviation at each point.
#' }
#'
#' @examples
#' # Simulated example data
#' set.seed(42)
#' data <- cumsum(rnorm(200))  # Random walk data
#' result <- ThresholdingAlgo(y = data, lag = 30, threshold = 3, influence = 0.5)
#'
#' # Plot results
#' plot(data, type = "l", main = "Peak and Valley Detection", col = "black")
#' points(which(result$signals == 1), data[result$signals == 1], col = "blue", pch = 19)  # Peaks
#' points(which(result$signals == -1), data[result$signals == -1], col = "red", pch = 19)  # Valleys
#'
ThresholdingAlgo <- function(y,
                             lag = 100,
                             threshold = 3,
                             influence = 0) {
  # Initialize the signal vector, filtered values, and moving statistics
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  # Compute initial average and standard deviation for the lag period
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  # Loop over the time series, starting from lag+1
  for (i in (lag+1):length(y)){
    # Check if the current observation is a significant deviation from the moving average
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;# Peak detected
      } else {
        signals[i] <- -1;# Valley detected
      }
      # Update the filtered signal using the influence factor
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0 # No signal detected
      filteredY[i] <- y[i] # No influence, use raw value
    }
    # Update the moving average and standard deviation
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

#' Replace Sequences of Zeros in a Vector with Neighboring Values
#'
#' This function replaces sequences of zeros in a vector with the value of the closest non-zero element, provided the sequence length is within a defined threshold.
#'
#' @param vec A numeric vector that may contain sequences of zeros.
#' @param threshold An integer specifying the maximum length of consecutive zeros to replace. If the length of a zero sequence is less than or equal to `threshold`, the zeros are replaced by the nearest non-zero value. Sequences longer than the threshold are not modified.
#'
#' @details
#' The function scans the input vector for sequences of consecutive zeros and replaces them with the nearest preceding or following non-zero value, if the length of the zero sequence is within the specified `threshold`. Sequences longer than `threshold` are left unchanged.
#'
#' @return A modified version of the input vector, where certain sequences of zeros have been replaced by nearby non-zero values.
#'
#' @examples
#' vec <- c(1, 2, 0, 0, 3, 4, 0, 0, 0, 5)
#' replace_0(vec, threshold = 2)
#' # Expected output: c(1, 2, 2, 2, 3, 4, 0, 0, 0, 5)
#'
replace_0 <- function(vec, threshold) {
  # Find the start and end positions of sequences of 0s
  zero_sequences <- rle(vec)
  
  # Initialize a vector to store the results
  result <- vec
  
  # Keep track of the cumulative position in the original vector
  position <- 1
  
  for (i in seq_along(zero_sequences$lengths)) {
    # Check if the current run is 0s and if it's within the threshold
    if (zero_sequences$values[i] == 0 && zero_sequences$lengths[i] <= threshold) {
      # Determine the value to replace 0s with
      replace_value <- NA
      if (i > 1 && zero_sequences$values[i - 1] != 0) {
        replace_value <- zero_sequences$values[i - 1]
      } else if (i < length(zero_sequences$values) && zero_sequences$values[i + 1] != 0) {
        replace_value <- zero_sequences$values[i + 1]
      }
      
      # Replace the 0s with the identified value
      if (!is.na(replace_value)) {
        result[position:(position + zero_sequences$lengths[i] - 1)] <- replace_value
      }
    }
    
    # Update the position to the start of the next run
    position <- position + zero_sequences$lengths[i]
  }
  
  return(result)
}

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

#' Apply Basic Cues Function to a Specific Site
#'
#' This function applies the `runing_basic_cues` method to analyze climate and seed data for a specific site. It identifies weather cues based on the provided results from a CSP analysis and climate data for the site.
#'
#' @param Results_CSPsub A data frame containing the results of the CSP method for the site, including slope and R-squared values.
#' @param siteforsub String. The unique identifier for the site (e.g., longitude and latitude).
#' @param Fagus.seed A data frame containing seed production data for Fagus (beech) trees, with site-specific information such as `plotname.lon.lat` and `Year`.
#' @param climate.beech.path String. The path to the directory containing the climate data for the site.
#' @param refday Integer. The reference day used for the analysis (default is 305).
#' @param lastdays Integer. The total number of days used for the analysis (default is 600).
#' @param rollwin Integer. The size of the rolling window for climate data smoothing (default is 1).
#'
#' @details
#' The function filters the seed data for the specified site and applies the `runing_basic_cues` function to analyze the relationship between climate variables and biological data. It identifies key weather cues by applying a thresholding algorithm to the results of a CSP analysis.
#'
#' @return
#' A data frame summarizing the fitted model for each identified weather cue window. If no signals are detected, the function will return a message indicating that no windows were identified.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' cues_summary <- basiccues_function_site(Results_CSPsub = Results_CSPsub,
#'                                         siteforsub = "longitude=-0.15_latitude=50.85",
#'                                         Fagus.seed = Fagus.seed,
#'                                         climate.beech.path = "path/to/climate/data",
#'                                         refday = 305,
#'                                         lastdays = 600,
#'                                         rollwin = 1)
#' }
#'
#' @seealso \code{\link{runing_basic_cues}}
#'
#' @import dplyr qs lubridate stringr purrr
#' @export
basiccues_function_site <- function(Results_CSPsub, 
                                    siteforsub, 
                                    seed.data, 
                                    climate.path, 
                                    refday, lastdays, rollwin) {
  
  # Filter the Fagus seed data by the site name
  data.sub <- seed.data %>%
    dplyr::filter(plotname.lon.lat == siteforsub)
  
  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub$plotname.lon.lat), 
    path = climate.path, 
    scale.climate = TRUE
  )
  
  # Run the CSP site analysis function
  runing_basic_cues(data = data.sub,
                    siteforsub = siteforsub,
                    lag  = 100,
                    threshold = 3,
                    influence = 0,
                    tolerancedays = 7,
                    refday = 305,
                    lastdays = 600,
                    rollwin = rollwin,
                    climate_csv = climate_data,
                    Results_CSPsub = Results_CSPsub,
                    variablemoving = 'TMEAN',
                    myform.fin = formula('log.seed~TMEAN'),
                    model_type = 'lm')
}

