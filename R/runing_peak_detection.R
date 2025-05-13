#' Identify Weather Cue Windows Using Threshold-Based Peak Detection
#'
#' This function identifies potential weather cue windows from climate data using a thresholding algorithm
#' applied to the product of slope and R² values obtained from a prior CSP (climate sensitivity profile) analysis.
#' Peaks are identified when this product exceeds a specified threshold of deviation from a moving average.
#'
#' @param lag Integer. The number of past observations used to compute the rolling mean and standard deviation
#' in the thresholding algorithm. Larger values smooth the baseline. Default is 100.
#' @param threshold Numeric. The number of standard deviations a value must exceed the rolling mean to be
#' considered a signal. Default is 3.
#' @param influence Numeric between 0 and 1. Determines the influence of detected signals on the rolling
#' statistics; 0 means no influence. Default is 0 (fully robust).
#' @param tolerancedays Integer. Number of days allowed as gaps within a signal to still consider it part of
#' a continuous window. Default is 7.
#' @param refday Integer. Day of year used as a reference (e.g., 305 = Nov 1). Default is 305.
#' @param lastdays Integer. Number of days before `refday` to include in the analysis. Default is 600.
#' @param rollwin Integer. Window size for rolling average in climate data. Default is 1.
#' @param siteforsub Character. Unique site identifier (e.g., "longitude=16.71_latitude=53.09").
#' @param climate_data Data frame. Daily climate data, must include a `year` column and the covariate(s) used
#' in modeling (e.g., `TMEAN`).
#' @param Results_days Data frame. Output from CSP analysis including slope (`estimate`) and R² (`r.squared`)
#' columns for each day.
#' @param bio_data Data frame. Biological data for the site (e.g., seed production), must include `plotname.lon.lat`,
#' `Year`, and modeled response variable (e.g., `log.seed`).
#' @param yearneed Integer. Number of years required for cross-validation or fitting. Default is 2.
#' @param formula_model Formula. The model to fit (e.g., `log.seed ~ TMEAN`). Default is `log.seed ~ TMEAN`.
#' @param model_type Character. Type of model to fit. Either `"lm"` or `"gam"`. Default is `"lm"`.
#'
#' @details
#' The function applies a z-score-based peak detection algorithm to identify the most informative time windows
#' based on signal strength (product of slope and R² from CSP). These windows are then evaluated using linear or
#' generalized additive models over historical climate data.
#'
#' @return A data frame with model summaries (e.g., slope, R², AIC) for each identified cue window. If no
#' significant windows are found, the function returns an empty object or message.
#'
#' @examples
#' \dontrun{
#' result <- runing_peak_detection_site(
#'   lag = 100,
#'   threshold = 3,
#'   influence = 0,
#'   siteforsub = "longitude=16.71_latitude=53.09",
#'   climate_data = climate_df,
#'   Results_days = csp_output,
#'   bio_data = seed_data
#' )
#' }
#'
#' @export
runing_peak_detection = function(
  lag = 100,
  threshold = 3,
  influence = 0,
  tolerancedays = 7,
  refday = 305,
  lastdays = 600,
  rollwin = 1,
  siteforsub = "longitude=-0.15_latitude=50.85",
  climate_data = climate_data,
  Results_days = Results_days,
  bio_data = bio_data,
  yearneed = 2,
  formula_model = formula('log.seed ~ TMEAN'),
  model_type = 'lm'
) {
  #mostly similar strucutre to CSP methods
  list_slope <- as.list(Results_days$estimate)
  list_rs <- as.list(Results_days$r.squared)
  day = seq(rollwin, (lastdays - 1), 1)
  slope = list_slope
  r_s = list_rs

  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), day = day, r_s = unlist(r_s))

  data_beforesign = temporary %>%
    mutate(combo_r2_slope = r_s * slope)

  y = data_beforesign$combo_r2_slope

  #use algo found on stakoverflow
  result <- Thresholding_algorithm(y, lag, threshold, influence)

  #now made them together
  formatzek = data.frame(
    day = seq(rollwin, (lastdays - 1), 1),
    signal = replace_0((result$signals), threshold = tolerancedays)
  ) %>%
    dplyr::filter(signal != 0)

  window.basic = as.vector(formatzek$day)

  #tolerance of 20 days gaps between 0
  sequences_days = extract_sequences_auto(
    window.basic,
    tolerance = tolerancedays
  )
  #i guess now will do the same shit as the other methods
  window_ranges_df <- save_window_ranges(sequences_days) %>%
    mutate(windows.sequences.number = 1:nrow(.))

  yearneed <- yearneed #2
  yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)

  # Apply the function across all years in yearperiod and combine results
  covariates.of.interest <- all.vars(formula_model)[-1] # Remove the response variable

  rolling.data <- map_dfr(
    yearperiod,
    reformat.climate.backtothepast,
    climate_data = climate_data,
    yearneed,
    refday = refday,
    lastdays = lastdays,
    rollwin = 1,
    covariates.of.interest = covariates.of.interest
  )

  output_fit_summary.temp.basic <- map_dfr(
    1:nrow(window_ranges_df),
    ~ reruning_windows_modelling(
      .,
      bio_data = bio_data,
      window_ranges_df = window_ranges_df,
      rolling.data = rolling.data,
      formula_model = formula_model,
      yearneed,
      model_type = model_type
    )
  )
}
