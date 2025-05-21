#' Detect Weather Cue Windows Using Peak Signal Detection
#'
#' This function identifies influential weather cue windows based on a peak signal detection algorithm applied
#' to the product of regression slope and R² values over time. It detects peaks where signal strength significantly
#' deviates from the local mean using a z-score thresholding method, then fits models (e.g., linear regression) over
#' these periods to evaluate biological-climate relationships.
#'
#' @param lag Integer. Number of past days used for calculating the rolling mean and standard deviation in the thresholding algorithm. Default is 100.
#' @param threshold Numeric. Z-score threshold. A signal is detected if the product of slope × R² exceeds the rolling mean by this many standard deviations. Default is 3.
#' @param influence Numeric between 0 and 1. Determines how strongly new signals influence the rolling statistics. `0` means no influence. Default is 0.
#' @param tolerancedays Integer. Maximum gap (in days) between signals to be considered as part of the same window. Default is 7.
#' @param refday Integer. Reference day-of-year (DOY) to count backwards from when building the climate window. Default is 305 (Nov 1).
#' @param lastdays Integer. Number of days to look back from `refday`. Default is 600.
#' @param rollwin Integer. Size of the rolling window for climate data smoothing. Default is 1.
#' @param siteforsub Character. Unique site identifier used for labeling output. Should match identifiers in `bio_data`.
#' @param climate_data Data frame of daily climate data. Must include `year`, `yday`, and the climate variable used in `formula_model`.
#' @param Results_days Data frame. Output from a daily sensitivity analysis (e.g., CSP), must contain `estimate` (slope) and `r.squared` columns.
#' @param bio_data Data frame of biological observations, including `year`, the response variable, and site identifiers such as `plotname.lon.lat`.
#' @param yearneed Integer. Number of years needed prior to each observation to compute rolling climate data. Default is 2.
#' @param formula_model Formula. The model formula to be used for evaluation (e.g., `log.seed ~ TMEAN`).
#' @param model_type Character. Type of model to fit. One of `"lm"` (linear regression) or `"betareg"` (beta regression). Default is `"lm"`.
#'
#' @details
#' The function calculates a signal vector using a rolling z-score approach on the product of slope and R² from daily regressions.
#' Significant signals are grouped into windows, which are then used to aggregate historical climate data and fit predictive models.
#' Model performance metrics are extracted using `reruning_windows_modelling()`.
#'
#' @return A data frame summarizing model fit statistics for each identified cue window, including:
#' \itemize{
#'   \item \code{window.open}, \code{window.close}
#'   \item Model coefficients and significance, and other metrics like R2 and AIC
#'   \item Site identifiers
#' }
#'
#' @seealso \code{\link{Thresholding_algorithm}}, \code{\link{reruning_windows_modelling}}, \code{\link{runing_csp}}, \code{\link{reformat_climate_backtothepast}}
#'
#' @examples
#' \dontrun{
#' result <- runing_peak_detection(
#'   lag = 30,
#'   threshold = 2.5,
#'   influence = 0.1,
#'   climate_data = climate_df,
#'   Results_days = csp_results,
#'   bio_data = seed_data,
#'   formula_model = formula(log.seed ~ TMEAN)
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
    reformat_climate_backtothepast,
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
