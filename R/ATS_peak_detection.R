#' Perform Peak signal detection Analysis Across all Time Series (ATS)
#'
#' This function performs peak signal detection on climate and seed production data for a specific site
#' using the CSP-derived slope and R² values. The function identifies strong signal periods (weather cue windows)
#' by detecting peaks in the product of slope and R² using a thresholding algorithm. It then evaluates those periods
#' using linear models.
#'
#' @param Results_days A data frame containing daily results from prior CSP analysis. Must include columns `estimate` (slope) and `r.squared`.
#' @param siteforsub Character. Site identifier, typically formatted as `longitude=..._latitude=...`, used to match `bio_data_all`.
#' @param bio_data_all A data frame of biological data across multiple sites. Must include `plotname.lon.lat`, `sitenewname`, `log.seed`, and `Year`.
#' @param climate.path Character. File path to the folder containing site-specific climate `.qs` files.
#' @param refday Integer. Day of year used as the reference point for backwards windowing. Default is 305 (November 1).
#' @param lastdays Integer. Number of days before `refday` to include in the analysis window. Default is 600.
#' @param rollwin Integer. Size of the rolling window to apply to climate data (in days). Default is 1 (no smoothing).
#' @param lag Integer. Number of past days used in rolling mean/std for the thresholding algorithm.
#' @param threshold Numeric. Threshold in standard deviations to classify a signal as a peak.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Filters the biological dataset for a specific site using `siteforsub`.
#'   \item Loads climate data for the matching site using `format_climate_data()`.
#'   \item Applies `runing_peak_detection()` to identify windows with strong signal (high slope × R² product).
#'   \item Evaluates each identified window using regression and returns performance summaries.
#' }
#'
#' The internal peak detection algorithm uses a moving z-score filter to identify signal periods that exceed
#' a given number of standard deviations above the mean.
#'
#' @return A data frame summarizing all identified weather cue windows and their corresponding model results (e.g., slope, R², RMSE). If no window is detected, a message is returned instead.
#'
#' @examples
#' \dontrun{
#' cues_summary <- ATS_peak_detection(
#'   Results_days = Results_CSP[["longitude=3.5_latitude=44.1"]],
#'   siteforsub = "longitude=3.5_latitude=44.1",
#'   bio_data_all = Fagus.seed,
#'   climate.path = "data/climate/",
#'   refday = 305,
#'   lastdays = 600,
#'   rollwin = 1,
#'   lag = 100,
#'   threshold = 3
#' )
#' }
#'
#' @seealso \code{\link{runing_peak_detection}}, \code{\link{Thresholding_algorithm}}, \code{\link{save_window_ranges}}
#' @import dplyr purrr lubridate stringr
#' @export

ATS_peak_detection <- function(
  Results_days,
  siteforsub,
  bio_data_all,
  climate.path,
  refday,
  lastdays,
  rollwin,
  lag,
  threshold
) {
  # Filter the Fagus seed data by the site name
  data.sub <- bio_data_all %>%
    dplyr::filter(plotname.lon.lat == siteforsub)

  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub$plotname.lon.lat),
    path = climate.path,
    scale.climate = TRUE
  )

  # Run the CSP site analysis function
  runing_peak_detection(
    bio_data = data.sub,
    siteforsub = siteforsub,
    lag = lag,
    threshold = threshold,
    influence = 0,
    tolerancedays = 7,
    refday = 305,
    lastdays = 600,
    rollwin = rollwin,
    climate_data = climate_data,
    Results_days = Results_days,
    formula_model = formula('log.seed~TMEAN'),
    model_type = 'lm'
  )
}
