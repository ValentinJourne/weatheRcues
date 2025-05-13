#note here it is additional function to extract best windows identified by CSP method
#' Run Climate Signal Processing for a Site
#'
#' This function processes climate signal data for a specific site by fitting a generalized additive model (GAM) and running window-based modeling.
#'
#' @param Results_days A subset of the `Results_CSP` data frame for a specific site.
#' @param bio_data A data frame containing site-level bio_data, typically including seed production and climate-related variables.
#' @param siteneame.forsub A string representing the name of the site to subset from the `Results_CSP` bio_data and climate data.
#' @param refday An integer specifying the reference day (default is 305).
#' @param lastdays An integer representing the last day of the range used for the analysis (default is the maximum of the range).
#' @param rollwin An integer specifying the size of the rolling window used for calculating temperature averages (default is 1).
#'
#' @return A data frame (`output_fit_summary.temp`) containing the results of the climate signal processing, including fitted model coefficients, model statistics, and window ranges.
#'
#' @details
#' The function processes climate signal data for a specific site by running the following steps:
#' \itemize{
#'   \item It extracts slopes and R-squared values from `Results_days` and constructs a temporary data frame for further analysis.
#'   \item It fits a generalized additive model (GAM) to the slopes and R-squared values using the `optimize_and_fit_gam` function.
#'   \item It identifies consecutive sequences of significant days (`extract_consecutive_sequences`) and generates window ranges.
#'   \item It loads the climate data for the site from the specified path, formats the dates, and scales the climate variables.
#'   \item It computes rolling temperature data for the selected time period.
#'   \item It applies the `reruning_windows_modelling` function across all identified window ranges to model the relationship between seed production and climate variables.
#' }
#'
#' @note The function assumes that the `bio_data` argument corresponds to the site-level data and contains a column named `plotname.lon.lat` for site identification.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' output <- runing_csp_site(Results_days = Results_CSP[i],
#'                           bio_data = bio_data,
#'                           siteneame.forsub = "site_123",
#'                           climate.path = "/path/to/climate/data")
#' }
#'
#' @import dplyr
#' @import qs
#' @import purrr
#' @import lubridate
#' @export
runing_csp = function(
  Results_days = Results_days,
  bio_data = bio_data,
  siteneame.forsub = siteneame.forsub,
  option.check.name = TRUE,
  climate_data = climate_data,
  refday = 305,
  lastdays = max(range),
  rollwin = 1,
  optim.k = F,
  formula_model = formula('seed~TMEAN'),
  model_type = 'lm',
  yearneed = 2
) {
  if (
    !is.data.frame(Results_days) &&
      !is_tibble(climate_data) ||
      ncol(Results_days) < 2
  ) {
    stop(
      "Results_days must be a data frame with at least two columns: estimate and r.squared."
    )
  }
  if (!is.data.frame(bio_data) && !is_tibble(bio_data) || nrow(bio_data) == 0) {
    stop("bio_data must be a non-empty data frame.")
  }
  if (!is.data.frame(climate_data) || nrow(climate_data) == 0) {
    stop("climate_data must be a non-empty data frame.")
  }
  if (!is.numeric(lastdays) || length(lastdays) != 1) {
    stop("lastdays must be a numeric value of length 1.")
  }
  if (!is.numeric(refday) || length(refday) != 1) {
    stop("refday must be a numeric value of length 1.")
  }
  if (!is.numeric(yearneed) || length(yearneed) != 1) {
    stop("yearneed must be a numeric value of length 1.")
  }
  if (!is.numeric(rollwin) || length(rollwin) != 1) {
    stop("rollwin must be a numeric value of length 1.")
  }
  if (!is.character(siteneame.forsub) || length(siteneame.forsub) != 1) {
    stop("siteneame.forsub must be a character string of length 1.")
  }
  if (!is.logical(option.check.name) || length(option.check.name) != 1) {
    stop("option.check.name must be a TRUE/FALSE.")
  }
  if (!is.logical(optim.k) || length(optim.k) != 1) {
    stop("optim.k must be a TRUE/FALSE.")
  }

  covariates.of.interest <- all.vars(formula_model)[-1] # Remove the response variable

  list_slope <- as.list(Results_days$estimate)
  list_rs <- as.list(Results_days$r.squared)

  day = seq(rollwin, (lastdays - 1), 1)
  slope = list_slope
  r_s = list_rs
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), day = day, r_s = unlist(r_s))

  #do not optimize, too long
  results <- optimize_and_fit_gam(
    temporary,
    optim.k = optim.k,
    plots = F,
    k = (nrow(bio_data) - 1)
  ) #(12+6) #I will specify just number of month

  days = get_predictions_windows(
    slope_gam = results[[1]],
    rs_gam = results[[2]],
    temporary
  )$days
  sequences_days = extract_consecutive_sequences(days, keep_all = TRUE)
  window_ranges_df <- save_window_ranges(sequences_days) %>%
    mutate(windows.sequences.number = 1:nrow(.))

  if (option.check.name == T) {
    if (
      (siteneame.forsub == unique(bio_data$sitenewname) |
        siteneame.forsub == unique(bio_data$plotname.lon.lat)) ==
        F
    ) {
      stop()
    }
  }

  if (option.check.name == T & is.na(siteneame.forsub) == T) {
    stop(paste("siteneame.forsub argument for subseting site is missing :O "))
  }

  # Define the year period
  yearneed <- yearneed #2
  yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)
  rolling.data <- purrr::map_dfr(
    yearperiod,
    reformat.climate.backtothepast,
    climate_data = climate_data,
    yearneed,
    refday = refday,
    lastdays = lastdays,
    rollwin = rollwin,
    covariates.of.interest = covariates.of.interest
  )
  output_fit_summary.temp <- purrr::map_dfr(
    1:nrow(window_ranges_df),
    ~ reruning_windows_modelling(
      .,
      bio_data = bio_data,
      window_ranges_df = window_ranges_df,
      rolling.data = rolling.data,
      formula_model = formula_model,
      model_type = model_type,
      yearneed
    )
  )
  return(output_fit_summary.temp)
}
