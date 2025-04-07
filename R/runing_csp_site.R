#note here it is additional function to extract best windows identified by CSP method
#' Run Climate Signal Processing for a Site
#'
#' This function processes climate signal data for a specific site by fitting a generalized additive model (GAM) and running window-based modeling.
#'
#' @param Results_CSPsub A subset of the `Results_CSP` data frame for a specific site. 
#' @param data A data frame containing site-level data, typically including seed production and climate-related variables.
#' @param siteneame.forsub A string representing the name of the site to subset from the `Results_CSP` data and climate data.
#' @param climate.path A string specifying the file path to the folder containing the climate data files.
#' @param refday An integer specifying the reference day (default is 305).
#' @param lastdays An integer representing the last day of the range used for the analysis (default is the maximum of the range).
#' @param rollwin An integer specifying the size of the rolling window used for calculating temperature averages (default is 1).
#'
#' @return A data frame (`output_fit_summary.temp`) containing the results of the climate signal processing, including fitted model coefficients, model statistics, and window ranges.
#'
#' @details 
#' The function processes climate signal data for a specific site by running the following steps:
#' \itemize{
#'   \item It extracts slopes and R-squared values from `Results_CSPsub` and constructs a temporary data frame for further analysis.
#'   \item It fits a generalized additive model (GAM) to the slopes and R-squared values using the `optimize_and_fit_gam` function.
#'   \item It identifies consecutive sequences of significant days (`extract_consecutive_sequences`) and generates window ranges.
#'   \item It loads the climate data for the site from the specified path, formats the dates, and scales the climate variables.
#'   \item It computes rolling temperature data for the selected time period.
#'   \item It applies the `reruning_windows_modelling` function across all identified window ranges to model the relationship between seed production and climate variables.
#' }
#' 
#' @note The function assumes that the `data` argument corresponds to the site-level data and contains a column named `plotname.lon.lat` for site identification.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' output <- runing_csp_site(Results_CSPsub = Results_CSP[i],
#'                           data = seed.data,
#'                           siteneame.forsub = "site_123",
#'                           climate.path = "/path/to/climate/data")
#' }
#'
#' @import dplyr
#' @import qs
#' @import purrr
#' @import lubridate
#' @export
runing_csp_site = function(Results_CSPsub = Results_CSPsub,
                           data = data,
                           siteneame.forsub = siteneame.forsub,
                           option.check.name = TRUE,
                           climate_csv = climate_csv,
                           refday = 305,
                           lastdays = max(range),
                           rollwin = 1,
                           optim.k = F,
                           variablemoving = 'TMEAN',
                           myform.fin = formula('seed~TMEAN'),
                           model_type = 'lm',
                           yearneed = 2){
  
  if (!is.data.frame(Results_CSPsub) && !is_tibble(climate_data) || ncol(Results_CSPsub) < 2) {
    stop("Results_CSPsub must be a data frame with at least two columns: estimate and r.squared.")
  }
  if (!is.data.frame(data) && !is_tibble(data) || nrow(data) == 0) {
    stop("data must be a non-empty data frame.")
  }
  if (!is.data.frame(climate_csv) || nrow(climate_csv) == 0) {
    stop("climate_csv must be a non-empty data frame.")
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
  if (!is.character(variablemoving) || length(variablemoving) != 1) {
    stop("variablemoving must be a character string of length 1.")
  }
  if (!is.logical(option.check.name) || length(option.check.name) != 1) {
    stop("option.check.name must be a TRUE/FALSE.")
  }
  if (!is.logical(optim.k) || length(optim.k) != 1) {
    stop("optim.k must be a TRUE/FALSE.")
  }
  
  list_slope <- as.list(Results_CSPsub$estimate)
  list_rs <- as.list(Results_CSPsub$r.squared)
  
  day = seq(rollwin,(lastdays-1),1)
  slope = list_slope
  r_s = list_rs
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), 
                          day = day, 
                          r_s = unlist(r_s))
  
  #do not optimize, too long 
  results <- optimize_and_fit_gam(temporary, optim.k = optim.k, plots = F, k = (nrow(data)-1) ) #(12+6) #I will specify just number of month 
  
  days = get_predictions_windows(slope_gam = results[[1]], 
                                 rs_gam = results[[2]], 
                                 temporary)$days
  sequences_days = extract_consecutive_sequences(days, keep_all = TRUE)
  window_ranges_df <- save_window_ranges(sequences_days) %>% 
    mutate(windows.sequences.number = 1:nrow(.))
  
  if(option.check.name == T){
    if((siteneame.forsub == unique(data$sitenewname)|siteneame.forsub == unique(data$plotname.lon.lat))==F){
      stop()
    }
  }
  
  if(option.check.name==T & is.na(siteneame.forsub) == T){
    stop(paste("siteneame.forsub argument for subseting site is missing :O "))
  }
  
  
  # Define the year period
  yearneed <- yearneed#2
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  rolling.data <- purrr::map_dfr(yearperiod, reformat.climate.backtothepast, 
                                 climate = climate_csv, 
                                 yearneed, 
                                 refday = refday, 
                                 lastdays = lastdays, 
                                 rollwin = rollwin, 
                                 variablemoving = variablemoving)
  output_fit_summary.temp <- purrr::map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = data, 
                                                                                                  window_ranges_df = window_ranges_df,
                                                                                                  rolling.data = rolling.data,
                                                                                                  myform.fin = myform.fin,
                                                                                                  model_type = model_type,
                                                                                                  yearneed))
  return(output_fit_summary.temp)
  
}