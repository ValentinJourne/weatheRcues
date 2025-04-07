#' This function applies the PSR (Penalized Spline Regression) method to identify weather cues for a specific site using biological data. The function is based on methods adapted from Simmonds et al. It filters data for a particular site and calls the `runing_psr_site` function to analyze the climate and biological data.
#'
#' @param site A string representing the site name, which should match the `plotname.lon.lat` in `seed.data`.
#' @param seed.data A data frame containing biological data (e.g., seed production) with columns such as `plotname.lon.lat`, `Year`, and `log.seed`.
#' @param lastdays Integer specifying the total number of days to include in the analysis (default is 600).
#' @param lastdays Integer specifying the final day for the analysis (usually the length of the time series).
#' @param matrice A numeric vector of length 2 indicating the penalties to apply in the smoothing function of the model (default is `c(3,1)`).
#' @param knots Integer specifying the number of knots for the GAM model. If `NULL` (default), the function will set the knots to the number of years minus one.
#'
#' @details
#' The function filters biological data for the specified site and applies the PSR method to analyze the effect of weather variables (e.g., temperature) on the biological data (e.g., seed production). It uses a rolling window to calculate average climate conditions over a specified period.
#'
#' @return
#' A data frame summarizing the significant weather cues for the site based on the PSR method. If no significant windows are found, the function returns a data frame with `NA` values.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' site <- "site_name"
#' seed.data <- your_fagus_data
#' result <- PSR_function_site(site = site,
#'                             seed.data = seed.data,
#'                             lastdays = 600,
#'                             lastdays = 600,
#'                             matrice = c(3,1),
#'                             knots = NULL)
#' }
#'
#' @export
PSR_function_site <- function(site, 
                              seed.data, 
                              lastdays,
                              refday,
                              climate.path, 
                              matrice = c(3,1),
                              knots = NULL,
                              tolerancedays = 7,
                              yearneed = 2) {
  
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- seed.data %>%
    dplyr::filter(plotname.lon.lat == site)
  
  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub.fagus$plotname.lon.lat), 
    path = climate.path, 
    scale.climate = TRUE
  )
  
  # Run the CSP site analysis function
  runing_psr_site(bio_data = data.sub.fagus,
                  site = site,
                  climate_csv = climate_data,
                  lastdays = 600,
                  refday = refday,
                  rollwin = 1,
                  covariates.of.interest = 'TMEAN',
                  matrice = matrice,
                  knots = knots,
                  yearneed = yearneed,
                  tolerancedays = tolerancedays#,
                  #plot = TRUE
  )
}