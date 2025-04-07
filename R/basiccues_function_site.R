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
                                    refday, lastdays, rollwin, 
                                    lag, threshold) {
  
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
                    lag  = lag,
                    threshold = threshold,
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