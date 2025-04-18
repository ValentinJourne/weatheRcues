#' Perform Climate Signal Processing (CSP) Site Analysis
#'
#' This function performs a CSP site analysis by filtering the provided seed data for a specific site
#' and then calling the `runing_csp_site` function to analyze climate and seed data over a specified time window.
#'
#' @param Results_CSPsub A subset of results for a specific site. Typically a data frame or list that contains climate and seed-related estimates.
#' @param siteneame.forsub A character string indicating the site name (longitude and latitude) used to filter the seed dataset.
#' @param seed.data A data frame containing seed data for multiple sites, including columns like `plotname.lon.lat` that are used to filter for a specific site.
#' @param climate.path A string specifying the file path to the climate data related to seed production in trees.
#' @param refday An integer representing the reference day of the year used to align the rolling climate window.
#' @param lastdays An integer representing the total number of days for the climate window analysis.
#' @param rollwin An integer specifying the size of the rolling window for the climate data (number of days to average).
#'
#' @details
#' This function filters the provided `seed.data` data by the site specified in `siteneame.forsub` and 
#' then runs the `runing_csp_site` function with the filtered seed data, climate data, and other parameters. 
#' The `runing_csp_site` function processes climate data within a rolling window to find relationships between 
#' climate signals and seed data.
#'
#' @return
#' A data frame or list, which is the output of the `runing_csp_site` function. It typically contains model estimates, 
#' climate window ranges, and other statistical summaries related to the climate-seed analysis.
#'
#' @examples
#' \dontrun{
#' # Example usage of the CSP_function_site function
#' CSP_function_site(
#'   Results_CSPsub = Results_CSP[[1]], 
#'   siteneame.forsub = "51.0_10.0", 
#'   seed.data = seed.data, 
#'   climate.path = "data/climate_beech", 
#'   refday = 305, 
#'   lastdays = 600, 
#'   rollwin = 1
#' )
#' }
#'
#' @export
CSP_function_site <- function(Results_CSPsub, 
                              siteneame.forsub, 
                              seed.data, 
                              climate.path, 
                              refday, lastdays, rollwin) {
  
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- seed.data %>%
    dplyr::filter(sitenewname == siteneame.forsub | plotname.lon.lat == siteneame.forsub)
  
  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub.fagus$plotname.lon.lat), 
    path = climate.path, 
    scale.climate = TRUE
  )
  
  # Run the CSP site analysis function
  runing_csp_site(Results_CSPsub = Results_CSPsub,
                  data = data.sub.fagus,
                  siteneame.forsub = siteneame.forsub,
                  climate_csv = climate_data,
                  refday = refday,
                  lastdays = lastdays,   # Change from max(range) to lastdays
                  rollwin = rollwin,
                  myform.fin = formula('log.seed~TMEAN'),
                  model_type = 'lm')
}