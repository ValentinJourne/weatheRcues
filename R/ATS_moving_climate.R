#' Main Function of correlation - regression and process accross all sites
#'
#' This function applies the moving climate analysis to all sites listed in the `seed.data` data frame. It iterates over each site, performs the moving window analysis, and combines the results into a single data frame.
#'
#' @param seed.data A data frame containing seed count and site-level information. It should include a column `plotname.lon.lat` to specify site names and a column `log.seed` for seed counts.
#' @param climate.path A string specifying the path to the directory containing climate data files. These files should be named according to the site names.
#' @param refday An integer representing the reference day for the rolling window analysis. Default is 305.
#' @param lastdays An integer specifying the number of days to include in the rolling window analysis. Default is the maximum of a specified range.
#' @param rollwin An integer specifying the rolling window size. Default is 1.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Identifies unique site names from the `seed.data` data frame.
#'   \item Applies the `site.moving.climate.analysis` function to each site using the unique site names.
#'   \item Combines the results from all sites into a single data frame.
#' }
#'
#' @return A data frame containing the results of the moving window analysis for all sites. Includes correlation coefficients, model coefficients, and other relevant statistics for each site.
#'
#' @examples
#' # Example usage:
#' Fagus_seed <- data.frame(
#'   plotname.lon.lat = rep(c('Site1', 'Site2'), each = 10),
#'   log.seed = rnorm(20),
#'   year = rep(2000:2009, 2)
#' )
#' climate_beech_path <- 'path/to/climate/data'
#' result <- FULL.moving.climate.analysis(
#'   seed.data = Fagus_seed,
#'   climate.path = climate_beech_path,
#'   refday = 305,
#'   lastdays = 30,
#'   rollwin = 1
#' )
#'
ATS_moving_climate <- function(
  bio_data_all,
  climate.path,
  refday = 305,
  lastdays = 700,
  rollwin = 1,
  formula_model = formula('log.seed ~ TMEAN'),
  model_type = 'lm'
) {
  # Get the list of unique site names from seed data
  al.sites <- unique(bio_data_all$sitenewname)

  # Iterate over each site and perform the moving climate analysis
  results.moving <- purrr::map_dfr(al.sites, function(site.name) {
    # Filter biological data for the current site
    bio_data <- bio_data_all %>%
      dplyr::filter(sitenewname == site.name) %>%
      as.data.frame()

    # format climate
    climate_data <- format_climate_data(
      site = unique(bio_data$plotname.lon.lat),
      path = climate.path,
      scale.climate = TRUE
    )

    # run csp site level
    runing_daily_relationship(
      bio_data = bio_data,
      climate_data = climate_data,
      lastdays = lastdays,
      refday = refday,
      rollwin = 1,
      formula_model = formula_model,
      model_type = model_type
    )
  })

  return(results.moving)
}
