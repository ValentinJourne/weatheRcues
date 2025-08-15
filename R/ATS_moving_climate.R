#' Perform Daily Relationship Analysis Across all Time Series (ATS)
#'
#' This function performs moving window correlation/regression analysis across all sites in a biological dataset.
#' For each site, it loads the corresponding climate data, applies rolling climate smoothing, and runs the `runing_daily_relationship()` function
#' to extract day-wise associations between climate variables and biological responses.
#'
#' @param bio_data_all A data frame containing biological data across multiple sites. Must include columns like `sitenewname`, `plotname.lon.lat`, `log.seed`, and `year`.
#' @param climate.path Character. File path to the folder containing site-level climate data files. File names should include the site identifiers.
#' @param refday Integer. Reference day of the year from which backward windows are computed (e.g., 305 = November 1). Default is 305.
#' @param lastdays Integer. Number of days before `refday` to consider in the rolling window. Default is 700.
#' @param rollwin Integer. Size of the rolling window used to smooth climate variables. Default is 1 (no smoothing).
#' @param formula_model Formula. A formula specifying the biological response and climate predictor (e.g., `log.seed ~ TMEAN`). Default is `log.seed ~ TMEAN`.
#' @param model_type Character. Either `"lm"` or `"betareg"` indicating which model to fit per day. Default is `"lm"`.
#'
#' @details
#' For each site, this function:
#' \itemize{
#'   \item Subsets the biological data for that site.
#'   \item Loads the corresponding climate file from `climate.path` using `format_climate_data()`.
#'   \item Applies `runing_daily_relationship()` to generate slope, R², and correlation time series.
#'   \item Aggregates results across all sites into a single tidy output.
#' }
#'
#' @return A data frame containing results from daily moving window models for each site, including slope estimates, correlation coefficients, and R² values for each day.
#'
#' @examples
#' \dontrun{
#' bio_data <- data.frame(
#'   sitenewname = rep(c("site1", "site2"), each = 10),
#'   plotname.lon.lat = rep(c("site1", "site2"), each = 10),
#'   log.seed = rnorm(20),
#'   year = rep(2000:2009, 2)
#' )
#'
#' results <- ATS_moving_climate(
#'   bio_data_all = bio_data,
#'   climate.path = "data/climate/",
#'   refday = 305,
#'   lastdays = 600,
#'   rollwin = 1,
#'   formula_model = log.seed ~ TMEAN
#' )
#' }
#'
#' @seealso \code{\link{runing_daily_relationship}}, \code{\link{format_climate_data}}
#' @export

ATS_moving_climate <- function(
  bio_data_all,
  climate.path,
  refday = 305,
  lastdays = 700,
  rollwin = 1,
  formula_model = formula('log.seed ~ TMEAN'),
  model_type = 'lm',
  ...
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
      scale.climate = TRUE,
      date_column = "DATEB"
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
