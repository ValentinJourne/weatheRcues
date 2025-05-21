#' Perform Moving Window Analysis for a Specific Site
#'
#' This function performs a moving window analysis on seed and climate data for a specific site. It involves loading and formatting climate data, applying rolling window functions or not if set rolling at 1 - as done in our study, and then running a moving window analysis to investigate the relationship between seed count and climate variables.
#'
#' @param site.name A string specifying the site name for which to perform the analysis. This should match entries in the `plotname.lon.lat` column of `bio_data`.
#' @param bio_data A data frame containing seed count and site-level information, including `plotname.lon.lat` and `log.seed`.
#' @param climate.path A string specifying the path to the directory containing climate data files. These files should include the site name in their names.
#' @param lastdays An integer specifying the number of days to consider in the rolling window analysis.
#' @param formula_model A formula specifying the model to fit for each moving window. Default is `formula('log.seed~rolling_avg_tmean')`.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Loads and formats the biological data for the specified site.
#'   \item Loads and formats the climate data, including scaling the climate variables.
#'   \item Defines the period for rolling climate data based on a specified number of years.
#'   \item Applies a rolling window function to the climate data for each year in the defined period.
#'   \item Runs a moving window analysis using the `runing.movingwin.analysis` function.
#' }
#'
#' @return A data frame containing the results of the moving window analysis for the specified site. Includes correlation coefficients, model coefficients, and other relevant statistics.
#'
#' @examples
#' # Example usage:
#' site_name <- 'Site1'
#' Fagus_seed <- data.frame(
#'   plotname.lon.lat = rep(c('Site1', 'Site2'), each = 10),
#'   log.seed = rnorm(20),
#'   year = rep(2000:2009, 2)
#' )
#' climate_beech_path <- 'path/to/climate/data'
#' last_days <- 30
#' my_form <- formula('log.seed~rolling_avg_tmean')
#'
#' result <- site.moving.climate.analysis(
#'   site.name = site_name,
#'   bio_data = Fagus_seed,
#'   climate.path = climate_beech_path,
#'   lastdays = last_days,
#'   formula_model = my_form
#' )
#'
#' @export
runing_daily_relationship <- function(
  bio_data,
  climate_data,
  lastdays,
  formula_model,
  refday = 305,
  yearneed = 2,
  model_type = 'lm',
  rollwin = 1
) {
  if (!is.data.frame(bio_data) && !is_tibble(bio_data)) {
    stop("bio_data must be a data frame or tibble.")
  }

  if (!is.data.frame(climate_data) && !is_tibble(climate_data)) {
    stop("climate_data must be a data frame or tibble.")
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

  if (
    !is.character(model_type) ||
      length(model_type) != 1 ||
      !model_type %in% c('lm', 'betareg')
  ) {
    stop("model_type must be either 'lm' or 'betareg'.")
  }

  if (!inherits(formula_model, "formula")) {
    stop(
      "formula_model must be a valid formula. Double check your formula, because both response and predictors are needed well defined"
    )
  }

  if (!as.character(formula_model)[2] %in% colnames(bio_data)) {
    stop(paste(
      "Column",
      as.character(formula_model)[2],
      "not found in bio_data."
    ))
  }
  if (!as.character(formula_model)[3] %in% colnames(climate_data)) {
    stop(paste(
      "Column",
      as.character(formula_model)[3],
      "not found in climate_data"
    ))
  }
  # Define the year period
  #yearneed <- 2
  yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)

  # Apply the function across all years in yearperiod and combine results
  rolling.data <- purrr::map_dfr(
    yearperiod,
    reformat_climate_backtothepast,
    climate = climate_data,
    yearneed = yearneed,
    refday = refday,
    lastdays = lastdays,
    rollwin = rollwin,
    covariates.of.interest = as.character(formula_model)[3]
  )

  results.moving.site = daily_parameters_rolling_climate(
    bio_data = bio_data,
    rolling.data = rolling.data,
    method = 'spearman',
    formula_model = formula_model,
    model_type = model_type
  )

  return(results.moving.site)
}
