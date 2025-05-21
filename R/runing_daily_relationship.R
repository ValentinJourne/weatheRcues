#' Run Daily Rolling Climate-Response Analysis
#'
#' This function performs a daily moving window analysis to estimate the relationship between a biological response (e.g., seed production)
#' and daily climate variables. For each day prior to a reference date, the function computes a rolling average climate value and fits a model
#' (e.g., linear or beta regression) to estimate effect sizes and correlation with the biological variable.
#'
#' @param bio_data A data frame containing the biological response data. Must include the response variable defined in `formula_model`, `year`, and a site identifier.
#' @param climate_data A data frame of daily climate values. Must include columns for `year`, `yday`, and the climate variable defined in `formula_model`.
#' @param lastdays Integer. Number of days before the reference day to include in the rolling analysis. For example, `lastdays = 365` evaluates the full previous year.
#' @param formula_model A formula specifying the model to fit (e.g., `log.seed ~ TMEAN`). The right-hand side is used to extract the relevant climate variable.
#' @param refday Integer. Reference day-of-year (DOY) from which to compute days backwards. Default is 305 (November 1st).
#' @param yearneed Integer. Minimum number of years to retain when evaluating historical climate data. Default is 2.
#' @param model_type Character. Type of model to fit: either `'lm'` (linear regression) or `'betareg'` (beta regression). Default is `'lm'`.
#' @param rollwin Integer. Size of the rolling average window applied to the climate variable. Default is 1 (i.e., daily resolution without smoothing).
#'
#' @details
#' This function aligns climate data with biological observations across years. It applies a rolling average to the climate variable over
#' the specified `lastdays`, centered on each day before `refday`. For each day, the function fits a univariate model using the
#' provided `formula_model`, returning slope estimates, correlation statistics, and model fit metrics.
#'
#' This is the preparatory step for methods like the Climate Sensitivity Profile (CSP) or peak signal detection.
#'
#' @return A data frame containing one row per day (reversed time), with columns:
#' \itemize{
#'   \item \code{days.reversed}: Days before the reference date.
#'   \item \code{estimate}, \code{p.value}, \code{r.squared}: Model coefficients and fit statistics.
#'   \item \code{correlation}, \code{pvalue.cor}: Spearman or Pearson correlation between the response and climate variable.
#' }
#'
#' @seealso \code{\link{runing_csp}}, \code{\link{runing_peak_detection}}, \code{\link{daily_parameters_rolling_climate}}
#'
#' @examples
#' \dontrun{
#' result <- runing_daily_relationship(
#'   bio_data = seed_data,
#'   climate_data = climate_data,
#'   lastdays = 365,
#'   refday = 305,
#'   formula_model = formula('log.seed ~ TMEAN')
#' )
#' }
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
