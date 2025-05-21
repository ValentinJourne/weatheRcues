#' Re-run Climate Window Model for a Specific Window Sequence
#'
#' This function re-fits a regression model for a given climate window (defined by index `z`) using previously calculated window bounds.
#' It aggregates the relevant climate covariate over the selected window, merges it with biological data, and fits a model (either linear or beta regression) to assess the effect.
#'
#' @param z Integer. The index (row) of the `window_ranges_df` data frame indicating which climate window to evaluate.
#' @param bio_data Data frame. Biological site-level data containing at least `year`, `sitenewname`, and `plotname.lon.lat`, as well as the response variable in `formula_model`.
#' @param window_ranges_df Data frame. Contains window definitions (e.g., `window.open`, `window.close`) and optional metadata like `windows.sequences.number`.
#' @param rolling.data Data frame. Rolling climate data with columns `days.reversed`, `year`, `LONGITUDE`, `LATITUDE`, and the climate covariate used in `formula_model`.
#' @param formula_model Formula. A model formula (e.g., `log.seed ~ mean.temperature`) specifying the response and climate covariate to be modeled.
#' @param model_type Character. The type of model to fit. Options are `"lm"` (default) for linear regression or `"betareg"` for beta regression.
#' @param refday Integer. The reference day (day of year) from which `days.reversed` were calculated. Included for metadata.
#' @param rollwin Integer. The size of the rolling window applied to the climate data. Included for metadata.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts the corresponding open/close bounds for window `z`.
#'   \item Aggregates the covariate over the selected time window.
#'   \item Joins the climate summary with biological observations by `year`.
#'   \item Fits the model and extracts relevant statistics (e.g., coefficients, R², AIC, etc.).
#' }
#' If `model_type = "betareg"`, ensure the response variable is bounded between 0 and 1.
#'
#' @return A data frame (1 row) summarizing the model results for window `z`, including:
#' \itemize{
#'   \item Site info: \code{sitenewname}, \code{plotname.lon.lat}
#'   \item Window info: \code{window.open}, \code{window.close}, \code{windows.size}, \code{nsequence.id}
#'   \item Model results: \code{intercept.estimate}, \code{slope.estimate}, \code{pvalue}, \code{r2}, \code{AIC}, \code{sigma}, \code{nobs}
#' }
#'
#' @examples
#' \dontrun{
#' result <- reruning_windows_modelling(
#'   z = 1,
#'   bio_data = bio_data,
#'   window_ranges_df = window_ranges_df,
#'   rolling.data = rolling_climate,
#'   formula_model = formula('log.seed ~ TMEAN'),
#'   model_type = 'lm'
#' )
#' }
#'
#' @seealso \code{\link{save_window_ranges}}, \code{\link{reformat_climate_backtothepast}}, \code{\link{DescTools::CCC}}, \code{\link{broom::tidy}}
#' @export

reruning_windows_modelling = function(
  z,
  bio_data = bio_data,
  window_ranges_df = window_ranges_df,
  rolling.data = rolling.data,
  formula_model = formula('log.seed ~ mean.temperature'),
  model_type = 'lm',
  refday = 305,
  rollwin = 1
) {
  if (!is.numeric(z) || length(z) != 1) {
    stop("`z` must be a numeric value of length 1.")
  }

  if (!is.data.frame(bio_data)) {
    stop("`bio_data` must be a data frame or tibble.")
  }

  if (!is.data.frame(window_ranges_df)) {
    stop("`window_ranges_df` must be a data frame or tibble.")
  }

  if (!is.data.frame(rolling.data)) {
    stop("`rolling.data` must be a data frame or tibble.")
  }

  if (!inherits(formula_model, "formula")) {
    stop("`formula_model` must be a formula.")
  }

  if (!model_type %in% c("lm", "betareg")) {
    stop("`model_type` must be either 'lm' or 'betareg'.")
  }

  if (!is.numeric(refday) || length(refday) != 1) {
    stop("`refday` must be a numeric value of length 1.")
  }
  if (!is.numeric(rollwin) || length(rollwin) != 1) {
    stop("`rollwin` must be a numeric value of length 1.")
  }

  covariates.of.interest = as.character(formula_model)[3]

  # Extract the window open and close for the current iteration
  window.open <- window_ranges_df$window.open[z]
  window.close <- window_ranges_df$window.close[z]
  window_number <- window_ranges_df$windows.sequences.number[z]

  # Filter the rolling temperature data according to the current window range
  climate_windows_best <- rolling.data %>%
    dplyr::filter(
      days.reversed <= window.open & days.reversed >= window.close
    ) %>%
    dplyr::group_by(LONGITUDE, LATITUDE, year) %>%
    dplyr::summarise(
      variableofinterest.aggregate = mean(
        !!sym(covariates.of.interest),
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(window_number = window_number) # Add the window number for identification

  #i am chaning col 4
  colnames(climate_windows_best)[4] <- covariates.of.interest

  for (col in colnames(bio_data)) {
    if (col == "Year") {
      colnames(bio_data)[colnames(bio_data) == "Year"] <- "year"
    }
  }

  # Merge temperature data with the site-level data and fit the model
  fit.best.mode.data <- bio_data %>%
    #mutate(year = Year) %>%
    dplyr::left_join(climate_windows_best, by = "year")

  if (model_type == 'lm') {
    fit.best.model = lm(formula_model, data = fit.best.mode.data)
  } else {
    # Ensure the dependent variable is between 0 and 1 for `betareg`
    if (model_type == 'betareg')
      fit.best.model = betareg::betareg(
        formula_model,
        data = fit.best.mode.data
      )
  }

  # Extract model coefficients and statistics
  tidy_model <- broom::tidy(fit.best.model)
  glance_model <- broom::glance(fit.best.model)

  if (model_type == 'lm') {
    r2 = glance_model$r.squared
    sigma = glance_model$sigma
  } else {
    # Ensure the dependent variable is between 0 and 1 for `betareg`
    if (model_type == 'betareg') r2 = glance_model$pseudo.r.squared
    sigma = NA
  }

  # Create a data frame for the results
  data.frame(
    sitenewname = unique(bio_data$sitenewname),
    plotname.lon.lat = unique(bio_data$plotname.lon.lat),
    reference.day = refday,
    windows.size = rollwin,
    window.open = window.open,
    window.close = window.close,
    intercept.estimate = tidy_model$estimate[1],
    intercept.std.error = tidy_model$std.error[1],
    slope.estimate = tidy_model$estimate[2],
    slope.std.error = tidy_model$std.error[2],
    pvalue = tidy_model$p.value[2],
    r2 = r2,
    AIC = glance_model$AIC,
    sigma = sigma,
    nobs = glance_model$nobs,
    nsequence.id = z
  )
}
