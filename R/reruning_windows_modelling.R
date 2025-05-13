#' Re-run climate windows modeling for a specified sequence
#'
#' This function re-runs climate windows modeling for a given index `z`, extracting
#' window ranges, filtering rolling temperature data, and merging it with biological
#' site-level data to fit a linear model (`lm`) using a predefined formula.
#'
#' @param z Integer. The index representing the row number in the `window_ranges_df` data frame, which contains the window ranges for analysis.
#' @param bio_data Data frame. The biological site-level data, which contains
#'   seed production and other relevant variables for modeling. Default is `bio_data`.
#' @param rolling.data Data frame. The climate data containing
#'   rolling temperature averages and other climate variables. Default is `rolling.data`.
#'
#' @return A data frame containing the model results, which include:
#'   - `sitenewname`: The name of the biological site.
#'   - `plotname.lon.lat`: The plot name and its longitude/latitude.
#'   - `reference.day`: The day used as the reference point for windowing.
#'   - `windows.size`: The size of the window in days.
#'   - `window.open`: Start of the window period (days before the reference day).
#'   - `window.close`: End of the window period (days before the reference day).
#'   - `intercept`: The intercept of the fitted linear model.
#'   - `intercept.se`: Standard error of the intercept.
#'   - `estimate`: The coefficient estimate for the temperature variable in the model.
#'   - `estimate.se`: Standard error of the coefficient estimate.
#'   - `pvalue`: P-value for the temperature coefficient.
#'   - `r2`: R-squared value of the model.
#'   - `AIC`: Akaike Information Criterion for the model.
#'   - `nobs`: Number of observations used in the model.
#'   - `nsequence.id`: The index corresponding to the window range used.
#'
#' @details
#' The function:
#' 1. Extracts the `window.open` and `window.close` values from the `window_ranges_df`
#'    for the specified index `z`.
#' 2. Filters `rolling.data` to retrieve climate data within the window range.
#' 3. Aggregates temperature data (mean temperature) for each combination of `LONGITUDE`,
#'    `LATITUDE`, and `year`.
#' 4. Joins the filtered climate data with biological site-level data (`bio_data`).
#' 5. Fits a linear model (`lm`) using the combined data and the predefined formula `formula_model`.
#' 6. Returns a data frame with the model's coefficients, standard errors, p-values, R-squared, AIC,
#'    number of observations, and sequence ID.
#'
#' @examples
#' # Assuming window_ranges_df, rolling.data, and bio_data are defined:
#' result <- reruning_windows_modelling(1, bio_data = bio_data, rolling.data = rolling.data)
#' print(result)
#'
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
