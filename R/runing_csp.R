#' Run Climate Sensitivity Profile (CSP) Window Identification for a Site
#'
#' Identifies the most influential weather cue windows at a given site using the Climate Sensitivity Profile (CSP) method.
#' This method smooths the relationship between climate covariates and a biological response over daily lags using GAMs, and extracts
#' periods of jointly strong slope and R² signals. These windows are then used to fit predictive models and evaluate their performance.
#'
#' @param Results_days A data frame with daily results from a previous CSP analysis. Must contain `estimate` (slope) and `r.squared` for each lag day.
#' @param bio_data A data frame of biological data (e.g., seed production). Must include identifiers such as `plotname.lon.lat`, `sitenewname`, `year`, and the response variable.
#' @param siteneame.forsub A character string identifying the target site. Used to validate data alignment if `option.check.name = TRUE`.
#' @param climate_data A data frame of daily climate data. Must include columns for `year`, `yday`, and the climate variable used in the model.
#' @param refday Integer. The reference day of year (e.g., 305 for Nov 1) from which to count backwards. Default is 305.
#' @param lastdays Integer. Number of days prior to `refday` to include in the climate window. Default is `max(range)`.
#' @param rollwin Integer. Size of the rolling window applied to the climate variable. Default is 1 (no smoothing).
#' @param optim.k Logical. Whether to optimize the GAM smoothness (`k`). Default is `FALSE`. Recommended to be `FALSE` for small time series.
#' @param formula_model A model formula (e.g., `log.seed ~ TMEAN`). Used for downstream model fitting on extracted windows.
#' @param model_type Model type for window-based modeling. Either `'lm'` (default) or `'betareg'`. For now only lm option is available (but we are planing for future option)
#' @param option.check.name Logical. Whether to check if `siteneame.forsub` matches identifiers in `bio_data`. Default is `TRUE`.
#' @param yearneed Integer. Minimum number of years of data to include before modeling. Used to define valid years. Default is 2.
#' @param k.provided Optional. An integer to manually set the number of knots (`k`) for GAM smoothing. If `NA` (default), it uses `nrow(bio_data) - 1`.
#'
#' @details
#' The Climate Sensitivity Profile (CSP) approach builds on daily linear models to estimate the strength of association between a climate variable
#' and a biological response. These daily slopes and R² values are smoothed with GAMs, and time windows are extracted where both metrics show consistent signal.
#'
#' The extracted windows are used to compute rolling means over historical climate data, which are then used to fit simple models on the biological data.
#'
#' @return A data frame summarizing the model fit for each identified weather cue window, including:
#' \itemize{
#'   \item \code{window.open}, \code{window.close}: Start and end of the cue window (in days before `refday`)
#'   \item \code{estimate}, \code{intercept}: Model coefficients from fitting the biological response to the cue
#'   \item \code{r2}, \code{AIC}, \code{nobs}: Model diagnostics and number of observations
#'   \item \code{sitenewname}, \code{plotname.lon.lat}: Site identifiers
#' }
#'
#' @seealso \code{\link{optimize_and_fit_gam}}, \code{\link{get_predictions_windows}},
#' \code{\link{extract_consecutive_sequences}}, \code{\link{reruning_windows_modelling}}
#'
#' @examples
#' \dontrun{
#' runing_csp(
#'   Results_days = my_daily_CSP_results,
#'   bio_data = bio_df,
#'   climate_data = climate_df,
#'   siteneame.forsub = "Site_001",
#'   refday = 305,
#'   lastdays = 600,
#'   formula_model = formula('log.seed ~ TMEAN')
#' )
#' }
#'
#' @export

runing_csp = function(
  Results_days = Results_days,
  bio_data = bio_data,
  siteneame.forsub = siteneame.forsub,
  option.check.name = TRUE,
  climate_data = climate_data,
  refday = 305,
  lastdays = max(range),
  rollwin = 1,
  optim.k = F,
  formula_model = formula('seed~TMEAN'),
  model_type = 'lm',
  yearneed = 2,
  k.provided = NA
) {
  if (
    !is.data.frame(Results_days) &&
      !is_tibble(climate_data) ||
      ncol(Results_days) < 2
  ) {
    stop(
      "Results_days must be a data frame with at least two columns: estimate and r.squared."
    )
  }
  if (!is.data.frame(bio_data) && !is_tibble(bio_data) || nrow(bio_data) == 0) {
    stop("bio_data must be a non-empty data frame.")
  }
  if (!is.data.frame(climate_data) || nrow(climate_data) == 0) {
    stop("climate_data must be a non-empty data frame.")
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
  if (!is.logical(option.check.name) || length(option.check.name) != 1) {
    stop("option.check.name must be a TRUE/FALSE.")
  }
  if (!is.logical(optim.k) || length(optim.k) != 1) {
    stop("optim.k must be a TRUE/FALSE.")
  }

  covariates.of.interest <- all.vars(formula_model)[-1] # Remove the response variable

  list_slope <- as.list(Results_days$estimate)
  list_rs <- as.list(Results_days$r.squared)

  day = seq(rollwin, (lastdays - 1), 1)
  slope = list_slope
  r_s = list_rs
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), day = day, r_s = unlist(r_s))

  #do not optimize, too long
  if (is.na(k.provided)) {
    k.provided = (nrow(bio_data) - 1)
  } else {
    k.provided = k.provided
  }
  results <- optimize_and_fit_gam(
    temporary,
    optim.k = optim.k,
    plots = F,
    k = k.provided
  ) #(12+6) #I will specify just number of month

  days = get_predictions_windows(
    slope_gam = results[[1]],
    rs_gam = results[[2]],
    temporary
  )$days
  sequences_days = extract_consecutive_sequences(days, keep_all = TRUE)
  window_ranges_df <- save_window_ranges(sequences_days) %>%
    mutate(windows.sequences.number = 1:nrow(.))

  if (option.check.name == T) {
    if (
      (siteneame.forsub == unique(bio_data$sitenewname) |
        siteneame.forsub == unique(bio_data$plotname.lon.lat)) ==
        F
    ) {
      stop()
    }
  }

  if (option.check.name == T & is.na(siteneame.forsub) == T) {
    stop(paste("siteneame.forsub argument for subseting site is missing :O "))
  }

  # Define the year period
  yearneed <- yearneed #2
  yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)
  rolling.data <- purrr::map_dfr(
    yearperiod,
    reformat_climate_backtothepast,
    climate_data = climate_data,
    yearneed,
    refday = refday,
    lastdays = lastdays,
    rollwin = rollwin,
    covariates.of.interest = covariates.of.interest
  )
  output_fit_summary.temp <- purrr::map_dfr(
    1:nrow(window_ranges_df),
    ~ reruning_windows_modelling(
      .,
      bio_data = bio_data,
      window_ranges_df = window_ranges_df,
      rolling.data = rolling.data,
      formula_model = formula_model,
      model_type = model_type,
      yearneed
    )
  )
  return(output_fit_summary.temp)
}
