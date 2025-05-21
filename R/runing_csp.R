#' Run Climate Sensitivity Profile (CSP) Window Identification for a Site
#'
#' This function identifies the most influential weather cue windows at a site using the Climate Sensitivity Profile (CSP) method. It fits generalized additive models (GAMs) to the relationship between slope and $R^2$ values over a series of daily windows and identifies periods where both metrics jointly indicate a strong climate signal. It then extracts these windows and fits linear models to the biological response using climate data smoothed with a rolling average.
#'
#' @param Results_days A data frame containing at least columns `estimate` (slope) and `r.squared` (model $R^2$) for each lagged day.
#' @param bio_data A data frame containing biological observations (e.g., seed production), with identifiers such as `plotname.lon.lat`, `sitenewname`, and the response variable (e.g., `log.seed`).
#' @param siteneame.forsub Character. The unique site name to be validated against `bio_data`.
#' @param climate_data A data frame containing daily climate values for the site, with `year`, `yday`, and the climate variable of interest.
#' @param refday Integer. The reference day of year (DOY) from which to compute backward days (default is 305, i.e., November 1st).
#' @param lastdays Integer. Number of days to look back from the reference day. Default is the maximum of the defined range.
#' @param rollwin Integer. Size of the rolling window for smoothing the climate variable. Default is 1 (no smoothing).
#' @param optim.k Logical. Whether to optimize the GAM knot complexity automatically. Default is FALSE (recommended to avoid overfitting in short time series).
#' @param formula_model A formula object defining the relationship to fit (e.g., `log.seed ~ TMEAN`).
#' @param model_type Character. Either `'lm'` or `'betareg'` for model fitting type. Default is `'lm'`.
#' @param option.check.name Logical. Whether to validate that `siteneame.forsub` matches site identifiers in `bio_data`. Default is TRUE.
#' @param yearneed Integer. Minimum number of years of historical data needed to calculate lagged climate windows. Default is 2.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Fits smoothed GAMs to the slopes and $R^2$ values from the CSP results.
#'   \item Detects sequences of days with high combined signal strength using thresholds derived from GAM predictions.
#'   \item Extracts and formats these windows, then uses `reruning_windows_modelling()` to fit models over the selected periods.
#'   \item Applies a rolling window to climate data across all years needed to support modeling.
#' }
#'
#' @return A data frame containing, for each selected window:
#' \itemize{
#'   \item Window opening and closing dates (`window.open`, `window.close`)
#'   \item Model performance metrics (e.g., `r2`, `AIC`, `estimate`, `intercept`)
#'   \item Number of observations used, window identifier, and reference day
#' }
#'
#' @note This function assumes daily resolution and that column names (e.g., `log.seed`, `TMEAN`) are consistent across inputs. Validate `bio_data` and `climate_data` for consistency before using.
#'
#' @examples
#' \dontrun{
#' result <- runing_csp(
#'   Results_days = CSP_output[["Site_1"]],
#'   bio_data = bio_data,
#'   climate_data = climate_data,
#'   siteneame.forsub = "Site_1",
#'   refday = 305,
#'   lastdays = 600
#' )
#' }
#'
#' @seealso \code{\link{optimize_and_fit_gam}}, \code{\link{reruning_windows_modelling}}, \code{\link{get_predictions_windows}}, \code{\link{extract_consecutive_sequences}}
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
