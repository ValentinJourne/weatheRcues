#' Cross-Validation Prediction Based on Weather Cue Windows
#'
#' This function performs prediction and model validation for a given weather cue window identified from prior analysis.
#' It calculates performance metrics such as MAE, RMSE, R², and Lin's Concordance Correlation Coefficient by applying the estimated window-specific model to the validation dataset.
#'
#' @param z Integer index specifying which row of `window_ranges_df` to use for extracting the window range and model parameters.
#' @param bio_data A data frame of biological observations (e.g., seed production), including at least `year`, `sitenewname`, and `plotname.lon.lat`.
#' @param window_ranges_df A data frame containing the best window parameters including `window.open`, `window.close`, `intercept.estimate`, and `slope.estimate`.
#' @param rolling.data A data frame of rolling climate data containing `days.reversed`, `year`, and the climate variable to be matched with the biological data.
#' @param refday Integer. Reference day of the year from which windows are computed. Default is 305 (i.e., November 1).
#' @param rollwin Integer specifying the rolling window size used in the climate smoothing. Default is 1.
#' @param formula_model A formula used for the model (e.g., `log.seed ~ temp`). The right-hand side defines the climate covariate to be extracted from `rolling.data`.
#' @param model_type Character string. Currently only `'lm'` is supported (for linear regression).
#'
#' @return A data frame (1 row) containing:
#' \itemize{
#'   \item \code{sitenewname}, \code{plotname.lon.lat}, and \code{reference.day}
#'   \item The \code{window.open} and \code{window.close} used
#'   \item Model performance metrics: \code{mae}, \code{scaled.mae}, \code{r2.validation}, \code{pbiais}, \code{rmse}, \code{scaled.rmse}, \code{rmse.normalized}
#'   \item Lin's Concordance Correlation Coefficient (\code{ccc.rhos.est}, \code{ccc.rhos.lowCI}, \code{ccc.rhos.highCI})
#'   \item \code{method.cues}: The name of the method used to derive the window
#' }
#'
#' @examples
#' \dontrun{
#' result <- cross_validation_outputs_windows_modelling(
#'   z = 1,
#'   bio_data = bio_train,
#'   window_ranges_df = window_params,
#'   rolling.data = rolling_climate,
#'   formula_model = formula('log.seed ~ temp')
#' )
#' }
#'
#' @export
cross_validation_outputs_windows_modelling = function(
  z,
  bio_data = bio_data,
  window_ranges_df = window_ranges_df,
  rolling.data = rolling.data,
  refday = 305,
  rollwin = 1,
  formula_model = formula('log.seed ~ mean.temperature'),
  model_type = 'lm'
) {
  response.var <- all.vars(formula_model)[1]
  covariates.of.interest <- all.vars(formula_model)[2]

  #en gros apres je veux predire avec les valeurs de parametres que jai eu
  window.open <- window_ranges_df$window.open[z]
  window.close <- window_ranges_df$window.close[z]
  window_method <- window_ranges_df$method[z]
  intercept = window_ranges_df$intercept.estimate[z]
  slope = window_ranges_df$slope.estimate[z]

  # Filter the rolling temperature data according to the current window range
  climate_windows_best <- rolling.data %>%
    filter(days.reversed <= window.open & days.reversed >= window.close) %>%
    group_by(LONGITUDE, LATITUDE, year) %>%
    dplyr::summarise(
      covariates.of.interest.aggregate = mean(
        !!sym(covariates.of.interest),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    mutate(window_method = window_method)

  colnames(climate_windows_best)[4] <- covariates.of.interest

  for (col in colnames(bio_data)) {
    if (col == "Year") {
      colnames(bio_data)[colnames(bio_data) == "Year"] <- "year"
    }
  }

  bio.with.window.climate <- bio_data %>%
    left_join(climate_windows_best, by = "year")

  #make predictions
  observed <- bio.with.window.climate[[response.var]]
  predicted <- intercept +
    slope * bio.with.window.climate[[covariates.of.interest]]

  predictions <- intercept +
    slope * bio.with.window.climate[, covariates.of.interest]

  #get MAE and other metrics
  absolute_errors <- abs(observed - predicted)
  mae <- mean(absolute_errors, na.rm = TRUE)
  mae.scaled <- mae / sd(observed, na.rm = TRUE)
  r2.p <- cor(predicted, observed, use = "complete.obs")^2
  pbiais <- Metrics::percent_bias(observed, predicted)
  rmse <- Metrics::rmse(observed, predicted)
  rmse.scaled <- rmse / sd(observed, na.rm = TRUE)
  rmse.normalized <- rmse / mean(observed, na.rm = TRUE)
  ccc <- DescTools::CCC(observed, predicted)

  out = data.frame(
    sitenewname = unique(bio_data$sitenewname),
    plotname.lon.lat = unique(bio_data$plotname.lon.lat),
    reference.day = refday,
    windows.size = rollwin,
    window.open = window.open,
    window.close = window.close,
    mae = mae,
    scaled.mae = mae.scaled,
    r2.validation = r2.p,
    pbiais = pbiais,
    rmse = rmse,
    scaled.rmse = rmse.scaled,
    rmse.normalized = rmse.normalized,
    ccc.rhos.est = unname(ccc$rho.c[1]),
    ccc.rhos.lowCI = unname(ccc$rho.c[2]),
    ccc.rhos.highCI = unname(ccc$rho.c[3]),
    method.cues = window_method
  )
}
