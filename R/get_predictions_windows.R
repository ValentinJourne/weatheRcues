#' Get Window Ranges from Predictions
#'
#' This function calculates prediction windows for given Generalized Additive Models (GAMs) based on slopes and R-squared values. It determines the range of days where the predicted slopes and R-squared values fall within specific quantiles, and then identifies the concurrent period with extreme values.
#'
#' @param slope_gam A GAM model object used for predicting slopes. This model should be fitted using `mgcv::gam` and provide predictions related to the response variable of interest.
#' @param rs_gam A GAM model object used for predicting R-squared values. This model should also be fitted using `mgcv::gam` and provide predictions related to the goodness-of-fit of the model.
#' @param temporary A placeholder parameter that is not used in the current function implementation. It is included for potential future use or consistency with other function signatures.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Predicts values and standard errors for the slope and R-squared GAM models.
#'   \item Identifies the lower and upper range of coefficients where the predictions fall within the 2.5th and 97.5th percentiles of the fitted values.
#'   \item Identifies the range of days where the R-squared values fall within the 97.5th percentile.
#'   \item Finds the intersection of these ranges to determine the concurrent periods with extreme values.
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pred_S}: A list with predicted values and standard errors for the slope GAM model.
#'   \item \code{pred_R}: A list with predicted values and standard errors for the R-squared GAM model.
#'   \item \code{days}: A numeric vector representing the days where the predictions meet the specified criteria.
#' }
#'
#' @examples
#' # Example usage
#' # Assuming slope_gam and rs_gam are previously defined GAM models
#' result <- get_predictions_windows(slope_gam = my_slope_gam, rs_gam = my_rs_gam, temporary = NULL)
#' 
#' # Accessing results
#' pred_S <- result$pred_S
#' pred_R <- result$pred_R
#' days <- result$days
#'
get_predictions_windows <- function(slope_gam, rs_gam, temporary) {
  pred_S <- predict(slope_gam, se.fit = TRUE, type = "response")
  pred_R <- predict(rs_gam, se.fit = TRUE, type = "response")
  
  coef_range_l <- which(((pred_S$fit - 1.96 * pred_S$se.fit) - 100) <= quantile(pred_S$fit - 100, 0.025))
  coef_range_u <- which(((pred_S$fit + 1.96 * pred_S$se.fit) - 100) >= quantile(pred_S$fit - 100, 0.975))
  r_sq_range <- which(((pred_R$fit + 1.96 * pred_R$se.fit) - 100) >= quantile(pred_R$fit - 100, 0.975))
  
  temp_window_l <- coef_range_l[is.element(coef_range_l, r_sq_range)]
  temp_window_u <- coef_range_u[is.element(coef_range_u, r_sq_range)]
  
  if (length(temp_window_l) == 0) temp_window_l <- coef_range_l
  if (length(temp_window_u) == 0) temp_window_u <- coef_range_u
  
  days <- find_concurrent_period(temp_window = c(temp_window_l, temp_window_u), pred_C = pred_S)
  list(pred_S = pred_S, pred_R = pred_R, days = days)
}