#' Identify Key Prediction Windows from GAM Outputs
#'
#' This function identifies time windows (days) where both the slope and R-squared predictions from two Generalized Additive Models (GAMs) exceed specified quantile thresholds. It is particularly useful for climate cue window detection where meaningful biological responses are expected during periods of strong and consistent climate signals.
#'
#' @param slope_gam A GAM object (from \code{mgcv::gam}) modeling the slope estimates across days (e.g., daily climate–response effects).
#' @param rs_gam A GAM object modeling R-squared values across days (indicating model fit quality).
#' @param temporary (Currently unused) A placeholder parameter kept for compatibility or future use.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Computes predicted values and standard errors for both slope and R-squared GAMs.
#'   \item Identifies days where slope predictions fall below the 2.5th percentile or above the 97.5th percentile.
#'   \item Identifies days where R-squared predictions exceed the 97.5th percentile.
#'   \item Finds the overlap (intersection) between high-signal slope and R-squared days.
#'   \item Returns the concurrent period as a set of biologically relevant days.
#' }
#'
#' This method is inspired by Simmonds et al. and is used for detecting biologically significant climate signal windows.
#'
#' @return A list with the following elements:
#' \item{pred_S}{Predicted values and SEs for the slope GAM model.}
#' \item{pred_R}{Predicted values and SEs for the R-squared GAM model.}
#' \item{days}{Vector of days where high signal overlap occurs (i.e., meaningful prediction window).}
#'
#' @examples
#' \dontrun{
#' result <- get_predictions_windows(slope_gam = my_slope_gam,
#'                                   rs_gam = my_rs_gam,
#'                                   temporary = NULL)
#' selected_days <- result$days
#' }
#'
#' @seealso \code{\link[mgcv]{gam}}, \code{\link{optimize_and_fit_gam}}, \code{\link{find_concurrent_period}}
#'
#' @export

get_predictions_windows <- function(slope_gam, rs_gam, temporary) {
  pred_S <- predict(slope_gam, se.fit = TRUE, type = "response")
  pred_R <- predict(rs_gam, se.fit = TRUE, type = "response")

  coef_range_l <- which(
    ((pred_S$fit - 1.96 * pred_S$se.fit) - 100) <=
      quantile(pred_S$fit - 100, 0.025)
  )
  coef_range_u <- which(
    ((pred_S$fit + 1.96 * pred_S$se.fit) - 100) >=
      quantile(pred_S$fit - 100, 0.975)
  )
  r_sq_range <- which(
    ((pred_R$fit + 1.96 * pred_R$se.fit) - 100) >=
      quantile(pred_R$fit - 100, 0.975)
  )

  temp_window_l <- coef_range_l[is.element(coef_range_l, r_sq_range)]
  temp_window_u <- coef_range_u[is.element(coef_range_u, r_sq_range)]

  if (length(temp_window_l) == 0) temp_window_l <- coef_range_l
  if (length(temp_window_u) == 0) temp_window_u <- coef_range_u

  days <- find_concurrent_period(
    temp_window = c(temp_window_l, temp_window_u),
    pred_C = pred_S
  )
  list(pred_S = pred_S, pred_R = pred_R, days = days)
}
