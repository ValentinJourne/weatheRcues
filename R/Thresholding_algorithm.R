#' Detect Peaks and Valleys Using a Robust Thresholding Algorithm
#'
#' Implements a real-time peak and valley detection algorithm based on a moving average and standard deviation filter.
#' This method is robust to signal autocorrelation and was adapted from a community solution on StackOverflow:
#' \url{https://stackoverflow.com/questions/72784873/conditional-peak-valley-signal-detection-in-realtime-timeseries-data-r}.
#'
#' @param y A numeric vector representing the time series to analyze.
#' @param lag Integer. The size of the moving window used to compute the rolling mean and standard deviation. Default is 100.
#' @param threshold Numeric. The number of standard deviations a new value must differ from the rolling mean to be classified as a signal. Default is 3.
#' @param influence Numeric (between 0 and 1). Determines how much influence a detected signal has on the recalculation of mean and standard deviation.
#' An influence of 0 means signals are completely excluded from future statistics (more robust); 1 means full influence (less robust). Default is 0.
#'
#' @details
#' The algorithm identifies peaks and valleys by comparing each new observation to a rolling mean and standard deviation.
#' If a point deviates from the mean by more than `threshold` times the rolling standard deviation, it is classified as a signal:
#' \itemize{
#'   \item \code{1} for a peak,
#'   \item \code{-1} for a valley,
#'   \item \code{0} for no signal.
#' }
#'
#' After each signal is detected, the moving mean and standard deviation are updated using the `influence` parameter to moderate the impact of the signal on future values.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{signals}{An integer vector of the same length as \code{y} containing 1 (peak), -1 (valley), or 0 (no signal).}
#'   \item{avgFilter}{A numeric vector with the rolling mean used at each step.}
#'   \item{stdFilter}{A numeric vector with the rolling standard deviation used at each step.}
#' }
#'
#' @examples
#' set.seed(42)
#' data <- cumsum(rnorm(200))  # Simulate time series
#' result <- Thresholding_algorithm(data, lag = 30, threshold = 3, influence = 0.5)
#'
#' plot(data, type = "l", main = "Peak and Valley Detection", col = "gray")
#' points(which(result$signals == 1), data[result$signals == 1], col = "blue", pch = 19)  # Peaks
#' points(which(result$signals == -1), data[result$signals == -1], col = "red", pch = 19)  # Valleys
#'
#' @export

Thresholding_algorithm <- function(y, lag = 100, threshold = 3, influence = 0) {
  # Initialize the signal vector, filtered values, and moving statistics
  signals <- rep(0, length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  # Compute initial average and standard deviation for the lag period
  avgFilter[lag] <- mean(y[0:lag], na.rm = TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm = TRUE)
  # Loop over the time series, starting from lag+1
  for (i in (lag + 1):length(y)) {
    # Check if the current observation is a significant deviation from the moving average
    if (abs(y[i] - avgFilter[i - 1]) > threshold * stdFilter[i - 1]) {
      if (y[i] > avgFilter[i - 1]) {
        signals[i] <- 1 # Peak detected
      } else {
        signals[i] <- -1 # Valley detected
      }
      # Update the filtered signal using the influence factor
      filteredY[i] <- influence * y[i] + (1 - influence) * filteredY[i - 1]
    } else {
      signals[i] <- 0 # No signal detected
      filteredY[i] <- y[i] # No influence, use raw value
    }
    # Update the moving average and standard deviation
    avgFilter[i] <- mean(filteredY[(i - lag):i], na.rm = TRUE)
    stdFilter[i] <- sd(filteredY[(i - lag):i], na.rm = TRUE)
  }
  return(list(
    "signals" = signals,
    "avgFilter" = avgFilter,
    "stdFilter" = stdFilter
  ))
}
