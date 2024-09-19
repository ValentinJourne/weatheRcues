#' Thresholding Algorithm for Peak and Valley Detection
#' obtained from Stackoverflow : https://stackoverflow.com/questions/72784873/conditional-peak-valley-signal-detection-in-realtime-timeseries-data-r
#' 
#' This function detects peaks and valleys in a time series using a threshold-based approach. It identifies significant deviations from a moving average using a threshold based on standard deviations and an influence factor to control the response to new signals.
#'
#' @param y A numeric vector representing the time series data.
#' @param lag An integer specifying the number of observations to use for the moving window average and standard deviation. Defaults to 100.
#' @param threshold A numeric value defining the number of standard deviations that a data point must deviate from the moving average to be considered a signal (either peak or valley). Defaults to 3.
#' @param influence A numeric value between 0 and 1 that controls how much influence a new signal has on the filtered time series. If `influence = 1`, the algorithm reacts fully to new signals, while `influence = 0` ignores new signals. Defaults to 0.
#'
#' @details
#' The algorithm works by:
#' \itemize{
#'   \item Calculating a moving average and standard deviation over a specified `lag`.
#'   \item Comparing new data points with the previous moving average.
#'   \item Flagging points that deviate more than `threshold` times the standard deviation from the average as a signal (1 for peak, -1 for valley).
#'   \item Updating the moving average and standard deviation based on the influence of the new signals.
#' }
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item `signals`: A vector where 1 represents a peak, -1 represents a valley, and 0 indicates no signal.
#'   \item `avgFilter`: The filtered moving average at each point.
#'   \item `stdFilter`: The filtered moving standard deviation at each point.
#' }
#'
#' @examples
#' # Simulated example data
#' set.seed(42)
#' data <- cumsum(rnorm(200))  # Random walk data
#' result <- ThresholdingAlgo(y = data, lag = 30, threshold = 3, influence = 0.5)
#'
#' # Plot results
#' plot(data, type = "l", main = "Peak and Valley Detection", col = "black")
#' points(which(result$signals == 1), data[result$signals == 1], col = "blue", pch = 19)  # Peaks
#' points(which(result$signals == -1), data[result$signals == -1], col = "red", pch = 19)  # Valleys
#'
ThresholdingAlgo <- function(y,
                             lag = 100,
                             threshold = 3,
                             influence = 0) {
  # Initialize the signal vector, filtered values, and moving statistics
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  # Compute initial average and standard deviation for the lag period
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  # Loop over the time series, starting from lag+1
  for (i in (lag+1):length(y)){
    # Check if the current observation is a significant deviation from the moving average
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;# Peak detected
      } else {
        signals[i] <- -1;# Valley detected
      }
      # Update the filtered signal using the influence factor
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0 # No signal detected
      filteredY[i] <- y[i] # No influence, use raw value
    }
    # Update the moving average and standard deviation
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

#' Replace Sequences of Zeros in a Vector with Neighboring Values
#'
#' This function replaces sequences of zeros in a vector with the value of the closest non-zero element, provided the sequence length is within a defined threshold.
#'
#' @param vec A numeric vector that may contain sequences of zeros.
#' @param threshold An integer specifying the maximum length of consecutive zeros to replace. If the length of a zero sequence is less than or equal to `threshold`, the zeros are replaced by the nearest non-zero value. Sequences longer than the threshold are not modified.
#'
#' @details
#' The function scans the input vector for sequences of consecutive zeros and replaces them with the nearest preceding or following non-zero value, if the length of the zero sequence is within the specified `threshold`. Sequences longer than `threshold` are left unchanged.
#'
#' @return A modified version of the input vector, where certain sequences of zeros have been replaced by nearby non-zero values.
#'
#' @examples
#' vec <- c(1, 2, 0, 0, 3, 4, 0, 0, 0, 5)
#' replace_0(vec, threshold = 2)
#' # Expected output: c(1, 2, 2, 2, 3, 4, 0, 0, 0, 5)
#'
replace_0 <- function(vec, threshold) {
  # Find the start and end positions of sequences of 0s
  zero_sequences <- rle(vec)
  
  # Initialize a vector to store the results
  result <- vec
  
  # Keep track of the cumulative position in the original vector
  position <- 1
  
  for (i in seq_along(zero_sequences$lengths)) {
    # Check if the current run is 0s and if it's within the threshold
    if (zero_sequences$values[i] == 0 && zero_sequences$lengths[i] <= threshold) {
      # Determine the value to replace 0s with
      replace_value <- NA
      if (i > 1 && zero_sequences$values[i - 1] != 0) {
        replace_value <- zero_sequences$values[i - 1]
      } else if (i < length(zero_sequences$values) && zero_sequences$values[i + 1] != 0) {
        replace_value <- zero_sequences$values[i + 1]
      }
      
      # Replace the 0s with the identified value
      if (!is.na(replace_value)) {
        result[position:(position + zero_sequences$lengths[i] - 1)] <- replace_value
      }
    }
    
    # Update the position to the start of the next run
    position <- position + zero_sequences$lengths[i]
  }
  
  return(result)
}
