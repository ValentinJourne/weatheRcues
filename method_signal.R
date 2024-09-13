#method for signal processing 
#obtained from Stackoverflow : https://stackoverflow.com/questions/72784873/conditional-peak-valley-signal-detection-in-realtime-timeseries-data-r
#the idea is to identify peak and valley, but I do not want to just find the peak, because sometimes peak might be too low to be true signal
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

#function to replace 0 
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
