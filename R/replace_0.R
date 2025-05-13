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
    if (
      zero_sequences$values[i] == 0 && zero_sequences$lengths[i] <= threshold
    ) {
      # Determine the value to replace 0s with
      replace_value <- NA
      if (i > 1 && zero_sequences$values[i - 1] != 0) {
        replace_value <- zero_sequences$values[i - 1]
      } else if (
        i < length(zero_sequences$values) && zero_sequences$values[i + 1] != 0
      ) {
        replace_value <- zero_sequences$values[i + 1]
      }

      # Replace the 0s with the identified value
      if (!is.na(replace_value)) {
        result[
          position:(position + zero_sequences$lengths[i] - 1)
        ] <- replace_value
      }
    }

    # Update the position to the start of the next run
    position <- position + zero_sequences$lengths[i]
  }

  return(result)
}
