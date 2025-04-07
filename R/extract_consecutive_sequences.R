#' Extract Consecutive Sequences from a Vector
#'
#' This function identifies and extracts sequences of consecutive numbers from a given vector. It can either return all sequences of consecutive numbers or only the longest sequence.
#'
#' @param values Numeric vector. A vector of numeric values from which to extract consecutive sequences.
#' @param keep_all Logical. If `TRUE`, the function returns all sequences of consecutive numbers found in the input vector. If `FALSE`, it returns only the longest sequence. Default is `FALSE`.
#'
#' @details
#' The function iterates through the numeric vector to find sequences where each value is consecutive (i.e., each value is exactly one greater than the previous value). It creates a list of these sequences. Depending on the `keep_all` parameter, it either returns all found sequences or just the longest one.
#'
#' @return
#' If `keep_all` is `FALSE`, returns a numeric vector of the longest sequence of consecutive numbers. If `keep_all` is `TRUE`, returns a list where each element is a numeric vector representing a sequence of consecutive numbers.
#'
#' @examples
#' # Example usage of extract_consecutive_sequences
#' values <- c(1, 2, 3, 5, 6, 7, 10, 11, 12)
#' longest_sequence <- extract_consecutive_sequences(values, keep_all = FALSE)
#' all_sequences <- extract_consecutive_sequences(values, keep_all = TRUE)
#' print(longest_sequence) # Prints: [1, 2, 3]
#' print(all_sequences)    # Prints: [[1, 2, 3], [5, 6, 7], [10, 11, 12]]
#'
extract_consecutive_sequences <- function(values, keep_all = FALSE) {
  # Initialize variables
  sequences <- list()
  current_sequence <- integer(0)
  
  # Iterate through the values
  #if does not match afte +1 then it will create new sequence length 
  for (i in seq_along(values)) {
    if (i == 1 || values[i] == values[i - 1] + 1) {
      current_sequence <- c(current_sequence, values[i])
    } else {
      sequences <- c(sequences, list(current_sequence))
      current_sequence <- values[i]  # Start new sequence
    }
  }
  
  # Add the last sequence
  sequences <- c(sequences, list(current_sequence))
  
  if (!keep_all) {
    # Find the longest sequence
    longest_sequence <- sequences[[which.max(lengths(sequences))]]
    return(longest_sequence)
  } else {
    return(sequences)
  }
}