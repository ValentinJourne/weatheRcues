#' Extract Sequences with Tolerance for Gaps
#'
#' This function extracts sequences of consecutive numbers from a numeric vector, allowing for gaps of up to a specified tolerance. It is useful when you want to group consecutive days or events with a small allowed gap.
#'
#' @param vec Numeric vector. A vector of numeric values to extract sequences from. The vector should be sorted if not using the function.
#' @param tolerance Numeric. The maximum allowed gap between consecutive values for them to be considered part of the same sequence. For example, a tolerance of 1 means that a gap of up to 1 day between consecutive values is allowed.
#'
#' @details
#' The function sorts the input vector and iterates through it to identify sequences of numbers where the gap between consecutive numbers does not exceed the specified tolerance. If the gap is greater than the tolerance, the current sequence is saved, and a new sequence is started. This is useful for analyzing time series data where short gaps might be considered part of the same event or sequence.
#'
#' @return
#' A list of numeric vectors, where each vector represents a sequence of consecutive numbers with gaps up to the specified tolerance. If there are no sequences found, an empty list is returned.
#'
#' @examples
#' # Example usage of extract_sequences_auto
#' vec <- c(1, 2, 3, 5, 6, 8, 9, 10, 15)
#' sequences <- extract_sequences_auto(vec, tolerance = 1)
#' print(sequences) # Prints: [[1, 2, 3], [5, 6], [8, 9, 10], [15]]
#'
# extract_sequences_auto <- function(vec, tolerance) {
#   # Sort the vector to handle sequences in ascending order - start - end 
#   vec <- sort(vec)
#   sequences <- list()
#   current_seq <- c(vec[1])
#   # Iterate over the vector to find folowing sequences
#   for (i in 2:length(vec)) {
#     if (vec[i] - current_seq[length(current_seq)] <= tolerance) {
#       current_seq <- c(current_seq, vec[i])
#     } else {
#       sequences[[length(sequences) + 1]] <- current_seq
#       current_seq <- c(vec[i])
#     }
#   }
#   sequences[[length(sequences) + 1]] <- current_seq
#   return(sequences)
# }

#UPDATED 
extract_sequences_auto <- function(vec, tolerance) {
  # Handle case where the vector has only one element
  if (length(vec) == 1) {
    return(list(vec))
  }
  
  # Sort the vector to handle sequences in ascending order
  vec <- sort(vec)
  sequences <- list()
  current_seq <- c(vec[1])
  
  # Iterate over the vector to find following sequences
  for (i in 2:length(vec)) {
    if (vec[i] - current_seq[length(current_seq)] <= tolerance) {
      current_seq <- c(current_seq, vec[i])
    } else {
      sequences[[length(sequences) + 1]] <- current_seq
      current_seq <- c(vec[i])
    }
  }
  
  # Add the last sequence
  sequences[[length(sequences) + 1]] <- current_seq
  
  return(sequences)
}