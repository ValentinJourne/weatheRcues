#' Extract Sequences of Numbers with Tolerated Gaps
#'
#' This function identifies sequences of nearly consecutive values from a numeric vector, allowing for small gaps between values (defined by a tolerance). It is particularly useful for detecting temporally clustered events (e.g., climate signal windows) with some missing or noisy data.
#'#UPDATED in May 2025
#'
#' @param vec Numeric vector. A vector of values (e.g., days or indices). The vector does not need to be pre-sorted.
#' @param tolerance Integer. The maximum allowed gap between values for them to be considered part of the same sequence. For instance, a `tolerance = 2` means values up to 2 units apart are grouped.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Sorts the input vector.
#'   \item Iterates through values, grouping them into sub-vectors if they are within `tolerance` of each other.
#'   \item Returns a list of these sub-vectors (i.e., sequences).
#' }
#' Useful in time-series applications for identifying signal windows or event clusters with small interruptions.
#'
#' @return A list of numeric vectors. Each vector represents a detected sequence where gaps between elements do not exceed the specified tolerance.
#'
#' @examples
#' vec <- c(1, 2, 3, 6, 7, 10, 11, 12)
#' extract_sequences_auto(vec, tolerance = 1)
#' # Returns: list(c(1, 2, 3), c(6, 7), c(10, 11, 12))
#'
#' vec2 <- c(3, 5, 6, 9)
#' extract_sequences_auto(vec2, tolerance = 2)
#' # Returns: list(c(3, 5, 6), c(9))
#'
#' @export
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
