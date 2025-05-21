#' Save Window Ranges from Sequences of Days
#'
#' This function processes a list of sequences of days and saves the ranges for each sequence. For each sequence, it identifies the start and end days of the window and compiles these into a data frame.
#'
#' @param sequences_days List. A list of sequences of days, where each sequence is a vector of day numbers. Each sequence represents a period of interest.
#'
#' @details
#' The function performs the following:
#' \itemize{
#'   \item Checks if there is only one sequence in the list and handles it accordingly.
#'   \item Loops through each sequence of days in the list.
#'   \item For each sequence, it extracts the first and last day to define the range.
#'   \item Converts the results into a data frame for easier handling.
#' }
#'
#' @return A data frame with two columns:
#' \item{window.close}{The start day of each sequence.}
#' \item{window.open}{The end day of each sequence.}
#'
#' @examples
#' # Example usage of the save_window_ranges function
#' sequences <- list(c(1, 2, 3, 4, 5), c(10, 11, 12, 13))
#' window_ranges <- save_window_ranges(sequences)
#' print(window_ranges)
#'@export
save_window_ranges <- function(sequences_days) {
  window_ranges <- list()

  # Check the length of sequences_days
  if (length(sequences_days) == 1) {
    # If only one sequence, save it directly
    window.close <- sequences_days[[1]][1]
    window.open <- tail(sequences_days[[1]], n = 1)

    # Store the result in the list
    window_ranges[[1]] <- c(window.close, window.open)
  } else {
    # Loop through each sequence and save the windows
    for (p in 1:length(sequences_days)) {
      window.close <- sequences_days[[p]][1]
      window.open <- tail(sequences_days[[p]], n = 1) #should come from utils base R

      # Store the result in the list
      window_ranges[[p]] <- c(window.close, window.open)
    }
  }

  # Convert the list to a data frame for easier handling
  window_ranges_df <- do.call(
    rbind,
    lapply(window_ranges, function(x) {
      data.frame(window.close = x[1], window.open = x[2])
    })
  )

  return(window_ranges_df)
}
