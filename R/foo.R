#' Adjust Date for a Specific Year Range
#'
#' This foo function adjusts the year of a date object to ensure that
#' the date falls within the range from 1950 to 2050. The function updates
#' the year of the input date based on the provided reference year.
#'
#' @param x A date object. The function adjusts the year of this date to
#'           ensure it falls between 1950 and 2050.
#' @param year A numeric value representing the reference year (default is 1968).
#'             The function uses this year to determine whether to assign the
#'             date to the 1900s or 2000s.
#'
#' @return A date object with the adjusted year. The year is set to fall
#'         within the range of 1950 to 2050 based on the provided reference year.
#'
#' @details
#' The function extracts the year component from the date object and compares
#' it with the reference year. If the year in the date object is greater than
#' the reference year, it is assigned to the 1900s; otherwise, it is assigned
#' to the 2000s. This ensures that the date falls within the specified range.
#'
#' @examples
#' # Example usage of the foo function
#' original_date <- as.Date("01/01/80", format="%m/%d/%y")
#' adjusted_date <- foo(original_date, year = 1968)
#' print(adjusted_date)
#'
#' @export
foo <- function(x, year = 1968) {
  m <- year(x) %% 100
  year(x) <- ifelse(m > year %% 100, 1900 + m, 2000 + m)
  x
}
