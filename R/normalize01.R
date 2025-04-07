#' Normalize Values to a Range of 0 to 1
#' 
#' This function normalizes a numeric vector so that its values are scaled to 
#' fall within the range [0, 1]. The normalization is performed by subtracting 
#' the minimum value of the vector from each element and then dividing by the 
#' range of the vector.
#' 
#' @param x A numeric vector. The function normalizes this vector to a range of [0, 1].
#' 
#' @return A numeric vector with the same length as the input vector `x`, where the values 
#'         are scaled to fall within the range [0, 1].
#' 
#' @details 
#' The function performs min-max normalization. The minimum value in the input vector 
#' is scaled to 0, and the maximum value is scaled to 1. All other values are scaled 
#' proportionally between 0 and 1.
#' 
#' @examples
#' # Example usage of the normalize01 function
#' original_vector <- c(10, 20, 30, 40, 50)
#' normalized_vector <- normalize01(original_vector)
#' print(normalized_vector)
#' 
normalize01 <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}