#' Get Minimum and Maximum of a Vector
#'
#' This function computes the minimum and maximum values of a numeric vector.
#'
#' @param vector Numeric vector. The input data from which the minimum and maximum values will be computed.
#'
#' @return A numeric vector of length two containing:
#' - The minimum value (first element).
#' - The maximum value (second element).
#'
#' @details This function returns the smallest and largest values from the input numeric vector.
#' If the vector contains missing values (`NA`), they will be ignored in the calculation.
#'
#'
#' @examples
#' vec <- c(1, 2, 3, 4, 5)
#' get.min.max(vec) # returns c(1, 5)
#'
#' vec_with_na <- c(NA, 2, 3, 4, 5)
#' get.min.max(vec_with_na) # returns c(2, 5)
#'
#' @export
get.min.max = function(vector) {
  o.min = min(vector, na.rm = T)
  o.max = max(vector, na.rm = T)
  c(o.min, o.max)
}
