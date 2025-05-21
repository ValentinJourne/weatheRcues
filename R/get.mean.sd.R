#' Get Mean and Standard Deviation of a Vector
#'
#' This function computes the mean and standard deviation of a numeric vector.
#'
#' @param vector Numeric vector. The input data from which the mean and standard deviation will be computed.
#'
#' @return A numeric vector of length two containing:
#' - The mean value (first element).
#' - The standard deviation (second element).
#'
#' @details This function returns the average (mean) and standard deviation of the input numeric vector.
#' If the vector contains missing values (`NA`), they will be ignored in the calculation.
#'
#' @examples
#' vec <- c(1, 2, 3, 4, 5)
#' get.mean.sd(vec) # returns c(3, 1.581139)
#'
#' vec_with_na <- c(NA, 2, 3, 4, 5)
#' get.mean.sd(vec_with_na) # returns c(3.5, 1.290994)
#'
#'@export
get.mean.sd = function(vector) {
  o.mean = mean(vector, na.rm = T)
  o.sd = sd(vector, na.rm = T)
  c(o.mean, o.sd)
}
