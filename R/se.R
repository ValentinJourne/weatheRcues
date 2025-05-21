#' Calculate Standard Error
#'
#' This function calculates the standard error of a numeric vector, excluding missing values (`NA`).
#'
#' @param x Numeric vector. The data for which the standard error is to be calculated.
#'
#' @details The function first removes missing values using `na.omit()`. It then calculates the sample size (`n`), and returns the standard deviation of the non-missing values divided by the square root of the sample size (i.e., the standard error).
#'
#' @return A numeric value representing the standard error of the input vector.
#'
#' @importFrom stats sd
#' @export
#'
#' @examples
#' # Example with a numeric vector
#' values <- c(1.5, 2.3, 3.1, NA, 4.0)
#' standard_error <- se(values)
#' print(standard_error)
#'@export
se <- function(x) {
  x2 <- na.omit(x)
  n <- length(x2)
  sd(x2) / sqrt(n)
}
