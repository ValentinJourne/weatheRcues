#' Calculate Quantiles and Return as Tibble
#' 
#' This function calculates specified quantiles of a numeric vector and returns 
#' the results as a tibble. The default quantiles are the 25th, 50th, and 75th percentiles.
#' 
#' @param x A numeric vector. The function computes quantiles for this vector.
#' @param q A numeric vector of quantiles to compute. Default is c(0.25, 0.5, 0.75). 
#'          Values should be between 0 and 1.
#' 
#' @return A tibble with two columns:
#' \describe{
#'   \item{\code{{{x}}}}{A column containing the quantile values of the input vector `x`.}
#'   \item{\code{{{x}}_q}}{A column containing the quantile values that were computed.}
#' }
#' 
#' @details 
#' The function calculates quantiles specified by the `q` parameter for the numeric vector `x`.
#' The result is a tibble where each row represents a quantile and its corresponding value.
#' 
#' @examples
#' # Example usage of the quibble2 function
#' data_vector <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' result <- quibble2(data_vector, q = c(0.1, 0.5, 0.9))
#' print(result)
#' 
quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}