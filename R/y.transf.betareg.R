#' Transform Response Variable for Beta Regression
#'
#' This function transforms a response variable \( y \) for use in beta regression
#' by scaling it to fit the range required by beta regression models. The transformation
#' is performed as follows:
#' 
#' \deqn{y' = \frac{y \cdot (n - 1) + 0.5}{n}} 
#'
#' where \( n \) is the number of non-missing observations in \( y \).
#'
#' @param y A numeric vector of response values that are to be transformed. The vector 
#'          can contain NA values which will be ignored in the calculation of the transformation.
#'
#' @return A numeric vector of the same length as `y`, containing the transformed response values. 
#'         Values are scaled to the range (0, 1), suitable for beta regression.
#'
#' @examples
#' # Example usage
#' y <- c(0.1, 0.5, NA, 0.3, 0.9)
#' transformed_y <- y.transf.betareg(y)
#' print(transformed_y)
#'
#' @export
y.transf.betareg <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}