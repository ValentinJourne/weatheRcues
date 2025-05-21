#' Calculate Parameter Ranges
#'
#' This function processes the input data for parameters alpha, beta, and sigma, and returns either the mean and standard deviation or the minimum and maximum values, based on the specified option.
#'
#' @param raw.data.param.alpha Numeric vector or data for the parameter alpha.
#' @param raw.data.param.beta Numeric vector or data for the parameter beta.
#' @param raw.data.param.sigma Numeric vector or data for the parameter sigma.
#' @param option Character. Either `"mean.sd"` to return the mean and standard deviation, or `"min.max"` to return the minimum and maximum values. Default is `"mean.sd"`.
#'
#' @details This function takes in raw data for parameters alpha, beta, and sigma. Depending on the `option` argument, it calculates either:
#' - the minimum and maximum values (`min.max`), or
#' - the mean and standard deviation (`mean.sd`).
#' The function uses helper functions `get.min.max()` and `get.mean.sd()` to perform these calculations.
#' If an invalid option is provided, the function will stop and return an error message.
#'
#' @return A list containing the processed values for:
#' - `alpha`: Processed values for the alpha parameter.
#' - `beta`: Processed values for the beta parameter.
#' - `sigma`: Processed values for the sigma parameter.
#'
#' @examples
#' # Example with mean and standard deviation option
#' alpha <- rnorm(100)
#' beta <- rnorm(100)
#' sigma <- rnorm(100)
#' param_ranges <- parameter.range(alpha, beta, sigma, option = "mean.sd")
#' print(param_ranges)
#'
#' # Example with min and max option
#' param_ranges_min_max <- parameter.range(alpha, beta, sigma, option = "min.max")
#' print(param_ranges_min_max)
#'
#' @export
parameter.range = function(
  raw.data.param.alpha,
  raw.data.param.beta,
  raw.data.param.sigma,
  option = 'mean.sd'
) {
  if (option == 'min.max') {
    raw.data.param.alpha = get.min.max(as_vector(raw.data.param.alpha))
    raw.data.param.beta = get.min.max(as_vector(raw.data.param.beta))
    raw.data.param.sigma = get.min.max(as_vector(raw.data.param.sigma))
  } else if (option == 'mean.sd') {
    raw.data.param.alpha = get.mean.sd(as_vector(raw.data.param.alpha))
    raw.data.param.beta = get.mean.sd(as_vector(raw.data.param.beta))
    raw.data.param.sigma = get.mean.sd(as_vector(raw.data.param.sigma))
  } else {
    stop('wrong option: must be either "mean.sd" or "min.max"')
  }

  list(
    alpha = raw.data.param.alpha,
    beta = raw.data.param.beta,
    sigma = raw.data.param.sigma
  )
}
