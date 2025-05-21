#' Optimize and Fit Generalized Additive Models (GAMs)
#'
#' This function fits two Generalized Additive Models (GAMs) to describe the relationship between a time variable (`day`)
#' and two response variables (`slope` and `r.squared`) typically derived from a Climate Sensitivity Profile (CSP).
#' It includes the option to automatically optimize the number of knots (`k`) for the smooth term using `mgcv::k.check`.
#'
#' @param temporary A data frame or tibble containing the input data. Must contain at least the columns: `day`, `slope`, and `r_s`.
#' @param optim.k Logical. If `TRUE`, the function attempts to optimize the number of knots `k` by checking significance
#' of smooth terms using `mgcv::k.check`. If `FALSE`, a fixed value for `k` is used. Default is `TRUE`.
#' @param plots Logical. If `TRUE`, plots of the fitted GAM models with partial effects and residuals are displayed using `gratia`. Default is `FALSE`.
#' @param k Integer. The number of knots to use for the GAM smooth term if `optim.k = FALSE`. Ignored if `optim.k = TRUE`. Default is `20`.
#'
#' @details
#' If `optim.k = TRUE`, the function iteratively fits GAMs across a range of `k` values (from 10 to 365)
#' and performs `k.check` until it finds a value of `k` with all significant smooth term p-values (p < 0.05).
#' If no suitable `k` is found, it defaults to `k = -1`.
#'
#' Two separate models are fit:
#' \itemize{
#'   \item \code{slope_gam}: GAM model of `slope ~ s(day)`
#'   \item \code{rs_gam}: GAM model of `r_s ~ s(day)`
#' }
#'
#' If `plots = TRUE`, both models are visualized using `gratia::draw()` and displayed via `cowplot::plot_grid()`.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{slope_gam}: Fitted GAM for the slope variable.
#'   \item \code{rs_gam}: Fitted GAM for the R² variable.
#'   \item \code{k}: The number of knots used in the GAM fitting.
#' }
#'
#' @examples
#' set.seed(42)
#' temp_df <- data.frame(
#'   day = 1:365,
#'   slope = rnorm(365, mean = 0, sd = 0.5),
#'   r_s = runif(365, min = 0, max = 1)
#' )
#'
#' # Automatic knot optimization
#' result <- optimize_and_fit_gam(temp_df, optim.k = TRUE, plots = TRUE)
#'
#' # Use fixed number of knots
#' result_fixed <- optimize_and_fit_gam(temp_df, optim.k = FALSE, k = 30, plots = TRUE)
#'
#' @importFrom mgcv gam k.check
#' @importFrom cowplot plot_grid
#' @importFrom gratia draw
#' @export

optimize_and_fit_gam <- function(temporary, optim.k = TRUE, plots = F, k = 20) {
  if (
    !is.data.frame(temporary) && !is_tibble(temporary) || ncol(temporary) < 2
  ) {
    stop(
      "You need data frame with at least two columns to fit the gam optimize."
    )
  }

  if (optim.k) {
    # Function to check significance in k.check
    is_significant <- function(check) {
      p_values <- check[, 'p-value']
      all(p_values < 0.05) # Returns TRUE if all p-values are significant
    }

    # Range of k values to try
    k_values <- seq(10, 365) # Adjust to one per day
    optimal_k <- NULL

    for (k in 1:length(k_values)) {
      # Fit the model
      slope_gam <- mgcv::gam(
        slope ~ s(day, k = k_values[k], bs = "cr"),
        data = temporary
      )

      # Perform k.check
      check <- k.check(slope_gam)
      pvalue <- check[, 'p-value']
      kfin <- check[, "k'"]

      # Check if all p-values are significant
      if (!is_significant(check)) {
        optimal_k <- k_values[k]
        break
      }
    }
    if (!is.null(optimal_k)) {
      cat("Optimal k found:", optimal_k, "\n")
      slope_gam <- mgcv::gam(
        slope ~ s(day, k = optimal_k, bs = "cr"),
        data = temporary
      )
      rs_gam <- mgcv::gam(
        r_s ~ s(day, k = optimal_k, bs = "cr"),
        data = temporary
      )
      k <- optimal_k
    } else {
      cat(
        "No optimal k found within the specified range. Using default k = -1.\n"
      )
      slope_gam <- mgcv::gam(
        slope ~ s(day, k = -1, bs = "cr"),
        data = temporary
      )
      rs_gam <- mgcv::gam(r_s ~ s(day, k = -1, bs = "cr"), data = temporary)
      k <- -1
    }
  } else {
    slope_gam <- mgcv::gam(slope ~ s(day, k = k, bs = "cr"), data = temporary)
    rs_gam <- mgcv::gam(r_s ~ s(day, k = k, bs = "cr"), data = temporary)
    k <- k
  }

  if (plots) {
    results <- cowplot::plot_grid(
      gratia::draw(slope_gam, residuals = TRUE) +
        ylab('Slope (partial effect)'),
      gratia::draw(rs_gam, residuals = TRUE) + ylab('R2 (partial effect)')
    )
    print(results)
  }

  list(slope_gam = slope_gam, rs_gam = rs_gam, k = k)
}
