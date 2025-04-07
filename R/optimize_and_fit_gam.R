#' Optimize and Fit Generalized Additive Models (GAMs)
#'
#' This function fits Generalized Additive Models (GAMs) to predict slope and R-squared values based on the provided data. It includes an option to optimize the number of knots used in the smooth term of the GAM. If optimization is enabled, it searches for the optimal number of knots that results in significant p-values for the model's smooth terms.
#'
#' @param temporary A data frame containing the data to fit the GAM models. It must include at least two columns: \code{slope} (the response variable for the slope model) and \code{r_s} (the response variable for the R-squared model), as well as \code{day} (the predictor variable).
#' @param optim.k A logical flag indicating whether to optimize the number of knots (\code{k}) for the smooth term. If \code{TRUE}, the function will search for the optimal number of knots. Default is \code{TRUE}.
#' @param plots A logical flag indicating whether to plot the fitted GAM models. If \code{TRUE}, plots of the fitted models and their residuals will be displayed. Default is \code{FALSE}.
#' @param k An integer specifying the number of knots to use if \code{optim.k} is \code{FALSE}. Default is \code{20}. If \code{optim.k} is \code{TRUE}, this parameter is ignored.
#'
#' @details
#' If \code{optim.k} is \code{TRUE}, the function iterates over a range of possible knot values (from 10 to 365). For each value, it fits a GAM model and checks the significance of the smooth term using \code{k.check}. The first value of \code{k} that results in significant p-values for all tests is chosen as the optimal number of knots. If no optimal value is found, the default value of \code{k = -1} is used.
#' 
#' The function fits two GAM models: one for the \code{slope} and one for \code{r_s}, both with the chosen number of knots. If \code{plots} is \code{TRUE}, it plots the fitted models and their residuals using \code{cowplot}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{slope_gam}: The fitted GAM model for the slope.
#'   \item \code{rs_gam}: The fitted GAM model for R-squared.
#'   \item \code{k}: The number of knots used in the models.
#' }
#'
#' @examples
#' # Example data
#' set.seed(123)
#' temporary <- data.frame(
#'   day = 1:365,
#'   slope = rnorm(365),
#'   r_s = rnorm(365)
#' )
#'
#' # Fit GAM models with automatic optimization of knots
#' result <- optimize_and_fit_gam(temporary, optim.k = TRUE, plots = TRUE)
#' 
#' # Fit GAM models with a fixed number of knots
#' result_fixed_k <- optimize_and_fit_gam(temporary, optim.k = FALSE, k = 30, plots = TRUE)
#'
optimize_and_fit_gam <- function(temporary, optim.k = TRUE, plots = F, k = 20) {
  if (!is.data.frame(temporary) && !is_tibble(temporary) || ncol(temporary) < 2) {
    stop("You need data frame with at least two columns to fit the gam optimize.")
  }
  
  if (optim.k) {
    # Function to check significance in k.check
    is_significant <- function(check) {
      p_values <- check[,'p-value']
      all(p_values < 0.05) # Returns TRUE if all p-values are significant
    }
    
    # Range of k values to try
    k_values <- seq(10, 365)  # Adjust to one per day 
    optimal_k <- NULL
    
    for (k in 1:length(k_values)) {
      # Fit the model
      slope_gam <- mgcv::gam(slope ~ s(day, k = k_values[k], bs = "cr"), data = temporary)
      
      # Perform k.check
      check <- k.check(slope_gam)
      pvalue <- check[,'p-value']
      kfin <- check[,"k'"]
      
      # Check if all p-values are significant
      if (!is_significant(check)) {
        optimal_k <- k_values[k]
        break
      }
    }
    if (!is.null(optimal_k)) {
      cat("Optimal k found:", optimal_k, "\n")
      slope_gam <- mgcv::gam(slope ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      rs_gam <- mgcv::gam(r_s ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      k <- optimal_k
    } else {
      cat("No optimal k found within the specified range. Using default k = -1.\n")
      slope_gam <- mgcv::gam(slope ~ s(day, k = -1, bs = "cr"), data = temporary)
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
      gratia::draw(slope_gam, residuals = TRUE) + ylab('Slope (partial effect)'), 
      gratia::draw(rs_gam, residuals = TRUE) + ylab('R2 (partial effect)')
    )
    print(results)
  }
  
  list(slope_gam = slope_gam, rs_gam = rs_gam, k = k)
  
}