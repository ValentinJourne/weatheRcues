#' Apply P-spline Regression (PSR) Across all Time Series (ATS)
#'
#' This function applies the PSR (Penalized Spline Regression) method to identify influential climate windows
#' for a specific site using biological (e.g., seed production) and climate data. It wraps the `runing_psr()` function
#' for streamlined application across multiple sites.
#'
#' @param site Character. The site identifier, typically matching the `plotname.lon.lat` in `bio_data_all`.
#' @param bio_data_all A data frame of biological data for all sites, including columns such as `plotname.lon.lat`, `log.seed`, and `Year`.
#' @param lastdays Integer. Number of days before `refday` to include in the backward time window. Default is 600.
#' @param refday Integer. Day of year (DOY) used as the reference date for analysis (e.g., 305 = November 1).
#' @param climate.path Character. File path to the folder containing site-level climate `.qs` files.
#' @param matrice Numeric vector of length 2. Controls the smoothness and penalty in the PSR model. Passed to `mgcv::s()` as `m`. Default is `c(3, 1)`.
#' @param knots Integer or NULL. Number of knots to use in the spline. If `NULL`, defaults to `n_years - 1`.
#' @param tolerancedays Integer. Gap tolerance (in days) to identify continuous signal windows. Default is 7.
#' @param yearneed Integer. Number of years of prior data required for climate-window modeling. Default is 2.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Filters the biological dataset for the selected site.
#'   \item Loads matching climate data using `format_climate_data()`.
#'   \item Applies the `runing_psr()` function to detect windows of significant climate influence on seed output.
#'   \item Returns modeled summaries (e.g., window bounds, slope, R², AIC).
#' }
#'
#' @return A data frame summarizing model fit for each identified window, or a placeholder with `NA` values if no significant signal is detected.
#'
#' @examples
#' \dontrun{
#' ATS_PSR(
#'   site = "longitude=4.6_latitude=46.5",
#'   bio_data_all = Fagus.seed,
#'   lastdays = 600,
#'   refday = 305,
#'   climate.path = "data/climate/",
#'   matrice = c(3, 1),
#'   knots = NULL
#' )
#' }
#'
#' @seealso \code{\link{runing_psr}}, \code{\link[mgcv]{gam}}, \code{\link{format_climate_data}}
#' @import dplyr purrr
#' @export
ATS_PSR <- function(
  site,
  bio_data_all,
  lastdays,
  refday,
  climate.path,
  matrice = c(3, 1),
  knots = NULL,
  tolerancedays = 7,
  yearneed = 2
) {
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- bio_data_all %>%
    dplyr::filter(plotname.lon.lat == site)

  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub.fagus$plotname.lon.lat),
    path = climate.path,
    scale.climate = TRUE,
    date_column = "DATEB"
  )

  # Run the CSP site analysis function
  runing_psr(
    bio_data = data.sub.fagus,
    site = site,
    climate_data = climate_data,
    lastdays = 600,
    refday = refday,
    rollwin = 1,
    formula_model = formula('log.seed ~ TMEAN'),
    matrice = matrice,
    knots = knots,
    yearneed = yearneed,
    tolerancedays = tolerancedays #,
    #plot = TRUE
  )
}
