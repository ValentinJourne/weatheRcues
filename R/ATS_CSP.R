#' Apply Climate Sensitivity Profile (CSP) Analysis Across all Time Series (ATS)
#'
#' This function performs Climate Sensitivity Profile (CSP) analysis for a specific site across its entire available time series (ATS).
#' It filters biological seed data by site, extracts the corresponding climate data, and applies the `runing_csp()` function
#' to identify influential climate windows using generalized additive models (GAMs).
#'
#' @param Results_days A data frame containing slope and R² values from a previous daily relationship analysis (e.g., output of `runing_daily_relationship()`).
#' @param siteneame.forsub Character string. Site identifier (e.g., `"longitude=5.12_latitude=45.21"` or `plotname.lon.lat`) used to filter biological and climate data.
#' @param bio_data_all A data frame with biological data from multiple sites, including `sitenewname`, `plotname.lon.lat`, and `log.seed`.
#' @param climate.path Character string. File path to the folder containing climate data files for each site.
#' @param refday Integer. Reference day of the year from which backwards climate windows are calculated (e.g., 305 = Nov 1).
#' @param lastdays Integer. Number of days to include before the reference day (i.e., total window length).
#' @param rollwin Integer. Size of the rolling window used to smooth climate data. Default is 1 (no smoothing).
#' @param optim.k Logical. If `TRUE`, will optimize the number of knots in the GAM used during CSP fitting. Default is `FALSE`.
#'
#' @details
#' This function simplifies the pipeline for applying CSP across full time series. It:
#' \itemize{
#'   \item Filters the biological dataset to extract site-specific seed data.
#'   \item Loads and scales the corresponding climate data using `format_climate_data()`.
#'   \item Runs `runing_csp()` to detect significant weather cue windows based on GAM-predicted slope and R² profiles.
#' }
#'
#' @return A data frame with model fit summaries for each identified cue window at the specified site. Includes window bounds, estimates, and performance metrics.
#'
#' @seealso \code{\link{runing_csp}}, \code{\link{runing_daily_relationship}}, \code{\link{format_climate_data}}
#'
#' @examples
#' \dontrun{
#' result <- ATS_CSP(
#'   Results_days = Results_CSP[[1]],
#'   siteneame.forsub = "site_123",
#'   bio_data_all = bio_data,
#'   climate.path = "data/climate_beech",
#'   refday = 305,
#'   lastdays = 600,
#'   rollwin = 1
#' )
#' }
#'
#' @export

ATS_CSP <- function(
  Results_days,
  siteneame.forsub,
  bio_data_all,
  climate.path,
  refday,
  lastdays,
  rollwin,
  optim.k = F
) {
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- bio_data_all %>%
    dplyr::filter(
      sitenewname == siteneame.forsub | plotname.lon.lat == siteneame.forsub
    )

  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub.fagus$plotname.lon.lat),
    path = climate.path,
    scale.climate = TRUE
  )

  # Run the CSP site analysis function
  runing_csp(
    Results_days = Results_days,
    bio_data = data.sub.fagus,
    siteneame.forsub = siteneame.forsub,
    climate_data = climate_data,
    refday = refday,
    optim.k = optim.k,
    lastdays = lastdays, # Change from max(range) to lastdays
    rollwin = rollwin,
    formula_model = formula('log.seed~TMEAN'),
    model_type = 'lm'
  )
}
