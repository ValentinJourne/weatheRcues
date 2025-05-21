#' Reformat Climate Data for Time-Lagged Window Analysis
#'
#' This function extracts and processes historical daily climate data relative to a reference date.
#' It constructs a "reversed day" timeline for backward-looking analyses (e.g., climate cue detection),
#' applies a rolling average (if specified), and prepares the data for use in time-window modeling.
#'
#' @param yearsref Integer. The reference year from which to anchor the backward window. Typically corresponds to a seed or biological event year.
#' @param climate_data Data frame. Must include `year`, `yday`, `LONGITUDE`, `LATITUDE`, `date`, and the target climate variable.
#' @param yearneed Integer. Number of years of prior data to include. Controls how far back the data will be gathered. Default is 2.
#' @param refday Integer. Day of year (DOY) representing the reference event (e.g., 305 for Nov 1st). Adjusted for leap years if applicable.
#' @param lastdays Integer. Number of days to include in the backward time window. Default is 1095 (≈ 3 years).
#' @param rollwin Integer. Rolling average window size (in days). Default is 1 (no smoothing).
#' @param covariates.of.interest Character. Name of the climate variable column to apply the rolling average to (e.g., `"TMEAN"`).
#' @param align.moving Character. Alignment of the rolling mean (`"right"`, `"left"`, or `"center"`). Passed to \code{zoo::rollmeanr}. Default is `"right"`.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Validates inputs and ensures required columns are present.
#'   \item Handles leap years by adjusting `refday` if needed.
#'   \item Filters the climate data for the target window (`yearref - yearneed` to `yearref`).
#'   \item Reverses the daily time vector for backward analysis (e.g., day -1, -2, ..., -1095).
#'   \item Computes a rolling average on the climate variable (optional).
#'   \item Outputs a clean data frame suitable for merging with biological observations.
#' }
#'
#' @return A data frame with columns:
#' \item{LONGITUDE, LATITUDE, year, date, yday}{Metadata and temporal structure.}
#' \item{days.reversed}{The backward time index from the reference day (1 = most recent).}
#' \item{<climate_variable>}{Smoothed (or raw) values of the climate variable specified.}
#'
#' @examples
#' \dontrun{
#' processed <- reformat_climate_backtothepast(
#'   yearsref = 2020,
#'   climate_data = my_climate_df,
#'   yearneed = 2,
#'   refday = 305,
#'   lastdays = 700,
#'   rollwin = 7,
#'   covariates.of.interest = "TMEAN"
#' )
#' head(processed)
#' }
#'
#' @seealso \code{\link[zoo]{rollmean}}, \code{\link[lubridate]{leap_year}}, \code{\link[dplyr]{mutate}}
#' @export

reformat_climate_backtothepast <- function(
  yearsref = 2000,
  climate_data = climate_data,
  yearneed = 2,
  refday = 274,
  lastdays = 1095,
  rollwin = 1,
  covariates.of.interest = 'temperature.degree',
  align.moving = 'right'
) {
  validate_inputs <- function() {
    # Check that yearsref and yearneed are integers
    if (!is.numeric(yearsref) || length(yearsref) != 1)
      stop("yearsref must be a single numeric value.")
    if (!is.numeric(yearneed) || length(yearneed) != 1)
      stop("yearneed must be a single numeric value.")

    # Check that climate_data is a data frame with required columns
    if (!is.data.frame(climate_data)) stop("climate_data must be a data frame.")

    required_columns <- c(
      "year",
      "yday",
      "LONGITUDE",
      "LATITUDE",
      "date",
      covariates.of.interest
    )
    missing_columns <- setdiff(required_columns, names(climate_data))
    if (length(missing_columns) > 0) {
      stop(paste(
        "The climate_data data frame is missing required columns:",
        paste(missing_columns, collapse = ", ")
      ))
    }

    # Ensure covariates.of.interest column exists and is numeric
    if (!covariates.of.interest %in% names(climate_data)) {
      stop(paste(
        "Column",
        covariates.of.interest,
        "not found in the climate_data dataset."
      ))
    } else if (!is.numeric(climate_data[[covariates.of.interest]])) {
      stop(paste(
        "The column",
        covariates.of.interest,
        "must be numeric for rolling mean calculations."
      ))
    }

    # Check that refday and lastdays are positive integers
    if (!is.numeric(refday) || length(refday) != 1 || refday <= 0)
      stop("refday must be a positive numeric value.")
    if (!is.numeric(lastdays) || length(lastdays) != 1 || lastdays <= 0)
      stop("lastdays must be a positive numeric value.")

    # Check that rollwin is a positive integer
    if (!is.numeric(rollwin) || length(rollwin) != 1 || rollwin <= 0)
      stop("rollwin must be a positive numeric value.")
  }

  # Run the validation tests
  validate_inputs()

  #need to adjust for leap years?
  if (lubridate::leap_year(yearsref) == TRUE) {
    refday = refday + 1
  }

  yearrefminusOne <- yearsref - yearneed
  tt <- climate_data %>%
    dplyr::filter(year <= yearsref & year >= yearrefminusOne) %>%
    dplyr::mutate(
      referenceFin = ifelse(
        year == yearsref & yday == refday,
        1,
        ifelse(year == yearsref & yday > refday, NA, 0)
      )
    ) %>%
    dplyr::filter(!is.na(referenceFin)) %>%
    as.data.frame()

  # Create sequence going back lastdays days before the reference day
  seqDays <- seq(1, nrow(tt), 1)
  newsequance <- rep(seqDays)

  ttup <- tt %>%
    dplyr::mutate(days.reversed = rev(newsequance)) %>%
    dplyr::filter(days.reversed < lastdays)

  #use !!sym; convert a string, here my variable name, to a symbol
  ttupfin <- ttup %>%
    dplyr::arrange(days.reversed) %>%
    dplyr::mutate(
      rolling_avg_clim = zoo::rollmeanr(
        !!sym(covariates.of.interest),
        k = rollwin,
        fill = NA,
        align = align.moving
      )
    ) %>%
    dplyr::mutate(year = max(year)) %>%
    dplyr::select(
      LONGITUDE,
      LATITUDE,
      year,
      date,
      yday,
      days.reversed,
      rolling_avg_clim
    )

  colnames(ttupfin)[7] = paste(covariates.of.interest)
  return(ttupfin)
}
