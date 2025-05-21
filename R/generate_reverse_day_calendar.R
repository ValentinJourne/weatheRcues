#' Generate Reverse Calendar Mapping for Days Before a Reference Date
#'
#' This function builds a calendar that maps days counted backward from a given reference date
#' (e.g., November 1st = DOY 305) to their corresponding calendar dates, months, years, and
#' day-of-year (DOY).
#'
#' For instance, a `days.reversed` value of 10 indicates the date 10 days before the specified `refday`.
#'
#' @param refday Integer. The reference day of year (DOY) to count backwards from. For example, November 1st is 305. Default is 274.
#' @param lastdays Integer. Number of days to count backward from the reference date. Default is 1095 (approximately 3 years).
#' @param yearback Integer. Number of full years of historical data needed to cover `lastdays`.
#' Must be ≥2 if `lastdays > 365` or ≥3 if `lastdays > 730`. Default is 3.
#' @param start_year Integer. The base year from which the calendar sequence begins. Default is 1946.
#' This is useful only for display or internal indexing.
#'
#' @details
#' The function internally creates a synthetic date range starting from `start_year` and spanning multiple years
#' (1946 to 1953 by default) to simulate a "calendar of the past." It then identifies reference days in that
#' sequence and maps the days preceding each reference day for a given number of years (`yearback`).
#' The result is a tidy data frame that shows the relationship between reverse days and actual calendar dates.
#'
#' Input consistency is checked, and an error is returned if `lastdays` is too long for the specified `yearback`.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{DATE}{The actual calendar date.}
#'   \item{YEAR}{The year of the date. Reset to most recent year used for plotting convenience.}
#'   \item{MONTHab}{Abbreviated month name (e.g., "Jan", "Feb").}
#'   \item{DOY}{Day of the year (1–366).}
#'   \item{days.reversed}{The number of days prior to the reference day (1 = the day before refday).}
#' }
#'
#' @examples
#' # Generate 600 days of lagged calendar data prior to November 1st
#' calendar_data <- generate_reverse_day_calendar(refday = 305, lastdays = 600, yearback = 3)
#' head(calendar_data)
#'
#' @export
generate_reverse_day_calendar <- function(
  refday = 274,
  lastdays = 1095,
  yearback = 3,
  start_year = 1946
) {
  library(dplyr)
  library(lubridate)

  if (lastdays > 365 & yearback < 2) stop("For >365 days, set yearback >= 2")
  if (lastdays > 730 & yearback < 3) stop("For >730 days, set yearback >= 3")

  end_year <- start_year + yearback + 1

  DATE <- seq(
    as.Date(paste0(start_year, "-01-01")),
    as.Date(paste0(end_year, "-01-01")),
    by = "day"
  )
  dfata <- data.frame(
    DATE = DATE,
    YEAR = year(DATE),
    MONTHab = month.abb[month(DATE)],
    DOY = yday(DATE)
  )

  yearperiod <- start_year:(end_year - 1)
  vectotemp <- vector("list", length = length(yearperiod) - yearback)

  for (k in seq_along(vectotemp)) {
    yearsref <- yearperiod[k + yearback - 1]
    yearrefminusOne <- yearsref - yearback

    tt <- dfata %>%
      filter(YEAR >= yearrefminusOne & YEAR <= yearsref) %>%
      mutate(
        referenceFin = case_when(
          YEAR == yearsref & DOY == refday ~ 1,
          YEAR == yearsref & DOY > refday ~ NA_real_,
          TRUE ~ 0
        )
      ) %>%
      filter(!is.na(referenceFin)) %>%
      arrange(desc(DATE)) %>%
      mutate(days.reversed = 1:n())

    vectotemp[[k]] <- tt %>%
      filter(days.reversed < lastdays) %>%
      arrange(days.reversed) %>%
      mutate(YEAR = yearsref)
  }

  bind_rows(vectotemp)
}
