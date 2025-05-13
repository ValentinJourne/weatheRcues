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
  print(paste0(
    "refday:",
    refday,
    ' and lastdays:',
    lastdays,
    ' and yearback:',
    yearback
  ))
  if (lastdays > 365 & yearback < 2)
    stop(
      "You need to adjust last day and yearback. If lastdays is higher than 365, you need 2 years "
    )
  if (lastdays > 365 * 2 & yearback < 3)
    stop(
      "You need to adjust last day and yearback. If lastdays is higher than 365*2, you need 3 years "
    )

  #Need for the next function
  monthstart = c('-01-01')
  DATE = seq(
    as.Date(paste0(start_year, monthstart)),
    as.Date(paste0("1953", monthstart)),
    by = "days"
  )
  MONTH = format(as.Date(DATE, format = "%Y-%m-%d"), "%m") %>% as.numeric()
  MONTHab = month.abb[MONTH]
  YEAR = format(as.Date(DATE, format = "%Y-%m-%d"), "%Y") %>% as.numeric()
  DOY = yday(DATE)
  dfata = data.frame(DATE, YEAR, MONTHab, DOY)
  yearperiod = start_year:1953
  sizevec = length(unique(YEAR)) - yearback
  refday = refday
  vectotemp = NULL
  for (k in 1:sizevec) {
    yearsref = yearperiod[k]
    yearrefminusOne <- yearsref - yearback
    tt <- dfata %>%
      dplyr::filter(YEAR <= yearsref & YEAR >= yearrefminusOne) %>%
      dplyr::mutate(
        referenceFin = ifelse(
          YEAR == yearsref & DOY == refday,
          1,
          ifelse(YEAR == yearsref & DOY > refday, NA, 0)
        )
      ) %>%
      dplyr::filter(!is.na(referenceFin)) %>%
      as.data.frame()
    #create sequence going back 365 month before
    seqDays <- seq(1, nrow(tt), 1)
    newsequance <- rep(seqDays)
    ttup <- tt %>%
      dplyr::mutate(days.reversed = rev(newsequance)) %>%
      dplyr::filter(days.reversed < lastdays)
    ttupfin = ttup %>%
      dplyr::arrange(days.reversed) %>%
      dplyr::mutate(YEAR = max(YEAR))
    vectotemp <- rbind(vectotemp, ttupfin)
  }
  return(vectotemp)
}
