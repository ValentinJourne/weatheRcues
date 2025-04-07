#' Format Climate Data for Moving Windows analysis
#' 
#' This function processes climate data to create a historical moving window of climate variables based on a reference day and a period going back in time. 
#' It adjusts for leap years, applies a rolling average to the specified climate variable, and prepares the data for analysis over past periods.
#' Here you could adjust the rolling window to more than 1day (for example if you want average temperature over 7 days, or more or less), but in our method comparison, we set the default to 1
#' 
#' @param yearsref Numeric. The reference year from which to start the calculations. Default is 2000.
#' @param climate Data frame. The climate dataset containing variables such as temperature and date information.
#' @param yearneed Numeric. The number of years needed to adjust for fruit maturation time. Default is 2. But the idea is that if you want to look more than 1095 days to the past, you need to adjust here
#' @param refday Numeric. The reference day of the year (DOY) from which the moving window calculation begins. Default is 274.
#' @param lastdays Numeric. The total number of days to go back from the reference day. Default is 1095 (3 years).
#' @param rollwin Numeric. The window size for the rolling average calculation. Default is 1.
#' @param variablemoving Character. The name of the climate variable to apply the rolling average to. Default is 'temperature.degree' because it is my variable temperature name.
#' 
#' @details
#' The function performs the following tasks:
#' \itemize{
#'   \item Checks if the specified variable exists in the climate dataset.
#'   \item Adjusts the reference day if the reference year is a leap year.
#'   \item Filters the climate data to include only the relevant years.
#'   \item Creates a sequence of days going back from the reference day.
#'   \item Applies a rolling average to the specified climate variable.
#'   \item Returns a data frame with the processed climate data, including the rolling average.
#' }
#' 
#' If the specified `variablemoving` is not found in the climate dataset, a warning is issued, and the function returns `NULL`.
#' 
#' @return A data frame containing:
#' \item{LONGITUDE}{Longitude of the location.}
#' \item{LATITUDE}{Latitude of the location.}
#' \item{year}{The year of each record.}
#' \item{date}{The date of each record.}
#' \item{yday}{Day of Year (DOY) for each date.}
#' \item{days.reversed}{Number of days going back from the reference day.}
#' \item{rolling_avg_tmean}{Rolling average of the specified climate variable.}
#' 
#' @examples
#' # Example usage of the reformat.climate.backtothepast function
#' processed_climate_data <- reformat.climate.backtothepast(yearsref = 2000, 
#'                                                          climate = climate_data, 
#'                                                          yearneed = 2, 
#'                                                          refday = 274, 
#'                                                          lastdays = 1095, 
#'                                                          rollwin = 7, 
#'                                                          variablemoving = 'temperature.degree')
#' head(processed_climate_data)
#'  
reformat.climate.backtothepast <- function(yearsref = 2000, 
                                           climate = climate, 
                                           yearneed = 2, 
                                           refday = 274, 
                                           lastdays = 1095, 
                                           rollwin = 1, 
                                           variablemoving = 'temperature.degree',
                                           align.moving = 'right') {
  
  validate_inputs <- function() {
    # Check that yearsref and yearneed are integers
    if (!is.numeric(yearsref) || length(yearsref) != 1) stop("yearsref must be a single numeric value.")
    if (!is.numeric(yearneed) || length(yearneed) != 1) stop("yearneed must be a single numeric value.")
    
    # Check that climate is a data frame with required columns
    if (!is.data.frame(climate)) stop("climate must be a data frame.")
    
    required_columns <- c("year", "yday", "LONGITUDE", "LATITUDE", "date", variablemoving)
    missing_columns <- setdiff(required_columns, names(climate))
    if (length(missing_columns) > 0) {
      stop(paste("The climate data frame is missing required columns:", paste(missing_columns, collapse = ", ")))
    }
    
    # Ensure variablemoving column exists and is numeric
    if (!variablemoving %in% names(climate)) {
      stop(paste("Column", variablemoving, "not found in the climate dataset."))
    } else if (!is.numeric(climate[[variablemoving]])) {
      stop(paste("The column", variablemoving, "must be numeric for rolling mean calculations."))
    }
    
    # Check that refday and lastdays are positive integers
    if (!is.numeric(refday) || length(refday) != 1 || refday <= 0) stop("refday must be a positive numeric value.")
    if (!is.numeric(lastdays) || length(lastdays) != 1 || lastdays <= 0) stop("lastdays must be a positive numeric value.")
    
    # Check that rollwin is a positive integer
    if (!is.numeric(rollwin) || length(rollwin) != 1 || rollwin <= 0) stop("rollwin must be a positive numeric value.")
  }
  
  # Run the validation tests
  validate_inputs()
  
  
  #need to adjust for leap years?
  if(lubridate::leap_year(yearsref)==TRUE){
    refday = refday+1
  }
  
  yearrefminusOne <- yearsref - yearneed
  tt <- climate %>%
    dplyr::filter(year <= yearsref & year >= yearrefminusOne) %>%
    dplyr::mutate(referenceFin = ifelse(year == yearsref & yday == refday, 1, ifelse(year == yearsref & yday > refday, NA, 0))) %>%
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
    dplyr::mutate(rolling_avg_tmean = zoo::rollmeanr(!!sym(variablemoving), k = rollwin, fill = NA, align = align.moving)) %>%
    dplyr::mutate(year = max(year)) %>%
    dplyr::select(LONGITUDE, LATITUDE, year, date, yday, days.reversed, rolling_avg_tmean)
  
  colnames(ttupfin)[7] = paste(variablemoving)
  return(ttupfin)
}