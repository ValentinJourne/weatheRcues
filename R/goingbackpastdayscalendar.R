#' Generate Calendar for "Historical" Days Reference
#' 
#' This function generates a calendar to map days to their corresponding Day of Year (DOY) 
#' for a specified number of days going back into the past from a reference day. It helps in 
#' creating a historical reference for time series analysis by aligning days from past years.
#' For example if a day reversed starting from 1st november is 10, it means it is 10 days before 
#' 1st of November but will provide here the DOY
#' The function is not best, but still working 
#' meaning that you might need to adjust this function, because here I want this mostly for my figure 
#' 
#' @param refday Numeric. The reference day of the year (DOY) from which to start the backward calculation. Default is 274.
#' @param lastdays Numeric. The total number of days to go back into the past from the reference day. Default is 1095 (1095days, >3 years).
#' @param yearback Numeric. The number of years to go back from the reference year. Default is 3 years.
#' 
#' @details
#' The function creates a calendar by generating dates from a base period (1946 to 1953) and maps these dates to their corresponding DOY. It then calculates the days going back from a specified reference day and creates a sequence of dates accordingly. The function ensures that the number of days to go back and the number of years specified are consistent.
#' 
#' The function will stop and return an error if `lastdays` is greater than 365 and `yearback` is less than 2, or if `lastdays` is greater than 365*2 and `yearback` is less than 3.
#' 
#' @return A data frame containing:
#' \item{DATE}{The date of each day in the sequence.}
#' \item{YEAR}{The year of each date.}
#' \item{MONTHab}{The abbreviated month name.}
#' \item{DOY}{The Day of Year (DOY) for each date.}
#' \item{days.reversed}{The number of days going back from the reference day.}
#' 
#' @examples
#' # Example usage of the goingbackpastdayscalendar function
#' calendar_data <- goingbackpastdayscalendar(refday = 274, lastdays = 1095, yearback = 3)
#' head(calendar_data)
#' 
goingbackpastdayscalendar <-function(refday = 274, 
                                     lastdays = 1095, 
                                     yearback = 3){
  print(paste0("refday:", refday, ' and lastdays:', lastdays, ' and yearback:', yearback))
  if(lastdays>365 & yearback <2)stop("You need to adjust last day and yearback. If lastdays is higher than 365, you need 2 years ")
  if(lastdays>365*2 & yearback <3)stop("You need to adjust last day and yearback. If lastdays is higher than 365*2, you need 3 years ")
  
  #Need for the next function
  monthstart = c('-01-01')
  DATE = seq(as.Date(paste0("1946", monthstart)), as.Date(paste0("1953", monthstart)), by="days")
  MONTH =  format(as.Date(DATE, format="%Y-%m-%d"),"%m") %>% as.numeric()
  MONTHab = month.abb[MONTH]
  YEAR =  format(as.Date(DATE, format="%Y-%m-%d"),"%Y") %>% as.numeric()
  DOY = yday(DATE)
  dfata = data.frame(DATE,YEAR,MONTHab, DOY)
  yearperiod = 1946:1953
  sizevec = length(unique(YEAR))-yearback
  refday = refday
  vectotemp = NULL
  for(k in 1:sizevec){
    yearsref = yearperiod[k]
    yearrefminusOne <- yearsref-yearback
    tt <- dfata %>% 
      dplyr::filter(YEAR <= yearsref & YEAR >= yearrefminusOne) %>% 
      dplyr::mutate(referenceFin = ifelse(YEAR == yearsref & DOY == refday, 1,
                                          ifelse(YEAR == yearsref & DOY > refday, NA, 0))) %>% 
      dplyr::filter(!is.na(referenceFin)) %>% 
      as.data.frame()
    #create sequence going back 365 month before 
    seqDays <- seq(1,nrow(tt),1)
    newsequance <- rep(seqDays)
    ttup <- tt %>% 
      dplyr::mutate(days.reversed = rev(newsequance))%>% 
      dplyr::filter(days.reversed< lastdays )
    ttupfin = ttup %>%
      dplyr::arrange(days.reversed)  %>% 
      dplyr::mutate(YEAR = max(YEAR))  
    vectotemp <- rbind(vectotemp, ttupfin) 
  }
  return(vectotemp)
}