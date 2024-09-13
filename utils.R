################################################################################################
#util function, either for data management, or other small stuff
#Note that I am maybe not everytime using them :O 
################################################################################################

#' Adjust Date for a Specific Year Range
#' 
#' This foo function adjusts the year of a date object to ensure that 
#' the date falls within the range from 1950 to 2050. The function updates 
#' the year of the input date based on the provided reference year.
#' 
#' @param x A date object. The function adjusts the year of this date to 
#'           ensure it falls between 1950 and 2050.
#' @param year A numeric value representing the reference year (default is 1968). 
#'             The function uses this year to determine whether to assign the 
#'             date to the 1900s or 2000s.
#' 
#' @return A date object with the adjusted year. The year is set to fall 
#'         within the range of 1950 to 2050 based on the provided reference year.
#' 
#' @details 
#' The function extracts the year component from the date object and compares 
#' it with the reference year. If the year in the date object is greater than 
#' the reference year, it is assigned to the 1900s; otherwise, it is assigned 
#' to the 2000s. This ensures that the date falls within the specified range.
#' 
#' @examples
#' # Example usage of the foo function
#' original_date <- as.Date("01/01/80", format="%m/%d/%y")
#' adjusted_date <- foo(original_date, year = 1968)
#' print(adjusted_date)
#' 
foo <- function(x, year=1968){
  m <- year(x) %% 100
  year(x) <- ifelse(m > year %% 100, 1900+m, 2000+m)
  x
}

#' Normalize Values to a Range of 0 to 1
#' 
#' This function normalizes a numeric vector so that its values are scaled to 
#' fall within the range [0, 1]. The normalization is performed by subtracting 
#' the minimum value of the vector from each element and then dividing by the 
#' range of the vector.
#' 
#' @param x A numeric vector. The function normalizes this vector to a range of [0, 1].
#' 
#' @return A numeric vector with the same length as the input vector `x`, where the values 
#'         are scaled to fall within the range [0, 1].
#' 
#' @details 
#' The function performs min-max normalization. The minimum value in the input vector 
#' is scaled to 0, and the maximum value is scaled to 1. All other values are scaled 
#' proportionally between 0 and 1.
#' 
#' @examples
#' # Example usage of the normalize01 function
#' original_vector <- c(10, 20, 30, 40, 50)
#' normalized_vector <- normalize01(original_vector)
#' print(normalized_vector)
#' 
normalize01 <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

#' Calculate Quantiles and Return as Tibble
#' 
#' This function calculates specified quantiles of a numeric vector and returns 
#' the results as a tibble. The default quantiles are the 25th, 50th, and 75th percentiles.
#' 
#' @param x A numeric vector. The function computes quantiles for this vector.
#' @param q A numeric vector of quantiles to compute. Default is c(0.25, 0.5, 0.75). 
#'          Values should be between 0 and 1.
#' 
#' @return A tibble with two columns:
#' \describe{
#'   \item{\code{{{x}}}}{A column containing the quantile values of the input vector `x`.}
#'   \item{\code{{{x}}_q}}{A column containing the quantile values that were computed.}
#' }
#' 
#' @details 
#' The function calculates quantiles specified by the `q` parameter for the numeric vector `x`.
#' The result is a tibble where each row represents a quantile and its corresponding value.
#' 
#' @examples
#' # Example usage of the quibble2 function
#' data_vector <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' result <- quibble2(data_vector, q = c(0.1, 0.5, 0.9))
#' print(result)
#' 
quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}

#' Calculate Standard Error of Spearman Correlation
#' 
#' This function computes the standard error (SE) of Spearman's rank correlation coefficient
#' based on the sample size and the correlation coefficient.
#' 
#' @param cor A numeric value representing the Spearman correlation coefficient.
#' @param n An integer representing the sample size.
#' 
#' @return A numeric value representing the standard error of the Spearman correlation coefficient.
#' 
#' @details 
#' The standard error is calculated using the formula:
#' \deqn{ \text{SE}_{\text{cor}} = \sqrt{\frac{(1 - \text{cor}^2)^2 (1 + \text{cor}^2 / 2)}{n - 3}} }
#' where \code{cor} is the Spearman correlation coefficient and \code{n} is the sample size.
#' This formula is used to estimate the variability of the Spearman correlation coefficient.
#' 
#' @examples
#' # Example usage of the correlation.spearman.se function
#' spearman_cor <- 0.8
#' sample_size <- 50
#' se <- correlation.spearman.se(spearman_cor, sample_size)
#' print(se)
#' 
correlation.spearman.se <- function (cor, n) 
{
  se.cor <- sqrt((1 - cor^2)^2 * (1 + cor^2/2)/(n - 3))
  return(se.cor)
}

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
      filter(YEAR <= yearsref & YEAR >= yearrefminusOne) %>% 
      mutate(referenceFin = ifelse(YEAR == yearsref & DOY == refday, 1,
                                   ifelse(YEAR == yearsref & DOY > refday, NA, 0))) %>% 
      filter(!is.na(referenceFin)) %>% 
      as.data.frame()
    #create sequence going back 365 month before 
    seqDays <- seq(1,nrow(tt),1)
    newsequance <- rep(seqDays)
    ttup <- tt %>% 
      mutate(days.reversed = rev(newsequance))%>% 
      filter(days.reversed< lastdays )
    ttupfin = ttup %>%
      arrange(days.reversed)  %>% 
      mutate(YEAR = max(YEAR))  
    vectotemp <- rbind(vectotemp, ttupfin) 
  }
  return(vectotemp)
}

#' Save Window Ranges from Sequences of Days
#'
#' This function processes a list of sequences of days and saves the ranges for each sequence. For each sequence, it identifies the start and end days of the window and compiles these into a data frame.
#'
#' @param sequences_days List. A list of sequences of days, where each sequence is a vector of day numbers. Each sequence represents a period of interest.
#'
#' @details
#' The function performs the following:
#' \itemize{
#'   \item Checks if there is only one sequence in the list and handles it accordingly.
#'   \item Loops through each sequence of days in the list.
#'   \item For each sequence, it extracts the first and last day to define the range.
#'   \item Converts the results into a data frame for easier handling.
#' }
#'
#' @return A data frame with two columns:
#' \item{windowsclose}{The start day of each sequence.}
#' \item{windowsopen}{The end day of each sequence.}
#'
#' @examples
#' # Example usage of the save_window_ranges function
#' sequences <- list(c(1, 2, 3, 4, 5), c(10, 11, 12, 13))
#' window_ranges <- save_window_ranges(sequences)
#' print(window_ranges)
#'
save_window_ranges <- function(sequences_days) {
  # Initialize a list to store the results
  window_ranges <- list()
  
  # Check the length of sequences_days
  if (length(sequences_days) == 1) {
    # If only one sequence, save it directly
    windowsclose <- sequences_days[[1]][1]
    windowsopen <- tail(sequences_days[[1]], n = 1)
    
    # Store the result in the list
    window_ranges[[1]] <- c(windowsclose, windowsopen)
  } else {
    # Loop through each sequence and save the windows
    for (p in 1:length(sequences_days)) {
      windowsclose <- sequences_days[[p]][1]
      windowsopen <- tail(sequences_days[[p]], n = 1)
      
      # Store the result in the list
      window_ranges[[p]] <- c(windowsclose, windowsopen)
    }
  }
  
  # Convert the list to a data frame for easier handling
  window_ranges_df <- do.call(rbind, lapply(window_ranges, function(x) {
    data.frame(windowsclose = x[1], 
               windowsopen = x[2])
  }))
  
  return(window_ranges_df)
}

#' Extract Consecutive Sequences from a Vector
#'
#' This function identifies and extracts sequences of consecutive numbers from a given vector. It can either return all sequences of consecutive numbers or only the longest sequence.
#'
#' @param values Numeric vector. A vector of numeric values from which to extract consecutive sequences.
#' @param keep_all Logical. If `TRUE`, the function returns all sequences of consecutive numbers found in the input vector. If `FALSE`, it returns only the longest sequence. Default is `FALSE`.
#'
#' @details
#' The function iterates through the numeric vector to find sequences where each value is consecutive (i.e., each value is exactly one greater than the previous value). It creates a list of these sequences. Depending on the `keep_all` parameter, it either returns all found sequences or just the longest one.
#'
#' @return
#' If `keep_all` is `FALSE`, returns a numeric vector of the longest sequence of consecutive numbers. If `keep_all` is `TRUE`, returns a list where each element is a numeric vector representing a sequence of consecutive numbers.
#'
#' @examples
#' # Example usage of extract_consecutive_sequences
#' values <- c(1, 2, 3, 5, 6, 7, 10, 11, 12)
#' longest_sequence <- extract_consecutive_sequences(values, keep_all = FALSE)
#' all_sequences <- extract_consecutive_sequences(values, keep_all = TRUE)
#' print(longest_sequence) # Prints: [1, 2, 3]
#' print(all_sequences)    # Prints: [[1, 2, 3], [5, 6, 7], [10, 11, 12]]
#'
extract_consecutive_sequences <- function(values, keep_all = FALSE) {
  # Initialize variables
  sequences <- list()
  current_sequence <- integer(0)
  
  # Iterate through the values
  #if does not match afte +1 then it will create new sequence length 
  for (i in seq_along(values)) {
    if (i == 1 || values[i] == values[i - 1] + 1) {
      current_sequence <- c(current_sequence, values[i])
    } else {
      sequences <- c(sequences, list(current_sequence))
      current_sequence <- values[i]  # Start new sequence
    }
  }
  
  # Add the last sequence
  sequences <- c(sequences, list(current_sequence))
  
  if (!keep_all) {
    # Find the longest sequence
    longest_sequence <- sequences[[which.max(lengths(sequences))]]
    return(longest_sequence)
  } else {
    return(sequences)
  }
}

#' Extract Sequences with Tolerance for Gaps
#'
#' This function extracts sequences of consecutive numbers from a numeric vector, allowing for gaps of up to a specified tolerance. It is useful when you want to group consecutive days or events with a small allowed gap.
#'
#' @param vec Numeric vector. A vector of numeric values to extract sequences from. The vector should be sorted if not using the function.
#' @param tolerance Numeric. The maximum allowed gap between consecutive values for them to be considered part of the same sequence. For example, a tolerance of 1 means that a gap of up to 1 day between consecutive values is allowed.
#'
#' @details
#' The function sorts the input vector and iterates through it to identify sequences of numbers where the gap between consecutive numbers does not exceed the specified tolerance. If the gap is greater than the tolerance, the current sequence is saved, and a new sequence is started. This is useful for analyzing time series data where short gaps might be considered part of the same event or sequence.
#'
#' @return
#' A list of numeric vectors, where each vector represents a sequence of consecutive numbers with gaps up to the specified tolerance. If there are no sequences found, an empty list is returned.
#'
#' @examples
#' # Example usage of extract_sequences_auto
#' vec <- c(1, 2, 3, 5, 6, 8, 9, 10, 15)
#' sequences <- extract_sequences_auto(vec, tolerance = 1)
#' print(sequences) # Prints: [[1, 2, 3], [5, 6], [8, 9, 10], [15]]
#'
extract_sequences_auto <- function(vec, tolerance) {
  # Sort the vector to handle sequences in ascending order - start - end 
  vec <- sort(vec)
  sequences <- list()
  current_seq <- c(vec[1])
  # Iterate over the vector to find folowing sequences
  for (i in 2:length(vec)) {
    if (vec[i] - current_seq[length(current_seq)] <= tolerance) {
      current_seq <- c(current_seq, vec[i])
    } else {
      sequences[[length(sequences) + 1]] <- current_seq
      current_seq <- c(vec[i])
    }
  }
  sequences[[length(sequences) + 1]] <- current_seq
  return(sequences)
}


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
                                           variablemoving = 'temperature.degree') {
  # Print parameter values - to see if no shit 
  print(paste("rollwin - size window:", rollwin))
  print(paste("refday - reference day:", refday))
  print(paste("lastdays - last day:", lastdays))
  print(paste("variablemoving - variable to move:", variablemoving))
  print(paste("yearneed - adjust to maturation fruit time :", yearneed, 'years'))
  
  if (!variablemoving %in% names(climate)) {
    warning(paste("Warning: Column", variablemoving, "not found in the dataset."))
    return(NULL)
  }
  
  #need to adjust for leap years?
  if(lubridate::leap_year(yearsref)==TRUE){
    refday = refday+1
  }
  
  yearrefminusOne <- yearsref - yearneed
  tt <- climate %>%
    filter(year <= yearsref & year >= yearrefminusOne) %>%
    mutate(referenceFin = ifelse(year == yearsref & yday == refday, 1, ifelse(year == yearsref & yday > refday, NA, 0))) %>%
    filter(!is.na(referenceFin)) %>%
    as.data.frame()
  
  # Create sequence going back lastdays days before the reference day
  seqDays <- seq(1, nrow(tt), 1)
  newsequance <- rep(seqDays)
  
  ttup <- tt %>%
    mutate(days.reversed = rev(newsequance)) %>%
    filter(days.reversed < lastdays)
  
  #use !!sym; convert a string, here my variable name, to a symbol
  ttupfin <- ttup %>%
    arrange(days.reversed) %>%
    mutate(rolling_avg_tmean = zoo::rollmeanr(!!sym(variablemoving), k = rollwin, fill = NA, align = 'right')) %>%
    mutate(year = max(year)) %>%
    dplyr::select(LONGITUDE, LATITUDE, year, date, yday, days.reversed, rolling_avg_tmean)
  
  return(ttupfin)
}


#' Combine Multiple Random Forest Models
#'copy from https://stackoverflow.com/questions/19170130/combining-random-forests-built-with-different-training-sets-in-r

#' This function combines multiple `randomForest` objects into a single `randomForest` object. The combination is done by merging the trees from each individual model, allowing for the aggregation of predictions from multiple random forests built with different training sets.
#'
#' @param ... Multiple `randomForest` objects. These should be objects of class `"randomForest"` created using the `randomForest` package. They will be combined into a single random forest model.
#'
#' @details
#' The function first checks that all input objects are of class `"randomForest"` and that they have the same set of predictor variables. It then combines the trees from each model, adjusting for differences in tree structure if necessary. This is achieved by padding trees and node information to ensure consistency across models.
#' 
#' The function supports both classification and regression random forests. For classification, the `cutoff` values are retained, while for regression, the mean squared error (MSE) and R-squared values are not included in the final combined model.
#'
#' @return
#' A `randomForest` object that combines the trees from all input models. This object has the same structure as a typical `randomForest` object, with merged tree information and adjusted for the combined number of trees.
#'
#' @examples
#' # Example usage of my_combine
#' library(randomForest)
#' rf1 <- randomForest(Species ~ ., data = iris, ntree = 50)
#' rf2 <- randomForest(Species ~ ., data = iris, ntree = 50)
#' combined_rf <- my_combine(rf1, rf2)
#' print(combined_rf)
#'
my_combine <- function (...) 
{
  pad0 <- function(x, len) c(x, rep(0, len - length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow = len - 
                                              nrow(x), ncol = ncol(x)))
  rflist <- list(...)
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) 
    stop("Argument must be a list of randomForest objects")
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
  vlist <- lapply(rflist, function(x) rownames(importance(x)))
  numvars <- sapply(vlist, length)
  if (!all(numvars[1] == numvars[-1])) 
    stop("Unequal number of predictor variables in the randomForest objects.")
  for (i in seq_along(vlist)) {
    if (!all(vlist[[i]] == vlist[[1]])) 
      stop("Predictor variables are different in the randomForest objects.")
  }
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
  if (all(haveForest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodestatus, nrnodes)))
    rf$forest$bestvar <- do.call("cbind", lapply(rflist, 
                                                 function(x) padm0(x$forest$bestvar, nrnodes)))
    rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$xbestsplit, nrnodes)))
    rf$forest$nodepred <- do.call("cbind", lapply(rflist, 
                                                  function(x) padm0(x$forest$nodepred, nrnodes)))
    tree.dim <- dim(rf$forest$treemap)
    if (classRF) {
      rf$forest$treemap <- array(unlist(lapply(rflist, 
                                               function(x) apply(x$forest$treemap, 2:3, pad0, 
                                                                 nrnodes))), c(nrnodes, 2, ntree))
    }
    else {
      rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, 
                                                        function(x) padm0(x$forest$leftDaughter, nrnodes)))
      rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, 
                                                         function(x) padm0(x$forest$rightDaughter, nrnodes)))
    }
    rf$forest$ntree <- ntree
    if (classRF) 
      rf$forest$cutoff <- rflist[[1]]$forest$cutoff
  }
  else {
    rf$forest <- NULL
  }
  #
  #Tons of stuff removed here...
  #
  if (classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if (haveTest) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  }
  else {
    rf$mse <- rf$rsq <- NULL
    if (haveTest) 
      rf$test$mse <- rf$test$rsq <- NULL
  }
  rf
}

