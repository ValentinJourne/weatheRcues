################################################################################################
#util function, either for data management, or other small stuff
#Note that I am maybe not everytime using them :O 
################################################################################################

#' Check Folder Existence and Contents
#'
#' This function checks whether a folder exists at the specified path and optionally verifies that specific files are present in the folder.
#'
#' @param folder_path A string specifying the path to the folder.
#' @param file_pattern An optional regular expression pattern to match file names in the folder (e.g., "*.csv" for all CSV files).
#' 
#' @return A list with two elements:
#' \item{folder_exists}{A logical indicating whether the folder exists.}
#' \item{files_present}{A character vector of files present in the folder (filtered by `required_files` or `file_pattern`, if provided).}
#' @examples
#' \dontrun{
#' # Check if folder exists and if "file1.csv" and "file2.csv" are in the folder
#' check_folder_and_contents("/path/to/folder", required_files = c("file1.csv", "file2.csv"))
#' 
#' # Check if folder exists and if any CSV files are in the folder
#' check_folder_and_contents("/path/to/folder", file_pattern = "*.csv")
#' }
#' @export
check_folder_and_contents <- function(folder_path, file_pattern = NULL) {
  # Check if the folder exists
  folder_exists <- dir.exists(folder_path)
  
  if (!folder_exists) {
    return(list(folder_exists = FALSE, files_present = NULL))
  }
  
  # Get the list of files in the folder
  folder_contents <- list.files(path = folder_path, full.names = FALSE)
  # If a pattern is provided, filter the files by pattern
   if (!is.null(file_pattern)) {
    files_present <- folder_contents[grepl(file_pattern, folder_contents)]
  } else {
    # Return all files if no specific requirement
    files_present <- folder_contents
  }
  
  return(list(folder_exists = TRUE, files_present = files_present))
}

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

#UPDATED 
extract_sequences_auto <- function(vec, tolerance) {
  # Handle case where the vector has only one element
  if (length(vec) == 1) {
    return(list(vec))
  }
  
  # Sort the vector to handle sequences in ascending order
  vec <- sort(vec)
  sequences <- list()
  current_seq <- c(vec[1])
  
  # Iterate over the vector to find following sequences
  for (i in 2:length(vec)) {
    if (vec[i] - current_seq[length(current_seq)] <= tolerance) {
      current_seq <- c(current_seq, vec[i])
    } else {
      sequences[[length(sequences) + 1]] <- current_seq
      current_seq <- c(vec[i])
    }
  }
  
  # Add the last sequence
  sequences[[length(sequences) + 1]] <- current_seq
  
  return(sequences)
}
#' Re-run climate windows modeling for a specified sequence
#'
#' This function re-runs climate windows modeling for a given index `z`, extracting 
#' window ranges, filtering rolling temperature data, and merging it with biological 
#' site-level data to fit a linear model (`lm`) using a predefined formula.
#'
#' @param z Integer. The index representing the row number in the `window_ranges_df` data frame, which contains the window ranges for analysis.
#' @param tible.sitelevel Data frame. The biological site-level data, which contains 
#'   seed production and other relevant variables for modeling. Default is `tible.sitelevel`.
#' @param rolling.temperature.data Data frame. The climate data containing 
#'   rolling temperature averages and other climate variables. Default is `rolling.temperature.data`.
#'
#' @return A data frame containing the model results, which include:
#'   - `sitenewname`: The name of the biological site.
#'   - `plotname.lon.lat`: The plot name and its longitude/latitude.
#'   - `reference.day`: The day used as the reference point for windowing.
#'   - `windows.size`: The size of the window in days.
#'   - `window.open`: Start of the window period (days before the reference day).
#'   - `window.close`: End of the window period (days before the reference day).
#'   - `intercept`: The intercept of the fitted linear model.
#'   - `intercept.se`: Standard error of the intercept.
#'   - `estimate`: The coefficient estimate for the temperature variable in the model.
#'   - `estimate.se`: Standard error of the coefficient estimate.
#'   - `pvalue`: P-value for the temperature coefficient.
#'   - `r2`: R-squared value of the model.
#'   - `AIC`: Akaike Information Criterion for the model.
#'   - `nobs`: Number of observations used in the model.
#'   - `nsequence.id`: The index corresponding to the window range used.
#'
#' @details
#' The function:
#' 1. Extracts the `windowsopen` and `windowsclose` values from the `window_ranges_df`
#'    for the specified index `z`.
#' 2. Filters `rolling.temperature.data` to retrieve climate data within the window range.
#' 3. Aggregates temperature data (mean temperature) for each combination of `LONGITUDE`, 
#'    `LATITUDE`, and `year`.
#' 4. Joins the filtered climate data with biological site-level data (`tible.sitelevel`).
#' 5. Fits a linear model (`lm`) using the combined data and the predefined formula `myform.fin`.
#' 6. Returns a data frame with the model's coefficients, standard errors, p-values, R-squared, AIC, 
#'    number of observations, and sequence ID.
#'
#' @examples
#' # Assuming window_ranges_df, rolling.temperature.data, and tible.sitelevel are defined:
#' result <- reruning_windows_modelling(1, tible.sitelevel = tible.sitelevel, rolling.temperature.data = rolling.temperature.data)
#' print(result)
#'
#' @export
reruning_windows_modelling = function(z, tible.sitelevel = tible.sitelevel, window_ranges_df = window_ranges_df, rolling.temperature.data = rolling.temperature.data, myform.fin = formula('log.seed ~ mean.temperature'), refday = 305, rollwin = 1) {
  
  # Extract the window open and close for the current iteration
  windowsopen <- window_ranges_df$windowsopen[z]
  windowsclose <- window_ranges_df$windowsclose[z]
  window_number <- window_ranges_df$windows.sequences.number[z]
  
  # Filter the rolling temperature data according to the current window range
  climate_windows_best <- rolling.temperature.data %>%
    filter(days.reversed <= windowsopen & days.reversed >= windowsclose) %>%
    group_by(LONGITUDE, LATITUDE, year) %>%
    summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(window_number = window_number)  # Add the window number for identification
  
  # Merge temperature data with the site-level data and fit the model
  fit.best.model <- tible.sitelevel %>%
    mutate(year = Year) %>%
    left_join(climate_windows_best, by = "year") %>%
    lm(myform.fin, data = .)
  
  # Extract model coefficients and statistics
  tidy_model <- tidy(fit.best.model)
  glance_model <- glance(fit.best.model)
  
  # Create a data frame for the results
  data.frame(
    sitenewname = unique(tible.sitelevel$sitenewname),
    plotname.lon.lat = unique(tible.sitelevel$plotname.lon.lat),
    reference.day = refday,
    windows.size = rollwin,
    #knots.number = results$k,
    window.open = windowsopen,
    window.close = windowsclose,
    intercept = tidy_model$estimate[1],
    intercept.se = tidy_model$std.error[1],
    estimate = tidy_model$estimate[2],
    estimate.se = tidy_model$std.error[2],
    pvalue = tidy_model$p.value[2],
    r2 = glance_model$r.squared,
    AIC = glance_model$AIC,
    nobs = glance_model$nobs,
    nsequence.id = z
  )
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


#' Perform Moving Window Analysis on Climate Data
#'
#' This function performs a moving window analysis to investigate the relationship between a response variable (e.g., seed count) and a rolling climate variable (e.g., temperature) over time. It calculates correlation coefficients, fits linear models, and extracts relevant statistics.
#'
#' @param data A data frame containing the seed count and other site-level information. It should include columns `Year` (or `year`), `sitenewname`, and `log.seed`.
#' @param rolling.temperature.data A data frame containing rolling climate data, including the rolling average temperature and a `days.reversed` column for the moving window analysis.
#' @param method A string specifying the correlation method to use. Default is 'spearman'. Other options include 'pearson' and 'kendall'.
#' @param covariates.of.interest A string specifying the name of the climate variable to be used as a covariate in the analysis. Default is 'rolling_avg_tmean'.
#' @param myform A formula specifying the model to fit for each moving window. Default is \code{formula('log.seed~rolling_avg_tmean')}.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Merges the site-level data with rolling climate data.
#'   \item Calculates correlation coefficients and p-values for each window.
#'   \item Fits linear models for each moving window and extracts model coefficients and statistics.
#'   \item Computes the standard error of the Spearman correlation coefficients.
#' }
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{sitenewname}: The name of the site.
#'   \item \code{plotname.lon.lat}: The plot's longitude and latitude (if available).
#'   \item \code{days.reversed}: The days reversed for the moving window.
#'   \item \code{term}: The term of the model.
#'   \item \code{estimate}: The estimated coefficient for the term.
#'   \item \code{std.error}: The standard error of the coefficient estimate.
#'   \item \code{p.value}: The p-value for the coefficient estimate.
#'   \item \code{r.squared}: The R-squared value of the model.
#'   \item \code{logLik}: The log-likelihood of the model.
#'   \item \code{correlation}: The correlation coefficient.
#'   \item \code{pvalue.cor}: The p-value of the correlation.
#'   \item \code{correlation.se}: The standard error of the correlation coefficient.
#' }
#'
#' @examples
#' # Example data (here only one obs that would not work, but that's the idea)
#' data <- data.frame(
#' Year = as.numeric(2000:2019),
#' sitenewname = rep('site1', 20),
#' log.seed = rnorm(20))
#' rolling.temperature.data <- data.frame(
#' days.reversed = 1:365,
#' rolling_avg_tmean = rnorm(365), year = 2017)
#'
#' # Run the moving window analysis
#' result <- running.movingwin.analysis(data = data, rolling.temperature.data = rolling.temperature.data)
#'
runing.movingwin.analysis = function(data = data,
                                     rolling.temperature.data = rolling.temperature.data,
                                     method = 'spearman',
                                     covariates.of.interest = 'rolling_avg_tmean',
                                     myform = formula('log.seed~rolling_avg_tmean')){
  
  #merge data seed to moving climate
  tible.sitelevel = data %>% #site = bio_data 
    rename(year = Year) %>% 
    left_join(rolling.temperature.data) %>% 
    drop_na(!!sym(covariates.of.interest))
  
  #define correlation - calculate correlation and se, extract also p value 
  n = tible.sitelevel %>% dplyr::select(year, sitenewname) %>% distinct() %>% nrow()
  
  correlation.all <- tible.sitelevel %>% 
    nest(data = -days.reversed) %>%
    mutate(correlation = purrr::map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$estimate)) %>% 
    mutate(pvalue.cor = purrr::map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$p.value))
  
  
  cortemp = correlation.all %>% 
    unnest(c(correlation, pvalue.cor)) %>% 
    dplyr::select(days.reversed, correlation, pvalue.cor) %>% 
    dplyr::mutate(correlation.se = correlation.spearman.se(.$correlation, n)) %>% 
    dplyr::mutate(sitenewname = unique(tible.sitelevel$sitenewname))
  
  #use purr for iteration 
  #here broom package used, because it is simple lm (and not glmmTMB, need to adjust then )
  fitted_models = tible.sitelevel %>%
    nest(data = -days.reversed) %>%
    mutate(model = purrr::map(data, ~lm(myform, data = ., na.action = na.omit)),
           tidied = purrr::map(model, broom::tidy),
           glanced = purrr::map(model, broom::glance),
           augmented = purrr::map(model, broom::augment))
  
  slope = fitted_models %>%
    unnest(tidied) %>% 
    filter(str_detect(term, as.character(myform)[3])) %>% 
    dplyr::select(days.reversed,term, estimate, std.error, p.value) %>% 
    mutate(sitenewname  = unique(tible.sitelevel$sitenewname)) %>% 
    left_join(fitted_models %>%
                unnest(glanced) %>% 
                dplyr::select(days.reversed,r.squared, logLik) %>% 
                mutate(sitenewname  = unique(tible.sitelevel$sitenewname),
                       plotname.lon.lat = unique(tible.sitelevel$plotname.lon.lat))) %>% 
    dplyr::select(sitenewname, plotname.lon.lat, days.reversed, everything()) %>% 
    left_join(cortemp)
  
  return(slope)
}

#' Perform Moving Window Analysis for a Specific Site
#'
#' This function performs a moving window analysis on seed and climate data for a specific site. It involves loading and formatting climate data, applying rolling window functions or not if set rolling at 1 - as done in our study, and then running a moving window analysis to investigate the relationship between seed count and climate variables.
#'
#' @param site.name A string specifying the site name for which to perform the analysis. This should match entries in the `plotname.lon.lat` column of `seed.data`.
#' @param seed.data A data frame containing seed count and site-level information, including `plotname.lon.lat` and `log.seed`.
#' @param climate.path A string specifying the path to the directory containing climate data files. These files should include the site name in their names.
#' @param lastdays An integer specifying the number of days to consider in the rolling window analysis.
#' @param myform A formula specifying the model to fit for each moving window. Default is `formula('log.seed~rolling_avg_tmean')`.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Loads and formats the biological data for the specified site.
#'   \item Loads and formats the climate data, including scaling the climate variables.
#'   \item Defines the period for rolling climate data based on a specified number of years.
#'   \item Applies a rolling window function to the climate data for each year in the defined period.
#'   \item Runs a moving window analysis using the `runing.movingwin.analysis` function.
#' }
#'
#' @return A data frame containing the results of the moving window analysis for the specified site. Includes correlation coefficients, model coefficients, and other relevant statistics.
#'
#' @examples
#' # Example usage:
#' site_name <- 'Site1'
#' Fagus_seed <- data.frame(
#'   plotname.lon.lat = rep(c('Site1', 'Site2'), each = 10),
#'   log.seed = rnorm(20),
#'   year = rep(2000:2009, 2)
#' )
#' climate_beech_path <- 'path/to/climate/data'
#' last_days <- 30
#' my_form <- formula('log.seed~rolling_avg_tmean')
#'
#' result <- site.moving.climate.analysis(
#'   site.name = site_name,
#'   seed.data = Fagus_seed,
#'   climate.path = climate_beech_path,
#'   lastdays = last_days,
#'   myform = my_form
#' )
#'
site.moving.climate.analysis <- function(bio_data, 
                                         climate.data, 
                                         lastdays, 
                                         myform) {
  # Define the year period
  yearneed <- 2
  yearperiod <- (min(climate.data$year) + yearneed):max(climate.data$year)
  
  # Apply the function across all years in yearperiod and combine results
  rolling.temperature.data <- purrr::map_dfr(yearperiod, reformat.climate.backtothepast, 
                                             climate = climate.data, 
                                             yearneed = yearneed, 
                                             refday = 305, 
                                             lastdays = lastdays, 
                                             rollwin = 1, 
                                             variablemoving = 'TMEAN')
  
  results.moving.site = runing.movingwin.analysis(data = bio_data,
                                                  rolling.temperature.data = rolling.temperature.data,
                                                  method = 'spearman',
                                                  covariates.of.interest = 'rolling_avg_tmean',
                                                  myform = myform)
  
  return(results.moving.site)
}

#' Main Function of correlation - regression and process accross all sites
#'
#' This function applies the moving climate analysis to all sites listed in the `seed.data` data frame. It iterates over each site, performs the moving window analysis, and combines the results into a single data frame.
#'
#' @param seed.data A data frame containing seed count and site-level information. It should include a column `plotname.lon.lat` to specify site names and a column `log.seed` for seed counts.
#' @param climate.path A string specifying the path to the directory containing climate data files. These files should be named according to the site names.
#' @param refday An integer representing the reference day for the rolling window analysis. Default is 305.
#' @param lastdays An integer specifying the number of days to include in the rolling window analysis. Default is the maximum of a specified range.
#' @param rollwin An integer specifying the rolling window size. Default is 1.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Identifies unique site names from the `seed.data` data frame.
#'   \item Applies the `site.moving.climate.analysis` function to each site using the unique site names.
#'   \item Combines the results from all sites into a single data frame.
#' }
#'
#' @return A data frame containing the results of the moving window analysis for all sites. Includes correlation coefficients, model coefficients, and other relevant statistics for each site.
#'
#' @examples
#' # Example usage:
#' Fagus_seed <- data.frame(
#'   plotname.lon.lat = rep(c('Site1', 'Site2'), each = 10),
#'   log.seed = rnorm(20),
#'   year = rep(2000:2009, 2)
#' )
#' climate_beech_path <- 'path/to/climate/data'
#' result <- FULL.moving.climate.analysis(
#'   seed.data = Fagus_seed,
#'   climate.path = climate_beech_path,
#'   refday = 305,
#'   lastdays = 30,
#'   rollwin = 1
#' )
#'
FULL.moving.climate.analysis <- function(seed.data.all = seed.data.all,
                                         climate.path = climate.path,
                                         refday = 305,
                                         lastdays = max(range),
                                         rollwin = 1) {
  
  # Get the list of unique site names from seed data
  al.sites <- unique(seed.data.all$sitenewname)
  
  # Iterate over each site and perform the moving climate analysis
  results.moving <- purrr::map_dfr(al.sites, function(site.name) {
    
    # Filter biological data for the current site
    bio_data <- seed.data.all %>%
      filter(sitenewname == site.name) %>%
      as.data.frame() 
    
    # format climate
    climate_data <- format_climate_data(
      site = unique(bio_data$plotname.lon.lat), 
      path = climate.path, 
      scale.climate = TRUE
    )
    
    # run csp site level
    site.moving.climate.analysis(
      bio_data = bio_data, 
      climate.data = climate_data, 
      lastdays = lastdays,
      myform = formula('log.seed ~ rolling_avg_tmean')
    )
  })
  
  return(results.moving)
}

#calcuate standard error 
se <- function(x){
  x2 <- na.omit(x)
  n <- length(x2)
  sd(x2)/sqrt(n)
}

format_climate_data <- function(site,
                                path,
                                scale.climate = TRUE) {
  climate_data <- qs::qread(list.files(path = path, full.names = TRUE, pattern = site)) %>%
    as.data.frame() %>%
    mutate(DATEB = as.Date(DATEB, format = "%m/%d/%y")) %>%
    mutate(date = foo(DATEB, 1949)) %>%  
    mutate(yday = yday(date),
           year = year(date)) %>%
    mutate(TMEAN = as.numeric(TMEAN),
           TMAX = as.numeric(TMAX),
           TMIN = as.numeric(TMIN),
           PRP = as.numeric(PRP))
  
  if (scale.climate) {
    climate_data <- climate_data %>%
      mutate(across(c(TMEAN, TMAX, TMIN, PRP), scale)) %>%
      mutate(across(c(TMEAN, TMAX, TMIN, PRP), as.vector))
  }
 
  return(climate_data)
}

parameter.range = function(raw.data.param.alpha,
                           raw.data.param.beta,
                           raw.data.param.sigma,
                           option = 'mean.sd'){
  
  if(option == 'min.max'){
    raw.data.param.alpha = get.min.max(as_vector(raw.data.param.alpha))
    raw.data.param.beta = get.min.max(as_vector(raw.data.param.beta))
    raw.data.param.sigma = get.min.max(as_vector(raw.data.param.sigma))
  } else if(option == 'mean.sd'){
    raw.data.param.alpha = get.mean.sd(as_vector(raw.data.param.alpha))
    raw.data.param.beta = get.mean.sd(as_vector(raw.data.param.beta))
    raw.data.param.sigma = get.mean.sd(as_vector(raw.data.param.sigma))
  } else {
    stop('wrong option: must be either "mean.sd" or "min.max"')
  }
  
  list(alpha = raw.data.param.alpha, beta = raw.data.param.beta, sigma = raw.data.param.sigma)
  
}

get.min.max = function(vector){
  o.min = min(vector)
  o.max = max(vector)
  c(o.min,o.max)
}
get.mean.sd = function(vector){
  o.mean = mean(vector)
  o.sd = sd(vector)
  c(o.mean,o.sd)
}
