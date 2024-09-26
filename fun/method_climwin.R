#' Perform Climate Window Analysis - based on climwin R package from Bailey et al
#'
#' This function performs a climate window analysis using the `climwin` R package to identify
#' the best temperature windows that correlate with seed production. It processes climate and
#' biological data, runs a sliding window analysis, and returns statistics and windows identified
#' NOTE that most parameters are already well described in climwin R package ;O  
#'You also need to specify your null model as reference, here 
#'baseline = lm(log.seed~1, data = data)
#'but double check you have the proper response name variable
#'
#' @param climate.beech.path A character string specifying the path to the directory containing
#' the climate data files. The function searches for files matching the `site.name` pattern in this directory.
#' @param data A data frame containing biological data including seed production and relevant metadata.
#' It should include columns for `Date2` and `sitenewname`.
#' @param site.name A character string specifying the name of the site for which the analysis is conducted.
#' This is used to filter the climate data specific to this site.
#' @param range A numeric vector of length 2 specifying the range of days over which to perform the sliding
#' window analysis. This defines the maximum and minimun range for window lengths in days. Defaults to `c(600, 0)`.
#' @param cinterval A character string specifying the interval at which climate data is aggregated. 
#' Options include 'day', 'week', etc. Defaults to `'day'` same as climwin.
#' @param refday A numeric vector of length 2 specifying the reference day for absolute windows. 
#' Defaults to `c(01, 11)` (e.g., November 1st) - meaning will search before 1st of November.
#' @param optionwindows A character string specifying the type of window to use. Options include 'absolute'
#' for fixed dates and 'relative' for relative to a specific event. Defaults to `'absolute'`. Relative might be useful
#' if you have a specific date time event, which is not the case with our seed production data
#' @param climate_var A character string specifying the name of the climate variable to analyze (e.g., temperature
#' mean, maximum, or minimum). Defaults to `'TMEAN'`. But it should basically match your same column in your climate file
#'
#' @details
#' The function filters and processes climate data, performs a sliding window analysis using the `climwin` package,
#' and extracts model summaries and window. The function assumes the climate data files are in a format compatible with the `qs`
#' package (qs file are just smaller size and can be read by using qs::qread) and that the biological data contains a valid response variable as specified by `formulanull`.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{sitenewname}: Site-specific information.
#'   \item \code{climate.file}: Name of the climate data file used.
#'   \item \code{climwin_output$combos}: Results from the sliding window analysis, the output from climwin::slidingwin.
#'   \item Performance metrics of the best model (`performance::model_performance`).
#'   \item Coefficients of the best model extracted using `broom::tidy`.
#' }
#'
#' @examples
#' # Example usage:
#' # result <- climwin_site_days(
#' #   climate.beech.path = "path/to/climate/data",
#' #   data = biological_data,
#' #   site.name = "site1"
#' # )
#'
#' @export
climwin_site_days <- function(climate.beech.path,
                              data,
                              site.name,
                              range = c(600, 0),
                              cinterval = 'day',
                              refday = c(01, 11),
                              optionwindows = 'absolute',
                              climate_var = 'TMEAN',
                              formulanull = formula(log.seed ~ 1)) {
  
  print(paste("range:", range))
  print(paste("cinterval:", cinterval))
  print(paste("refday:", refday))
  print(paste("variable of interest:", climate_var))
  
  
  # Filter climate data for the site
  climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = site.name)
  if(length(climate_beech_unique) == 0) stop("You missed the file matching site name.")
  
  # Reformat the climate data, here it is a qs file, but you can change if you have different type files
  #climate_csv <- read_csv(climate_beech_unique, col_types = cols(...1 = col_skip())) %>%
  climate_data <- qs::qread(climate_beech_unique) %>%
    as.data.frame() %>%
    mutate(DATEB = as.Date(DATEB, format = "%m/%d/%y")) %>%
    mutate(date = foo(DATEB, 1949)) %>%
    mutate(yday = lubridate::yday(date)) %>%
    mutate(year = as.numeric(str_sub(as.character(date),1, 4)) ) %>%
    mutate(TMEAN = as.numeric(TMEAN),
           TMAX = as.numeric(TMAX),
           TMIN = as.numeric(TMIN),
           PRP = as.numeric(PRP)) %>%
    mutate(across(c(TMEAN, TMAX, TMIN, PRP), scale))%>% 
    mutate(across(c(TMEAN, TMAX, TMIN, PRP), as.vector))
  
  if(nrow(climate_data)==0)stop("Your file seems to be not in the good shape. Please grab a cup of tea or coffee")
  
  # Load the biological data for the site
  data <- data  %>% 
    filter(!is.na(eval(parse(text = strsplit(as.character(formulanull)[3], " ~ ")[[1]][1]))))
  
  # Run the climwin analysis
  climwin_output <- climwin::slidingwin(
    xvar = list(temperature.degree = climate_data[[climate_var]]),
    cdate = climate_data$date,
    bdate = data$Date2,
    baseline = lm(log.seed~1, data = data),#i Needed to specify the formula here, if not it is not working properly
    cinterval = cinterval,
    range = range,
    type = optionwindows,
    refday = refday,
    stat = "mean",
    cmissing = 'method2',
    func = "lin"
  )
  
  # Extract the summary for the best model
  broom_summary <- broom::tidy(climwin_output[[1]]$BestModel) %>%
    filter(term == 'climate') %>%
    select(-term)
  
  # Extract performance statistics and combine with site information
  statistics <- bind_cols(
    sitenewname = unique(data$sitenewname),
    climate.file = site.name,
    climwin_output$combos, # Extract statistics for all variants
    performance::model_performance(climwin_output[[1]]$BestModel) %>% as.data.frame(),
    broom_summary
  )
  
  return(statistics)
}
