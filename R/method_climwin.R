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
#' @param climate.path A character string specifying the path to the directory containing
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
#' #   climate.path = "path/to/climate/data",
#' #   data = biological_data,
#' #   site.name = "site1"
#' # )
#'
#' @export
climwin_site_days <- function(climate_data,
                              data,
                              site.name, #character will be added in the final form file
                              range = c(600, 0),
                              cinterval = 'day',
                              refday = c(01, 11),
                              optionwindows = 'absolute',
                              climate_var = 'TMEAN',
                              stat.aggregate = 'mean',
                              formulanull = formula(log.seed ~ 1)) {
  
  print(paste("range:", range))
  print(paste("cinterval:", cinterval))
  print(paste("refday:", refday))
  print(paste("stat aggregate:", stat.aggregate))
  print(paste("windows option:", optionwindows))
  print(paste("variable of interest:", climate_var))
  
  
  # Filter climate data for the site
  #if(length(climate_beech_unique) == 0) stop("You missed the file matching site name.")
  if(nrow(climate_data)==0)stop("Your climate file seems to be not in the good shape. Please double check")
  if(nrow(data)==0)stop("Your bio data file seems to be not in the good shape. Please double check")

  # Load the biological data for the site
  data <- data  %>% 
    dplyr::filter(!is.na(eval(parse(text = strsplit(as.character(formulanull)[3], " ~ ")[[1]][1]))))
  
  # Run the climwin analysis
  climwin_output <- climwin::slidingwin(
    xvar = list(temperature.degree = climate_data[[climate_var]]),
    cdate = climate_data$date,
    bdate = data$Date2,
    baseline = lm(log.seed~1, data = data),#i Needed to specify the formula here, if not it is not working properly (with future_map)
    cinterval = cinterval,
    range = range,
    type = optionwindows,
    refday = refday,
    stat = stat.aggregate,
    cmissing = 'method2',
    func = "lin"
  )
  
  # Extract the summary for the best model
  broom_summary_slope <- broom::tidy(climwin_output[[1]]$BestModel) %>%
    dplyr::filter(term == 'climate') %>%
    dplyr::select(-term)%>% 
    rename_with(.cols = everything(), function(x){paste0("slope.", x)})
  broom_summary_intercept <- broom::tidy(climwin_output[[1]]$BestModel) %>%
    dplyr::filter(term != 'climate') %>%
    dplyr::select(-term) %>% 
    dplyr::rename_with(.cols = everything(), function(x){paste0("intercept.", x)})
  sigma.model = sigma(climwin_output[[1]]$BestModel)
  
  
  # Extract performance statistics and combine with site information
  statistics <- dplyr::bind_cols(
    sitenewname = unique(data$sitenewname),
    climate.file = site.name,
    climwin_output$combos, # Extract statistics for all variants
    performance::model_performance(climwin_output[[1]]$BestModel) %>% as.data.frame(),
    broom_summary_slope,
    broom_summary_intercept,
    sigma = sigma.model
  ) %>% 
    dplyr::rename(window.open = WindowOpen,
                  window.close = WindowClose) %>% 
    dplyr::mutate(window.open = window.open, #a checker, mais je crois que climwin start vector at 0, and me at 1, so maybe need to add +1 
           window.close = window.close)
  
  return(statistics)
}
