#method to identify weather cues based on PSR method - Roberts et al studies 
#the code is a version adapted from Simmonds et al

#' Run PSR Method for Identifying Weather Cues at a Specific Site
#'
#' This function performs the PSR (Penalized Spline Regression) method to identify weather cues for biological data, based on the methodology from Roberts et al., and adapted from Simmonds et al. 
#' The function uses climate data and site-specific biological data to model the effect of temperature (or other climate variables) on a response variable (e.g., seed production).
#'
#' @param bio_data A data frame containing biological data with columns such as `Year`, `plotname.lon.lat`, and `log.seed`.
#' @param site A string representing the site name, which should match the `plotname.lon.lat` in `bio_data`.
#' @param climate.beech.path A string representing the path to the climate data files (e.g., daily temperature data).
#' @param tot_days Integer specifying the total number of days to include in the analysis (default is 600).
#' @param refday Integer representing the reference day (default is 305). This is the day from which the climate data is tracked backward.
#' @param rollwin Integer representing the size of the rolling window for calculating the rolling averages (default is 1).
#' @param covariates.of.interest A string indicating the climate variable of interest (e.g., 'TMEAN', 'TMAX', 'TMIN', 'PRP') (default is 'TMEAN').
#' @param matrice A numeric vector of length 2 indicating the penalties to apply in the smoothing function of the model (default is `c(3,1)`).
#' @param knots Integer specifying the number of knots for the GAM model. If `NULL` (default), the function will set the knots to `ny - 1`, where `ny` is the number of years of data.
#' @param tolerancedays Integer specifying the tolerance (in days) for identifying consecutive periods (default is 7).
#' @param plot Logical indicating whether to plot the partial effect of the model (default is `TRUE`).
#'
#' @details 
#' The function processes climate data, calculates rolling averages of climate variables, and fits a generalized additive model (GAM) to the biological data. It identifies periods of significant climate effects based on the model's output.
#' 
#' The function checks for significant temperature cues by calculating upper and lower limits (mean ± 1.96 * SD) and identifies time windows where the model exceeds these thresholds.
#'
#' @return A data frame summarizing the significant windows identified by the model. If no significant windows are found, the function returns a data frame with `NA` values.
#' 
#' @examples
#' \dontrun{
#' # Example usage:
#' bio_data <- your_biological_data
#' site <- "site_name"
#' climate_path <- "path_to_climate_data/"
#' result <- runing_psr_site(bio_data = bio_data, 
#'                           site = site, 
#'                           climate.beech.path = climate_path, 
#'                           tot_days = 600, 
#'                           refday = 305, 
#'                           rollwin = 1, 
#'                           covariates.of.interest = 'TMEAN',
#'                           matrice = c(3,1), 
#'                           knots = NULL, 
#'                           tolerancedays = 7, 
#'                           plot = TRUE)
#' }
#' 
#' @seealso \code{\link[mgcv]{gam}}
#' @export
runing_psr_site = function(bio_data = bio_data,
                           site = site,
                           climate.beech.path = here('climate_dailyEOBS/'),
                           tot_days = 600,
                           refday = 305,
                           rollwin = 1,
                           covariates.of.interest = 'TMEAN',
                           matrice = c(3,1),
                           knots = NULL,
                           tolerancedays = 7,
                           plot = TRUE){
  
  #just checking 
  if((site == unique(bio_data$plotname.lon.lat))==F){
    stop()
  }
  #path clim
  climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern =site)
  climate_csv <- qs::qread(climate_beech_unique) %>%
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
  
  #now psr method, based from Simmonds et al
  # need climate data to be arranged with year as row
  # need to reduce climate dataframe to only year, yday and temp
  climate2 <- data.frame(year = climate_csv$year, yday = climate_csv$yday, temp = climate_csv[,covariates.of.interest])
  tempmat <- climate2 %>% spread(, key = yday, value = temp)
  tempmat <- tempmat[,-1]
  #number years monitoring seeds 
  ny<-length(bio_data$Year)
  nt<-tot_days-1
  ## Formatting data
  index.matrix=matrix(1:nt,ny,nt,byrow=TRUE)
  # Define the year period
  yearneed <- 2
  #will fiter the year period here to the year needeed 
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  # Apply the function across all years in yearperiod and combine results
  rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate_csv, 
                                      yearneed = yearneed, 
                                      refday = refday, 
                                      lastdays = tot_days, 
                                      rollwin = rollwin, 
                                      variablemoving = covariates.of.interest)
  #merge data seed to moving climate
  tible.sitelevel = bio_data %>% #site = bio_data 
    rename(year = Year) %>% 
    left_join(rolling.temperature.data) %>% 
    drop_na(!!sym('rolling_avg_tmean'))
  
  climate2 <- data.frame(year = tible.sitelevel$year, 
                         yday = tible.sitelevel$days.reversed, 
                         temp = tible.sitelevel$rolling_avg_tmean)
  covariate.matrix = climate2 %>%
    spread(key = yday, value = temp) %>%
    dplyr::select(-year) %>%                    
    as.matrix()    
  covariate.matrix <- unname(as.matrix(covariate.matrix))
  #second order, but from https://link.springer.com/article/10.1007/s00484-007-0141-4
  #it is possible to adjust between 1-3 (instead of 2) like c(1,1) instead of c(2,1)
  #here I specify to 0 instead of 1 as in SImmonds et al, meaning that I have no penaties on slope
  #if using 1 instead of 0, the model is not able to identify a cues 
  #which make sense when looking https://link.springer.com/article/10.1007/s00484-011-0472-z 
  if(is.null(knots) ){
    K = ny-1
    model<-mgcv::gam(log.seed~s(index.matrix,k=K,m=matrice,bs="ps",by=covariate.matrix), 
                     data = bio_data, method="GCV.Cp")
    summary(model)
  }else{
    model<-gam(log.seed~s(index.matrix,k=K,m=c(2,1),bs="ps",by=covariate.matrix), 
               data = bio_data, method="GCV.Cp")
    summary(model)}
  
  if(plot==TRUE){
    plotted <- plot(model, ylab = c('partial effect'), xlab = c('days prior response'))
  }
  coefs <- data.frame(fit = plotted[[1]]$fit)
  #wihtout rounding now, will provide different output, and adjust by what is in the main study . 
  upper_limit <- mean(coefs$fit) + (1.96 * sd(coefs$fit))
  lower_limit <- mean(coefs$fit) - (1.96 * sd(coefs$fit))
  # Find the markers (days) where the fit exceeds the threshold
  marker <- which(coefs$fit > upper_limit | coefs$fit < lower_limit)
  
  if(length(marker) == 0){
    output_fit_summary.psr.temp = data.frame(sitenewname = unique(bio_data$sitenewname),
                                             plotname.lon.lat = unique(bio_data$plotname.lon.lat),
                                             reference.day = refday,
                                             windows.size = NA, 
                                             window.open = NA, window.close = NA, intercept = NA,
                                             intercept.se = NA, estimate = NA, estimate.se = NA,
                                             pvalue = NA, r2 = NA, AIC = NA, nobs = ny,
                                             nsequence.id = NA)
    print('no windows identified')
  }else{
    window <- round(plotted[[1]]$x[find_concurrent_period(marker, coefs)])
    
    #here does not work because not consecutive 
    #extract_consecutive_sequences(window, keep_all = T) 
    sequences_days = extract_sequences_auto(window, tolerance = tolerancedays) 
    #i guess now will do the same shit as the other methods 
    window_ranges_df <- save_window_ranges(sequences_days) %>% 
      mutate(windows.sequences.number = 1:nrow(.))
    
    output_fit_summary.psr.temp <- map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = bio_data, 
                                                                                                 window_ranges_df = window_ranges_df,
                                                                                                 rolling.temperature.data = rolling.temperature.data,
                                                                                                 myform.fin = formula('log.seed ~ mean.temperature')))
    
  }
  return(output_fit_summary.psr.temp)
}

#' This function applies the PSR (Penalized Spline Regression) method to identify weather cues for a specific site using biological data. The function is based on methods adapted from Simmonds et al. It filters data for a particular site and calls the `runing_psr_site` function to analyze the climate and biological data.
#'
#' @param site A string representing the site name, which should match the `plotname.lon.lat` in `Fagus.seed`.
#' @param Fagus.seed A data frame containing biological data (e.g., seed production) with columns such as `plotname.lon.lat`, `Year`, and `log.seed`.
#' @param tot_days Integer specifying the total number of days to include in the analysis (default is 600).
#' @param lastdays Integer specifying the final day for the analysis (usually the length of the time series).
#' @param matrice A numeric vector of length 2 indicating the penalties to apply in the smoothing function of the model (default is `c(3,1)`).
#' @param knots Integer specifying the number of knots for the GAM model. If `NULL` (default), the function will set the knots to the number of years minus one.
#'
#' @details
#' The function filters biological data for the specified site and applies the PSR method to analyze the effect of weather variables (e.g., temperature) on the biological data (e.g., seed production). It uses a rolling window to calculate average climate conditions over a specified period.
#'
#' @return
#' A data frame summarizing the significant weather cues for the site based on the PSR method. If no significant windows are found, the function returns a data frame with `NA` values.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' site <- "site_name"
#' Fagus.seed <- your_fagus_data
#' result <- PSR_function_site(site = site,
#'                             Fagus.seed = Fagus.seed,
#'                             tot_days = 600,
#'                             lastdays = 600,
#'                             matrice = c(3,1),
#'                             knots = NULL)
#' }
#'
#' @export
PSR_function_site <- function(site, 
                              Fagus.seed, 
                              tot_days,
                              lastdays, 
                              matrice = c(3,1),
                              knots = NULL) {
  
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- Fagus.seed %>%
    dplyr::filter(plotname.lon.lat == site)
  
  # Run the CSP site analysis function
  runing_psr_site(bio_data = data.sub.fagus,
                  site = site,
                  climate.beech.path = here('climate_dailyEOBS/'),
                  tot_days = 600,
                  refday = 305,
                  rollwin = 1,
                  covariates.of.interest = 'TMEAN',
                  matrice = matrice,
                  knots = knots,
                  tolerancedays = 7,
                  plot = TRUE)
}
