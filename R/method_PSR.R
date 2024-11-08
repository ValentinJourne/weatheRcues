#method to identify weather cues based on PSR method - from Roberts et al studies 
#' Run PSR Method for Identifying Weather Cues at a Specific Site
#'
#' This function performs the PSR (Penalized Spline Regression) method to identify weather cues for biological data, based on the methodology from Roberts et al., and adapted from Simmonds et al. 
#' The function uses climate data and site-specific biological data to model the effect of temperature (or other climate variables) on a response variable (e.g., seed production).
#'
#' @param bio_data A data frame containing biological data with columns such as `Year`, `plotname.lon.lat`, and `log.seed`.
#' @param site A string representing the site name, which should match the `plotname.lon.lat` in `bio_data`.
#' @param climate.path A string representing the path to the climate data files (e.g., daily temperature data).
#' @param lastdays Integer specifying the total number of days to include in the analysis (default is 600).
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
#'                           climate.path = climate_path, 
#'                           lastdays = 600, 
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
                           climate_csv = climate_csv,
                           lastdays = 600,
                           refday = 305,
                           rollwin = 1,
                           covariates.of.interest = 'TMEAN',
                           myform.fin = formula('log.seed ~ TMEAN'),
                           matrice = c(3,1),
                           knots = NULL,
                           tolerancedays = 7,
                           #plot = TRUE,
                           yearneed = 2){
  
  if (!is.data.frame(bio_data)) {
    stop("bio_data must be a data frame or tibble.")
  }
  if (!is.data.frame(climate_csv)) {
    stop("climate_csv must be a data frame or tibble.")
  }
  if (!is.numeric(lastdays) || length(lastdays) != 1) {
    stop("lastdays must be a numeric value of length 1.")
  }
  if (!is.numeric(refday) || length(refday) != 1) {
    stop("refday must be a numeric value of length 1.")
  }
  if (!is.numeric(rollwin) || length(rollwin) != 1) {
    stop("rollwin must be a numeric value of length 1.")
  }
  if (!is.character(covariates.of.interest) || length(covariates.of.interest) != 1) {
    stop("covariates.of.interest must be a single character string.")
  }
  if (!is.numeric(matrice) || length(matrice) != 2) {
    stop("matrice must be a numeric vector of length 2.")
  }
  if (!is.null(knots) && (!is.numeric(knots) || length(knots) == 0)) {
    stop("knots must be NULL or a numeric vector with at least one element.")
  }
  if (!is.numeric(tolerancedays) || length(tolerancedays) != 1) {
    stop("tolerancedays must be a numeric value of length 1. This arg is about how many days you will tolerate for window identification")
  }
  # if (!is.logical(plot) || length(plot) != 1) {
  #   stop("plot must be a logical value of length 1 (TRUE or FALSE). This will plot the gam predictions")
  # }
  if (!is.numeric(yearneed) || length(yearneed) != 1) {
    stop("yearneed must be a numeric value of length 1.")
  }
  #just checking 
  if((site == unique(bio_data$plotname.lon.lat))==F){
    stop()
  }
  for (col in colnames(bio_data)) {
    if (col == "Year") {
      colnames(bio_data)[colnames(bio_data) == "Year"] <- "year"
    }
  }
  for (col in colnames(climate_csv)) {
    if (col == "Year") {
      colnames(climate_csv)[colnames(climate_csv) == "Year"] <- "year"
    }
  }
  
  
  #now psr method, based from Simmonds et al
  # need climate data to be arranged with year as row
  # need to reduce climate dataframe to only year, yday and temp
  climate2 <- data.frame(year = climate_csv$year, yday = climate_csv$yday, temp = climate_csv[,covariates.of.interest])
  tempmat <- climate2 %>% spread(, key = yday, value = temp)
  tempmat <- tempmat[,-1]
  #number years monitoring seeds 
  ny<-length(bio_data$year)
  nt<-lastdays-1
  ## Formatting data
  index.matrix=matrix(1:nt,ny,nt,byrow=TRUE)
  # Define the year period
  yearneed <- yearneed
  #will fiter the year period here to the year needeed 
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  # Apply the function across all years in yearperiod and combine results
  rolling.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate_csv, 
                                      yearneed = yearneed, 
                                      refday = refday, 
                                      lastdays = lastdays, 
                                      rollwin = rollwin, 
                                      variablemoving = covariates.of.interest)
  #merge data seed to moving climate
  tible.sitelevel = bio_data %>% #site = bio_data 
    #rename(year = Year) %>% 
    left_join(rolling.data) %>% 
    tidyr::drop_na(!!sym(covariates.of.interest))
  
  climate2 <- data.frame(year = tible.sitelevel$year, 
                         yday = tible.sitelevel$days.reversed, 
                         varOfInterst = tible.sitelevel[,covariates.of.interest])
  covariate.matrix = climate2 %>%
    spread(key = yday, value = varOfInterst) %>%
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
  
  #if(plot==TRUE){
    plotted <- plot(model, ylab = c('partial effect'), xlab = c('days prior response'))
  #}
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
    
    #it will use the data from roll data climate 
    output_fit_summary.psr.temp <- map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = bio_data, 
                                                                                                 window_ranges_df = window_ranges_df,
                                                                                                 rolling.data  = rolling.data,
                                                                                                 myform.fin = myform.fin,
                                                                                                 refday = refday,
                                                                                                 model_type = 'lm'))
    
  }
  return(output_fit_summary.psr.temp)
}

#' This function applies the PSR (Penalized Spline Regression) method to identify weather cues for a specific site using biological data. The function is based on methods adapted from Simmonds et al. It filters data for a particular site and calls the `runing_psr_site` function to analyze the climate and biological data.
#'
#' @param site A string representing the site name, which should match the `plotname.lon.lat` in `seed.data`.
#' @param seed.data A data frame containing biological data (e.g., seed production) with columns such as `plotname.lon.lat`, `Year`, and `log.seed`.
#' @param lastdays Integer specifying the total number of days to include in the analysis (default is 600).
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
#' seed.data <- your_fagus_data
#' result <- PSR_function_site(site = site,
#'                             seed.data = seed.data,
#'                             lastdays = 600,
#'                             lastdays = 600,
#'                             matrice = c(3,1),
#'                             knots = NULL)
#' }
#'
#' @export
PSR_function_site <- function(site, 
                              seed.data, 
                              lastdays,
                              refday,
                              climate.path, 
                              matrice = c(3,1),
                              knots = NULL,
                              tolerancedays = 7,
                              yearneed = 2) {
  
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- seed.data %>%
    dplyr::filter(plotname.lon.lat == site)
  
  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub.fagus$plotname.lon.lat), 
    path = climate.path, 
    scale.climate = TRUE
  )
  
  # Run the CSP site analysis function
  runing_psr_site(bio_data = data.sub.fagus,
                  site = site,
                  climate_csv = climate_data,
                  lastdays = 600,
                  refday = refday,
                  rollwin = 1,
                  covariates.of.interest = 'TMEAN',
                  matrice = matrice,
                  knots = knots,
                  yearneed = yearneed,
                  tolerancedays = tolerancedays,
                  plot = TRUE)
}
