#' Perform Climate Window Analysis - based on Climate sensitivity profil from Thackeray et al
#' note that some function have been obtained or adjusted from this paper OR also from Simmonds et al
#' 
#' 
#' 
#' Find Concurrent Periods with Extreme Averages
#'
#' This function identifies periods within a given temperature window where the average of predicted coefficients is most extreme. The periods are determined based on gaps in sequential dates and are selected based on the maximum average of the absolute values of the predicted coefficients.
#'
#' @param temp_window A numeric vector representing the temperature window or sequence of dates. This vector should be ordered.
#' @param pred_C A data frame or matrix containing the predicted coefficients. It must include a column named `fit` that contains the predicted values corresponding to the dates in `temp_window`.
#'
#' @details
#' The function first calculates the differences between consecutive dates in `temp_window` and identifies gaps greater than 1. It then segments the temperature window into periods based on these gaps. Each segment is assigned a unique identifier. The function then computes the mean absolute value of predicted coefficients for each segment and selects the period with the maximum average. This period is returned as the most extreme.
#'
#' @return A numeric vector representing the temperature window where the average of the absolute predicted coefficients is the highest.
#'
#' @examples
#' # Example data
#' temp_window <- 1:10
#' pred_C <- data.frame(fit = c(2, 3, 1, 5, 7, 8, 6, 4, 9, 2))
#' 
#' # Find the concurrent period with the most extreme average
#' find_concurrent_period(temp_window, pred_C)
#'
find_concurrent_period=function(temp_window,pred_C){
  
  #create dummy index
  idx=1
  #find differences between submitted dates and note which are greater than 1
  diff_id = c(1,which(diff(temp_window)>1))
  
  if(length(diff_id)>1){if((which(diff(temp_window)>1))[1]==1){diff_id=diff_id[-1]}}
  if(length(temp_window)>1){
    #Find all the series of sequential dates and give unique id
    for(rv in 1:length(diff_id)){
      
      if(rv==length(diff_id)){
        idx=c(idx,rep(rv,length((diff_id[rv]+1):length(temp_window))))
      }else{
        idx=c(idx,rep(rv,length((diff_id[rv]+1):(diff_id[rv+1]))))
      }
      
    }
  }
  
  #from the estimated coefficients and the concurrent window indices (idx) 
  #find which period has the most extreme average. call that the period to use
  mx_coef = which.max(tapply(abs(pred_C$fit[temp_window]),idx,mean))
  
  return(temp_window[idx==mx_coef])
  
}

#' Get Window Ranges from Predictions
#'
#' This function calculates prediction windows for given Generalized Additive Models (GAMs) based on slopes and R-squared values. It determines the range of days where the predicted slopes and R-squared values fall within specific quantiles, and then identifies the concurrent period with extreme values.
#'
#' @param slope_gam A GAM model object used for predicting slopes. This model should be fitted using `mgcv::gam` and provide predictions related to the response variable of interest.
#' @param rs_gam A GAM model object used for predicting R-squared values. This model should also be fitted using `mgcv::gam` and provide predictions related to the goodness-of-fit of the model.
#' @param temporary A placeholder parameter that is not used in the current function implementation. It is included for potential future use or consistency with other function signatures.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Predicts values and standard errors for the slope and R-squared GAM models.
#'   \item Identifies the lower and upper range of coefficients where the predictions fall within the 2.5th and 97.5th percentiles of the fitted values.
#'   \item Identifies the range of days where the R-squared values fall within the 97.5th percentile.
#'   \item Finds the intersection of these ranges to determine the concurrent periods with extreme values.
#' }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pred_S}: A list with predicted values and standard errors for the slope GAM model.
#'   \item \code{pred_R}: A list with predicted values and standard errors for the R-squared GAM model.
#'   \item \code{days}: A numeric vector representing the days where the predictions meet the specified criteria.
#' }
#'
#' @examples
#' # Example usage
#' # Assuming slope_gam and rs_gam are previously defined GAM models
#' result <- get_predictions_windows(slope_gam = my_slope_gam, rs_gam = my_rs_gam, temporary = NULL)
#' 
#' # Accessing results
#' pred_S <- result$pred_S
#' pred_R <- result$pred_R
#' days <- result$days
#'
get_predictions_windows <- function(slope_gam, rs_gam, temporary) {
  pred_S <- predict(slope_gam, se.fit = TRUE, type = "response")
  pred_R <- predict(rs_gam, se.fit = TRUE, type = "response")
  
  coef_range_l <- which(((pred_S$fit - 1.96 * pred_S$se.fit) - 100) <= quantile(pred_S$fit - 100, 0.025))
  coef_range_u <- which(((pred_S$fit + 1.96 * pred_S$se.fit) - 100) >= quantile(pred_S$fit - 100, 0.975))
  r_sq_range <- which(((pred_R$fit + 1.96 * pred_R$se.fit) - 100) >= quantile(pred_R$fit - 100, 0.975))
  
  temp_window_l <- coef_range_l[is.element(coef_range_l, r_sq_range)]
  temp_window_u <- coef_range_u[is.element(coef_range_u, r_sq_range)]
  
  if (length(temp_window_l) == 0) temp_window_l <- coef_range_l
  if (length(temp_window_u) == 0) temp_window_u <- coef_range_u
  
  days <- find_concurrent_period(temp_window = c(temp_window_l, temp_window_u), pred_C = pred_S)
  list(pred_S = pred_S, pred_R = pred_R, days = days)
}

#' Optimize and Fit Generalized Additive Models (GAMs)
#'
#' This function fits Generalized Additive Models (GAMs) to predict slope and R-squared values based on the provided data. It includes an option to optimize the number of knots used in the smooth term of the GAM. If optimization is enabled, it searches for the optimal number of knots that results in significant p-values for the model's smooth terms.
#'
#' @param temporary A data frame containing the data to fit the GAM models. It must include at least two columns: \code{slope} (the response variable for the slope model) and \code{r_s} (the response variable for the R-squared model), as well as \code{day} (the predictor variable).
#' @param optim.k A logical flag indicating whether to optimize the number of knots (\code{k}) for the smooth term. If \code{TRUE}, the function will search for the optimal number of knots. Default is \code{TRUE}.
#' @param plots A logical flag indicating whether to plot the fitted GAM models. If \code{TRUE}, plots of the fitted models and their residuals will be displayed. Default is \code{FALSE}.
#' @param k An integer specifying the number of knots to use if \code{optim.k} is \code{FALSE}. Default is \code{20}. If \code{optim.k} is \code{TRUE}, this parameter is ignored.
#'
#' @details
#' If \code{optim.k} is \code{TRUE}, the function iterates over a range of possible knot values (from 10 to 365). For each value, it fits a GAM model and checks the significance of the smooth term using \code{k.check}. The first value of \code{k} that results in significant p-values for all tests is chosen as the optimal number of knots. If no optimal value is found, the default value of \code{k = -1} is used.
#' 
#' The function fits two GAM models: one for the \code{slope} and one for \code{r_s}, both with the chosen number of knots. If \code{plots} is \code{TRUE}, it plots the fitted models and their residuals using \code{cowplot}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{slope_gam}: The fitted GAM model for the slope.
#'   \item \code{rs_gam}: The fitted GAM model for R-squared.
#'   \item \code{k}: The number of knots used in the models.
#' }
#'
#' @examples
#' # Example data
#' set.seed(123)
#' temporary <- data.frame(
#'   day = 1:365,
#'   slope = rnorm(365),
#'   r_s = rnorm(365)
#' )
#'
#' # Fit GAM models with automatic optimization of knots
#' result <- optimize_and_fit_gam(temporary, optim.k = TRUE, plots = TRUE)
#' 
#' # Fit GAM models with a fixed number of knots
#' result_fixed_k <- optimize_and_fit_gam(temporary, optim.k = FALSE, k = 30, plots = TRUE)
#'
optimize_and_fit_gam <- function(temporary, optim.k = TRUE, plots = F, k = 20) {
  if (optim.k) {
    # Function to check significance in k.check
    is_significant <- function(check) {
      p_values <- check[,'p-value']
      all(p_values < 0.05) # Returns TRUE if all p-values are significant
    }
    
    # Range of k values to try
    k_values <- seq(10, 365)  # Adjust to one per day 
    optimal_k <- NULL
    
    for (k in 1:length(k_values)) {
      # Fit the model
      slope_gam <- gam(slope ~ s(day, k = k_values[k], bs = "cr"), data = temporary)
      
      # Perform k.check
      check <- k.check(slope_gam)
      pvalue <- check[,'p-value']
      kfin <- check[,"k'"]
      
      # Check if all p-values are significant
      if (!is_significant(check)) {
        optimal_k <- k_values[k]
        break
      }
    }
    if (!is.null(optimal_k)) {
      cat("Optimal k found:", optimal_k, "\n")
      slope_gam <- gam(slope ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      rs_gam <- gam(r_s ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      k <- optimal_k
    } else {
      cat("No optimal k found within the specified range. Using default k = -1.\n")
      slope_gam <- gam(slope ~ s(day, k = -1, bs = "cr"), data = temporary)
      rs_gam <- gam(r_s ~ s(day, k = -1, bs = "cr"), data = temporary)
      k <- -1
    }
  } else {
    slope_gam <- gam(slope ~ s(day, k = k, bs = "cr"), data = temporary)
    rs_gam <- gam(r_s ~ s(day, k = k, bs = "cr"), data = temporary)
    k <- k
  }
  
  if (plots) {
    results <- cowplot::plot_grid(
      draw(slope_gam, residuals = TRUE) + ylab('Slope (partial effect)'), 
      draw(rs_gam, residuals = TRUE) + ylab('R2 (partial effect)')
    )
    print(results)
  }
  
  list(slope_gam = slope_gam, rs_gam = rs_gam, k = k)
  
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
    mutate(correlation = map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$estimate)) %>% 
    mutate(pvalue.cor = map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$p.value))
  
  
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
           tidied = purrr::map(model, tidy),
           glanced = purrr::map(model, glance),
           augmented = purrr::map(model, augment))
  
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
#' @param site.name A string specifying the site name for which to perform the analysis. This should match entries in the `plotname.lon.lat` column of `Fagus.seed`.
#' @param Fagus.seed A data frame containing seed count and site-level information, including `plotname.lon.lat` and `log.seed`.
#' @param climate.beech.path A string specifying the path to the directory containing climate data files. These files should include the site name in their names.
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
#'   Fagus.seed = Fagus_seed,
#'   climate.beech.path = climate_beech_path,
#'   lastdays = last_days,
#'   myform = my_form
#' )
#'
site.moving.climate.analysis <- function(site.name, 
                                         Fagus.seed, 
                                         climate.beech.path, 
                                         lastdays, 
                                         myform) {
  # Load the biological data for the site
  bio_data <- Fagus.seed %>%
    filter(plotname.lon.lat == site.name) %>%
    as.data.frame() 
  
  # Climate load and format
  climate_beech_unique <- list.files(path = climate.beech.path, full.names = TRUE, pattern = site.name)
  
  
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
  
  # Define the year period
  yearneed <- 2
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  
  # Apply the function across all years in yearperiod and combine results
  rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate_csv, 
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
#' This function applies the moving climate analysis to all sites listed in the `Fagus.seed` data frame. It iterates over each site, performs the moving window analysis, and combines the results into a single data frame.
#'
#' @param Fagus.seed A data frame containing seed count and site-level information. It should include a column `plotname.lon.lat` to specify site names and a column `log.seed` for seed counts.
#' @param climate.beech.path A string specifying the path to the directory containing climate data files. These files should be named according to the site names.
#' @param refday An integer representing the reference day for the rolling window analysis. Default is 305.
#' @param lastdays An integer specifying the number of days to include in the rolling window analysis. Default is the maximum of a specified range.
#' @param rollwin An integer specifying the rolling window size. Default is 1.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Identifies unique site names from the `Fagus.seed` data frame.
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
#'   Fagus.seed = Fagus_seed,
#'   climate.beech.path = climate_beech_path,
#'   refday = 305,
#'   lastdays = 30,
#'   rollwin = 1
#' )
#'
FULL.moving.climate.analysis <- function(Fagus.seed = Fagus.seed,
                                         climate.beech.path = climate.beech.path,
                                         refday = 305,
                                         lastdays = max(range),
                                         rollwin = 1) {
  al.sites <- unique(Fagus.seed$plotname.lon.lat)
  
  results.moving <- map_dfr(al.sites, 
                            ~site.moving.climate.analysis(site.name = .x, 
                                                          Fagus.seed = Fagus.seed, 
                                                          climate.beech.path = climate.beech.path,                            
                                                          lastdays = lastdays,
                                                          myform = formula('log.seed~rolling_avg_tmean')))
  
  return(results.moving)
}
