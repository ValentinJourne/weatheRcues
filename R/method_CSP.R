#' Perform Climate Window Analysis - based on Climate sensitivity profil from Thackeray et al
#' note that some function have been obtained or adjusted from Thackeray et al paper OR also from Simmonds et al
#' Find Concurrent Periods with Extreme Averages
#'
#' @description This function identifies periods within a given temperature window where the average of predicted coefficients is most extreme. The periods are determined based on gaps in sequential dates and are selected based on the maximum average of the absolute values of the predicted coefficients.
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
#' @export
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
      slope_gam <- mgcv::gam(slope ~ s(day, k = k_values[k], bs = "cr"), data = temporary)
      
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
      slope_gam <- mgcv::gam(slope ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      rs_gam <- mgcv::gam(r_s ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      k <- optimal_k
    } else {
      cat("No optimal k found within the specified range. Using default k = -1.\n")
      slope_gam <- mgcv::gam(slope ~ s(day, k = -1, bs = "cr"), data = temporary)
      rs_gam <- mgcv::gam(r_s ~ s(day, k = -1, bs = "cr"), data = temporary)
      k <- -1
    }
  } else {
    slope_gam <- mgcv::gam(slope ~ s(day, k = k, bs = "cr"), data = temporary)
    rs_gam <- mgcv::gam(r_s ~ s(day, k = k, bs = "cr"), data = temporary)
    k <- k
  }
  
  if (plots) {
    results <- cowplot::plot_grid(
      gratia::draw(slope_gam, residuals = TRUE) + ylab('Slope (partial effect)'), 
      gratia::draw(rs_gam, residuals = TRUE) + ylab('R2 (partial effect)')
    )
    print(results)
  }
  
  list(slope_gam = slope_gam, rs_gam = rs_gam, k = k)
  
}

#note here it is additional function to extract best windows identified by CSP method
#' Run Climate Signal Processing for a Site
#'
#' This function processes climate signal data for a specific site by fitting a generalized additive model (GAM) and running window-based modeling.
#'
#' @param Results_CSPsub A subset of the `Results_CSP` data frame for a specific site. 
#' @param data A data frame containing site-level data, typically including seed production and climate-related variables.
#' @param siteneame.forsub A string representing the name of the site to subset from the `Results_CSP` data and climate data.
#' @param climate.path A string specifying the file path to the folder containing the climate data files.
#' @param refday An integer specifying the reference day (default is 305).
#' @param lastdays An integer representing the last day of the range used for the analysis (default is the maximum of the range).
#' @param rollwin An integer specifying the size of the rolling window used for calculating temperature averages (default is 1).
#'
#' @return A data frame (`output_fit_summary.temp`) containing the results of the climate signal processing, including fitted model coefficients, model statistics, and window ranges.
#'
#' @details 
#' The function processes climate signal data for a specific site by running the following steps:
#' \itemize{
#'   \item It extracts slopes and R-squared values from `Results_CSPsub` and constructs a temporary data frame for further analysis.
#'   \item It fits a generalized additive model (GAM) to the slopes and R-squared values using the `optimize_and_fit_gam` function.
#'   \item It identifies consecutive sequences of significant days (`extract_consecutive_sequences`) and generates window ranges.
#'   \item It loads the climate data for the site from the specified path, formats the dates, and scales the climate variables.
#'   \item It computes rolling temperature data for the selected time period.
#'   \item It applies the `reruning_windows_modelling` function across all identified window ranges to model the relationship between seed production and climate variables.
#' }
#' 
#' @note The function assumes that the `data` argument corresponds to the site-level data and contains a column named `plotname.lon.lat` for site identification.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' output <- runing_csp_site(Results_CSPsub = Results_CSP[i],
#'                           data = seed.data,
#'                           siteneame.forsub = "site_123",
#'                           climate.path = "/path/to/climate/data")
#' }
#'
#' @import dplyr
#' @import qs
#' @import purrr
#' @import lubridate
#' @export
runing_csp_site = function(Results_CSPsub = Results_CSPsub,
                           data = data,
                           siteneame.forsub = siteneame.forsub,
                           option.check.name = TRUE,
                           climate_csv = climate_csv,
                           refday = 305,
                           lastdays = max(range),
                           rollwin = 1,
                           optim.k = F,
                           variablemoving = 'TMEAN'){
  
  
  list_slope <- as.list(Results_CSPsub$estimate)
  list_rs <- as.list(Results_CSPsub$r.squared)
  
  day = seq(rollwin,(lastdays-1),1)
  slope = list_slope
  r_s = list_rs
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), 
                          day = day, 
                          r_s = unlist(r_s))
  
  #do not optimize, too long 
  results <- optimize_and_fit_gam(temporary, optim.k = optim.k, plots = F, k = (nrow(data)-1) ) #(12+6) #I will specify just number of month 
  
  days = get_predictions_windows(slope_gam = results[[1]], 
                                 rs_gam = results[[2]], 
                                 temporary)$days
  sequences_days = extract_consecutive_sequences(days, keep_all = TRUE)
  window_ranges_df <- save_window_ranges(sequences_days) %>% 
    mutate(windows.sequences.number = 1:nrow(.))
  
  if(option.check.name == T){
    if((siteneame.forsub == unique(data$sitenewname)|siteneame.forsub == unique(data$plotname.lon.lat))==F){
      stop()
    }
  }
  
  if(option.check.name==T & is.na(siteneame.forsub) == T){
    stop(paste("siteneame.forsub argument for subseting site is missing :O "))
  }

  
  # Define the year period
  yearneed <- 2
  yearperiod <- (min(climate_csv$year) + yearneed):max(climate_csv$year)
  rolling.temperature.data <- purrr::map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate_csv, 
                                      yearneed = yearneed, 
                                      refday = refday, 
                                      lastdays = lastdays, 
                                      rollwin = rollwin, 
                                      variablemoving = variablemoving)
  output_fit_summary.temp <- purrr::map_dfr(1:nrow(window_ranges_df), ~reruning_windows_modelling(.,tible.sitelevel = data, 
                                                                                           window_ranges_df = window_ranges_df,
                                                                                           rolling.temperature.data = rolling.temperature.data,
                                                                                           myform.fin = formula('log.seed ~ mean.temperature')))
  return(output_fit_summary.temp)
  
}

#' Perform Climate Signal Processing (CSP) Site Analysis
#'
#' This function performs a CSP site analysis by filtering the provided seed data for a specific site
#' and then calling the `runing_csp_site` function to analyze climate and seed data over a specified time window.
#'
#' @param Results_CSPsub A subset of results for a specific site. Typically a data frame or list that contains climate and seed-related estimates.
#' @param siteneame.forsub A character string indicating the site name (longitude and latitude) used to filter the seed dataset.
#' @param seed.data A data frame containing seed data for multiple sites, including columns like `plotname.lon.lat` that are used to filter for a specific site.
#' @param climate.path A string specifying the file path to the climate data related to seed production in trees.
#' @param refday An integer representing the reference day of the year used to align the rolling climate window.
#' @param lastdays An integer representing the total number of days for the climate window analysis.
#' @param rollwin An integer specifying the size of the rolling window for the climate data (number of days to average).
#'
#' @details
#' This function filters the provided `seed.data` data by the site specified in `siteneame.forsub` and 
#' then runs the `runing_csp_site` function with the filtered seed data, climate data, and other parameters. 
#' The `runing_csp_site` function processes climate data within a rolling window to find relationships between 
#' climate signals and seed data.
#'
#' @return
#' A data frame or list, which is the output of the `runing_csp_site` function. It typically contains model estimates, 
#' climate window ranges, and other statistical summaries related to the climate-seed analysis.
#'
#' @examples
#' \dontrun{
#' # Example usage of the CSP_function_site function
#' CSP_function_site(
#'   Results_CSPsub = Results_CSP[[1]], 
#'   siteneame.forsub = "51.0_10.0", 
#'   seed.data = seed.data, 
#'   climate.path = "data/climate_beech", 
#'   refday = 305, 
#'   lastdays = 600, 
#'   rollwin = 1
#' )
#' }
#'
#' @export
CSP_function_site <- function(Results_CSPsub, 
                              siteneame.forsub, 
                              seed.data, 
                              climate.path, 
                              refday, lastdays, rollwin) {
  
  # Filter the Fagus seed data by the site name
  data.sub.fagus <- seed.data %>%
    dplyr::filter(sitenewname == siteneame.forsub | plotname.lon.lat == siteneame.forsub)
  
  #extract climate matching site
  climate_data <- format_climate_data(
    site = unique(data.sub.fagus$plotname.lon.lat), 
    path = climate.path, 
    scale.climate = TRUE
  )
  
  # Run the CSP site analysis function
  runing_csp_site(Results_CSPsub = Results_CSPsub,
                  data = data.sub.fagus,
                  siteneame.forsub = siteneame.forsub,
                  climate_csv = climate_data,
                  refday = refday,
                  lastdays = lastdays,   # Change from max(range) to lastdays
                  rollwin = rollwin)
}
