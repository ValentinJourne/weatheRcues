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