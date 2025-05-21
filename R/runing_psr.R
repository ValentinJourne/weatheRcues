#' Run P-spline Regression (PSR) to Identify Weather Cues
#'
#' This function applies Penalized Spline Regression (PSR) to identify temporal windows during which climate variables (e.g., temperature) significantly influence a biological response variable (e.g., seed production) at a given site. The method is adapted from Roberts et al. and Simmonds et al., and is suitable for exploring time-lagged climate–biology relationships using generalized additive models (GAMs) with a penalized spline term.
#'
#' @param bio_data A data frame containing biological data, including columns such as `year`, `plotname.lon.lat`, `sitenewname`, and the response variable (e.g., `log.seed`).
#' @param site Character. The site identifier, which must match values in the `plotname.lon.lat` column of `bio_data`.
#' @param climate_data A data frame containing climate data, with columns for `year`, `yday`, and the climate variable of interest (e.g., `TMEAN`).
#' @param lastdays Integer. The number of days before the reference date to include in the analysis. Default is 600.
#' @param refday Integer. The reference day of year (DOY) from which to track backward. Default is 305 (i.e., November 1st).
#' @param rollwin Integer. Size of the rolling window to smooth climate data. Default is 1 (i.e., no smoothing).
#' @param formula_model A formula describing the relationship between the biological response and the climate predictor. Default is \code{log.seed ~ TMEAN}.
#' @param matrice Numeric vector of length 2. Controls the order of the B-spline and the penalty difference (see \code{mgcv::s}). Default is \code{c(3, 1)}.
#' @param knots Integer or NULL. Number of knots for the smooth spline. If NULL (default), it is set to \code{number of years - 1}.
#' @param tolerancedays Integer. Tolerance in days for identifying contiguous periods of significance. Default is 7.
#' @param yearneed Integer. Number of years of data required for modeling each subset. Default is 2.
#'
#' @details
#' The function reformats the climate data into a matrix where each row corresponds to a year and columns to daily values. A GAM is fit using penalized splines to assess the relationship between the response and predictor across time. Significant periods are identified as those where partial effects exceed ±1.96 standard deviations from the mean effect.
#'
#' If no significant period is found, a data frame with \code{NA} values is returned. If significant periods are found, the function returns a summary data frame for each detected window, including the model's fit statistics and effect estimates.
#'
#' @return A data frame summarizing:
#' \itemize{
#'   \item \code{window.open} and \code{window.close} – Identified window bounds.
#'   \item \code{estimate}, \code{intercept}, \code{r2}, \code{AIC} – Model fit statistics.
#'   \item \code{nobs}, \code{nsequence.id} – Number of observations and identified window sequences.
#' }
#'
#' @examples
#' \dontrun{
#' result <- runing_psr(
#'   bio_data = bio_df,
#'   site = "site_001",
#'   climate_data = climate_df,
#'   lastdays = 600,
#'   refday = 305,
#'   formula_model = log.seed ~ TMEAN
#' )
#' }
#'
#' @seealso \code{\link[mgcv]{gam}}, \code{\link{extract_sequences_auto}}, \code{\link{reruning_windows_modelling}}
#' @export

runing_psr = function(
  bio_data = bio_data,
  site = site,
  climate_data = climate_data,
  lastdays = 600,
  refday = 305,
  rollwin = 1,
  formula_model = formula('log.seed ~ TMEAN'),
  matrice = c(3, 1),
  knots = NULL,
  tolerancedays = 7,
  #plot = TRUE,
  yearneed = 2
) {
  if (!is.data.frame(bio_data)) {
    stop("bio_data must be a data frame or tibble.")
  }
  if (!is.data.frame(climate_data)) {
    stop("climate_data must be a data frame or tibble.")
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
  if (!is.numeric(matrice) || length(matrice) != 2) {
    stop("matrice must be a numeric vector of length 2.")
  }
  if (!is.null(knots) && (!is.numeric(knots) || length(knots) == 0)) {
    stop("knots must be NULL or a numeric vector with at least one element.")
  }
  if (!is.numeric(tolerancedays) || length(tolerancedays) != 1) {
    stop(
      "tolerancedays must be a numeric value of length 1. This arg is about how many days you will tolerate for window identification"
    )
  }
  # if (!is.logical(plot) || length(plot) != 1) {
  #   stop("plot must be a logical value of length 1 (TRUE or FALSE). This will plot the gam predictions")
  # }
  if (!is.numeric(yearneed) || length(yearneed) != 1) {
    stop("yearneed must be a numeric value of length 1.")
  }
  #just checking
  if ((site == unique(bio_data$plotname.lon.lat)) == F) {
    stop()
  }
  for (col in colnames(bio_data)) {
    if (col == "Year") {
      colnames(bio_data)[colnames(bio_data) == "Year"] <- "year"
    }
  }
  for (col in colnames(climate_data)) {
    if (col == "Year") {
      colnames(climate_data)[colnames(climate_data) == "Year"] <- "year"
    }
  }

  covariates.of.interest = as.character(formula_model)[3]
  #now psr method, based from Simmonds et al
  # need climate data to be arranged with year as row
  # need to reduce climate dataframe to only year, yday and temp
  climate2 <- data.frame(
    year = climate_data$year,
    yday = climate_data$yday,
    temp = climate_data[, covariates.of.interest]
  )
  tempmat <- climate2 %>% spread(, key = yday, value = temp)
  tempmat <- tempmat[, -1]
  #number years monitoring seeds
  ny <- length(bio_data$year)
  nt <- lastdays - 1
  ## Formatting data
  index.matrix = matrix(1:nt, ny, nt, byrow = TRUE)
  # Define the year period
  yearneed <- yearneed
  #will fiter the year period here to the year needeed
  yearperiod <- (min(climate_data$year) + yearneed):max(climate_data$year)
  # Apply the function across all years in yearperiod and combine results
  rolling.data <- map_dfr(
    yearperiod,
    reformat_climate_backtothepast,
    climate_data = climate_data,
    yearneed = yearneed,
    refday = refday,
    lastdays = lastdays,
    rollwin = rollwin,
    covariates.of.interest = covariates.of.interest
  )
  #merge data seed to moving climate
  tible.sitelevel = bio_data %>% #site = bio_data
    #rename(year = Year) %>%
    left_join(rolling.data) %>%
    tidyr::drop_na(!!sym(covariates.of.interest))

  climate2 <- data.frame(
    year = tible.sitelevel$year,
    yday = tible.sitelevel$days.reversed,
    covariates.of.interest = tible.sitelevel[, covariates.of.interest]
  )
  covariate.matrix = climate2 %>%
    spread(key = yday, value = covariates.of.interest) %>%
    dplyr::select(-year) %>%
    as.matrix()
  covariate.matrix <- unname(as.matrix(covariate.matrix))
  #second order, but from https://link.springer.com/article/10.1007/s00484-007-0141-4
  #it is possible to adjust between 1-3 (instead of 2) like c(1,1) instead of c(2,1)
  #here I specify to 0 instead of 1 as in SImmonds et al, meaning that I have no penaties on slope
  #if using 1 instead of 0, the model is not able to identify a cues
  #which make sense when looking https://link.springer.com/article/10.1007/s00484-011-0472-z
  if (is.null(knots)) {
    K = ny - 1
    model <- mgcv::gam(
      log.seed ~
        s(index.matrix, k = K, m = matrice, bs = "ps", by = covariate.matrix),
      data = bio_data,
      method = "GCV.Cp",
      control = gam.control(trace = TRUE)
    ) #try GCV.Cp or REML, different MacOS/Window
    summary(model)
  } else {
    model <- gam(
      log.seed ~
        s(index.matrix, k = K, m = c(2, 1), bs = "ps", by = covariate.matrix),
      data = bio_data,
      method = "GCV.Cp",
      control = gam.control(trace = TRUE)
    ) #try GCV.Cp or REML, different MacOS/Window
    summary(model)
  }

  #if(plot==TRUE){
  plotted <- plot(
    model,
    ylab = c('partial effect'),
    xlab = c('days prior response')
  )
  #}
  coefs <- data.frame(fit = plotted[[1]]$fit)
  #wihtout rounding now, will provide different output, and adjust by what is in the main study .
  upper_limit <- mean(coefs$fit) + (1.96 * sd(coefs$fit))
  lower_limit <- mean(coefs$fit) - (1.96 * sd(coefs$fit))
  # Find the markers (days) where the fit exceeds the threshold
  marker <- which(coefs$fit > upper_limit | coefs$fit < lower_limit)

  if (length(marker) == 0) {
    output_fit_summary.psr.temp = data.frame(
      sitenewname = unique(bio_data$sitenewname),
      plotname.lon.lat = unique(bio_data$plotname.lon.lat),
      reference.day = refday,
      windows.size = NA,
      window.open = NA,
      window.close = NA,
      intercept = NA,
      intercept.se = NA,
      estimate = NA,
      estimate.se = NA,
      pvalue = NA,
      r2 = NA,
      AIC = NA,
      nobs = ny,
      nsequence.id = NA
    )
    print('no windows identified')
  } else {
    window <- round(plotted[[1]]$x[find_concurrent_period(marker, coefs)])

    #here does not work because not consecutive
    #extract_consecutive_sequences(window, keep_all = T)
    sequences_days = extract_sequences_auto(window, tolerance = tolerancedays)
    #i guess now will do the same shit as the other methods
    window_ranges_df <- save_window_ranges(sequences_days) %>%
      mutate(windows.sequences.number = 1:nrow(.))

    #it will use the data from roll data climate
    output_fit_summary.psr.temp <- purrr::map_dfr(
      1:nrow(window_ranges_df),
      ~ reruning_windows_modelling(
        .,
        bio_data = bio_data,
        window_ranges_df = window_ranges_df,
        rolling.data = rolling.data,
        formula_model = formula_model,
        refday = refday,
        model_type = 'lm'
      )
    )
  }
  return(output_fit_summary.psr.temp)
}
