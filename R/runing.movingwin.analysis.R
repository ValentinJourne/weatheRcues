#' Perform Moving Window Analysis on Climate Data
#'
#' This function performs a moving window analysis to investigate the relationship between a response variable (e.g., seed count) and a rolling climate variable (e.g., temperature) over time. It calculates correlation coefficients, fits linear models, and extracts relevant statistics.
#'
#' @param data A data frame containing the seed count and other site-level information. It should include columns `Year` (or `year`), `sitenewname`, and `log.seed`.
#' @param rolling.data A data frame containing rolling climate data, including the rolling average temperature and a `days.reversed` column for the moving window analysis.
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
#' rolling.data <- data.frame(
#' days.reversed = 1:365,
#' rolling_avg_tmean = rnorm(365), year = 2017)
#'
#' # Run the moving window analysis
#' result <- running.movingwin.analysis(data = data, rolling.data = rolling.data)
#'
runing.movingwin.analysis = function(data = data,
                                     rolling.data = rolling.data,
                                     method = 'spearman',
                                     myform = formula('log.seed~temperature'),
                                     model_type = 'lm'){
  
  
  if (!is.data.frame(data) && !is_tibble(data)) {
    stop("data must be a data frame or tibble.")
  }
  
  if (!is.data.frame(rolling.data) && !is_tibble(rolling.data)) {
    stop("rolling.data must be a data frame or tibble.")
  }
  
  if (!model_type %in% c('lm', 'betareg')) {
    stop("model_type must be either 'lm' or 'betareg'")
  }
  
  if (!inherits(myform, "formula")) {
    stop("myform must be a valid formula.")
  }
  
  covariates.of.interest = as.character(myform)[3]
  if (!covariates.of.interest %in% colnames(rolling.data)) {
    stop(paste("Column", covariates.of.interest, "not found in rolling.data"))
  }
  
  #for me because sometimes not same name Year or year
  for (col in colnames(data)) {
    if (col == "Year") {
      colnames(data)[colnames(data) == "Year"] <- "year"
    }
  }
  
  #merge data seed to moving climate
  tible.sitelevel = data %>% #site = bio_data 
    #rename(year = Year) %>% 
    dplyr::left_join(rolling.data) %>% 
    tidyr::drop_na(!!sym(covariates.of.interest))
  
  #define correlation - calculate correlation and se, extract also p value 
  n = tible.sitelevel %>% dplyr::select(year, sitenewname) %>% distinct() %>% nrow()
  
  correlation.all <- tible.sitelevel %>% 
    tidyr::nest(data = -days.reversed) %>%
    dplyr::mutate(correlation = purrr::map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$estimate)) %>% 
    dplyr::mutate(pvalue.cor = purrr::map(data, ~cor.test(y=.$log.seed, x=.[[covariates.of.interest]], method = method)$p.value))
  
  
  cortemp = correlation.all %>% 
    tidyr::unnest(c(correlation, pvalue.cor)) %>% 
    dplyr::select(days.reversed, correlation, pvalue.cor) %>% 
    dplyr::mutate(correlation.se = correlation.spearman.se(.$correlation, n)) %>% 
    dplyr::mutate(sitenewname = unique(tible.sitelevel$sitenewname))
  
  #use purr for iteration 
  #here broom package used, because it is simple lm (and not glmmTMB, need to adjust then )
  
  fitted_models <- tible.sitelevel %>%
    tidyr::nest(data = -days.reversed) %>%
    mutate(model = purrr::map(data, ~{
      if (model_type == 'lm') {
        lm(myform, data = ., na.action = na.omit)
      } else {
        # Ensure the dependent variable is between 0 and 1 for `betareg`
        if (model_type == 'betareg')
          betareg::betareg(myform, data = .)
      }
    }),
    tidied = purrr::map(model, broom::tidy),
    glanced = purrr::map(model, broom::glance),
    augmented = purrr::map(model, broom::augment))
  
  
  modelr2 = fitted_models %>%
    tidyr::unnest(glanced) %>%
    dplyr::select(days.reversed, dplyr::contains('r.squared'), logLik) 
  
  octopus = fitted_models %>%
    tidyr::unnest(tidied) %>% 
    dplyr::filter(str_detect(term, as.character(myform)[3])) %>% 
    dplyr::select(days.reversed,term, estimate, std.error, p.value) %>% 
    dplyr::mutate(sitenewname  = unique(tible.sitelevel$sitenewname)) %>% 
    dplyr::left_join(modelr2) %>% 
    dplyr::mutate(sitenewname = unique(tible.sitelevel$sitenewname),
                  plotname.lon.lat = unique(tible.sitelevel$plotname.lon.lat)) %>% 
    dplyr::select(sitenewname, plotname.lon.lat, days.reversed, everything()) %>% 
    dplyr::left_join(cortemp)
  
  return(octopus)
}