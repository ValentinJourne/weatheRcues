#' Run a daily relationship between weather and biological variable
#'
#' This function performs a day-by-day analysis of the relationship between a rolling climate variable and a biological response (e.g., seed production). For each day (reversed day format), it fits a regression model and computes correlation metrics to identify lagged climate cues.
#'
#' @param bio_data A data frame containing biological data, including at minimum the columns `year`, `log.seed`, `sitenewname`, and `plotname.lon.lat`.
#' @param rolling.data A data frame containing rolling climate data, including `days.reversed`, `year`, and the rolling climate covariate used in modeling.
#' @param method Character. Correlation method used to estimate correlation coefficients. Options: `"spearman"` (default), `"pearson"`, or `"kendall"`.
#' @param formula_model A formula indicating the response and predictor (e.g., `log.seed ~ TMEAN_rolling`). The response must exist in `bio_data`, and the predictor in `rolling.data`.
#' @param model_type Character. Type of model to fit. Either `"lm"` (default) for linear regression or `"betareg"` for beta regression. The latter assumes the response is bounded between 0 and 1.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Merges `bio_data` and `rolling.data` by year.
#'   \item Fits a linear or beta regression model for each unique value of `days.reversed`.
#'   \item Extracts model slope, standard error, p-value, R², and log-likelihood.
#'   \item Computes correlation and correlation standard error for each time point.
#' }
#'
#' This function is intended to help detect temporal weather cues influencing biological processes across time.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{sitenewname}, \code{plotname.lon.lat}, \code{days.reversed}
#'   \item Model estimates: \code{estimate}, \code{std.error}, \code{p.value}
#'   \item Performance metrics: \code{r.squared}, \code{logLik}
#'   \item Correlation results: \code{correlation}, \code{pvalue.cor}, \code{correlation.se}
#' }
#'
#' @examples
#' \dontrun{
#' result <- daily_parameters_rolling_climate(
#'   bio_data = my_bio_df,
#'   rolling.data = my_climate_rolling,
#'   method = "spearman",
#'   formula_model = log.seed ~ TMEAN_rolling
#' )
#' }
#'
#' @seealso \code{\link{runing_daily_relationship}}, \code{\link{cor.test}}, \code{\link[betareg]{betareg}}, \code{\link[broom]{tidy}}
#' @export
daily_parameters_rolling_climate = function(
  bio_data = bio_data,
  rolling.data = rolling.data,
  method = 'spearman',
  formula_model = formula('log.seed~temperature'),
  model_type = 'lm'
) {
  if (!is.data.frame(bio_data) && !is_tibble(bio_data)) {
    stop("data must be a data frame or tibble.")
  }

  if (!is.data.frame(rolling.data) && !is_tibble(rolling.data)) {
    stop("rolling.data must be a data frame or tibble.")
  }

  if (!model_type %in% c('lm', 'betareg')) {
    stop("model_type must be either 'lm' or 'betareg'")
  }

  if (!inherits(formula_model, "formula")) {
    stop("formula_model must be a valid formula.")
  }

  covariates.of.interest = as.character(formula_model)[3]
  response_var <- all.vars(formula_model)[1]

  if (!covariates.of.interest %in% colnames(rolling.data)) {
    stop(paste("Column", covariates.of.interest, "not found in rolling.data"))
  }

  #for me because sometimes not same name Year or year
  for (col in colnames(bio_data)) {
    if (col == "Year") {
      colnames(bio_data)[colnames(bio_data) == "Year"] <- "year"
    }
  }

  #merge data seed to moving climate
  tible.sitelevel = bio_data %>% #site = bio_data
    #rename(year = Year) %>%
    dplyr::left_join(rolling.data) %>%
    tidyr::drop_na(!!sym(covariates.of.interest))

  #define correlation - calculate correlation and se, extract also p value
  n = tible.sitelevel %>%
    dplyr::select(year, sitenewname) %>%
    distinct() %>%
    nrow()

  correlation.all <- tible.sitelevel %>%
    tidyr::nest(data = -days.reversed) %>%
    dplyr::mutate(
      correlation = purrr::map(
        data,
        ~ cor.test(
          y = .[[response_var]],
          x = .[[covariates.of.interest]],
          method = method
        )$estimate
      )
    ) %>%
    dplyr::mutate(
      pvalue.cor = purrr::map(
        data,
        ~ cor.test(
          y = .[[response_var]],
          x = .[[covariates.of.interest]],
          method = method
        )$p.value
      )
    )

  cortemp = correlation.all %>%
    tidyr::unnest(c(correlation, pvalue.cor)) %>%
    dplyr::select(days.reversed, correlation, pvalue.cor) %>%
    dplyr::mutate(
      correlation.se = correlation.spearman.se(.$correlation, n)
    ) %>%
    dplyr::mutate(sitenewname = unique(tible.sitelevel$sitenewname))

  #use purr for iteration
  #here broom package used, because it is simple lm (and not glmmTMB, need to adjust then )

  fitted_models <- tible.sitelevel %>%
    tidyr::nest(data = -days.reversed) %>%
    mutate(
      model = purrr::map(
        data,
        ~ {
          if (model_type == 'lm') {
            lm(formula_model, data = ., na.action = na.omit)
          } else {
            # Ensure the dependent variable is between 0 and 1 for `betareg`
            if (model_type == 'betareg')
              betareg::betareg(formula_model, data = .)
          }
        }
      ),
      tidied = purrr::map(model, broom::tidy),
      glanced = purrr::map(model, broom::glance),
      augmented = purrr::map(model, broom::augment)
    )

  modelr2 = fitted_models %>%
    tidyr::unnest(glanced) %>%
    dplyr::select(days.reversed, dplyr::contains('r.squared'), logLik)

  octopus = fitted_models %>%
    tidyr::unnest(tidied) %>%
    dplyr::filter(str_detect(term, as.character(formula_model)[3])) %>%
    dplyr::select(days.reversed, term, estimate, std.error, p.value) %>%
    dplyr::mutate(sitenewname = unique(tible.sitelevel$sitenewname)) %>%
    dplyr::left_join(modelr2) %>%
    dplyr::mutate(
      sitenewname = unique(tible.sitelevel$sitenewname),
      plotname.lon.lat = unique(tible.sitelevel$plotname.lon.lat)
    ) %>%
    dplyr::select(
      sitenewname,
      plotname.lon.lat,
      days.reversed,
      everything()
    ) %>%
    dplyr::left_join(cortemp)

  return(octopus)
}
