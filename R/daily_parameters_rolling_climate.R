#' Run Moving Window Regression and Correlation Analysis
#'
#' This function performs a moving window analysis to explore how a climate covariate (e.g., temperature)
#' influences a biological response (e.g., seed production) over time. It combines regression modeling and
#' correlation testing within rolling windows to identify periods of strong climate-response relationships.
#'
#' @param bio_data A data frame or tibble containing biological data. Must include columns for `year`,
#' `sitenewname`, `log.seed`, and optionally `plotname.lon.lat`.
#' @param rolling.data A data frame or tibble with rolling climate variables. Must include a `days.reversed`
#' column and a climate covariate (e.g., `rolling_avg_tmean`).
#' @param method Character. The correlation method to use (`"spearman"`, `"pearson"`, or `"kendall"`).
#' Default is `"spearman"`.
#' @param formula_model A formula defining the relationship to be modeled (e.g., `log.seed ~ rolling_avg_tmean`).
#' The response variable must exist in `bio_data`, and the predictor in `rolling.data`.
#' @param model_type Character. Type of regression model to use. Options are `"lm"` (linear model, default) or `"betareg"`
#' (beta regression for bounded outcomes between 0 and 1).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Merges biological data with rolling climate data.
#'   \item Computes correlation coefficients and their p-values for each moving window.
#'   \item Fits regression models for each time window and extracts slopes, standard errors, R², and log-likelihood.
#'   \item Calculates standard error for the correlation using the sample size.
#' }
#'
#' It is particularly useful for identifying lagged weather cues in time series datasets.
#'
#' @return A data frame including:
#' \itemize{
#'   \item `sitenewname`, `plotname.lon.lat`: Site identifiers.
#'   \item `days.reversed`: Time window in reversed day units.
#'   \item `term`, `estimate`, `std.error`, `p.value`: Model slope summary.
#'   \item `r.squared`, `logLik`: Model performance metrics.
#'   \item `correlation`, `pvalue.cor`, `correlation.se`: Correlation stats.
#' }
#'
#' @examples
#' # Simulate small example
#' data <- data.frame(
#'   year = 2000:2019,
#'   sitenewname = rep("site1", 20),
#'   log.seed = rnorm(20),
#'   plotname.lon.lat = rep("1_1", 20)
#' )
#' rolling.data <- data.frame(
#'   days.reversed = 1:365,
#'   rolling_avg_tmean = rnorm(365),
#'   year = 2017
#' )
#' formula_used <- formula(log.seed ~ rolling_avg_tmean)
#'
#' result <- runing_moving_window_analysis(
#'   bio_data = data,
#'   rolling.data = rolling.data,
#'   method = "spearman",
#'   formula_model = formula_used
#' )
#'
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
          y = .$log.seed,
          x = .[[covariates.of.interest]],
          method = method
        )$estimate
      )
    ) %>%
    dplyr::mutate(
      pvalue.cor = purrr::map(
        data,
        ~ cor.test(
          y = .$log.seed,
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
