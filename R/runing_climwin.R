#' Run Climate Window Analysis Using Sliding Window Approach
#'
#' This function wraps the `climwin::slidingwin()` function to detect climate cue windows that best explain biological responses
#' (e.g., seed production) using a sliding window approach. The function fits models across various climate time windows
#' relative to a specified reference date, identifies the best-fitting window based on model performance, and optionally summarizes results.
#'
#' @param climate_data A data frame containing daily climate data. Must include a `date` column (class `Date`) and one column for the climate variable (e.g., `TMEAN`).
#' @param bio_data A data frame with biological data. Must include a `Date2` column (class `Date`) and a response variable that matches the left-hand side of `formulanull`.
#' @param site.name A character string used to label the output, e.g., the site or population name.
#' @param range A numeric vector of length 2 indicating the search window, e.g., `c(600, 0)` to test up to 600 days before the reference date.
#' @param cinterval Time unit used for aggregating the climate variable: `"day"` (default), `"week"`, etc.
#' @param refday A vector specifying the day and month of the reference event (e.g., `c(1, 11)` for November 1).
#' @param optionwindows One of `"absolute"` (fixed calendar date) or `"relative"` (relative to biological event). Default is `"absolute"`.
#' @param climate_var A string specifying the name of the climate variable column (e.g., `"TMEAN"`). Must exist in `climate_data`.
#' @param stat.aggregate Aggregation function to apply over the window: `"mean"` (default), `"sum"`, `"min"`, or `"max"`.
#' @param formulanull A formula (e.g., `log.seed ~ 1`) specifying the null model to which the climate models are compared.
#' @param fun Functional form of the climate effect. Use `"lin"` (default) for linear, `"quad"` for quadratic, etc.
#' @param cmissing Method to handle missing values in climate data. See \code{climwin::slidingwin} for options. Default is `"method2"`.
#' @param give.clean Logical. If `TRUE` (default), returns a cleaned summary of the best model results. If `FALSE`, returns the full climwin object.
#'
#' @return If `give.clean = TRUE`, a data frame with:
#' \itemize{
#'   \item \code{sitenewname}: The site identifier
#'   \item \code{climate.file}: Name of the climate file or site
#'   \item \code{window.open}, \code{window.close}: The best-fit window bounds (in days before reference)
#'   \item \code{r2}, \code{AIC}, \code{sigma}: Model performance metrics
#'   \item \code{slope.estimate}, \code{intercept.estimate}: Coefficients from the best model
#' }
#' If `give.clean = FALSE`, the full list returned by `climwin::slidingwin()` is returned instead.
#'
#' @details
#' The function evaluates the effect of climate on a biological variable by systematically testing a series of time windows
#' before a reference date. It compares the performance of models fitted on climate summaries from each window against a null model.
#' The `climwin` package is used internally.
#'
#' This wrapper adds input checking, flexible formula input, and a standardized output format compatible with batch analyses across multiple sites.
#'
#' @seealso \code{\link[climwin]{slidingwin}}
#'
#' @examples
#' \dontrun{
#' runing_climwin(
#'   climate_data = daily_temp,
#'   bio_data = seed_data,
#'   site.name = "MySite",
#'   range = c(365, 0),
#'   refday = c(1, 11),
#'   climate_var = "TMEAN",
#'   formulanull = log.seed ~ 1
#' )
#' }
#'
#' @importFrom dplyr %>% left_join mutate rename_with filter select bind_cols
#' @export

runing_climwin <- function(
  ...,
  climate_data,
  bio_data,
  site.name, #character will be added in the final form file
  range = c(600, 0),
  cinterval = 'day',
  refday = c(01, 11),
  optionwindows = 'absolute',
  climate_var = 'TMEAN',
  stat.aggregate = 'mean',
  formulanull = stats::formula('log.seed ~ 1'),
  fun = "lin",
  cmissing = 'method2',
  give.clean = TRUE
) {
  # Test function at the beginning for input validation
  test_inputs <- function() {
    # Check that climate_data and data are data frames
    if (
      !(inherits(climate_data, "data.frame") ||
        inherits(climate_data, "tbl_df"))
    ) {
      stop("climate_data must be a data frame or a tibble")
    }
    if (!(inherits(bio_data, "data.frame") || inherits(bio_data, "tbl_df"))) {
      stop("bio_data (your bio dataset) must be a data frame or a tibble")
    }
    # Check that site.name is a character
    if (!is.character(site.name) || length(site.name) != 1) {
      stop("site.name must be a single character string")
    }
    # Check that climate_data has required columns
    required_columns <- c("date", climate_var)
    missing_columns <- setdiff(required_columns, names(climate_data))
    if (length(missing_columns) > 0) {
      stop(paste(
        "climate_data is missing the required columns of interest:",
        paste(missing_columns, collapse = ", ")
      ))
    }

    # Check that data has required columns
    if (!"Date2" %in% names(bio_data))
      stop("bio_data must contain a 'Date2' column")

    # Check that range, refday are correctly specified
    if (!is.numeric(range) || length(range) != 2)
      stop("range must be a numeric vector of length 2")
    if (!is.numeric(refday) || length(refday) != 2)
      stop("refday must be a numeric vector of length 2")
    if (!optionwindows %in% c('absolute', 'relative')) {
      stop(
        "optionwindows must be either 'absolute' or 'relative'. More details can be find in climwin vignette (https://cran.r-project.org/web/packages/climwin/vignettes/climwin.html, sep 2024)"
      )
    }
    if (!stat.aggregate %in% c('mean', 'sum', 'min', 'max')) {
      stop(
        "stat.aggregate must be one of 'mean', 'sum', 'min', or 'max'. More details can be find in climwin vignette (https://cran.r-project.org/web/packages/climwin/vignettes/climwin.html, sep 2024)"
      )
    }
    if (!inherits(try(as.formula(formulanull), silent = TRUE), "formula")) {
      stop(
        "formulanull must be a valid formula in the form of 'response ~ predictor'. More details can be find in climwin vignette (https://cran.r-project.org/web/packages/climwin/vignettes/climwin.html, sep 2024)"
      )
    }
    if (nrow(climate_data) == 0)
      stop(
        "Your climate file seems to be not in the good shape. Please double check"
      )
    if (nrow(bio_data) == 0)
      stop(
        "Your bio_data file seems to be not in the good shape. Please double check"
      )
  }

  # Run the test function
  test_inputs()

  # Extract the response variable from the formula
  response.drop.na <- strsplit(
    as.character(as.formula(formulanull))[2],
    " ~ "
  )[[1]][1]

  # Debugging: Print the response variable and check if it exists in the data
  print(paste("Response variable:", response.drop.na))
  if (!response.drop.na %in% names(bio_data)) {
    stop(paste(
      "Response variable",
      response.drop.na,
      "not found in the bio_data"
    ))
  }

  # Filter the data to remove rows with NA in the response variable
  # data.ana <- data %>%
  #   tidyr::drop_na(!!sym(response.drop.na))
  #
  # # Debugging: Print the filtered data
  # print(paste("Number of rows in data.ana:", nrow(data.ana)))
  # if (nrow(data.ana) == 0) {
  #   stop("No rows remaining in data.ana after filtering. Check for missing values in the response variable.")
  # }

  # Debugging: Print the first few rows of data.ana
  #print("First few rows of data.ana:")
  #print(head(data.ana))

  #print val
  print(paste("range:", range))
  print(paste("cinterval:", cinterval))
  print(paste("ref day starting:", refday[1], 'with month:', refday[2]))
  print(paste("stat aggregate:", stat.aggregate))
  print(paste("windows option:", optionwindows))
  print(paste("variable of interest:", climate_var))

  #formula <- as.formula(formulanull) #this one is not working, need dynamic one by using eval and subsitute
  model_formula <- formulanull
  #model_base <- eval(substitute(
  #  lm(formula.here, data = bio_data),
  #  list(formula.here = model_formula)
  #))

  model_base <- eval(
    substitute(
      lm(formula.here, data = data_input),
      list(formula.here = model_formula, data_input = bio_data)
    )
  )

  # Run the climwin analysis
  climwin_output <- climwin::slidingwin(
    xvar = list(temperature.degree = climate_data[[climate_var]]),
    cdate = climate_data$date,
    bdate = bio_data$Date2,
    baseline = model_base, #eval(formulanull); i Needed to specify the formula here, if not it is not working properly (with future_map)
    cinterval = cinterval,
    range = range,
    type = optionwindows,
    refday = refday,
    stat = stat.aggregate,
    cmissing = cmissing,
    func = fun
  )

  # Extract the summary for the best model
  if (give.clean == TRUE) {
    broom_summary_slope <- broom::tidy(climwin_output[[1]]$BestModel) %>%
      dplyr::filter(term == 'climate') %>%
      dplyr::select(-term) %>%
      rename_with(.cols = everything(), function(x) {
        paste0("slope.", x)
      })
    broom_summary_intercept <- broom::tidy(climwin_output[[1]]$BestModel) %>%
      dplyr::filter(term != 'climate') %>%
      dplyr::select(-term) %>%
      dplyr::rename_with(.cols = everything(), function(x) {
        paste0("intercept.", x)
      })
    sigma.model = sigma(climwin_output[[1]]$BestModel)

    # Extract performance statistics and combine with site information
    statistics <- dplyr::bind_cols(
      sitenewname = unique(bio_data$sitenewname),
      climate.file = site.name,
      climwin_output$combos, # Extract statistics for all variants
      performance::model_performance(climwin_output[[1]]$BestModel) %>%
        as.data.frame(),
      broom_summary_slope,
      broom_summary_intercept,
      sigma = sigma.model
    ) %>%
      dplyr::rename(
        window.open = WindowOpen,
        window.close = WindowClose,
        r2 = R2
      ) #%>%
    #dplyr::mutate(
    #  window.open = window.open, #a checker, mais je crois que climwin start vector at 0, and me at 1, so maybe need to add +1
    #  window.close = window.close
    #)
  } else {
    #give output from climwin directly
    statistics = climwin_output
  }

  return(statistics)
}
