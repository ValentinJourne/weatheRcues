#' Run Climate Window Analysis (based on `climwin`)
#'
#' This function performs a climate window analysis using the `climwin` package to identify
#' the optimal climate window (i.e., period of influence) that best explains variations in a biological response variable,
#' typically seed production. The analysis uses a sliding window approach to correlate climate variables with biological observations.
#'
#' @param climate_data A data frame or tibble containing climate data. Must include a `date` column and a column matching the name provided in `climate_var` (e.g., `"TMEAN"`).
#' @param bio_data A data frame or tibble containing biological data, including a `Date2` column (date of biological observation) and the response variable used in `formulanull`.
#' @param site.name A character string used to label the output (e.g., the name or ID of the site being analyzed).
#' @param range A numeric vector of length two indicating the maximum and minimum time before the reference date (in days) over which to search for climate windows. Default is `c(600, 0)`.
#' @param cinterval The interval over which to aggregate the climate variable. Options include `"day"`, `"week"`, etc. Default is `"day"`.
#' @param refday Either a numeric DOY (e.g., `305` for November 1) or a vector of day and month (e.g., `c(1, 11)` for November 1). This sets the reference date for `absolute` windows.
#' @param optionwindows One of `"absolute"` or `"relative"`. `"absolute"` windows are fixed in time (e.g., fixed season), while `"relative"` windows move with the biological event. Default is `"absolute"`.
#' @param climate_var Name of the column in `climate_data` containing the climate variable of interest (e.g., `"TMEAN"`). Default is `"TMEAN"`.
#' @param stat.aggregate Aggregation function to apply within the window. One of `"mean"`, `"sum"`, `"min"`, or `"max"`. Default is `"mean"`.
#' @param formulanull A formula specifying the null model (e.g., `log.seed ~ 1`). Must match a column in `bio_data`. Default is `log.seed ~ 1`.
#' @param fun Functional form used to test the climate effect. Common values are `"lin"` (linear), `"quad"` (quadratic), etc. Default is `"lin"`.
#' @param cmissing Method for handling missing climate data. See `climwin::slidingwin()` documentation. Default is `"method2"`.
#' @param give.clean Logical. If `TRUE`, returns a cleaned summary of results including slope, intercept, R², AIC, and window bounds. If `FALSE`, returns the full output from `climwin::slidingwin()`. Default is `TRUE`.
#'
#' @return A data frame containing the best-fit window statistics:
#' \itemize{
#'   \item \code{sitenewname}: Site label from `bio_data`
#'   \item \code{climate.file}: Site name provided to `site.name`
#'   \item \code{window.open}, \code{window.close}: Days before the reference date indicating the climate window
#'   \item \code{r2}, \code{AIC}: Performance metrics for the best-fit model
#'   \item \code{slope.estimate}, \code{intercept.estimate}: Coefficient summaries from the best-fit model
#'   \item \code{sigma}: Model residual standard deviation
#' }
#'
#' @details
#' This function is a wrapper around `climwin::slidingwin()` with additional input checks and a standardized output format.
#' It is intended for automated analysis across multiple sites or datasets. The `refday` input supports both numeric DOY and day-month vector formats.
#'
#' @examples
#' \dontrun{
#' result <- runing_climwin_site(
#'   climate_data = my_climate,
#'   bio_data = my_seeds,
#'   site.name = "site1"
#' )
#' }
#'
#' @importFrom dplyr %>%
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
  model_base <- eval(substitute(
    lm(formula.here, data = bio_data),
    list(formula.here = model_formula)
  ))

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
