#' Perform Climate Window Analysis - based on climwin R package from Bailey et al.
#'
#' This function performs a climate window analysis using the `climwin` R package to identify
#' the best temperature windows that correlate with seed production. It processes climate and
#' biological data, runs a sliding window analysis, and returns statistics for the best windows identified.
#'
#' Note that most parameters are already well described in the climwin R package; however, this
#' function also requires the specification of a null model (baseline). For example:
#' \code{baseline = lm(log.seed ~ 1, data = data)}.
#' Double-check that the response variable is named correctly in the formula.
#'
#' @param climate_data A data frame or tibble containing climate data, with at least two columns:
#' \code{date} and the climate variable specified in \code{climate_var}. This is used for climate window analysis.
#' @param data A data frame or tibble containing biological data including seed production and relevant metadata.
#' It should include a column for \code{Date2}, which represents the date of observation.
#' @param site.name A character string specifying the name of the site for which the analysis is conducted.
#' This is used to filter the climate data specific to this site.
#' @param range A numeric vector of length 2 specifying the range of days over which to perform the sliding
#' window analysis. Defines the maximum and minimum range for window lengths in days. Defaults to \code{c(600, 0)}.
#' @param cinterval A character string specifying the interval at which climate data is aggregated.
#' Options include \code{'day'}, \code{'week'}, etc. Defaults to \code{'day'} (same as climwin).
#' @param refday A numeric vector of length 2 specifying the reference day for absolute windows.
#' Defaults to \code{c(01, 11)} (e.g., November 1st). This means the analysis will search for windows before
#' the 1st of November.
#' @param optionwindows A character string specifying the type of window to use. Options include \code{'absolute'}
#' for fixed dates and \code{'relative'} for windows relative to a specific event. Defaults to \code{'absolute'}.
#' Relative windows are useful if you have a specific event date for seed production.
#' @param climate_var A character string specifying the name of the climate variable to analyze (e.g., temperature
#' mean, maximum, or minimum). Defaults to \code{'TMEAN'} but should match the corresponding column name in the climate data.
#' @param stat.aggregate A character string specifying the aggregation method for the climate variable.
#' Options include \code{'mean'}, \code{'sum'}, \code{'min'}, or \code{'max'}. Defaults to \code{'mean'}.
#' @param formulanull A formula specifying the null model for the analysis. Defaults to \code{formula('log.seed ~ 1')},
#' but this should reflect the response variable in your biological data.
#' @param fun A character string specifying the functional form used for the sliding window analysis. Defaults to \code{"lin"}.
#' @param cmissing A character string specifying the method for handling missing climate data.
#' Defaults to \code{"method2"}.
#'
#' @details
#' This function processes climate data and biological data, runs a sliding window analysis using the `climwin` package,
#' and extracts statistics for the best climate window. The climate data files are assumed to be in a format compatible
#' with the `qs` package (small file size, which can be read using `qs::qread`).
#'
#' The function returns a set of performance metrics for the best-fitting model and the window that corresponds to the
#' best climate variables for seed production.
#'
#' @return A data frame containing the following:
#' \itemize{
#'   \item \code{sitenewname}: Site-specific information from the biological data.
#'   \item \code{climate.file}: Name of the climate data file used.
#'   \item \code{climwin_output$combos}: Results from the sliding window analysis (\code{climwin::slidingwin}).
#'   \item Performance metrics of the best model (from \code{performance::model_performance}).
#'   \item Coefficients of the best model extracted using \code{broom::tidy}.
#' }
#'
#' @examples
#' # Example usage:
#' # result <- climwin_site_days(
#' #   climate_data = climate_data,
#' #   data = biological_data,
#' #   site.name = "site1"
#' # )
#'
#' @export
climwin_site_days <- function(
  ...,
  climate_data,
  data,
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
  give.clean = T
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
    if (!(inherits(data, "data.frame") || inherits(data, "tbl_df"))) {
      stop("data (your bio dataset) must be a data frame or a tibble")
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
    if (!"Date2" %in% names(data)) stop("data must contain a 'Date2' column")

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
    if (nrow(data) == 0)
      stop(
        "Your bio data file seems to be not in the good shape. Please double check"
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
  if (!response.drop.na %in% names(data)) {
    stop(paste("Response variable", response.drop.na, "not found in the data"))
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
    lm(formula.here, data = data),
    list(formula.here = model_formula)
  ))

  # Run the climwin analysis
  climwin_output <- climwin::slidingwin(
    xvar = list(temperature.degree = climate_data[[climate_var]]),
    cdate = climate_data$date,
    bdate = data$Date2,
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
  if (give.clean == T) {
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
      sitenewname = unique(data$sitenewname),
      climate.file = site.name,
      climwin_output$combos, # Extract statistics for all variants
      performance::model_performance(climwin_output[[1]]$BestModel) %>%
        as.data.frame(),
      broom_summary_slope,
      broom_summary_intercept,
      sigma = sigma.model
    ) %>%
      dplyr::rename(window.open = WindowOpen, window.close = WindowClose) %>%
      dplyr::mutate(
        window.open = window.open, #a checker, mais je crois que climwin start vector at 0, and me at 1, so maybe need to add +1
        window.close = window.close
      )
  } else {
    #give output from climwin directly
    statistics = climwin_output
  }

  return(statistics)
}
