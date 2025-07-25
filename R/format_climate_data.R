#' Format Climate Data for a Specified Site
#'
#' This function reads climate data from a specified directory for a given site, applies necessary transformations,
#' and optionally scales the climate variables. It supports both `.qs` and `.csv` formats and can flexibly handle different
#' column names for climate variables and dates.
#'
#' @param site Character. The site identifier used to match climate data files (e.g., `"site_name"`). Must be a single string.
#' @param path Character. The directory path where the climate data files are stored. Must be a valid path as a single string.
#' @param scale.climate Logical. Should the climate variables (`TMEAN`, `TMAX`, `TMIN`, `PRP`) be standardized? Defaults to `TRUE`.
#' @param file_type Character. The file type of the climate data. Must be either `"qs"` or `"csv"`. Defaults to `"qs"`.
#' @param origin_year Numeric. The origin year used in the `foo()` date transformation function. Defaults to `1949`.
#' @param date_column Character or `NULL`. Name of the column containing date values. If `NULL`, the function attempts to auto-detect a column of class `Date`.
#' @param variable_mapping Named list. A mapping between standard variable names (`TMEAN`, `TMAX`, `TMIN`, `PRP`) and the actual column names in the input file.
#'                         Defaults to `list(TMEAN = "TMEAN", TMAX = "TMAX", TMIN = "TMIN", PRP = "PRP")`.
#'
#' @details The function performs the following steps:
#' - Locates and loads the file matching the specified site name in the provided directory (`path`).
#' - Reads the file using `qs::qread()` for `.qs` or `readr::read_csv()` for `.csv`. I did not manage for now to include more format
#' - Identifies and parses the date column (using `date_column` or automatic detection).
#' - Renames the climate variables (`TMEAN`, `TMAX`, `TMIN`, `PRP`) according to the user-provided `variable_mapping`.
#' - Applies a transformation to the date using the `foo()` function (e.g., adjusting to a common origin year).
#' - Extracts the year and day of year (`yday`) from the transformed date.
#' - Converts climate variables to numeric format.
#'
#' If `scale.climate = TRUE`, the climate variables are standardized using `scale()` and converted to vectors.
#'
#' @return A `data.frame` containing the formatted climate data with the following columns:
#' - `date`: Transformed date (via `foo()`).
#' - `yday`: Day of the year.
#' - `year`: Extracted year from the transformed date.
#' - `TMEAN`, `TMAX`, `TMIN`, `PRP`: Standardized or raw numeric climate variables.
#'
#' @importFrom dplyr mutate across
#' @importFrom qs qread
#' @importFrom readr read_csv
#' @importFrom lubridate yday year
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with default column names in qs format
#' formatted_data <- format_climate_data(
#'   site = "example_site",
#'   path = "/data/climate",
#'   scale.climate = TRUE,
#'   file_type = "qs"
#' )
#'
#' # Example with custom column names in CSV format
#' formatted_data <- format_climate_data(
#'   site = "example_site",
#'   path = "/data/climate",
#'   file_type = "csv",
#'   variable_mapping = list(
#'     TMEAN = "temperature_C_mean",
#'     TMAX = "temperature_C_max",
#'     TMIN = "temperature_C_min",
#'     PRP = "precipitation_mm"
#'   )
#' )
#' }
#' @export
format_climate_data <- function(
  site = 'sitename',
  path = 'D/mypath',
  scale.climate = TRUE,
  file_type = "qs",
  origin_year = 1949,
  date_column = NULL,
  variable_mapping = list(
    TMEAN = "TMEAN",
    TMAX = "TMAX",
    TMIN = "TMIN",
    PRP = "PRP"
  )
) {
  # Input validation
  if (!is.character(site) || length(site) != 1) stop("site must be a string.")
  if (!is.character(path) || length(path) != 1 || !dir.exists(path))
    stop("path must be a valid directory.")
  if (!is.logical(scale.climate) || length(scale.climate) != 1)
    stop("scale.climate must be TRUE or FALSE.")
  if (!file_type %in% c("qs", "csv")) stop("file_type must be 'qs' or 'csv'.")

  # Find file
  file_path <- list.files(path = path, full.names = TRUE, pattern = site)
  if (length(file_path) == 0) stop("No file found matching site name in path.")

  # Load file, qs or csv only sorry
  if (file_type == "qs") {
    climate_data <- qs::qread(file_path) %>% as.data.frame()
  } else {
    climate_data <- readr::read_csv(file_path, show_col_types = FALSE) %>%
      as.data.frame()
  }

  # Detect or assign date column
  if (is.null(date_column)) {
    date_column_candidates <- names(climate_data)[sapply(
      climate_data,
      inherits,
      "Date"
    )]
    if (length(date_column_candidates) == 0)
      stop("No Date column found. Specify date_column manually.")
    date_column <- date_column_candidates[1]
  }

  # Parse date if not Date class
  if (!inherits(climate_data[[date_column]], "Date")) {
    try_date <- suppressWarnings(as.Date(
      climate_data[[date_column]],
      format = "%Y-%m-%d"
    ))
    if (all(!is.na(try_date))) {
      climate_data[[date_column]] <- try_date
    } else {
      climate_data[[date_column]] <- as.Date(
        climate_data[[date_column]],
        format = "%m/%d/%y"
      )
    }
  }

  # Rename variables to standard names (TMEAN, TMAX, TMIN, PRP) that I need for after
  for (std_name in names(variable_mapping)) {
    actual_name <- variable_mapping[[std_name]]
    if (!(actual_name %in% names(climate_data))) {
      stop(paste("Missing variable in dataset:", actual_name))
    }
    names(climate_data)[names(climate_data) == actual_name] <- std_name
  }

  # Apply transformation
  climate_data <- climate_data %>%
    dplyr::mutate(
      date = foo(.data[[date_column]], origin_year),
      yday = yday(date),
      year = year(date),
      TMEAN = as.numeric(TMEAN),
      TMAX = as.numeric(TMAX),
      TMIN = as.numeric(TMIN),
      PRP = as.numeric(PRP)
    )

  # Optional scaling
  if (scale.climate) {
    climate_data <- climate_data %>%
      dplyr::mutate(across(c(TMEAN, TMAX, TMIN, PRP), ~ as.vector(scale(.))))
  }

  return(climate_data)
}
