#' Format Climate Data for a Specified Site
#'
#' This function reads climate data from a specified directory for a given site, applies necessary transformations,
#' and optionally scales the climate variables. The function can read files in either `.qs` or `.csv` format.
#'
#' @param site Character. The site identifier used to match climate data files (e.g., `"site_name"`). Must be a single string.
#' @param path Character. The directory path where the climate data files are stored. Must be a valid path as a single string.
#' @param scale.climate Logical. Should the climate variables (`TMEAN`, `TMAX`, `TMIN`, `PRP`) be standardized? Defaults to `TRUE`.
#' @param file_type Character. The file type of the climate data. Must be either `"qs"` or `"csv"`. Defaults to `"qs"`.
#'
#' @details The function performs the following steps:
#' - Reads the climate data file for the specified site using `qs::qread()` for `.qs` files or `readr::read_csv()` for `.csv` files.
#' - Converts the `DATEB` field to a proper date format.
#' - Applies a transformation using the `foo()` function (assumed to adjust the year based on an origin date of 1949).
#' - Extracts the year and day of the year (`yday`) from the transformed date.
#' - Converts temperature (`TMEAN`, `TMAX`, `TMIN`) and precipitation (`PRP`) variables to numeric format.
#'
#' If `scale.climate` is `TRUE`, the climate variables (`TMEAN`, `TMAX`, `TMIN`, `PRP`) are standardized using `scale()`.
#' The standardized values are then converted to vectors using `as.vector()`.
#'
#' @return A `data.frame` containing the formatted climate data with columns:
#' - `DATEB`: Original date (before transformation).
#' - `date`: Transformed date.
#' - `yday`: Day of the year.
#' - `year`: Year extracted from the date.
#' - `TMEAN`, `TMAX`, `TMIN`, `PRP`: Numeric values for mean temperature, maximum temperature, minimum temperature, and precipitation, respectively.
#'
#' If `scale.climate` is `TRUE`, these variables are standardized.
#'
#' @importFrom dplyr mutate across
#' @importFrom qs qread
#' @importFrom readr read_csv
#' @importFrom lubridate yday year
#' @export
#'
#' @examples
#' \dontrun{
#' # Format climate data for site "example_site" located at "/data/climate"
#' formatted_data <- format_climate_data(
#'   site = "example_site",
#'   path = "/data/climate",
#'   scale.climate = TRUE,
#'   file_type = "qs"
#' )
#' }
#' @export
format_climate_data <- function(
  site = 'sitename',
  path = 'D/mypath',
  scale.climate = TRUE,
  file_type = "qs"
) {
  # Input validation checks
  if (!is.character(site) || length(site) != 1) {
    stop("Error: 'site' should be a single string representing the site name.")
  }

  if (!is.character(path) || length(path) != 1 || !dir.exists(path)) {
    stop("Error: 'path' should be a valid directory path in character format.")
  }

  if (!is.logical(scale.climate) || length(scale.climate) != 1) {
    stop(
      "Error: 'scale.climate' should be a single logical value (TRUE or FALSE)."
    )
  }
  if (!file_type %in% c("qs", "csv")) {
    stop("Error: 'file_type' should be either 'qs' or 'csv'.")
  }

  #get the fil path matching site name
  file_path <- list.files(path = path, full.names = TRUE, pattern = site)

  # Load climate data based on file type
  if (file_type == "qs") {
    climate_data <- qs::qread(file_path) %>%
      as.data.frame()
  } else if (file_type == "csv") {
    climate_data <- readr::read_csv(file_path) %>%
      as.data.frame()
  }

  climate_data <- climate_data %>%
    dplyr::mutate(
      DATEB = as.Date(DATEB, format = "%m/%d/%y"),
      date = foo(DATEB, 1949),
      yday = yday(date),
      year = year(date),
      TMEAN = as.numeric(TMEAN),
      TMAX = as.numeric(TMAX),
      TMIN = as.numeric(TMIN),
      PRP = as.numeric(PRP)
    )

  if (scale.climate) {
    climate_data <- climate_data %>%
      dplyr::mutate(across(c(TMEAN, TMAX, TMIN, PRP), scale)) %>%
      dplyr::mutate(across(c(TMEAN, TMAX, TMIN, PRP), as.vector))
  }

  return(climate_data)
}
