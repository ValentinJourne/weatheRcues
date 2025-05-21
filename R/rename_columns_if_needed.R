#' Rename Climwin Output Columns If Present
#'
#' This helper function checks whether the columns `WindowOpen` and `WindowClose` exist in a data frame
#' and renames them to `window.open` and `window.close`, respectively. It is useful for standardizing
#' output from the `climwin` package.
#'
#' @param df A data frame potentially containing the columns `WindowOpen` and `WindowClose`.
#'
#' @return A data frame with renamed columns if applicable.
#'
#' @examples
#' df <- data.frame(WindowOpen = 100, WindowClose = 200)
#' rename_columns_if_needed(df)
#'
#' @export
rename_columns_if_needed <- function(df) {
  if ("WindowOpen" %in% colnames(df)) {
    df <- df %>%
      dplyr::rename(window.open = WindowOpen)
  }
  if ("WindowClose" %in% colnames(df)) {
    df <- df %>%
      dplyr::rename(window.close = WindowClose)
  }
  return(df)
}
