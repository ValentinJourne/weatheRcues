#for me to change the name from climwin output
rename_columns_if_needed <- function(df) {
  if ("WindowOpen" %in% colnames(df)) {
    df <- df %>%
      rename(window.open = WindowOpen)
  }
  if ("WindowClose" %in% colnames(df)) {
    df <- df %>%
      rename(window.close = WindowClose)
  }
  return(df)
}
