#' Check Folder Existence and Contents
#'
#' This function checks whether a folder exists at the specified path and optionally verifies that specific files are present in the folder.
#'
#' @param folder_path A string specifying the path to the folder.
#' @param file_pattern An optional regular expression pattern to match file names in the folder (e.g., "*.csv" for all CSV files).
#' 
#' @return A list with two elements:
#' \item{folder_exists}{A logical indicating whether the folder exists.}
#' \item{files_present}{A character vector of files present in the folder (filtered by `required_files` or `file_pattern`, if provided).}
#' @examples
#' \dontrun{
#' # Check if folder exists and if "file1.csv" and "file2.csv" are in the folder
#' check_folder_and_contents("/path/to/folder", required_files = c("file1.csv", "file2.csv"))
#' 
#' # Check if folder exists and if any CSV files are in the folder
#' check_folder_and_contents("/path/to/folder", file_pattern = "*.csv")
#' }
#' @export
check_folder_and_contents <- function(folder_path, file_pattern = NULL) {
  # Check if the folder exists
  folder_exists <- dir.exists(folder_path)
  
  if (!folder_exists) {
    return(list(folder_exists = FALSE, files_present = NULL))
  }
  
  # Get the list of files in the folder
  folder_contents <- list.files(path = folder_path, full.names = FALSE)
  # If a pattern is provided, filter the files by pattern
  if (!is.null(file_pattern)) {
    files_present <- folder_contents[grepl(file_pattern, folder_contents)]
  } else {
    # Return all files if no specific requirement
    files_present <- folder_contents
  }
  
  return(list(folder_exists = TRUE, files_present = files_present))
}