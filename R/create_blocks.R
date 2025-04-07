#' Create Non-Overlapping Data Blocks
#' 
#' This function divides a given dataset into a specified number of non-overlapping blocks. 
#' Each block contains a contiguous subset of rows from the dataset. The number of rows in 
#' each block is determined by the total number of rows divided by the number of blocks.
#' 
#' @param data A data frame or tibble to be divided into blocks. The data must be 
#'              arranged such that the order of the rows reflects the desired 
#'              sequential structure (e.g., time series data).
#' @param num_blocks An integer specifying the number of non-overlapping blocks to create 
#'                   from the dataset. Must be a positive integer.
#' 
#' @return A list containing the specified number of blocks, where each block is a 
#'         data frame. If the specified number of blocks cannot be created from the 
#'         dataset, an error is raised.
#' 
#' @examples
#' # Example dataset
#' example_data <- data.frame(year = 2000:2029, value = rnorm(30))
#' 
#' # Create 5 non-overlapping blocks from the dataset
#' blocks <- create_blocks(example_data, num_blocks = 5)
#' 
#' @export
create_blocks <- function(data, num_blocks) {
  n <- nrow(data)
  block_size <- n %/% num_blocks  # Calculate rows per block
  blocks <- list()
  
  for (i in seq(1, n, by = block_size)) {
    blocks[[length(blocks) + 1]] <- data[i:min(i + block_size - 1, n), ]
  }
  
  # Ensure we only return the specified number of blocks
  if (length(blocks) < num_blocks) {
    stop("Not enough data to create the specified number of blocks.")
  }
  
  return(blocks[1:num_blocks])  # Return only the requested number of blocks
}