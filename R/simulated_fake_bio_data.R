#' Simulate Biological Data Based on Climate Windows
#'
#' This function generates simulated biological datasets (e.g., seed production) using specified
#' alpha, beta, and sigma values to create controlled relationships with climate predictors.
#' You can save the simulated results to a specified folder using `.qs` format. Feel free to adjust this code for your own ;D
#'
#' @param fakeclimatewindow A data frame of simulated climate data (must include `year` and `TMEAN`).
#' @param setup.param A named list with numeric ranges for parameters: `alpha`, `beta`, and `sigma`.
#' @param num_simulations Integer. Number of simulated biological datasets to generate.
#' @param save_fake_data Logical. Whether to save the simulations to disk. Default is `TRUE`.
#' @param overwrite Logical. Whether to overwrite existing `.qs` file if it exists. Default is `FALSE`.
#' @param save_path Character. Folder path where the file should be saved. Default is `'outputs'`.
#'
#' @return A list of data frames containing simulated biological datasets with metadata.
#' @export
simulated_fake_bio_data <- function(
  fakeclimatewindow,
  setup.param,
  num_simulations,
  save_fake_data = TRUE,
  overwrite = FALSE,
  save_path = here("Application_MASTREE/outputs")
) {
  results_list <- vector("list", num_simulations)

  for (i in seq_len(num_simulations)) {
    alpha.random.value <- runif(1, setup.param$alpha[1], setup.param$alpha[2])
    beta.random.value <- runif(1, setup.param$beta[1], setup.param$beta[2])
    sigma.random.value <- runif(1, setup.param$sigma[1], setup.param$sigma[2])

    seed_production_simulated <- fakeclimatewindow %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(TMEAN = mean(TMEAN)) %>%
      dplyr::mutate(
        log.seed = alpha.random.value +
          beta.random.value * TMEAN +
          rnorm(n(), mean = 0, sd = sigma.random.value)
      ) %>%
      dplyr::filter(year > startyear) %>%
      dplyr::mutate(
        year = as.numeric(as.character(year)),
        Date = paste0("15/06/", year),
        Date2 = strptime(as.character(Date), format = "%d/%m/%Y"),
        sitenewname = as.character(i),
        plotname.lon.lat = as.character(i)
      )

    r2sim.before.climwin <- summary(
      lm(log.seed ~ TMEAN, data = seed_production_simulated)
    )$adj.r.squared

    result <- seed_production_simulated %>%
      dplyr::mutate(
        r2.before.climwin = r2sim.before.climwin,
        alpha.random.value.setup = alpha.random.value,
        beta.random.value.setup = beta.random.value,
        sigma.random.value.setup = sigma.random.value
      )

    results_list[[i]] <- result
  }

  if (save_fake_data) {
    filename <- paste0("simulated_bio_data_nsim_", num_simulations, ".qs")
    full_path <- file.path(save_path, filename)

    if (file.exists(full_path)) {
      if (overwrite) {
        message("File exists. Overwriting at: ", full_path)
        qs::qsave(results_list, full_path)
      } else {
        message("File exists. Skipping overwrite as 'overwrite = FALSE'.")
      }
    } else {
      message("Saving new file at: ", full_path)
      qs::qsave(results_list, full_path)
    }
  }

  return(results_list)
}
