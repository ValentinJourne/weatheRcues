% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSR_function_site.R
\name{PSR_function_site}
\alias{PSR_function_site}
\title{This function applies the PSR (Penalized Spline Regression) method to identify weather cues for a specific site using biological data. The function is based on methods adapted from Simmonds et al. It filters data for a particular site and calls the \code{runing_psr_site} function to analyze the climate and biological data.}
\usage{
PSR_function_site(
  site,
  seed.data,
  lastdays,
  refday,
  climate.path,
  matrice = c(3, 1),
  knots = NULL,
  tolerancedays = 7,
  yearneed = 2
)
}
\arguments{
\item{site}{A string representing the site name, which should match the \code{plotname.lon.lat} in \code{seed.data}.}

\item{seed.data}{A data frame containing biological data (e.g., seed production) with columns such as \code{plotname.lon.lat}, \code{Year}, and \code{log.seed}.}

\item{lastdays}{Integer specifying the final day for the analysis (usually the length of the time series).}

\item{matrice}{A numeric vector of length 2 indicating the penalties to apply in the smoothing function of the model (default is \code{c(3,1)}).}

\item{knots}{Integer specifying the number of knots for the GAM model. If \code{NULL} (default), the function will set the knots to the number of years minus one.}
}
\value{
A data frame summarizing the significant weather cues for the site based on the PSR method. If no significant windows are found, the function returns a data frame with \code{NA} values.
}
\description{
This function applies the PSR (Penalized Spline Regression) method to identify weather cues for a specific site using biological data. The function is based on methods adapted from Simmonds et al. It filters data for a particular site and calls the \code{runing_psr_site} function to analyze the climate and biological data.
}
\details{
The function filters biological data for the specified site and applies the PSR method to analyze the effect of weather variables (e.g., temperature) on the biological data (e.g., seed production). It uses a rolling window to calculate average climate conditions over a specified period.
}
\examples{
\dontrun{
# Example usage:
site <- "site_name"
seed.data <- your_fagus_data
result <- PSR_function_site(site = site,
                            seed.data = seed.data,
                            lastdays = 600,
                            lastdays = 600,
                            matrice = c(3,1),
                            knots = NULL)
}

}
