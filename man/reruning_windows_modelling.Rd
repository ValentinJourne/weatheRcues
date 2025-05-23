% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reruning_windows_modelling.R
\name{reruning_windows_modelling}
\alias{reruning_windows_modelling}
\title{Re-run Climate Window Model for a Specific Window Sequence}
\usage{
reruning_windows_modelling(
  z,
  bio_data = bio_data,
  window_ranges_df = window_ranges_df,
  rolling.data = rolling.data,
  formula_model = formula("log.seed ~ mean.temperature"),
  model_type = "lm",
  refday = 305,
  rollwin = 1
)
}
\arguments{
\item{z}{Integer. The index (row) of the \code{window_ranges_df} data frame indicating which climate window to evaluate.}

\item{bio_data}{Data frame. Biological site-level data containing at least \code{year}, \code{sitenewname}, and \code{plotname.lon.lat}, as well as the response variable in \code{formula_model}.}

\item{window_ranges_df}{Data frame. Contains window definitions (e.g., \code{window.open}, \code{window.close}) and optional metadata like \code{windows.sequences.number}.}

\item{rolling.data}{Data frame. Rolling climate data with columns \code{days.reversed}, \code{year}, \code{LONGITUDE}, \code{LATITUDE}, and the climate covariate used in \code{formula_model}.}

\item{formula_model}{Formula. A model formula (e.g., \code{log.seed ~ mean.temperature}) specifying the response and climate covariate to be modeled.}

\item{model_type}{Character. The type of model to fit. Options are \code{"lm"} (default) for linear regression or \code{"betareg"} for beta regression.}

\item{refday}{Integer. The reference day (day of year) from which \code{days.reversed} were calculated. Included for metadata.}

\item{rollwin}{Integer. The size of the rolling window applied to the climate data. Included for metadata.}
}
\value{
A data frame (1 row) summarizing the model results for window \code{z}, including:
\itemize{
\item Site info: \code{sitenewname}, \code{plotname.lon.lat}
\item Window info: \code{window.open}, \code{window.close}, \code{windows.size}, \code{nsequence.id}
\item Model results: \code{intercept.estimate}, \code{slope.estimate}, \code{pvalue}, \code{r2}, \code{AIC}, \code{sigma}, \code{nobs}
}
}
\description{
This function re-fits a regression model for a given climate window (defined by index \code{z}) using previously calculated window bounds.
It aggregates the relevant climate covariate over the selected window, merges it with biological data, and fits a model (either linear or beta regression) to assess the effect.
}
\details{
The function:
\itemize{
\item Extracts the corresponding open/close bounds for window \code{z}.
\item Aggregates the covariate over the selected time window.
\item Joins the climate summary with biological observations by \code{year}.
\item Fits the model and extracts relevant statistics (e.g., coefficients, R², AIC, etc.).
}
If \code{model_type = "betareg"}, ensure the response variable is bounded between 0 and 1.
}
\examples{
\dontrun{
result <- reruning_windows_modelling(
  z = 1,
  bio_data = bio_data,
  window_ranges_df = window_ranges_df,
  rolling.data = rolling_climate,
  formula_model = formula('log.seed ~ TMEAN'),
  model_type = 'lm'
)
}

}
\seealso{
\code{\link{save_window_ranges}}, \code{\link{reformat_climate_backtothepast}}, \code{\link{DescTools::CCC}}, \code{\link{broom::tidy}}
}
