% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runing_peak_detection.R
\name{runing_peak_detection}
\alias{runing_peak_detection}
\title{Detect Weather Cue Windows Using Peak Signal Detection}
\usage{
runing_peak_detection(
  lag = 100,
  threshold = 3,
  influence = 0,
  tolerancedays = 7,
  refday = 305,
  lastdays = 600,
  rollwin = 1,
  siteforsub = "longitude=-0.15_latitude=50.85",
  climate_data = climate_data,
  Results_days = Results_days,
  bio_data = bio_data,
  yearneed = 2,
  formula_model = formula("log.seed ~ TMEAN"),
  model_type = "lm"
)
}
\arguments{
\item{lag}{Integer. Number of past days used for calculating the rolling mean and standard deviation in the thresholding algorithm. Default is 100.}

\item{threshold}{Numeric. Z-score threshold. A signal is detected if the product of slope × R² exceeds the rolling mean by this many standard deviations. Default is 3.}

\item{influence}{Numeric between 0 and 1. Determines how strongly new signals influence the rolling statistics. \code{0} means no influence. Default is 0.}

\item{tolerancedays}{Integer. Maximum gap (in days) between signals to be considered as part of the same window. Default is 7.}

\item{refday}{Integer. Reference day-of-year (DOY) to count backwards from when building the climate window. Default is 305 (Nov 1).}

\item{lastdays}{Integer. Number of days to look back from \code{refday}. Default is 600.}

\item{rollwin}{Integer. Size of the rolling window for climate data smoothing. Default is 1.}

\item{siteforsub}{Character. Unique site identifier used for labeling output. Should match identifiers in \code{bio_data}.}

\item{climate_data}{Data frame of daily climate data. Must include \code{year}, \code{yday}, and the climate variable used in \code{formula_model}.}

\item{Results_days}{Data frame. Output from a daily sensitivity analysis (e.g., CSP), must contain \code{estimate} (slope) and \code{r.squared} columns.}

\item{bio_data}{Data frame of biological observations, including \code{year}, the response variable, and site identifiers such as \code{plotname.lon.lat}.}

\item{yearneed}{Integer. Number of years needed prior to each observation to compute rolling climate data. Default is 2.}

\item{formula_model}{Formula. The model formula to be used for evaluation (e.g., \code{log.seed ~ TMEAN}).}

\item{model_type}{Character. Type of model to fit. One of \code{"lm"} (linear regression) or \code{"betareg"} (beta regression). Default is \code{"lm"}.}
}
\value{
A data frame summarizing model fit statistics for each identified cue window, including:
\itemize{
\item \code{window.open}, \code{window.close}
\item Model coefficients and significance, and other metrics like R2 and AIC
\item Site identifiers
}
}
\description{
This function identifies influential weather cue windows based on a peak signal detection algorithm applied
to the product of regression slope and R² values over time. It detects peaks where signal strength significantly
deviates from the local mean using a z-score thresholding method, then fits models (e.g., linear regression) over
these periods to evaluate biological-climate relationships.
}
\details{
The function calculates a signal vector using a rolling z-score approach on the product of slope and R² from daily regressions.
Significant signals are grouped into windows, which are then used to aggregate historical climate data and fit predictive models.
Model performance metrics are extracted using \code{reruning_windows_modelling()}.
}
\examples{
\dontrun{
result <- runing_peak_detection(
  lag = 30,
  threshold = 2.5,
  influence = 0.1,
  climate_data = climate_df,
  Results_days = csp_results,
  bio_data = seed_data,
  formula_model = formula(log.seed ~ TMEAN)
)
}

}
\seealso{
\code{\link{Thresholding_algorithm}}, \code{\link{reruning_windows_modelling}}, \code{\link{runing_csp}}, \code{\link{reformat_climate_backtothepast}}
}
