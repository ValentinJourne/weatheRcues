% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runing_climwin.R
\name{runing_climwin}
\alias{runing_climwin}
\title{Run Climate Window Analysis Using Sliding Window Approach}
\usage{
runing_climwin(
  ...,
  climate_data,
  bio_data,
  site.name,
  range = c(600, 0),
  cinterval = "day",
  refday = c(1, 11),
  optionwindows = "absolute",
  climate_var = "TMEAN",
  stat.aggregate = "mean",
  formulanull = stats::formula("log.seed ~ 1"),
  fun = "lin",
  cmissing = "method2",
  give.clean = TRUE
)
}
\arguments{
\item{climate_data}{A data frame containing daily climate data. Must include a \code{date} column (class \code{Date}) and one column for the climate variable (e.g., \code{TMEAN}).}

\item{bio_data}{A data frame with biological data. Must include a \code{Date2} column (class \code{Date}) and a response variable that matches the left-hand side of \code{formulanull}.}

\item{site.name}{A character string used to label the output, e.g., the site or population name.}

\item{range}{A numeric vector of length 2 indicating the search window, e.g., \code{c(600, 0)} to test up to 600 days before the reference date.}

\item{cinterval}{Time unit used for aggregating the climate variable: \code{"day"} (default), \code{"week"}, etc.}

\item{refday}{A vector specifying the day and month of the reference event (e.g., \code{c(1, 11)} for November 1).}

\item{optionwindows}{One of \code{"absolute"} (fixed calendar date) or \code{"relative"} (relative to biological event). Default is \code{"absolute"}.}

\item{climate_var}{A string specifying the name of the climate variable column (e.g., \code{"TMEAN"}). Must exist in \code{climate_data}.}

\item{stat.aggregate}{Aggregation function to apply over the window: \code{"mean"} (default), \code{"sum"}, \code{"min"}, or \code{"max"}.}

\item{formulanull}{A formula (e.g., \code{log.seed ~ 1}) specifying the null model to which the climate models are compared.}

\item{fun}{Functional form of the climate effect. Use \code{"lin"} (default) for linear, \code{"quad"} for quadratic, etc.}

\item{cmissing}{Method to handle missing values in climate data. See \code{climwin::slidingwin} for options. Default is \code{"method2"}.}

\item{give.clean}{Logical. If \code{TRUE} (default), returns a cleaned summary of the best model results. If \code{FALSE}, returns the full climwin object.}
}
\value{
If \code{give.clean = TRUE}, a data frame with:
\itemize{
\item \code{sitenewname}: The site identifier
\item \code{climate.file}: Name of the climate file or site
\item \code{window.open}, \code{window.close}: The best-fit window bounds (in days before reference)
\item \code{r2}, \code{AIC}, \code{sigma}: Model performance metrics
\item \code{slope.estimate}, \code{intercept.estimate}: Coefficients from the best model
}
If \code{give.clean = FALSE}, the full list returned by \code{climwin::slidingwin()} is returned instead.
}
\description{
This function wraps the \code{climwin::slidingwin()} function to detect climate cue windows that best explain biological responses
(e.g., seed production) using a sliding window approach. The function fits models across various climate time windows
relative to a specified reference date, identifies the best-fitting window based on model performance, and optionally summarizes results.
}
\details{
The function evaluates the effect of climate on a biological variable by systematically testing a series of time windows
before a reference date. It compares the performance of models fitted on climate summaries from each window against a null model.
The \code{climwin} package is used internally.

This wrapper adds input checking, flexible formula input, and a standardized output format compatible with batch analyses across multiple sites.
}
\examples{
\dontrun{
runing_climwin(
  climate_data = daily_temp,
  bio_data = seed_data,
  site.name = "MySite",
  range = c(365, 0),
  refday = c(1, 11),
  climate_var = "TMEAN",
  formulanull = log.seed ~ 1
)
}

}
\seealso{
\code{\link[climwin]{slidingwin}}
}
