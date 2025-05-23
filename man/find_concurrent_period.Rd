% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_concurrent_period.R
\name{find_concurrent_period}
\alias{find_concurrent_period}
\title{Perform Climate Window Analysis - based on Climate sensitivity profil from Thackeray et al}
\usage{
find_concurrent_period(temp_window, pred_C)
}
\arguments{
\item{temp_window}{A numeric vector representing the temperature window or sequence of dates. This vector should be ordered.}

\item{pred_C}{A data frame or matrix containing the predicted coefficients. It must include a column named \code{fit} that contains the predicted values corresponding to the dates in \code{temp_window}.}
}
\value{
A numeric vector representing the temperature window where the average of the absolute predicted coefficients is the highest.
}
\description{
This function identifies periods within a given temperature window where the average of predicted coefficients is most extreme. The periods are determined based on gaps in sequential dates and are selected based on the maximum average of the absolute values of the predicted coefficients.
note that some function have been obtained or adjusted from Thackeray et al paper OR also from Simmonds et al.
}
\details{
The function first calculates the differences between consecutive dates in \code{temp_window} and identifies gaps greater than 1. It then segments the temperature window into periods based on these gaps. Each segment is assigned a unique identifier. The function then computes the mean absolute value of predicted coefficients for each segment and selects the period with the maximum average. This period is returned as the most extreme.
}
\examples{
# Example data
temp_window <- 1:10
pred_C <- data.frame(fit = c(2, 3, 1, 5, 7, 8, 6, 4, 9, 2))

# Find the concurrent period with the most extreme average
find_concurrent_period(temp_window, pred_C)

}
