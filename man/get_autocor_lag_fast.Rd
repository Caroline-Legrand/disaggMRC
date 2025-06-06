% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_autocor_lag_fast.R
\name{get_autocor_lag_fast}
\alias{get_autocor_lag_fast}
\title{get_autocor_lag_fast}
\usage{
get_autocor_lag_fast(x, k, na.rm = T, all = T)
}
\arguments{
\item{x}{A vector of precipitation, in mm}

\item{k}{A numerical value, the lag at which to calculate autocorrelation}

\item{na.rm}{Logical. If T, missing values (NA) will be ignored before calculation}

\item{all}{Logical. If F, the lag-k autocorrelation will be calculated from positive precipitation only}
}
\value{
A numerical value, the lag-k autocorrelation
}
\description{
A function that calculates the lag-k autocorrelation of a precipitation time series
}
\references{
Maloku, K., Hingray, B., and Evin, G. (2023).
Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
}
\author{
Kaltrina Maloku
}
