% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_params_MRC.R
\name{estimate_params_MRC}
\alias{estimate_params_MRC}
\title{estimate_params_MRC}
\usage{
estimate_params_MRC(
  vecPrecip,
  vecDates,
  resVecPrecip,
  aggLevels,
  by_season,
  I_start_class,
  threshold_int
)
}
\arguments{
\item{vecPrecip}{A vector of observed precipitation, in mm}

\item{vecDates}{The corresponding vector of dates}

\item{resVecPrecip}{Resolution of the time series, in minutes}

\item{aggLevels}{A vector of temporal aggregation levels in the cascade estimation procedure, in minutes}

\item{by_season}{Logical. Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis}

\item{I_start_class}{A numerical value, in mm/h. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1}

\item{threshold_int}{A numerical value, in mm/h. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1}
}
\value{
A list of data frames with different MRC parameters estimated from observed precipitation data
}
\description{
A function to get the MRC parameters from the observed weights
}
\references{
Maloku, K., Hingray, B., and Evin, G. (2023).
Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
}
\author{
Kaltrina Maloku
}
