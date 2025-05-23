% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_scaling_params_Intensity_aggLevel.R
\name{get_scaling_params_Intensity_aggLevel}
\alias{get_scaling_params_Intensity_aggLevel}
\title{get_scaling_params_Intensity_aggLevel}
\usage{
get_scaling_params_Intensity_aggLevel(
  vecPrecip,
  vecDates,
  resVecPrecip,
  aggLevels,
  by_season,
  listOptions
)
}
\arguments{
\item{vecPrecip}{A vector of observed precipitations, in mm}

\item{vecDates}{The corresponding vector of dates}

\item{resVecPrecip}{Resolution of the time series, in minutes}

\item{aggLevels}{A vector of temporal aggregation levels in the cascade estimation procedure, in minutes}

\item{by_season}{Logical. Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis}

\item{listOptions}{List with the following fields:
\itemize{
  \item \strong{I_min_fix}: A numerical value, in mm/h. For intensities above this value, the scaling model of alpha is constant. Recommended 0.01 mm/h. One might consider values between 0.001 and 0.1 mm/h. It corresponds to the I_zero of Eq. (2.13) Kaltrina Maloku's thesis
  \item \strong{I_max_fix}: A numerical value, in mm/h. For intensities below this value, the scaling model of alpha is constant. Recommended 7 mm/h. One might consider values between 5 and 10 mm/h. It corresponds to the I_one of Eq. (2.13) Kaltrina Maloku's thesis
  \item \strong{I_start_class}: A numerical value, in mm/h. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1 mm/h
  \item \strong{threshold_int}: A numerical value, in mm/h. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1 mm/h
}}
}
\value{
A list, a data frame of scaling parameters and a list of plots for showing empirical estimates and fitted scaling models
}
\description{
A function that estimates the scaling parameters of a MRC model, where the cascade weights are considered to depend on the temporal aggregation level and the intensity of the precipitation to be disaggregated.
In the first step, the observed weights are calculated for each temporal aggregation level. Also at this step, the \eqn{z} index is calculated for each observed precipitation intensity.
In the second step, the MRC parameters are estimated by considering a dependency on the temporal aggregation level, the precipitation intensity class and the asymmetry \eqn{z} index class.
In the last step, the scaling models are fitted to the MRC parameters estimated by the linear least squares method.
Plots showing the estimated parameters and their scaling models are also returned
}
\references{
Maloku, K., Hingray, B., and Evin, G. (2023).
Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
}
\author{
Kaltrina Maloku
}
