% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_Px_Intensity.R
\name{get_Px_Intensity}
\alias{get_Px_Intensity}
\title{get_Px_Intensity}
\usage{
get_Px_Intensity(Intensity, mu, Sigma)
}
\arguments{
\item{Intensity}{A numerical value, the precipitation intensity, in mm/h}

\item{mu, Sigma}{Model parameters}
}
\value{
A numerical value, the probability to divide the precipitation intensity
}
\description{
Intermittency model. 
It returns the probability to divide the precipitation intensity given the model parameters \eqn{\mu} and \eqn{\sigma}
}
\references{
Maloku, K., Hingray, B., and Evin, G. (2023).
Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
}
\author{
Kaltrina Maloku
}
