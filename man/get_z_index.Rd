% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_z_index.R
\name{get_z_index}
\alias{get_z_index}
\title{get_z_index}
\usage{
get_z_index(x)
}
\arguments{
\item{x}{A vector of precipitation amounts, in mm}
}
\value{
A vector of the same length as x, indicating the asymmetry \eqn{z} index. 
It is calculated as: \deqn{z_t =\frac{x_{t-1}+0.5 x_t}{x_{t-1}+x_t+x_{t+1}}}
}
\description{
A function that calculates the asymmetry \eqn{z} index for each precipitation amount
}
\references{
Maloku, K., Hingray, B., and Evin, G. (2023).
Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
}
\author{
Kaltrina Maloku
}
