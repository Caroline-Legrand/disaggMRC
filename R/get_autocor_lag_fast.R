#' get_autocor_lag_fast
#' 
#' @description A function that calculates the lag-k autocorrelation of a precipitation time series
#' 
#' @param x A vector of precipitation, in mm
#' @param k A numerical value, the lag at which to calculate autocorrelation
#' @param na.rm Logical. If T, missing values (NA) will be ignored before calculation
#' @param all Logical. If F, the lag-k autocorrelation will be calculated from positive precipitation only
#' 
#' @return A numerical value, the lag-k autocorrelation
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_autocor_lag_fast = function(x, k, na.rm = T, all = T){
  n = length(x)
  if(!all){ # For autocorrelation of positive precipitation only
    x[x == 0] = NA
  }
  # Estimate mean and variance
  mu_x = mean(x, na.rm = na.rm)
  sigma_x = stats::sd(x, na.rm = na.rm)
  
  x1 = x[(k+1):n]
  x2 = x[1:(n-k)]
  
  # Estimate autocorrelation 
  rho_k = mean((x1-mu_x)*(x2-mu_x), na.rm = T)*((n-k)/n)/(sigma_x^2)
  
  return(rho_k)
}