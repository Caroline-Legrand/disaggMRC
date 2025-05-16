#' get_alpha1
#' 
#' @description A function that estimates the parameter \eqn{\alpha_{1}} of the Beta distribution given the parameter \eqn{\alpha} of the symmetric Beta distribution and its mean
#' 
#' @param alpha A numerical value, the parameter \eqn{\alpha} of the symmetric Beta distribution
#' @param m A numerical value, the mean of the Beta distribution
#' 
#' @return A numerical value, the parameter \eqn{\alpha_{1}} of the Beta distribution
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#'
#' @export
get_alpha1 = function(alpha, m){
  # For the known alpha parameter of the symmetric Beta distribution, find the variance of the distribution
  var_weights = 1/(8*alpha+4)
  
  # The mean of the distribution
  mean_weights = m
  
  # The alpha1 parameter of the non-symmetric Beta distribution
  shape1 = mean_weights*((mean_weights*(1-mean_weights)/var_weights)-1)
  
  return(shape1)
}