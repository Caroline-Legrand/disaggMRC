#' get_alpha_intensity_1par
#' 
#' @description Model of the parameter \eqn{\alpha} of the Beta distribution as a function of precipitation intensity.
#' Its single parameter \eqn{K} explains the degree of dependency to intensity.
#' The model is considered to be constant for intensities smaller than \code{I_min} and higher than \code{I_max}
#' 
#' @param Intensity A numerical value, the precipitation intensity, in mm/h
#' @param I_min For intensities smaller than this value, the model is constant. Default 0.1 mm/h
#' @param I_max For intensities higher than this value, the model is constant. Default 10 mm/h
#' @param K Model parameter 
#' 
#' @return A numerical value, the parameter \eqn{\alpha} of the Beta distribution
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_alpha_intensity_1par = function(Intensity, I_min, I_max, K){
  I = log(Intensity)
  I1 = log(I_min)
  I2 = log(I_max)
  
  alpha = 1*(I<=I1) + (exp(K*(I-I1)^2))*(I>I1 & I<=I2) + (exp(K*(I2-I1)^2))*(I>I2)
  
  return(alpha)
}