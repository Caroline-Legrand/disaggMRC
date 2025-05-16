#' alpha_star_Intensity
#' 
#' @description Quadratic model of the scaled \eqn{\alpha} as a function of the precipitation intensity.
#' It returns the value of the scaled \eqn{\alpha} given the model parameters \eqn{c_0, c_1} and \eqn{c_2}
#' 
#' @param Intensity A numerical value, the precipitation intensity, in mm/h
#' @param c0,c1,c2 Model parameters
#' 
#' @return A numerical value, the value of the scaled \eqn{\alpha}
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
alpha_star_Intensity = function(Intensity, c0, c1, c2){
  return(exp(c0+c1*log(Intensity)+c2*log(Intensity)*log(Intensity)))
}