#' get_Px_Intensity
#' 
#' @description Intermittency model. 
#' It returns the probability to divide the precipitation intensity given the model parameters \eqn{\mu} and \eqn{\sigma}
#' 
#' @param Intensity A numerical value, the precipitation intensity, in mm/h
#' @param mu,Sigma Model parameters
#' 
#' @return A numerical value, the probability to divide the precipitation intensity
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_Px_Intensity = function(Intensity, mu, Sigma){
  E = pracma::erf((log(Intensity)-mu)/(sqrt(2)*Sigma))
  return(0.5+0.5*E)
}