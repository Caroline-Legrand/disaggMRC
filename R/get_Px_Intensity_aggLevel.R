#' get_Px_Intensity_aggLevel
#' 
#' @description Intermittency model as a function of the temporal aggregation level.
#' It returns the probability to divide the precipitation intensity given the model parameters \eqn{a_\mu, b_\mu, a_\sigma} and \eqn{b_\sigma}
#' 
#' @param Intensity A numerical value, the precipitation intensity, in mm/h
#' @param res_aggLevel A numerical value, the current aggregation level, in minutes
#' @param a_mu,b_mu,a_sigma,b_sigma Model parameters
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
get_Px_Intensity_aggLevel = function(Intensity, res_aggLevel, a_mu, b_mu, a_sigma, b_sigma){
  mu = a_mu*log(res_aggLevel)+b_mu
  Sigma = a_sigma*log(res_aggLevel)+b_sigma
  E = pracma::erf((log(Intensity)-mu)/(sqrt(2)*Sigma))
  return(0.5+0.5*E)
}