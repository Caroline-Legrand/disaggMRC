#' get_phi_z
#' 
#' @description Probability asymmetry ratio model.
#' It returns the ratio \eqn{\varphi = \frac{p_{01}}{p_{10}+p_{01}}} given the asymmetry \eqn{z} index and the model parameter \eqn{\nu}
#' 
#' @param z A vector of \eqn{z} index  
#' @param nu Model parameter
#' 
#' @return A vector of the same length as \eqn{z}
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_phi_z = function(z, nu){
  E = pracma::erf((0.5-z)/(sqrt(2)*nu))
  return(0.5+0.5*E)
}