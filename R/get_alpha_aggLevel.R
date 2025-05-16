#' get_alpha_aggLevel
#' 
#' @description Alpha model as a function of the temporal aggregation level. 
#' It returns the value of the parameter \eqn{\alpha} of the Beta distribution given the model parameters \eqn{\alpha_{0}} and \eqn{H}
#' 
#' @param res_aggLevel A numerical value, the current temporal aggregation level, in minutes
#' @param alpha0,H Model parameters
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
get_alpha_aggLevel = function(alpha0, H, res_aggLevel){
  return(alpha0*res_aggLevel^(H))
}