#' get_mean_z
#' 
#' @description Model of the mean of the cascade weights.
#' It returns the modelled mean of the distribution given the asymmetry \eqn{z} index and the model parameter \eqn{\lambda}
#'
#' @param z A vector of \eqn{z} index  
#' @param lambda Model parameter
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
get_mean_z = function(z, lambda){
  return(lambda*(z-0.5)+0.5)
}