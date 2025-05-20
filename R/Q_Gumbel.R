#' Q_Gumbel
#' 
#' @description A function that calculates the inverse of the cumulative distribution function (CDF) for a Gumbel distribution. It transforms a probability into a specific quantile of the Gumbel distribution
#' 
#' @param CDF A numerical value, the cumulative distribution function, i.e. the probability associated with a certain quantile
#' 
#' @return The numerical value associated with the Gumbel quantile for the given probability
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#'
#' @export
Q_Gumbel = function(CDF){
  return(-log(-log(CDF)))
}