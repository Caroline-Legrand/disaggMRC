#' prWet
#' 
#' @description A function that calculates the proportion of wet steps in a time series
#' 
#' @param x A vector of precipitation, in mm
#' 
#' @return A numerical value, the proportion of wet steps
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
prWet = function(x){
  x = na.omit(x)
  pr = mean(x > 0)
  return(pr)
}