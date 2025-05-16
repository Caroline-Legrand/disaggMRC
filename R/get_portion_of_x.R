#' get_portion_of_x
#' 
#' @description A function that finds the portion of the cascade weights that is strictly greater than 0 and strictly smaller than 1
#' 
#' @param x A vector of weights
#' 
#' @return The portion of weights that is strictly greater than 0 and strictly smaller than 1
#' 
#' @details If sample size is less than 10, return NA
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_portion_of_x = function(x){
  if(length(x) < 10){
    pr_x = NA
  } else {
    x = x[!is.na(x)]
    pr_x = mean(x > 0 & x < 1)
  }
  return(pr_x)
}