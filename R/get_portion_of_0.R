#' get_portion_of_0
#' 
#' @description A function that finds the portion of the cascade weights that is equal to 0 
#' 
#' @param x A vector of weights
#' 
#' @return The portion of weights equal to 0
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
get_portion_of_0 = function(x){
  if(length(x) < 10){
    pr_0 = NA
  } else {
    x = x[!is.na(x)]
    pr_0 = mean(x == 0)
  }
  return(pr_0)
}