#' get_indices_center_1280
#' 
#' @description A function that extracts "days" of 1280 minutes from a series of days of 1440 minutes
#' 
#' @param nb_days A numerical value, the number of days in the time series
#' 
#' @return A vector of indices
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_indices_center_1280 = function(nb_days){
  i = 1:nb_days
  ind_1280 = as.vector(sapply(i, function(x){
    x = (144*(x-1)+8+1):(144*(x)-8)
  }))
  return(ind_1280)
}