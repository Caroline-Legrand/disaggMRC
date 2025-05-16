#' change_res_22.5_to_60
#' 
#' @description Convert a precipitation time series at a resolution of 22.5 minutes into a 60-minute time series, assuming a uniform distribution in the 22.5-minute interval
#' 
#' @param x A vector of numerical values of precipitation at a resolution of 22.5 minutes
#' 
#' @return A vector of numerical values of precipitation at a resolution of 60 minutes
#' 
#' @author Kaltrina Maloku
#' 
#' @export
change_res_22.5_to_60 = function(x){
  # Distributed uniformly from 22.5-minute resolution to 30-minute resolution
  y = change_res_22.5_to_30(x)
  
  # To 60-minute resolution
  y = RcppRoll::roll_sum(x = y, n = 2, by = 2, fill = 0, align = "left", na.rm = F)
  
  y = y[seq(1, to = length(y), by = 2)]
  
  return(y)
}