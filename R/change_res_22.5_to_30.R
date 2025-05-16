#' change_res_22.5_to_30
#' 
#' @description Convert a precipitation time series at a resolution of 22.5 minutes into a 30-minute time series, assuming a uniform distribution in the 22.5-minute interval
#' 
#' @param x A vector of numerical values of precipitation at a resolution of 22.5 minutes
#' 
#' @return A vector of numerical values of precipitation at a resolution of 30 minutes
#' 
#' @author Kaltrina Maloku
#' 
#' @export
change_res_22.5_to_30 = function(x){
  y = rep(NA, length(x)*22.5/30)
  
  y1 = seq(1, by = 3, length.out = length(y)/3)
  
  y[y1] = x[y1+0:(length(y1)-1)]+(7.5/22.5)*x[y1+1:length(y1)]
  
  y[y1+1] = (15/22.5)*(x[y1+1:length(y1)]+x[y1+2:(length(y1)+1)])
  
  y[y1+2] = (7.5/22.5)*x[y1+2:(length(y1)+1)]+x[y1+3:(length(y1)+2)]
  
  y[is.na(y)] = 0
  
  return(y)
}