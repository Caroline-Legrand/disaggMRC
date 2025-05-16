#' change_res_11.25_to_30
#' 
#' @description Convert a precipitation time series at a resolution of 11.25 minutes into a 30-minute time series, assuming a uniform distribution in the 11.25-minute interval
#' 
#' @param x A vector of numerical values of precipitation at a resolution of 11.25 minutes
#' 
#' @return A vector of numerical values of precipitation at a resolution of 30 minutes
#' 
#' @author Kaltrina Maloku
#' 
#' @export
change_res_11.25_to_30 = function(x){
  # The length of the new time series
  n.y = length(x)*11.25/30
  
  # Define the vector that will contain the time series 
  y = rep(NA, n.y)
  
  # Aggregate to 22.5 minutes
  x_22.5_left = RcppRoll::roll_sum(x = x, n = 2, by = 2, fill = 0, align = "left", na.rm = F)
  x_22.5_left = x_22.5_left[seq(1, to = length(x_22.5_left), by = 2)]
  
  # Aggregate to 22.5 minutes
  x_22.5_right = RcppRoll::roll_sum(x = x[-1], n = 2, by = 2, fill = 0, align = "left", na.rm = F)
  x_22.5_right = x_22.5_right[seq(1, to = length(x_22.5_right), by = 2)]
  
  # Assign precipitation on the first 30 minutes
  y[seq(1, by = 3, length.out = n.y/3)] = x_22.5_left[seq(from = 1, by = 4, length.out = n.y/3)] + (7.5/11.25)*x[seq(from = 3, by = 8, length.out = n.y/3)]
  
  # Second 30 minutes
  y[seq(2, by = 3, length.out = n.y/3)] = (1-7.5/11.25)*x[seq(from = 3, by = 8, length.out = n.y/3)] + x_22.5_right[seq(from = 2, by = 4, length.out = n.y/3)] + (1-7.5/11.25)*x[seq(from = 6, by = 8, length.out = n.y/3)]
  
  # Third 30 minutes
  y[seq(3, by = 3, length.out = n.y/3)] = (7.5/11.25)*x[seq(from = 6, by = 8, length.out = n.y/3)] + x_22.5_left[seq(from = 4, by = 4, length.out = n.y/3)]
  
  y[is.na(y)] = 0
  
  return(y)
}