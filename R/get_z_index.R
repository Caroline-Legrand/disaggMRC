#' get_z_index
#' 
#' @description A function that calculates the asymmetry \eqn{z} index for each precipitation amount
#' 
#' @param x A vector of precipitation amounts, in mm
#' 
#' @return A vector of the same length as x, indicating the asymmetry \eqn{z} index. 
#' It is calculated as: \deqn{z_t =\frac{x_{t-1}+0.5 x_t}{x_{t-1}+x_t+x_{t+1}}}
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_z_index = function(x){
  # The length of the precipitation vector
  length_x = length(x)
  
  # R_{t-1}
  x_t_minus = x[1:(length_x-2)]
  
  # R_{t}
  x_t = x[2:(length_x-1)]
  
  # R_{t+1}
  x_t_plus = x[3:length_x]
  
  # Calculate the asymmetry z index for each precipitation amount
  z = (x_t_minus+0.5*x_t)/(x_t_minus+x_t+x_t_plus)
  
  # Add the first and last element of the vector
  z = c(0.5, z, 0.5)
  
  # If nearby values are unknown, assign 0.5. This happens in the case of missing data at times t-1 and t+1
  z[!is.na(x) & is.na(z)] = 0.5
  z[x == 0] = NA
  
  # Keep NA
  z[is.na(x)] = NA
  
  return(z)
}