#' sum_fixed_window
#' 
#' @description A function that calculates a sum over a fixed window
#' 
#' @param x A vector of precipitation amounts, in mm
#' @param k Length of the fixed window
#' 
#' @return A vector of sums over the window of length k, in mm
#' 
#' @details Sums are assigned to the leftmost element in the window
#' 
#' @author Kaltrina Maloku
#' 
#' @export
sum_fixed_window = function(x, k){
  # Check if the vector and window size are correct
  if(!class(x) %in% c("numeric", "integer")) stop("x must be a numeric vector")
  if(length(x) < k) stop("The length of the fixed window must be smaller than that of the precipitation vector")
  
  # Fixed sums
  x_sum = RcppRoll::roll_sum(x = x, n = k, by = k, fill = NULL, align = "left", na.rm = F) 
  
  return(x_sum)
}