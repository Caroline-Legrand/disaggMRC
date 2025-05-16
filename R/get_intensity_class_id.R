#' get_intensity_class_id
#' 
#' @description A function that, for a given precipitation amount, returns its intensity class id
#' 
#' @param x A vector of precipitation amounts, in mm
#' @param I_min A numerical value, the minimum precipitation observed, in mm/h. Default 0.001 mm/h for mean areal precipitation. For rain gauge data with a tipping volume of 0.1 mm, we suggest a value of 0.1*60/1440
#' @param res_aggLevel A numerical value, the current aggregation level, in minutes
#' 
#' @return A vector of the same length as x, indicating the intensity class id of each precipitation amount
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_intensity_class_id = function(x, I_min = 0.001, res_aggLevel){
  
  # Convert to mm/h 
  x = x*(60/res_aggLevel)
  
  # If the precipitation amount is 0, put NA 
  x[x <= 0] = NA
  
  # Find n such that I_min*2^(n-1) <= x < I_min*2^n
  n_int_class = floor(log(x/I_min)/log(2))+1
  
  # Set the maximal number of intensity classes 
  max_int_class = 30
  
  # Find the intensity value that represents each class at hourly resolution 
  int_bounds_intervals = I_min*2^(1:max_int_class-1)
  
  int_class_id = rep(NA, max_int_class-1)
  for (i in 1:(length(int_bounds_intervals)-1)){
    # See Rupp et al., 2009 for more details on the definition of classes
    int_class_id[i] = sqrt(int_bounds_intervals[i]*int_bounds_intervals[i+1])
  }
  
  # Find the intensity class id
  int_id = as.numeric(as.character(factor(n_int_class, levels = 1:(max_int_class-1), labels = int_class_id)))
  
  return(int_id)
}