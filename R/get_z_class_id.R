#' get_z_class_id
#' 
#' @description A function that, for a given \eqn{z} index, returns the id of its class
#' 
#' @param z A vector of \eqn{z} index
#' @param class_nb Number of classes
#' 
#' @return A vector of the same length as \eqn{z}, indicating the class id of each \eqn{z} index
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_z_class_id = function(z, class_nb = 10){
  # Create 10 classes
  classes_intervals = seq(0, 1, length.out = class_nb+1)
  
  # Get the id of each class
  classes_id = seq(from = 0.5*classes_intervals[2], by = classes_intervals[2]-classes_intervals[1], length.out = class_nb)
  
  # Determine the class of each z index
  z_classes = findInterval(z, classes_intervals, left.open = T, all.inside = T)
  
  # Determine the class id for each z index
  z_classes_id = as.numeric(as.character(factor(z_classes, levels = 1:class_nb, labels = classes_id)))
  
  return(z_classes_id)
}