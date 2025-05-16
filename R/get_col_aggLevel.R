#' get_col_aggLevel
#' 
#' @description A function for plots
#' 
#' @param aggLevel A vector providing the temporal aggregation levels, in minutes
#' 
#' @return A vector of colors, of the same length as \code{aggLevel}
#' 
#' @author Kaltrina Maloku
#' 
#' @export
get_col_aggLevel = function(aggLevel){
  col_aggLevel = viridis::viridis(length(aggLevel), direction = -1)
  return(col_aggLevel)
}