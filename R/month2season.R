#' month2season
#'
#' @description A function that transforms a vector of months to seasons
#'
#' @param vecMonth A vector of months given as integers 1:12
#' 
#' @author Guillaume Evin
#'
#' @export
month2season = function(vecMonth){
  iSeason = c(1,1,2,2,2,3,3,3,4,4,4,1)
  return(iSeason[vecMonth])
}