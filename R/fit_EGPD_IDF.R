#' fit_EGPD_IDF
#' 
#' @description A function that estimates the parameters of the EGPD
#' 
#' @param vec.precip A vector of observed hourly precipitation
#' @param vec.dates The corresponding vector of dates
#' @param durations Durations in hours (for hourly data)
#' @param declustering_duration Declustering to reduce serial dependence description
#' 
#' @return A list containing the EGPD parameters by month. Each element of the list is a matrix of size nStation x 3
#' 
#' @author Kaltrina Maloku
#' 
#' @export
fit_EGPD_IDF = function(vec.precip, 
                        vec.dates, 
                        durations = c(1,2,3,6,10,12,16,18,24,48,72),
                        declustering_duration = c(1,2,3,6,10,12,16,18,24,48,72)){
  
  # Vector of integers indicating the months (1,1,1,...)
  vec.month = as.numeric(strftime(vec.dates, "%m"))
  
  # Prepare a list of parameters
  parMargin = list()
  
  for(iMonth in 1:12){
    # Dates concerned in the month
    isMonth = vec.month==iMonth
    
    # Aggregate data to specified durations
    station_data = egpdIDF::aggregate_data(sample_data = data.frame(date = vec.dates[isMonth], Prec = vec.precip[isMonth]), st_code = "Prec", durations)
    
    # IDF fit
    fitted_idf = fit_IDF_model(station_data = station_data,
                               durations = durations,
                               declustering_duration = declustering_duration)
    
    # Get only parameters corresponding to 24-hour duration
    # Multiply the scale by the duration of interest to get the value in mm per day
    kappa.month = fitted_idf$kappa_param[which(durations == 24)]
    scale.month = 24*fitted_idf$scale_param[which(durations == 24)]
    xi.month = fitted_idf$shape_param[which(durations == 24)]
    
    parMargin[[iMonth]] = matrix(c(kappa.month, scale.month, xi.month), ncol = 3, dimnames = list(NULL, c("kappa", "sigma", "xi")))
  }
  
  # Return the list of parameters 
  return(parMargin)
}