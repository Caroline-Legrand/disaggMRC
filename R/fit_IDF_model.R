#' fit_IDF_model
#' 
#' @param station_data A matrix of precipitation data aggregated over different durations corresponding to one column
#' @param durations Durations in hours (for hourly data)
#' @param declustering_duration Declustering to reduce serial dependence description
# 
#' @return A list of parameters
#' 
#' @author Abubakar Haruna
#' 
#' @export
fit_IDF_model = function(station_data, durations, declustering_duration){
  # Local fit - EGPD fitting for each duration separately to get initial values for the IDF function
  initial_params = egpdIDF::egpd_idf_init(station_data = station_data,
                                          durations = durations,
                                          left_censoring_value = 0.5/durations,
                                          declustering_duration = declustering_duration,
                                          use_r_optim = T, 
                                          nrmse_tol = 0.05, 
                                          max_xi = 0.25)
  
  # IDF model fit
  fitted_idf = egpdIDF::fit_egpd_idf_data_driven(station_data = station_data,
                                                 durations = durations, 
                                                 left_censoring_value = initial_params$fits$lower_threshold,
                                                 declustering_duration = declustering_duration, 
                                                 initial_params = initial_params,
                                                 use_profile_likelihood = T, 
                                                 max_xi = 0.25, 
                                                 optim_algo = "Nelder-Mead")
  
  return(fitted_idf)
}