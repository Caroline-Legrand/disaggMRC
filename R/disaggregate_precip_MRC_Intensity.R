#' disaggregate_precip_MRC_Intensity
#' 
#' @description A function that disaggregates a coarse resolution time series into a fine resolution time series according to the MRC model, 
#' where the dependence on precipitation intensity is taken into account
#' 
#' @param vecPrecip_target A vector of observed precipitation at coarse resolution, in mm
#' @param vecDates_target The corresponding vector of dates
#' @param params_scaling A data frame of parameters needed for disaggregation
#' @param by_season Logical. Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
#' @param res_coarse_aggLevel Resolution of the time series to be disaggregated, in minutes
#' @param res_fine_aggLevel Target resolution of the disaggregated time series, in minutes
#' @param asymmetry_option Logical. Should the disaggregation be dependent on the asymmetry model?
#' 
#' @return A vector of high resolution precipitation
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
disaggregate_precip_MRC_Intensity = function(vecPrecip_target,
                                             vecDates_target,
                                             params_scaling,
                                             by_season,
                                             res_coarse_aggLevel,
                                             res_fine_aggLevel,
                                             asymmetry_option){
  
  # Check that the length of the date and precipitation vectors are the same
  if(length(vecPrecip_target) != length(vecDates_target))
    stop("The vector of dates must have the same length as the vector of precipitation")
  
  # Define if the estimation should be done by month or by season
  if(by_season){
    vecSeason_aggLevel = month2season(lubridate::month(vecDates_target))
  } else {
    vecSeason_aggLevel = lubridate::month(vecDates_target)
  }
  
  # Define the last level of aggregation before processing to achieve the target resolution
  if (res_fine_aggLevel == 30){
    res_stop_disagg = 11.25
  } else if (res_fine_aggLevel == 60){
    res_stop_disagg = 22.5
  } else {
    res_stop_disagg = res_fine_aggLevel
  }
  
  # Start from the coarsest level of aggregation
  i_aggLevel = res_coarse_aggLevel
  
  while (i_aggLevel > res_stop_disagg){
    # Get precipitation and seasonal vectors at the current temporal aggregation level
    if (i_aggLevel == res_coarse_aggLevel){
      # If the first level takes target days
      vecPrecip_aggLevel_current = vecPrecip_target
      vecSeason_aggLevel_current = vecSeason_aggLevel
    }
    
    # Length of current time series
    n_current = length(vecPrecip_aggLevel_current)
    
    # Prepare data frame
    # Convert precipitation amounts into precipitation intensities, in mm/h
    data_current = data.frame(Intensity = vecPrecip_aggLevel_current*(60/i_aggLevel),
                              zIndex = get_z_index(vecPrecip_aggLevel_current),
                              Season = vecSeason_aggLevel_current)
    
    # Add parameters to intensities and seasons, each season having its own parameters
    data_current_params = left_join(data_current, params_scaling, by = "Season")
    
    # Get MRC parameters from intensities and scaling parameters
    # rnd_unif is a random vector used for disaggregation
    data_current_params = data_current_params %>%
      mutate(Px = get_Px_Intensity(Intensity = Intensity, mu = mu, Sigma = Sigma),
             alpha = get_alpha_intensity_1par(Intensity = Intensity, I_min = 0.1, I_max = 10, K = K),
             phi = get_phi_z(z = zIndex, nu = nu),
             m = get_mean_z(z = zIndex, lambda = lambda),
             rnd_unif = runif(nrow(data_current_params)),
             bdc = NA)
    
    if (asymmetry_option){
      # If asymmetry option is true
      # First determine 0/1 weights
      data_current_params = data_current_params %>%
        mutate(P01 = phi*(1-Px), # Probability that all precipitation goes in the second half
               P10 = (1-phi)*(1-Px)) %>% # Probability that all precipitation goes in the first half
        mutate(zeros = 1*!(P01-rnd_unif >= 0), # Compare the probability to random deviates
               ones = 1*(1-P10-rnd_unif < 0)) %>% # And decide if 0, 1, or in between
        mutate(bdc = replace(bdc, zeros == 0,0)) %>%
        mutate(bdc = replace(bdc, ones == 1,1)) %>%
        dplyr::select(-c("zeros","ones","rnd_unif"))
      
      # From alpha and the mean as modelled by the asymmetry model, find alpha1 and alpha2
      data_current_params = data_current_params %>% 
        mutate(alpha1 = get_alpha1(alpha, m), alpha2 = get_alpha2(alpha, m))
      
      # Find indices where bdc is not yet defined, i.e. neither one nor zero
      indices_weights = is.na(data_current_params$bdc) & !is.na(data_current_params$Intensity) & !is.na(data_current_params$zIndex)
      
      # Draw random deviates from the Beta distribution only when the weight is not yet defined
      positive_weights = mapply(function(x, y) stats::rbeta(1, shape1 = x, shape2 = y),
                                data_current_params$alpha1[indices_weights],
                                data_current_params$alpha2[indices_weights])
      
    } else {
      # If asymmetry option is not true
      # First determine 0/1 weights
      data_current_params = data_current_params %>%
        mutate(zeros = 1*!(0.5*(1-Px)-rnd_unif >= 0),
               ones = 1*(0.5*(1+Px)-rnd_unif < 0)) %>%
        mutate(bdc = replace(bdc, zeros == 0,0)) %>%
        mutate(bdc = replace(bdc, ones == 1,1)) %>%
        dplyr::select(-c("zeros","ones","rnd_unif"))
      
      # Find indices where bdc is not yet defined, i.e. neither one nor zero
      indices_weights = is.na(data_current_params$bdc) & !is.na(data_current_params$Intensity) & !is.na(data_current_params$zIndex)
      
      # Draw random deviates from the Beta distribution only when the weight is not yet defined
      positive_weights = sapply(data_current_params$alpha[indices_weights], function(x) rbeta(1, shape1 = x, shape2 = x))
    }
    
    # Replace missing weights
    data_current_params = data_current_params %>%
      mutate(bdc = replace(bdc, is.na(bdc) & !is.na(zIndex) & !is.na(Intensity), positive_weights)) %>%
      mutate(bdc = replace(bdc, is.na(bdc), 0)) # Assign zero weight if there is no precipitation to disaggregate
    
    # Prepare a new precipitation vector with double the resolution
    vecPrecip_aggLevel_previous = rep(NA, 2*n_current)
    
    # The first half
    vecPrecip_aggLevel_previous[seq(1, by = 2, length.out = n_current)] = vecPrecip_aggLevel_current*data_current_params$bdc
    
    # To the second half
    vecPrecip_aggLevel_previous[seq(2, by = 2, length.out = n_current)] = vecPrecip_aggLevel_current*(1-data_current_params$bdc)
    
    # Update aggregation level
    i_aggLevel = i_aggLevel/2
    
    # Update precipitation and season vectors
    vecPrecip_aggLevel_current = vecPrecip_aggLevel_previous
    vecSeason_aggLevel_current = rep(vecSeason_aggLevel_current, each = 2)
  }
  
  if(res_fine_aggLevel == 30){ # If the resolution is 30 minutes, distribute uniformly to 30 min
    vecPrecip_aggLevel_current = change_res_11.25_to_30(vecPrecip_aggLevel_current)
  } else if (res_fine_aggLevel == 60){ # If the resolution is 60 minutes, distribute uniformly to 60 min 
    vecPrecip_aggLevel_current = change_res_22.5_to_60(vecPrecip_aggLevel_current)
  }
  
  # Return simulated scenario
  return(vecPrecip_aggLevel_current)
}