#' estimate_params_MRC
#' 
#' @description A function to get the MRC parameters from the observed weights
#' 
#' @param vecPrecip A vector of observed precipitation, in mm
#' @param vecDates The corresponding vector of dates
#' @param resVecPrecip Resolution of the time series, in minutes
#' @param aggLevels A vector of temporal aggregation levels in the cascade estimation procedure, in minutes
#' @param by_season Logical. Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
#' @param I_start_class A numerical value, in mm/h. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1
#' @param threshold_int A numerical value, in mm/h. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1
#' 
#' @return A list of data frames with different MRC parameters estimated from observed precipitation data
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
estimate_params_MRC = function(vecPrecip,
                               vecDates,
                               resVecPrecip,
                               aggLevels,
                               by_season,
                               I_start_class,
                               threshold_int){

  # ======================== Preliminary checks ========================
  
  # Preliminary checks on vector types and lengths
  if(length(vecPrecip) != length(vecDates)) stop("The length of the date and the precipitation vectors must be identical")
  if(!any(c("POSIXlt","POSIXct","Date") %in% class(vecDates))) stop("vecDates must be a vector of class date. Accepted formats are POSIXlt, POSIXct or Date")
  if(!class(vecPrecip) %in% c("numeric","integer")) stop("vecPrecip must be a numeric vector")
  
  # Resolution of input time series for parameter estimation 
  if(is.null(resVecPrecip)){ # If the resolution is not specified, find it in the date vector
    resVecPrecip = difftime(vecDates[2], vecDates[1], units = "mins")
    resVecPrecip = as.numeric(resVecPrecip)
  } else { # If not, check that it matches what has been declared
    if(resVecPrecip != difftime(vecDates[2],vecDates[1], units = "mins"))
      stop("The temporal resolution of the time series does not match what was declared")
  }
  
  # If it is a continuous time series with a resolution of 10 minutes, remove the first 80 minutes and the last 80 minutes of each day
  diff_time = difftime(vecDates[2:length(vecDates)], vecDates[1:(length(vecDates)-1)])
  diff_time = as.numeric(diff_time)
  
  if(all(diff_time == resVecPrecip) & resVecPrecip == 10){
    # Number of days
    nb_days = round(length(vecDates)/(24*60/resVecPrecip))
    # Get the indices for the 1280 minutes at the center of the 1440-minute days   
    indices_1280min = get_indices_center_1280(nb_days)
    vecDates = vecDates[indices_1280min]
    vecPrecip = vecPrecip[indices_1280min]
  } 
  
  # Check if the aggregation levels are correct
  if(!all(round(aggLevels/(2*resVecPrecip)) == aggLevels/(2*resVecPrecip))) stop("aggLevels must be multiples of 2*resVecPrecip")
  
  # Check if parameter estimation is done on a monthly or seasonal basis 
  if(!is.logical(by_season)){
    stop("by_season argument must be T (for seasonal estimation) or F (for monthly estimation)")
  }

  # Define if the estimation should be done by month or by season
  if(by_season){
    vecSeason = month2season(lubridate::month(vecDates))
  } else {
    vecSeason = lubridate::month(vecDates)
  }
  
  # ======================== Start aggregating ========================

  DF_weights_obs = NULL
  
  for(i_aggLevel in aggLevels){
    # Fixed sum of the previous aggregation level
    vecPrecip_aggLevel_previous = sum_fixed_window(x = vecPrecip, k = 0.5*i_aggLevel/resVecPrecip)
    
    # Fixed sum of the current aggregation level
    vecPrecip_aggLevel_current = sum_fixed_window(x = vecPrecip_aggLevel_previous, k = 2)
    n_current = length(vecPrecip_aggLevel_current)
    
    # Precipitation amount over the first half
    vecPrecip_aggLevel_previous_1st_half = vecPrecip_aggLevel_previous[seq(from = 1, by = 2, length.out = n_current)]
    
    # Observed weights
    weights_current = vecPrecip_aggLevel_previous_1st_half/vecPrecip_aggLevel_current
    
    # Vector of season for the current level
    vecSeason_current = vecSeason[seq(from = 1, by = i_aggLevel/resVecPrecip, length.out = n_current)]
    
    # Find z index at the current aggregation level
    vec_z_index = get_z_index(vecPrecip_aggLevel_current)
    
    # Find z class id
    vec_z_index_id = get_z_class_id(vec_z_index, class_nb = 10)
    
    # Find intensity class id at the current aggregation level
    vec_intensity_id = get_intensity_class_id(x = vecPrecip_aggLevel_current, I_min = I_start_class, res_aggLevel = i_aggLevel)
    
    # Create a data frame that holds all the necessary information
    DF_aggLevel_current = data.frame(Season = vecSeason_current,
                                     Precip = vecPrecip_aggLevel_current,
                                     Weight = weights_current,
                                     Intensity = vec_intensity_id,
                                     zIndex = vec_z_index_id,
                                     aggLevel = i_aggLevel) %>% filter(!is.na(Weight) & Precip > 0) %>% dplyr::select(-"Precip")
    
    # Merge with previous aggregation levels
    DF_weights_obs = rbind(DF_weights_obs, DF_aggLevel_current)
  }
  
  # If a vector Season element is NA, remove it
  DF_weights_obs = DF_weights_obs %>% dplyr::filter(!is.na(Season))
  
  # Remove the weights calculated for each first class of the aggregation level
  DF_weights_obs = left_join(DF_weights_obs, DF_weights_obs %>% group_by(aggLevel) %>% 
                             summarise_at("Intensity", function(x) sort(unique(x))[2]) %>% # Find the minimum intensity for each aggregation level
                             rename(Intensity_min = Intensity)) %>% # Assign it as minimum intensity 
                   dplyr::filter(Intensity > Intensity_min) # Select only the weights calculated on the second and higher classes for each aggregation level
  
  # Disregard weights calculated below a given threshold
  DF_weights_obs = DF_weights_obs %>% dplyr::filter(Intensity > threshold_int)
  
  # Select only positive weights and less than 1
  DF_weights_obs_pos = DF_weights_obs %>% dplyr::filter(Weight != 0 & Weight != 1)
  
  # ========================= Alpha parameter =========================
  
  # Fit a symmetric Beta distribution to the observed weights
  # depending on the season
  # depending on the temporal aggregation level
  alpha_estimates_by_aggLevel = DF_weights_obs_pos %>%
    group_by(Season, aggLevel) %>%
    dplyr::summarise_at("Weight", function(x) fit_symm_beta_dist(x)) %>% stats::na.omit() %>% 
    mutate(param = "alpha", zIndex = NA, Intensity = NA) %>% dplyr::rename(value = Weight) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # Fit a symmetric Beta distribution to the observed weights
  # depending on the season 
  # depending on the temporal aggregation level
  # depending on the intensity class
  alpha_estimates_by_aggLevel_Intensity = DF_weights_obs_pos %>% 
    group_by(Season, aggLevel, Intensity) %>% 
    dplyr::summarise_at("Weight", function(x) fit_symm_beta_dist(x)) %>% stats::na.omit() %>% 
    mutate(param = "alpha", zIndex = NA) %>% dplyr::rename(value = Weight) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # Divide each estimated alpha for the temporal aggregation level and intensity class by the estimated alpha for the temporal scale when intensity classes are ignored
  alpha_star_estimates_by_aggLevel_Intensity = dplyr::left_join(alpha_estimates_by_aggLevel_Intensity, alpha_estimates_by_aggLevel %>%
    dplyr::rename(alpha_t = value) %>% 
    dplyr::select(Season, aggLevel, alpha_t), by = c("Season", "aggLevel")) %>%
    mutate(value = value/alpha_t) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # ===================== Intermittency parameters =====================
  
  # Estimate Px, the probability to divide precipitation into two sub-steps 
  # by season, intensity and temporal aggregation level
  px_estimates_by_aggLevel_Intensity = DF_weights_obs %>% 
    group_by(Season, aggLevel, Intensity) %>% 
    dplyr::summarise_at("Weight", function(x) get_portion_of_x(x)) %>% stats::na.omit() %>% 
    mutate(param = "px", zIndex = NA) %>% dplyr::rename(value = Weight) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)
  
  # ======================= Asymmetry parameters =======================
  
  # Estimate the mean of the positive weight distribution, 
  # for the asymmetry model of weight distribution, 
  # by season and z index class
  mean_estimates_by_z_index = DF_weights_obs_pos %>% 
    group_by(Season, zIndex) %>% 
    summarise_at("Weight", function(x) mean(x, na.rm = T)) %>% 
    rename(value = Weight) %>%
    mutate(param = "emp_mean", aggLevel = NA, Intensity = NA) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex) 
  
  # Probability that all precipitation is distributed to the second half
  # by season and z index class
  p01_estimates_by_z_index = DF_weights_obs %>% 
    group_by(Season, zIndex) %>% 
    summarise_at("Weight", function(x) get_portion_of_0(x)) %>%
    rename(p01 = Weight)
  
  # Probability that all precipitation is distributed to the first half
  # by season and z index class
  p10_estimates_by_z_index = DF_weights_obs %>% 
    group_by(Season, zIndex) %>% 
    summarise_at("Weight", function(x) get_portion_of_1(x)) %>%
    rename(p10 = Weight)

  # Estimate phi, the asymmetry parameter for the intermittency model
  # by season and z index class
  phi_estimates_by_z_index = left_join(p01_estimates_by_z_index, p10_estimates_by_z_index, by = c("Season", "zIndex")) %>%
    mutate(value = p01/(p01+p10), param = "phi", aggLevel = NA, Intensity = NA) %>%
    dplyr::select(Season, param, value, aggLevel, Intensity, zIndex)

  # List of all estimated parameters 
  list_params_MRC = list(alpha_aggLevel = alpha_estimates_by_aggLevel,
                         alpha_aggLevel_Intensity = alpha_estimates_by_aggLevel_Intensity,
                         alpha_star_aggLevel_Intensity = alpha_star_estimates_by_aggLevel_Intensity,
                         px_aggLevel_Intensity = px_estimates_by_aggLevel_Intensity,
                         emp_mean_z = mean_estimates_by_z_index,
                         phi_z = phi_estimates_by_z_index)  
  
  return(list_params_MRC)
}