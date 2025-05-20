#' get_standard_stats_extremes
#' 
#' @description A function that calculates statistics from precipitation time series scenarios
#' 
#' @param df_data A data frame containing observed and disaggregated precipitation time series scenarios at a certain temporal resolution
#' 
#' @return A data frame of statistics
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#'
#' @export
get_standard_stats_extremes = function(df_data){
  # Variable names for all scenarios
  names_scenarios = colnames(df_data)[which(grepl("result.", colnames(df_data)))]
  
  # Add obs to variable names
  agg_by_names = c("obs", names_scenarios)

  # One day maxima
  df_max = df_data %>% group_by(Season, Year) %>% 
    summarise_at(agg_by_names, function(x) if(all(is.na(x))){return(NA)} else {return(max(x, na.rm = T))})
  df_max[df_max == -Inf] = NA
  
  # 50-year return period
  rt_period = 50; p = 1-1/rt_period
  rt_level_50 = df_max  %>% group_by(Season) %>% 
    summarise_at(agg_by_names, function(x) qq.emp(x = x, p = p)) %>% na.omit 
  rt_level_50$stat = "rt_50"
  
  # 20-year return period
  rt_period = 20; p = 1-1/rt_period
  rt_level_20 = df_max  %>% group_by(Season) %>% 
    summarise_at(agg_by_names, function(x) qq.emp(x = x, p = p)) %>% na.omit
  rt_level_20$stat = "rt_20"
  
  # 10-year return period
  rt_period = 10; p = 1-1/rt_period
  rt_level_10 = df_max  %>% group_by(Season) %>% 
    summarise_at(agg_by_names, function(x) qq.emp(x = x, p = p)) %>% na.omit
  rt_level_10$stat = "rt_10"
  
  # 5-year return period
  rt_period = 5; p = 1-1/rt_period
  rt_level_5 = df_max  %>% group_by(Season) %>% 
    summarise_at(agg_by_names, function(x) qq.emp(x = x, p = p)) %>% na.omit
  rt_level_5$stat = "rt_5"
  
  # 2-year return period
  rt_period = 2; p = 1-1/rt_period
  rt_level_2 = df_max  %>% group_by(Season) %>% 
    summarise_at(agg_by_names, function(x) qq.emp(x = x, p = p)) %>% na.omit
  rt_level_2$stat = "rt_2"
  
  # Precipitation maxima
  max_precip = df_data %>% group_by(Season) %>% 
    summarise_at(agg_by_names, function(x) if(all(is.na(x))){return(NA)} else {return(max(x, na.rm = T))})
  max_precip[max_precip == -Inf] = NA
  max_precip$stat = "maxP"
  
  # Data frame containing return levels and precipitation maxima
  TAB_return_levels = rbind(max_precip, rt_level_2, rt_level_5, rt_level_10, rt_level_20, rt_level_50)
  
  # Standard statistics
  if(T){
    TAB = NULL
    
    # Mean
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) mean(x, na.rm = T)) %>%
      add_column(stat = "MEAN")
    TAB = rbind(TAB, res_stat)
    
    # Standard deviation
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) sd(x, na.rm = T)) %>%
      add_column(stat = "SD")
    TAB = rbind(TAB, res_stat)
    
    # Coefficient of variation
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) sd(x, na.rm = T)/mean(x, na.rm = T)) %>%
      add_column(stat = "CV")
    TAB = rbind(TAB, res_stat)
    
    # Proportion of wet steps
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) prWet(x)) %>%
      add_column(stat = "prWet")
    TAB = rbind(TAB, res_stat)
    
    # Lag-1 autocorrelation
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) get_autocor_lag_fast(x, k = 1, na.rm = T, all = T)) %>%
      add_column(stat = "lag1")
    TAB = rbind(TAB, res_stat)
    
    # Lag-2 autocorrelation
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) get_autocor_lag_fast(x, k = 2, na.rm = T, all = T)) %>%
      add_column(stat = "lag2")
    TAB = rbind(TAB, res_stat)
    
    # Lag-3 autocorrelation
    res_stat = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) get_autocor_lag_fast(x, k = 3, na.rm = T, all = T)) %>%
      add_column(stat = "lag3")
    TAB = rbind(TAB, res_stat)
    
    # Mean and standard deviation of the length of dry and wet spells
    data_wet_dry = df_data  %>% group_by(Season) %>% 
      summarise_at(agg_by_names, function(x) get_length_spell_mean_sd(x)) %>%
      mutate(.before = "obs")
    
    TAB_mean_sd = NULL
    for(i.s in sort(unique(df_data$Season))){
      tab = NULL
      for(i in 1:4){
        tab = rbind(tab, sapply(data_wet_dry[data_wet_dry$Season == i.s,], "[[",i))
      }
      tab = as.data.frame(tab)
      tab$stat = c("event.mean.dry", "event.mean.wet", "event.sd.dry", "event.sd.wet")
      TAB_mean_sd = rbind.data.frame(tab, TAB_mean_sd)
    }
  }
  
  # Data frame containing all statistics
  TAB_complete  = rbind.data.frame(TAB, TAB_mean_sd, TAB_return_levels)
  
  # Find which columns contain the generated scenarios
  result_cols = TAB_complete %>% dplyr::select(contains("result"))
  
  # Find the mean, q05, q50, q95, standard deviation and score as assigned by CASE framework among scenarios for each statistic
  TAB_complete$mean = apply(result_cols, 1, function(row) mean(row, na.rm = TRUE))
  TAB_complete$q05 = apply(result_cols, 1, function(row) quantile(row, probs = 0.05, na.rm = TRUE))
  TAB_complete$q50 = apply(result_cols, 1, function(row) quantile(row, probs = 0.50, na.rm = TRUE))
  TAB_complete$q95 = apply(result_cols, 1, function(row) quantile(row, probs = 0.95, na.rm = TRUE))
  TAB_complete$sd = apply(result_cols, 1, function(row) sd(row, na.rm = TRUE))
  TAB_complete$score = sapply(1:nrow(TAB_complete), function(x) compute.score(sim.mean = TAB_complete$mean[x],
                                                                              sim.sd = TAB_complete$sd[x],
                                                                              q05 = TAB_complete$q05[x],
                                                                              q95 = TAB_complete$q95[x],
                                                                              stat.obs = TAB_complete$obs[x]))
  
  # Calculate MAE and MAPE
  COPY_stat.matrix = as.data.frame(TAB_complete)
  for(i_scen in which(grepl("result.", colnames(TAB_complete)))){
    COPY_stat.matrix[, i_scen] = abs(COPY_stat.matrix[, i_scen] - COPY_stat.matrix$obs)
  }
  COPY_stat.matrix = COPY_stat.matrix %>%
    mutate(MAE = rowMeans(as.matrix(COPY_stat.matrix[, grepl("result.", colnames(TAB_complete))]), na.rm = TRUE)) %>%
    mutate(MAPE = 100*MAE/abs(obs))
  
  # Add MAE and MAPE to the statistics data frame
  TAB_complete$MAE = COPY_stat.matrix$MAE
  TAB_complete$MAPE = COPY_stat.matrix$MAPE

  return(TAB_complete)
}