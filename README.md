
<!-- README.md is generated from README.Rmd. Please edit that file -->

# disaggMRC

The disaggMRC package provides R routines for fitting a stochastic
weather generator and generating sub-daily precipitation time series
scenarios. This manual consists of two parts, which can be run
independently of each other.

In the first part, we estimate the parameters of the 4 MRC models A, A+,
B and B+ presented by Maloku et al., 2023, generate sub-daily
precipitation time series scenarios and evaluate the 4 models by
calculating precipitation statistics.

In the second part, we provide an example of the generation of long
precipitation time series scenarios at sub-daily resolution. The
stochastic weather generator is composed of the GWEX model (Evin et al.,
2018), which generates daily precipitation time series scenarios, and
the MRC disaggregation model B+.

## 0. Installation

R, RStudio and Rtools have to be installed. The Rtools version has to
correspond to the R version (version 4.4.1 works, while version 4.5
seems not to). Launch RStudio and start by cleaning the entire
workspace.

``` r
rm(list=ls())
```

Then install and load the required packages.

``` r
# Packages should be installed only once
if(F){
  install.packages("devtools")
  install.packages("lubridate")
  install.packages("dplyr")
  install.packages("stats")
  install.packages("pracma")
  install.packages("RcppRoll")
  install.packages("tibble")
  install.packages("tidyr")
  install.packages("mev")
  install.packages("purrr")
  install.packages("ggplot2")
  
  devtools::install_github("guillaumeevin/GWEX")
  devtools::install_github("drabuharuna/egpdIDF")
  devtools::install_github("Caroline-Legrand/disaggMRC")
}

library(devtools)
library(lubridate)
library(dplyr)
library(stats)
library(pracma)
library(RcppRoll)
library(tibble)
library(tidyr)
library(mev)
library(purrr)
library(ggplot2)
library(GWEX)
library(egpdIDF)
library(disaggMRC)
```

Please note that the data used in the following were downloaded when the
disaggMRC package was installed, but are not visible in the RStudio
environment.

## 1.1. Parameter estimation of 4 MRC models

In this part, we estimate the parameters of the 4 MRC models presented
by Maloku et al., 2023. The data required for this is a time series of
precipitation at a resolution of 10 minutes. As an example, a 40-year
precipitation time series observed at a station in Switzerland at a
resolution of 10 minutes is provided.

``` r
# Show the first lines of the 10-minute precipitation dataset
head(PrecipData10min)
```

Define some options required to fit MRC models. We propose tested
values. If wanted, these values can be modified, but this should be done
with caution.

``` r
listOptionsMRC = list(I_min_fix = 0.01, # [mm/h]. For intensities above this value, the scaling model of alpha is constant. Recommended 0.01 mm/h. One might consider values between 0.001 and 0.1 mm/h. It corresponds to the I_zero of Eq. (2.13) Kaltrina Maloku's thesis
                      I_max_fix = 7, # [mm/h]. For intensities below this value, the scaling model of alpha is constant. Recommended 7 mm/h. One might consider values between 5 and 10 mm/h. It corresponds to the I_one of Eq. (2.13) Kaltrina Maloku's thesis
                      I_start_class = 0.001, # [mm/h]. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1 mm/h
                      threshold_int = 0.002) # [mm/h]. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1 mm/h
```

Estimate scaling parameters for models A and A+ for each season. The
cascade weights are considered to depend on the temporal aggregation
level and the intensity of the precipitation to be disaggregated. In
model A, the precipitation asymmetry is not taken into account, whereas
it is in model A+.

``` r
Model_A = get_scaling_params_Intensity_aggLevel(vecPrecip = PrecipData10min$obs, # The high-resolution observed data needed for parameter estimation
                                                vecDates = PrecipData10min$date, # The corresponding vector of dates
                                                resVecPrecip = 10, # The temporal resolution in minutes
                                                aggLevels = c(80,160,320,640,1280), # The aggregation levels at which MRC parameters are estimated
                                                by_season = T, # Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                listOptions = listOptionsMRC) # List of options for the MRC model
```

The function “get_scaling_params_Intensity_aggLevel” returns a data
frame of scaling parameters for each season and a list of plots for
showing empirical estimates and fitted scaling models.

``` r
# Define the directory where to save files
dir = "./param_4_models_MRC"

# Create this directory
dir.create(dir)

# Save scaling parameters
saveRDS(Model_A$params, paste0(dir,"/params_model_A.RData"))

# Save the plot of the non-zero subdivision probability Px
ggsave(paste0(dir,"/Px_model_A.png"), Model_A$fig_plots$Px_aggLevel_intensity)

# Save the plots of the two components h and g of the alpha parameter of the Beta distribution: alpha(tau, I) = h(tau)*g(I)
ggsave(paste0(dir,"/h_model_A.png"), Model_A$fig_plots$Alpha_aggLevel)
ggsave(paste0(dir,"/g_model_A.png"), Model_A$fig_plots$Alpha_star_aggLevel_intensity)

# Save the plot of the probability asymmetry ratio phi
ggsave(paste0(dir,"/phi_model_A.png"), Model_A$fig_plots$Asymm_ratio)

# Save the plot of the mean of the cascade weights
ggsave(paste0(dir,"/mean_model_A.png"), Model_A$fig_plots$Asymm_mean)
```

Estimate scaling parameters for models B and B+ for each season. The
cascade weights are considered to depend only on the intensity of the
precipitation to be disaggregated. In model B, the precipitation
asymmetry is not taken into account, whereas it is in model B+.

``` r
Model_B = get_scaling_params_Intensity(vecPrecip = PrecipData10min$obs, # The high-resolution observed data needed for parameter estimation
                                       vecDates = PrecipData10min$date, # vecDates = PrecipData10min$date, # The corresponding vector of dates
                                       resVecPrecip = 10, # The temporal resolution in minutes
                                       aggLevels = c(80,160,320,640,1280), # The aggregation levels at which MRC parameters are estimated
                                       by_season = T, # Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                       listOption = listOptionsMRC) # List of options for the MRC model
```

The function “get_scaling_params_Intensity” returns a data frame of
scaling parameters for each seasonn and a list of plots for showing
empirical estimates and fitted scaling models.

``` r
# Save scaling parameters
saveRDS(Model_B$params, paste0(dir,"/params_model_B.RData"))

# Save the plot of the non-zero subdivision probability Px
ggsave(paste0(dir,"/Px_model_B.png"), Model_B$fig_plots$Px_intensity)

# Save the plot of the alpha parameter of the Beta distribution
ggsave(paste0(dir,"/alpha_model_B.png"), Model_B$fig_plots$Alpha_intensity)

# Save the plot of the probability asymmetry ratio phi
ggsave(paste0(dir,"/phi_model_B.png"), Model_B$fig_plots$Asymm_ratio)

# Save the plot of the mean of the cascade weights
ggsave(paste0(dir,"/mean_model_B.png"), Model_B$fig_plots$Asymm_mean)
```

## 1.2. Disaggregation and evaluation of the 4 MRC models

In this part, we disaggregate a quasi-daily precipitation time series
(1280-minute resolution) to 40-minute resolution with the 4 MRC models
presented by Maloku et al., 2023. 30 precipitation time series scenarios
are generated by each model. The precipitation time series at
1280-minute resolution was obtained by aggregating the precipitation
time series at 10-minute resolution provided in part 1.1.

The precipitation time series at 10-minute resolution was also
aggregated to 40-minute resolution for model evaluation. Precipitation
statistics (e.g. standard deviation of precipitation amounts, proportion
of wet steps, lag-1 autocorrelation, mean length of wet steps, 5- and
20- year return levels) are calculated for each season at 40- and
160-minute resolutions from observed and disaggregated precipitation
time series scenarios.

``` r
# Show the first lines of the 40-minute aggregated precipitation dataset
head(PrecipData40min)

# Show the first lines of the 1280-minute aggregated precipitation dataset
head(PrecipData1280min)

# Define the number of scenarios to be generated
nb_scenarios_mrc = 30

# Define the directory where to save files
dir = "./disag_eval_4_models_MRC"

# Create this directory
dir.create(dir)
```

Disaggregation and evaluation of model A

``` r
# Define a data frame that will hold all the statistics calculated from observed and disaggregated precipitation time series scenarios
TAB_stats_model_A = NULL

# Define a data frame that will hold observed and disaggregated precipitation time series scenarios at 40-minute resolution
df_precip_40min_model_A = data.frame(date = PrecipData40min$date, obs = PrecipData40min$obs)

# Disaggregation of resolution 1280 minutes to 40 minutes
for(i_scen_mrc in 1:nb_scenarios_mrc){
  # Set a seed for random generation to be able to reproduce the results
  set.seed(i_scen_mrc)
  
  # Generation of one scenario
  one_scen_mrc_model_A = disaggregate_precip_MRC_Intensity_aggLevel(vecPrecip_target = PrecipData1280min$obs, # The vector of quasi-daily precipitations to be disaggregated
                                                                    vecDates_target = PrecipData1280min$date, # The corresponding vector of dates
                                                                    params_scaling = Model_A$params, # Parameters of the MRC model A
                                                                    by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                                    res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                                    res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                                    asymmetry_option = F) # In the MRC model A, the disaggregation does not depend on the asymmetry model
  
  # Add this scenario to the data frame
  df_precip_40min_model_A[, paste0("result.",i_scen_mrc)] = one_scen_mrc_model_A
}

# Save observed and disaggregated precipitation time series scenarios at 40-minute resolution for model A
saveRDS(df_precip_40min_model_A, paste0(dir,"/scenarios_prec_40min_model_A.RData"))

# Add a column for the season and year, which are necessary for evaluation on a seasonal basis
df_precip_40min_model_A = df_precip_40min_model_A %>% mutate(Season = month2season(lubridate::month(date)), Year = year(date))

# Calculate different precipitation statistics at 40-minutes resolution
TAB_stats_40min_model_A = get_standard_stats_extremes(df_data = df_precip_40min_model_A %>% dplyr::filter(!is.na(Season))) %>%
  mutate(resolution = 40, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Aggregate observed and disaggregated precipitation time series scenarios to 160-minute resolution
df_precip_160min_model_A = RcppRoll::roll_sum(x = df_precip_40min_model_A[, c("obs", paste0("result.", 1:nb_scenarios_mrc))] %>% as.matrix(), n = 4, by = 4, fill = NA, align = "left", na.rm = F) %>% as.data.frame

# Add date, season and year vectors
df_precip_160min_model_A = cbind.data.frame(date = df_precip_40min_model_A[ , "date"],
                                            Season = df_precip_40min_model_A[ , "Season"],
                                            Year = df_precip_40min_model_A[ , "Year"], 
                                            df_precip_160min_model_A)

# Select only values that are not NA
df_precip_160min_model_A = df_precip_160min_model_A[seq(1, to = nrow(df_precip_160min_model_A), by = 4), ]
df_precip_160min_model_A = df_precip_160min_model_A %>% dplyr::filter(!is.na(Season))

# Calculate different precipitation statistics at 160-minute resolution
TAB_stats_160min_model_A = get_standard_stats_extremes(df_data = df_precip_160min_model_A) %>%
  mutate(resolution = 160, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Combine precipitation statistics for 40- and 160-minute resolutions
TAB_stats_model_A = rbind.data.frame(TAB_stats_40min_model_A, TAB_stats_160min_model_A)

# Save precipitation statistics for model A
saveRDS(TAB_stats_model_A, paste0(dir,"/stats_eval_model_A.RData"))
```

Disaggregation and evaluation of model A+

``` r
# Define a data frame that will hold all the statistics calculated from observed and disaggregated precipitation time series scenarios
TAB_stats_model_A_plus = NULL

# Define a data frame that will hold observed and disaggregated precipitation time series scenarios at 40-minute resolution
df_precip_40min_model_A_plus = data.frame(date = PrecipData40min$date, obs = PrecipData40min$obs)

# Disaggregation of resolution 1280 minutes to 40 minutes
for(i_scen_mrc in 1:nb_scenarios_mrc){
  # Set a seed for random generation to be able to reproduce the results
  set.seed(i_scen_mrc)
  
  # Generation of one scenario
  one_scen_mrc_model_A_plus = disaggregate_precip_MRC_Intensity_aggLevel(vecPrecip_target = PrecipData1280min$obs, # The vector of quasi-daily precipitations to be disaggregated
                                                                         vecDates_target = PrecipData1280min$date, # The corresponding vector of dates
                                                                         params_scaling = Model_A$params, # Parameters of the MRC model A
                                                                         by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                                         res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                                         res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                                         asymmetry_option = T) # In the MRC model A+, the disaggregation depends on the asymmetry model
  
  # Add this scenario to the data frame
  df_precip_40min_model_A_plus[, paste0("result.",i_scen_mrc)] = one_scen_mrc_model_A_plus
}

# Save observed and disaggregated precipitation time series scenarios at 40-minute resolution for model A+
saveRDS(df_precip_40min_model_A_plus, paste0(dir,"/scenarios_prec_40min_model_A_plus.RData"))

# Add a column for the season and year, which are necessary for evaluation on a seasonal basis
df_precip_40min_model_A_plus = df_precip_40min_model_A_plus %>% mutate(Season = month2season(lubridate::month(date)), Year = year(date))

# Calculate different precipitation statistics at 40-minutes resolution
TAB_stats_40min_model_A_plus = get_standard_stats_extremes(df_data = df_precip_40min_model_A_plus %>% dplyr::filter(!is.na(Season))) %>%
  mutate(resolution = 40, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Aggregate observed and disaggregated precipitation time series scenarios to 160-minute resolution
df_precip_160min_model_A_plus = RcppRoll::roll_sum(x = df_precip_40min_model_A_plus[, c("obs", paste0("result.", 1:nb_scenarios_mrc))] %>% as.matrix(), n = 4, by = 4, fill = NA, align = "left", na.rm = F) %>% as.data.frame

# Add date, season and year vectors
df_precip_160min_model_A_plus = cbind.data.frame(date = df_precip_40min_model_A_plus[ , "date"],
                                                 Season = df_precip_40min_model_A_plus[ , "Season"],
                                                 Year = df_precip_40min_model_A_plus[ , "Year"], 
                                                 df_precip_160min_model_A_plus)

# Select only values that are not NA
df_precip_160min_model_A_plus = df_precip_160min_model_A_plus[seq(1, to = nrow(df_precip_160min_model_A_plus), by = 4), ]
df_precip_160min_model_A_plus = df_precip_160min_model_A_plus %>% dplyr::filter(!is.na(Season))

# Calculate different precipitation statistics at 160-minute resolution
TAB_stats_160min_model_A_plus = get_standard_stats_extremes(df_data = df_precip_160min_model_A_plus) %>%
  mutate(resolution = 160, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Combine precipitation statistics for 40- and 160-minute resolutions
TAB_stats_model_A_plus = rbind.data.frame(TAB_stats_40min_model_A_plus, TAB_stats_160min_model_A_plus)

# Save precipitation statistics for model A+
saveRDS(TAB_stats_model_A_plus, paste0(dir,"/stats_eval_model_A_plus.RData"))
```

Disaggregation and evaluation of model B

``` r
# Define a data frame that will hold all the statistics calculated from observed and disaggregated precipitation time series scenarios
TAB_stats_model_B = NULL

# Define a data frame that will hold observed and disaggregated precipitation time series scenarios at 40-minute resolution
df_precip_40min_model_B = data.frame(date = PrecipData40min$date, obs = PrecipData40min$obs)

# Disaggregation of resolution 1280 minutes to 40 minutes
for(i_scen_mrc in 1:nb_scenarios_mrc){
  # Set a seed for random generation to be able to reproduce the results
  set.seed(i_scen_mrc)
  
  # Generation of one scenario
  one_scen_mrc_model_B = disaggregate_precip_MRC_Intensity(vecPrecip_target = PrecipData1280min$obs, # The vector of quasi-daily precipitations to be disaggregated
                                                           vecDates_target = PrecipData1280min$date, # The corresponding vector of dates
                                                           params_scaling = Model_B$params, # Parameters of the MRC model B
                                                           by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                           res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                           res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                           asymmetry_option = F) # In the MRC model B, the disaggregation does not depend on the asymmetry model

  # Add this scenario to the data frame
  df_precip_40min_model_B[, paste0("result.",i_scen_mrc)] = one_scen_mrc_model_B
}

# Save observed and disaggregated precipitation time series scenarios at 40-minute resolution for model B
saveRDS(df_precip_40min_model_B, paste0(dir,"/scenarios_prec_40min_model_B.RData"))

# Add a column for the season and year, which are necessary for evaluation on a seasonal basis
df_precip_40min_model_B = df_precip_40min_model_B %>% mutate(Season = month2season(lubridate::month(date)), Year = year(date))

# Calculate different precipitation statistics at 40-minutes resolution
TAB_stats_40min_model_B = get_standard_stats_extremes(df_data = df_precip_40min_model_B %>% dplyr::filter(!is.na(Season))) %>%
  mutate(resolution = 40, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Aggregate observed and disaggregated precipitation time series scenarios to 160-minute resolution
df_precip_160min_model_B = RcppRoll::roll_sum(x = df_precip_40min_model_B[, c("obs", paste0("result.", 1:nb_scenarios_mrc))] %>% as.matrix(), n = 4, by = 4, fill = NA, align = "left", na.rm = F) %>% as.data.frame

# Add date, season and year vectors
df_precip_160min_model_B = cbind.data.frame(date = df_precip_40min_model_B[ , "date"],
                                            Season = df_precip_40min_model_B[ , "Season"],
                                            Year = df_precip_40min_model_B[ , "Year"], 
                                            df_precip_160min_model_B)

# Select only values that are not NA
df_precip_160min_model_B = df_precip_160min_model_B[seq(1, to = nrow(df_precip_160min_model_B), by = 4), ]
df_precip_160min_model_B = df_precip_160min_model_B %>% dplyr::filter(!is.na(Season))

# Calculate different precipitation statistics at 160-minute resolution
TAB_stats_160min_model_B = get_standard_stats_extremes(df_data = df_precip_160min_model_B) %>%
  mutate(resolution = 160, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Combine precipitation statistics for 40- and 160-minute resolutions
TAB_stats_model_B = rbind.data.frame(TAB_stats_40min_model_B, TAB_stats_160min_model_B)

# Save precipitation statistics for model B
saveRDS(TAB_stats_model_B, paste0(dir,"/stats_eval_model_B.RData"))
```

Disaggregation and evaluation of model B+

``` r
# Define a data frame that will hold all the statistics calculated from observed and disaggregated precipitation time series scenarios
TAB_stats_model_B_plus = NULL

# Define a data frame that will hold observed and disaggregated precipitation time series scenarios at 40-minute resolution
df_precip_40min_model_B_plus = data.frame(date = PrecipData40min$date, obs = PrecipData40min$obs)

# Disaggregation of resolution 1280 minutes to 40 minutes
for(i_scen_mrc in 1:nb_scenarios_mrc){
  # Set a seed for random generation to be able to reproduce the results
  set.seed(i_scen_mrc)
  
  # Generation of one scenario
  one_scen_mrc_model_B_plus = disaggregate_precip_MRC_Intensity(vecPrecip_target = PrecipData1280min$obs, # The vector of quasi-daily precipitations to be disaggregated
                                                                vecDates_target = PrecipData1280min$date, # The corresponding vector of dates
                                                                params_scaling = Model_B$params, # Parameters of the MRC model B
                                                                by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                                res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                                res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                                asymmetry_option = T) # In the MRC model B+, the disaggregation depends on the asymmetry model
  
  # Add this scenario to the data frame
  df_precip_40min_model_B_plus[, paste0("result.",i_scen_mrc)] = one_scen_mrc_model_B_plus
}

# Save observed and disaggregated precipitation time series scenarios at 40-minute resolution for model B+
saveRDS(df_precip_40min_model_B_plus, paste0(dir,"/scenarios_prec_40min_model_B_plus.RData"))

# Add a column for the season and year, which are necessary for evaluation on a seasonal basis
df_precip_40min_model_B_plus = df_precip_40min_model_B_plus %>% mutate(Season = month2season(lubridate::month(date)), Year = year(date))

# Calculate different precipitation statistics at 40-minutes resolution
TAB_stats_40min_model_B_plus = get_standard_stats_extremes(df_data = df_precip_40min_model_B_plus %>% dplyr::filter(!is.na(Season))) %>%
  mutate(resolution = 40, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Aggregate observed and disaggregated precipitation time series scenarios to 160-minute resolution
df_precip_160min_model_B_plus = RcppRoll::roll_sum(x = df_precip_40min_model_B_plus[, c("obs", paste0("result.", 1:nb_scenarios_mrc))] %>% as.matrix(), n = 4, by = 4, fill = NA, align = "left", na.rm = F) %>% as.data.frame

# Add date, season and year vectors
df_precip_160min_model_B_plus = cbind.data.frame(date = df_precip_40min_model_B_plus[ , "date"],
                                                 Season = df_precip_40min_model_B_plus[ , "Season"],
                                                 Year = df_precip_40min_model_B_plus[ , "Year"], 
                                                 df_precip_160min_model_B_plus)

# Select only values that are not NA
df_precip_160min_model_B_plus = df_precip_160min_model_B_plus[seq(1, to = nrow(df_precip_160min_model_B_plus), by = 4), ]
df_precip_160min_model_B_plus = df_precip_160min_model_B_plus %>% dplyr::filter(!is.na(Season))

# Calculate different precipitation statistics at 160-minute resolution
TAB_stats_160min_model_B_plus = get_standard_stats_extremes(df_data = df_precip_160min_model_B_plus) %>%
  mutate(resolution = 160, .after = "stat") %>% dplyr::filter(!is.na(Season))
  
# Combine precipitation statistics for 40- and 160-minute resolutions
TAB_stats_model_B_plus = rbind.data.frame(TAB_stats_40min_model_B_plus, TAB_stats_160min_model_B_plus)

# Save precipitation statistics for model B+
saveRDS(TAB_stats_model_B_plus, paste0(dir,"/stats_eval_model_B_plus.RData"))
```

## 2. Long scenarios generation

In this part, we provide an example of the generation of 5 precipitation
time series scenarios with a duration of 100 years and a resolution of
30 minutes. The stochastic weather generator is made up of the GWEX
model (Evin et al., 2018), which generates daily precipitation time
series scenarios, and the MRC disaggregation model B+.

The dataset used to fit the MRC model B+ is a 17-year hourly time series
of mean areal precipitation observed in Switzerland. This hourly
precipitation time series was also aggregated to a daily resolution to
fit the GWEX model.

``` r
# Show the first lines of the hourly (60 minutes) and daily (1440 minutes) mean areal precipitation datasets
head(MeanArealPrecipData60min)
head(MeanArealPrecipData1440min)

# Define the number of scenarios to be generated
nb_scenarios_gwex = 1
nb_scenarios_mrc = 5

# Define the vector of dates for generated daily scenarios
vecDatesSIM_GWEX = seq(from = as.Date("01/01/2000", format = "%d/%m/%Y"), to = as.Date("31/12/2099", format = "%d/%m/%Y"), by = "day")

# Vector of months for each day
vecDatesSIM_GWEX_month = lubridate::month(vecDatesSIM_GWEX)

# Starting date of the generated time series
d.start = vecDatesSIM_GWEX[1]

# Ending date of the generated time series
d.end = rev(vecDatesSIM_GWEX)[1]

# Define the directory where to save files
dir = "./long_scenarios"

# Create this directory
dir.create(dir)
```

Parameter estimation - MRC model B+

``` r
# Define some options required to fit the MRC model B+
# We propose tested values. If wanted, these values can be modified, but this should be done with caution
listOptionsMRC = list(I_min_fix = 0.01, # [mm/h]. For intensities above this value, the scaling model of alpha is constant. Recommended 0.01 mm/h. One might consider values between 0.001 and 0.1 mm/h. It corresponds to the I_zero of Eq. (2.13) Kaltrina Maloku's thesis
                      I_max_fix = 7, # [mm/h]. For intensities below this value, the scaling model of alpha is constant. Recommended 7 mm/h. One might consider values between 5 and 10 mm/h. It corresponds to the I_one of Eq. (2.13) Kaltrina Maloku's thesis
                      I_start_class = 0.001, # [mm/h]. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1 mm/h
                      threshold_int = 0.002) # [mm/h]. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1 mm/h

# Estimate scaling parameters for each month for the MRC model B+
params_scaling_MRC = get_scaling_params_Intensity(vecPrecip = MeanArealPrecipData60min$obs, # The high-resolution observed data needed for parameter estimation
                                                  vecDates = MeanArealPrecipData60min$date, # The corresponding vector of dates
                                                  resVecPrecip = 60, # The temporal resolution in minutes
                                                  aggLevels = 60*c(2,6,12,24), # The aggregation levels at which MRC parameters are estimated before fitting scaling models. Other pair values between 2 and 24 are allowed and may be considered in case of estimation problems
                                                  by_season = F, # Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis. For a more robust estimation, but slightly lower performance when evaluated on a monthly basis, consider an estimation by season
                                                  listOptions = listOptionsMRC) # List of options for the MRC model
                                                  
# Scaling parameters
params_scaling_MRC = params_scaling_MRC$params
```

Parameter estimation - GWEX model

``` r
# Estimate the parameters of the EGDP for daily data using the IDF model of Haruna et al., 2023 estimated from hourly data
params_egpd = fit_EGPD_IDF(vec.precip = MeanArealPrecipData60min$obs, # The vector of observed hourly precipitation
                           vec.dates = MeanArealPrecipData60min$date) # The corresponding vector of dates

# Create a GWEX object required to fit the GWEX model
myObsPrec = GwexObs(variable = "Prec", # The variable name
                    date = as.Date(MeanArealPrecipData1440min$date), # The vector of daily dates
                    obs = as.matrix(MeanArealPrecipData1440min$obs, ncol = 1)) # The daily observations as a one-column matrix

# Define some options
listOptionPrec = list(nLag = 2, # Order of the Markov chain for precipitation occurrences. Suggested value: 2. Other accepted values are 3 and 4. The length of wet/dry sequences might be better reproduced at 24h when 3 or 4, but be careful as parameter estimation may be less robust
                      th = 0, # Wet/dry days threshold for parameter estimation. Suggested value: 0 
                      isMAR = T, # Should we consider temporal autocorrelation? Suggested value: T
                      typeMargin = "EGPD") # The marginal distribution for positive daily precipitation amounts

# Fit the GWEX model
myParPrec = fitGwexModel(objGwexObs = myObsPrec, # The object of observation data
                         parMargin = params_egpd, # Parameters of the EGPD previously estimated using the IDF model
                         listOption = listOptionPrec) # List of options for parameter estimation
```

Generation of long precipitation time series scenarios at 30-minute
resolution

``` r
# Set a seed for random generation to be able to reproduce the results
set.seed(2024)

for (i_scen_gwex in 1:nb_scenarios_gwex){
  # Generate a daily precipitation scenario using the GWEX model
  mySimPrec = GWEX::simGwexModel(objGwexFit = myParPrec, # GWEX parameters
                                 nb.rep = 1, # Number of generated scenarios. This script is adapted for a value of 1
                                 d.start = d.start, # Starting date 
                                 d.end = d.end, # Ending date 
                                 objGwexObs = NULL, # Do not modify 
                                 prob.class = NULL, # Do not modify 
                                 objGwexSim = NULL) # Do not modify 
 
  # Disaggregate the daily scenario generated by the GWEX model 5 times using the MRC model B+
  for(i_scen_mrc in 1:nb_scenarios_mrc){
    one_scen_mrc = disaggregate_precip_MRC_Intensity(vecPrecip_target = mySimPrec@sim[,1,], # The vector of daily precipitations generated by the GWEX model
                                                     vecDates_target = vecDatesSIM_GWEX, # The corresponding vector of dates
                                                     params_scaling = params_scaling_MRC, # Parameters of the MRC model B+
                                                     by_season = F, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                     res_coarse_aggLevel = 1440, # The temporal resolution in minutes of the data to be disaggregated
                                                     res_fine_aggLevel = 30, # The target temporal resolution in minutes of the high-resolution data
                                                     asymmetry_option = T) # In the MRC model B+, the disaggregation depends on the asymmetry model
    
    # Save each 30-minute precipitation time series scenario
    saveRDS(one_scen_mrc, paste0(dir,"/scenarios_prec_30min_GWEX_",i_scen_gwex,"_MRC_B_plus_",i_scen_mrc,".RData"))
  }
}
```

## References

Evin, G., Favre, A.-C., and Hingray, B. (2018). Stochastic generation of
multi-site daily precipitation focusing on extreme events. Hydrol. Earth
Syst. Sci., 22, 655–672, <https://doi.org/10.5194/hess-22-655-2018>

Haruna, A., Blanchet, J., and Favre, A.-C. (2023). Modeling
intensity-duration-frequency curves for the whole range of non-zero
precipitation: A comparison of models. Water Resour. Res., 59,
e2022WR033362, <https://doi.org/10.1029/2022WR033362>

Maloku, K., Hingray, B., and Evin, G. (2023). Accounting for
precipitation asymmetry in a multiplicative random cascade
disaggregation model. Hydrol. Earth Syst. Sci., 27, 3643–3661,
<https://doi.org/10.5194/hess-27-3643-2023>
