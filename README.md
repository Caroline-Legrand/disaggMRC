
<!-- README.md is generated from README.Rmd. Please edit that file -->

# disaggMRC

The disaggMRC package provides R scripts for fitting a stochastic
weather generator and generating sub-daily precipitation time series
scenarios.

First, the parameter estimation of the 4 MRC models A, A+, B and B+
presented by Maloku et al., 2023 is provided. The data set used for the
parameter estimation is a 10-minute resolution precipitation data set
observed at a station in Switzerland. An example of disaggregation by
the 4 MRC models of the corresponding quasi-daily (1280 minutes)
observed precipitation time series at the target resolution of 40
minutes is provided.

Next, an example of the generation of 5 precipitation time series
scenarios with a duration of 100 years and a time step of 30 minutes is
provided. The stochastic weather generator is made up of the GWEX model
(Evin et al., 2018), which generates daily time series scenarios, and
the MRC disaggregation model B+. The data set used to fit the models is
an hourly dataset of mean areal precipitation observed in Switzerland.

## Installation

Clear the entire workspace

``` r
rm(list=ls())
```

Install and load the required packages

``` r
# Packages should be installed only once
if(F){
  devtools::install_github("guillaumeevin/GWEX")
  devtools::install_github("drabuharuna/egpdIDF")
  devtools::install_github("Caroline-Legrand/disaggMRC")
}

library(lubridate)
library(dplyr)
library(stats)
library(pracma)
library(RcppRoll)
library(tibble)
library(tidyr)
library(mev)
library(purrr)
library(GWEX)
library(egpdIDF)
library(disaggMRC)
```

## Loading and preparation of example data sets

Load 10-minute and hourly time series of observed precipitation

``` r
# Load an example of a 30-year precipitation time series observed at a station in Switzerland with a resolution of 10 minutes
load("~/disaggMRC/data/PrecipData10min.rda")

# Load an example of a 17-year hourly time series of mean areal precipitation observed in Switzerland
load("~/disaggMRC/data/PrecipData60min.rda")
```

Aggregate hourly observations to daily observations for GWEX model
fitting

``` r
# Create a vector of daily precipitation observations
vec_obs_daily_precip = RcppRoll::roll_sum(x = PrecipData60min$obs, n = 24, by = 24, fill = NA, align = "left", na.rm = F)[seq(1, by = 24, to = length(PrecipData60min$obs))]

# Create a vector of daily dates from the vector of hourly dates
vec_obs_daily_dates = as.Date(PrecipData60min$date[seq(1, by = 24, to = length(PrecipData60min$date))])
```

## Parameter estimation of the 4 MRC models presented by Maloku et al., 2023 and example of disaggregation by the 4 MRC models

Parameter estimation

``` r
# List of some options needed for fitting the MRC model
# We propose tested values. If wanted, these values can be modified, but this should be done with caution

listOptionsMRC = list(I_min_fix = 0.01, # [mm/h]. For intensities above this value, the scaling model of alpha is constant. Recommended 0.01 mm/h. One might consider values between 0.001 and 0.1 mm/h. It corresponds to the I_zero of Eq. (2.13) Kaltrina Maloku's thesis
                      I_max_fix = 7, # [mm/h]. For intensities below this value, the scaling model of alpha is constant. Recommended 7 mm/h. One might consider values between 5 and 10 mm/h. It corresponds to the I_one of Eq. (2.13) Kaltrina Maloku's thesis
                      I_start_class = 0.001, # [mm/h]. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1 mm/h
                      threshold_int = 0.002) # [mm/h]. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1 mm/h
```

Estimation of scaling parameters for models A and A+

``` r
# The cascade weights are considered to depend on the temporal aggregation level and the intensity of the precipitation to be disaggregated
# In model A, the precipitation asymmetry is not taken into account, whereas it is in model A+

Model_A = get_scaling_params_Intensity_aggLevel(vecPrecip = PrecipData10min$obs, # The high-resolution observed data needed for parameter estimation
                                                vecDates = PrecipData10min$date, # The corresponding vector of dates
                                                resVecPrecip = 10, # The temporal resolution in minutes
                                                aggLevels = c(80,160,320,640,1280), # The aggregation levels at which MRC parameters are estimated
                                                by_season = T, # Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                listOptions = listOptionsMRC) # List of options for the MRC model

# The function "get_scaling_params_Intensity_aggLevel" returns a data frame of scaling parameters and a list of plots for showing empirical estimates and fitted scaling models

# Scaling parameters. Seasons "1", "2", "3", "4" correspond to "DJF", "MAM", "JJA", "SON"
Model_A$params

# The non-zero subdivision probability Px
Model_A$fig_plots$Px_aggLevel_intensity

# The two components h and g of the alpha parameter of the Beta distribution: alpha(tau, I) = h(tau)*g(I)
Model_A$fig_plots$Alpha_aggLevel
Model_A$fig_plots$Alpha_star_aggLevel_intensity

# The probability asymmetry ratio phi
Model_A$fig_plots$Asymm_ratio

# The mean of the cascade weights
Model_A$fig_plots$Asymm_mean
```

Estimation of scaling parameters for models B and B+

``` r
# The cascade weights are considered to depend on the intensity of the precipitation to be disaggregated
# In model B, the precipitation asymmetry is not taken into account, whereas it is in model B+

Model_B = get_scaling_params_Intensity(vecPrecip = PrecipData10min$obs, # The high-resolution observed data needed for parameter estimation
                                       vecDates = PrecipData10min$date, # vecDates = PrecipData10min$date, # The corresponding vector of dates
                                       resVecPrecip = 10, # The temporal resolution in minutes
                                       aggLevels = c(80,160,320,640,1280), # The aggregation levels at which MRC parameters are estimated
                                       by_season = T, # Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                       listOption = listOptionsMRC) # List of options for the MRC model
                                      
# The function "get_scaling_params_Intensity" returns a data frame of scaling parameters and a list of plots for showing empirical estimates and fitted scaling models

# Scaling parameters. Seasons "1", "2", "3", "4" correspond to "DJF", "MAM", "JJA", "SON"
Model_B$params

# The non-zero subdivision probability Px 
Model_B$fig_plots$Px_intensity

# The alpha parameter of the Beta distribution
Model_B$fig_plots$Alpha_intensity

# The probability asymmetry ratio phi
Model_B$fig_plots$Asymm_ratio

# The mean of the cascade weights
Model_B$fig_plots$Asymm_mean
```

Get a coarse resolution precipitation time series (1280 minutes) that
will be used as time series to be disaggregated

``` r
# Vector of dates, one date per day
vecDates_target = unique(as.Date(PrecipData10min$date))

# Get 1280-minute indices
indices_1280min = get_indices_center_1280(length(vecDates_target))

# Get an accumulated precipitation time series of 1280 minutes from the 10-minute time series 
vecPrecip_target = sum_fixed_window(x = PrecipData10min$obs[indices_1280min], k = 128)
```

Disaggregation

``` r
# Define the number of generated scenarios
nb_scenarios_mrc = 10

# Define the name of the generation version
i_version_ensemble = "v1"

# Define the directory where to save files
dir = "./scenarios_4_models_MRC"

# Create this directory
dir.create(dir)

# Define the directory name for generated precipitation scenarios
dir_precip_scenarios = paste0(dir,"/",i_version_ensemble,"_prec_40min")
```

Disaggregation using the MRC model A

``` r
# Set a seed to be able to reproduce the results
set.seed(2024)

for(i_scen_mrc in 1:nb_scenarios_mrc){
  one_scen_mrc_model_A = disaggregate_precip_MRC_Intensity_aggLevel(vecPrecip_target = vecPrecip_target, # The vector of quasi-daily precipitations to be disaggregated
                                                                    vecDates_target = vecDates_target, # The corresponding vector of dates
                                                                    params_scaling = Model_A$params, # Parameters of the MRC model A
                                                                    by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                                    res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                                    res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                                    asymmetry_option = F) # In the MRC model A, the disaggregation does not depend on the asymmetry model
  
  # Save the 40-minute precipitation time series scenarios
  saveRDS(one_scen_mrc_model_A, paste0(dir_precip_scenarios,"_MRC_model_A_scen_",i_scen_mrc,".RData"))
    
  # Remove from memory
  rm(one_scen_mrc_model_A); gc()
}
```

Disaggregation using the MRC model A+

``` r
# Set a seed to be able to reproduce the results
set.seed(2024)

for(i_scen_mrc in 1:nb_scenarios_mrc){
  one_scen_mrc_model_A_plus = disaggregate_precip_MRC_Intensity_aggLevel(vecPrecip_target = vecPrecip_target, # The vector of quasi-daily precipitations to be disaggregated
                                                                         vecDates_target = vecDates_target, # The corresponding vector of dates
                                                                         params_scaling = Model_A$params, # Parameters of the MRC model A
                                                                         by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                                         res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                                         res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                                         asymmetry_option = T) # In the MRC model A+, the disaggregation depends on the asymmetry model
  
  # Save the 40-minute precipitation time series scenarios
  saveRDS(one_scen_mrc_model_A_plus, paste0(dir_precip_scenarios,"_MRC_model_A_plus_scen_",i_scen_mrc,".RData"))
    
  # Remove from memory
  rm(one_scen_mrc_model_A_plus); gc()
}
```

Disaggregation using the MRC model B

``` r
# Set a seed to be able to reproduce the results
set.seed(2024)

for(i_scen_mrc in 1:nb_scenarios_mrc){
  one_scen_mrc_model_B = disaggregate_precip_MRC_Intensity(vecPrecip_target = vecPrecip_target, # The vector of quasi-daily precipitations to be disaggregated
                                                           vecDates_target = vecDates_target, # The corresponding vector of dates
                                                           params_scaling = Model_B$params, # Parameters of the MRC model B
                                                           by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                           res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                           res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                           asymmetry_option = F) # In the MRC model B, the disaggregation does not depend on the asymmetry model
  
  # Save the 40-minute precipitation time series scenarios
  saveRDS(one_scen_mrc_model_B, paste0(dir_precip_scenarios,"_MRC_model_B_scen_",i_scen_mrc,".RData"))
    
  # Remove from memory
  rm(one_scen_mrc_model_B); gc()
}
```

Disaggregation using the MRC model B+

``` r
# Set a seed to be able to reproduce the results
set.seed(2024)

for(i_scen_mrc in 1:nb_scenarios_mrc){
  one_scen_mrc_model_B_plus = disaggregate_precip_MRC_Intensity(vecPrecip_target = vecPrecip_target, # The vector of quasi-daily precipitations to be disaggregated
                                                                vecDates_target = vecDates_target, # The corresponding vector of dates
                                                                params_scaling = Model_B$params, # Parameters of the MRC model B
                                                                by_season = T, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                                res_coarse_aggLevel = 1280, # The temporal resolution in minutes of the data to be disaggregated
                                                                res_fine_aggLevel = 40, # The target temporal resolution in minutes of the high-resolution data
                                                                asymmetry_option = T) # In the MRC model B+, the disaggregation depends on the asymmetry model
                                                                
  # Save the 40-minute precipitation time series scenarios
  saveRDS(one_scen_mrc_model_B_plus, paste0(dir_precip_scenarios,"_MRC_model_B_plus_scen_",i_scen_mrc,".RData"))
    
  # Remove from memory
  rm(one_scen_mrc_model_B_plus); gc()
}
```

## Example of the generation of 5 precipitation time series scenarios with a duration of 100 years and a time step of 30 minutes

Define some preliminary objects

``` r
# Define the number and length of generated scenarios
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
```

Define the directory and names of saved scenarios

``` r
# Define the name of the generation version
i_version_ensemble = "v1"
i_version_prec_1H = "mrc_v1"
i_version_prec_24H = "gwex_v1"

# Define the directory where to save files
dir = "./scenarios_GWEX_MRC"

# Create this directory
dir.create(dir)

# Define the directory name for generated precipitation scenarios
dir_precip_scenarios = paste0(dir,"/",i_version_ensemble,"_prec_30min")
```

Parameter estimation - MRC model B+

``` r
# List of some options needed for fitting the MRC model
# We propose tested values. If wanted, these values can be modified, but this should be done with caution

listOptionsMRC = list(I_min_fix = 0.01, # [mm/h]. For intensities above this value, the scaling model of alpha is constant. Recommended 0.01 mm/h. One might consider values between 0.001 and 0.1 mm/h. It corresponds to the I_zero of Eq. (2.13) Kaltrina Maloku's thesis
                      I_max_fix = 7, # [mm/h]. For intensities below this value, the scaling model of alpha is constant. Recommended 7 mm/h. One might consider values between 5 and 10 mm/h. It corresponds to the I_one of Eq. (2.13) Kaltrina Maloku's thesis
                      I_start_class = 0.001, # [mm/h]. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1 mm/h
                      threshold_int = 0.002) # [mm/h]. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1 mm/h

params_scaling_MRC = get_scaling_params_Intensity(vecPrecip = PrecipData60min$obs, # The high-resolution observed data needed for parameter estimation
                                                  vecDates = PrecipData60min$date, # The corresponding vector of dates
                                                  resVecPrecip = 60, # The temporal resolution in minutes
                                                  aggLevels = 60*c(2,6,12,24), # The aggregation levels at which MRC parameters are estimated before fitting scaling models. Other pair values between 2 and 24 are allowed and may be considered in case of estimation problems
                                                  by_season = F, # Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis. For a more robust estimation, but slightly lower performance when evaluated on a monthly basis, consider an estimation by season
                                                  listOptions = listOptionsMRC) # List of options for the MRC model
                                                  
# Scaling parameters. Numbers 1 to 12 of the variable Season correspond to the 12 months from January to December
params_scaling_MRC = params_scaling_MRC$params
params_scaling_MRC
```

Parameter estimation - GWEX model

``` r
# Estimate the parameters of the EGDP using the IDF model of Haruna et al., 2023
params_egpd = fit_EGPD_IDF(vec.precip = PrecipData60min$obs, # The vector of observed hourly precipitation
                           vec.dates = PrecipData60min$date) # The corresponding vector of dates

# Create GWEX object for precipitation
myObsPrec = GwexObs(variable = "Prec", # The variable name
                    date = vec_obs_daily_dates, # The vector of daily dates
                    obs = as.matrix(vec_obs_daily_precip, ncol = 1)) # The daily observations as a one-column matrix

# List of options for precipitation 
listOptionPrec = list(nLag = 2, # Order of the Markov chain for precipitation occurrences. Suggested value: 2. Other accepted values are 3 and 4. The length of wet/dry sequences might be better reproduced at 24h when 3 or 4, but be careful as parameter estimation may be less robust
                      th = 0, # Wet/dry days threshold for parameter estimation. Suggested value: 0 
                      isMAR = T, # Should we consider temporal autocorrelation? Suggested value: T
                      typeMargin = "EGPD") # The marginal distribution for positive daily precipitation amounts

# Fit the GWEX model for precipitation
myParPrec = fitGwexModel(objGwexObs = myObsPrec, # The object of observation data
                         parMargin = params_egpd, # Parameters of the EGPD previously estimated using the IDF model
                         listOption = listOptionPrec) # List of options for parameter estimation
```

Generation of precipitation time series scenarios at 30-minute
resolution

``` r
# Set a seed to be able to reproduce the results
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
 
  # Disaggregate the daily scenario generated 5 times using the MRC model B+
  for(i_scen_mrc in 1:nb_scenarios_mrc){
    one_scen_mrc = disaggregate_precip_MRC_Intensity(vecPrecip_target = mySimPrec@sim[,1,], # The vector of daily precipitations generated by the GWEX model
                                                     vecDates_target = vecDatesSIM_GWEX, # The corresponding vector of dates
                                                     params_scaling = params_scaling_MRC, # Parameters of the MRC model B+
                                                     by_season = F, # Are the parameters estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
                                                     res_coarse_aggLevel = 1440, # The temporal resolution in minutes of the data to be disaggregated
                                                     res_fine_aggLevel = 30, # The target temporal resolution in minutes of the high-resolution data
                                                     asymmetry_option = T) # In the MRC model B+, the disaggregation depends on the asymmetry model
    
    # Save the 30-minute precipitation time series scenarios
    saveRDS(one_scen_mrc, paste0(dir_precip_scenarios,"_scen_GWEX_",i_scen_gwex,"_MRC_",i_scen_mrc,".RData"))
    
    # Remove from memory
    rm(one_scen_mrc); gc()
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
