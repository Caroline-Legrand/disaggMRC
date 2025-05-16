#' get_scaling_params_Intensity
#' 
#' @description A function that estimates the scaling parameters of a MRC model, where the cascade weights are considered to depend on the intensity of the precipitation to be disaggregated.
#' In the first step, the observed weights are calculated for each temporal aggregation level. Also at this step, the \eqn{z} index is calculated for each observed precipitation intensity.
#' In the second step, the MRC parameters are estimated by considering a dependency on the temporal aggregation level, the precipitation intensity class and the asymmetry \eqn{z} index class.
#' In the last step, the scaling models are fitted to the MRC parameters estimated by the linear least squares method.
#' Plots showing the estimated parameters and their scaling models are also returned
#' 
#' @param vecPrecip A vector of observed precipitations, in mm
#' @param vecDates The corresponding vector of dates
#' @param resVecPrecip Resolution of the time series, in minutes
#' @param aggLevels A vector of temporal aggregation levels in the cascade estimation procedure, in minutes
#' @param by_season Logical. Should the parameters be estimated on a seasonal basis? If F, the parameters are estimated on a monthly basis, if T on a seasonal basis
#' @param listOptions List with the following fields:
#' \itemize{
#'   \item \strong{I_min_fix}: A numerical value, in mm/h. For intensities above this value, the scaling model of alpha is constant. Recommended 0.01 mm/h. One might consider values between 0.001 and 0.1 mm/h. It corresponds to the I_zero of Eq. (2.13) Kaltrina Maloku's thesis
#'   \item \strong{I_max_fix}: A numerical value, in mm/h. For intensities below this value, the scaling model of alpha is constant. Recommended 7 mm/h. One might consider values between 5 and 10 mm/h. It corresponds to the I_one of Eq. (2.13) Kaltrina Maloku's thesis
#'   \item \strong{I_start_class}: A numerical value, in mm/h. This value is used to create intensity classes. Recommended 0.001 mm/h. If modified, do so with care. It should be positive and not exceed 0.1 mm/h
#'   \item \strong{threshold_int}: A numerical value, in mm/h. A precipitation intensity threshold above which weights are ignored for the scaling models. Recommended 0.002 mm/h. Other accepted values may vary between 0 (no weights are ignored) and 0.1 mm/h
#' }
#' 
#' @return A list, a data frame of scaling parameters and a list of plots for showing empirical estimates and fitted scaling models
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
get_scaling_params_Intensity = function(vecPrecip,
                                        vecDates,
                                        resVecPrecip,
                                        aggLevels,
                                        by_season,
                                        listOptions){
  
  # ========================== List of options =========================
  
  if (is.null(listOptions)){
    listOptions = list()
  } else {
    if (!is.list(listOptions)) 
      stop("listOptions must be a list")
  }
  
  if ("I_min_fix" %in% names(listOptions)) {
    I_min_fix = listOptions[["I_min_fix"]]
    if (!is.numeric(I_min_fix)) 
      stop("I_min_fix must be numeric")
  } else {
    I_min_fix = 0.01
    listOptions[["I_min_fix"]] = I_min_fix
  }
  
  if ("I_max_fix" %in% names(listOptions)) {
    I_max_fix = listOptions[["I_max_fix"]]
    if (!is.numeric(I_max_fix)) 
      stop("I_max_fix must be numeric")
  } else {
    I_max_fix = 7
    listOptions[["I_max_fix"]] = I_max_fix
  }
  
  if ("I_start_class" %in% names(listOptions)) {
    I_start_class = listOptions[["I_start_class"]]
    if (!is.numeric(I_start_class))
      stop("I_start_class must be numeric")
  } else {
    I_start_class = 0.001
    listOptions[["I_start_class"]] = I_start_class
  }
  
  if ("threshold_int" %in% names(listOptions)) {
    threshold_int = listOptions[["threshold_int"]]
    if (!is.numeric(threshold_int)) 
      stop("threshold_int must be numeric")
  } else {
    threshold_int = 0.002
    listOptions[["threshold_int"]] = threshold_int
  }
  
  # ================ Estimate the parameters of the MRC ================
  
  params_MRC = estimate_params_MRC(vecPrecip = vecPrecip, 
                                   vecDates = vecDates,
                                   resVecPrecip = resVecPrecip,
                                   aggLevels = aggLevels,
                                   by_season = by_season,
                                   I_start_class = I_start_class,
                                   threshold_int = threshold_int)
  
  # Some initial parameters for model fitting 
  mu_start = -0.5
  sigma_start = 1
  a_start = -1 
  b_start = -1
  K_start = 0.1
  lambda_start = 0.5
  nu_start = 0.5
  
  # Prepare results 
  fit_intermittency_mu_sigma = NULL
  fit_alpha_intensity_K = NULL
  fit_mean_weights = NULL
  fit_intermittency_ratio_phi = NULL
  
  # Where to predict values 
  vecIntensity = exp(seq(from = log(0.003), to = log(150), length.out = 100))
  
  # ============================= Px model =============================
  
  # Estimate mu and sigma parameters 
  while("try-error" %in% class(fit_intermittency_mu_sigma) | is.null(fit_intermittency_mu_sigma)){
    fit_intermittency_mu_sigma = try(params_MRC$px_aggLevel_Intensity %>% group_by(Season) %>% 
                                       dplyr::do(fit = robustbase::nlrob(value ~ get_Px_Intensity(Intensity = Intensity, mu = b_mu, Sigma = b_sigma), 
                                                                         control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                                                         start = list(b_mu = mu_start, b_sigma = sigma_start), maxit = 2000, data = .data)) %>% 
                                       mutate(b_mu = fit$coefficients[1], b_sigma = fit$coefficients[2]) %>% dplyr::select(-fit), silent = T)
    
    mu_start = sample(x = seq(from = -0.9, to = -0.05, by = 0.2), size = 1)
    sigma_start = sample(x = seq(from = 0, to = 4, by = 0.2), size = 1)
  }
  
  # =========================== Alpha model ===========================

  # Estimate K parameter
  while("try-error" %in% class(fit_alpha_intensity_K) | is.null(fit_alpha_intensity_K)){
    fit_alpha_intensity_K = try(params_MRC$alpha_aggLevel_Intensity %>%
                                  mutate(I_min = I_min_fix, I_max = I_max_fix) %>% group_by(Season) %>%
                                  dplyr::do(fit = robustbase::nlrob(value ~ get_alpha_intensity_1par(Intensity = Intensity, I_min = 0.1, I_max = 10, K = K),
                                                                    control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                                                    start = list(K = K_start), maxit = 2000, data = .data)) %>%
                                  mutate(K = fit$coefficients[1], I_min = I_min_fix, I_max = I_max_fix) %>% dplyr::select(-fit), silent = T)
    
    K_start = sample(x = seq(from = -0.1, to = 0.2, by = 0.2), size = 1)
  }
  
  # ========================= Asymmetry model =========================
  
  # Fit linear model to estimated means for each z index class
  while("try-error" %in% class(fit_mean_weights) | is.null(fit_mean_weights)){
    fit_mean_weights = try(params_MRC$emp_mean_z %>% dplyr::filter(!is.na(zIndex) & !is.na(value)) %>% group_by(Season) %>% 
                             dplyr::do(fit = robustbase::nlrob(value ~ get_mean_z(z = zIndex, lambda = lambda), 
                                                               control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                                               start = list(lambda = lambda_start), maxit = 2000, data = .data)) %>%
                             mutate(lambda = fit$coefficients[1]) %>% dplyr::select(-fit), silent = T)
    
    lambda_start = sample(x = seq(from = 0.1,to = 0.9,by = 0.1), size = 1)
  }
  
  # Fit erf model to estimated ratios for each z index class
  while("try-error" %in% class(fit_intermittency_ratio_phi) | is.null(fit_intermittency_ratio_phi)){
    fit_intermittency_ratio_phi = try(params_MRC$phi_z %>% dplyr::filter(!is.na(zIndex) & !is.na(value)) %>% group_by(Season) %>% 
                                        dplyr::do(fit = robustbase::nlrob(value ~ get_phi_z(z = zIndex, nu = nu),
                                                                          control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                                                          start = list(nu = nu_start), maxit = 2000, data = .data)) %>%
                                        mutate(nu = fit$coefficients[1]) %>% dplyr::select(-fit), silent = T)
    
    nu_start = sample(x = seq(from = 0.1, to = 0.9, by = 0.01), size = 1)
  }
  
  params_asymm = left_join(fit_mean_weights, fit_intermittency_ratio_phi, by = "Season")
  
  # ========================= Plot: Px model ==========================
  
  list_plots = list()
  
  # Predict values based on model parameters
  predict_Px = left_join(expand.grid(Intensity = vecIntensity, Season = fit_intermittency_mu_sigma$Season, aggLevel = aggLevels), fit_intermittency_mu_sigma, by = "Season") %>%
    mutate(Px = get_Px_Intensity(Intensity, mu = b_mu, Sigma = b_sigma))
  
  # Px obs  
  ggplot_px_intensity_time = ggplot(data = params_MRC$px_aggLevel_Intensity)+
    geom_line(data = predict_Px, aes(x = Intensity, y = Px))+
    geom_point(aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_x_continuous(trans = 'log10',
                       limits = range(vecIntensity),
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100"))+
    scale_y_continuous(limits = c(0, 1))+ 
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]")+ 
    labs(x = latex2exp::TeX(r'(Intensity [mm h$^{-1}$])'),
         y = latex2exp::TeX(r'($P_x$)'),
         title = "Px as a function of intensity")+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+
    theme_bw()
  
  list_plots[["Px_intensity"]] = ggplot_px_intensity_time
  
  # ======================== Plot: Alpha model ========================
  
  # Predict alpha based on intensity
  predict_alpha_int = left_join(expand.grid(Intensity = vecIntensity, Season = fit_alpha_intensity_K$Season), fit_alpha_intensity_K, by = "Season") %>%
    mutate(alpha = get_alpha_intensity_1par(Intensity = Intensity, I_min = 0.1, I_max = 10, K = K))
  
  ggplot_alpha_int = ggplot()+
    geom_line(data = predict_alpha_int, aes(x = Intensity, y = alpha))+
    geom_point(data = params_MRC$alpha_aggLevel_Intensity, aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_x_continuous(trans = 'log10',
                       limits = range(vecIntensity),
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100"))+    
    scale_y_continuous(trans = 'log10', limits = c(0.2, 12))+ 
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]") +
    labs(x = latex2exp::TeX(r'(Intensity [mm h$^{-1}$])'),
         y = latex2exp::TeX(r'($\alpha$)'),
         title = "Alpha as a function of intensity")+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+
    theme_bw()
  
  list_plots[["Alpha_intensity"]] = ggplot_alpha_int
  
  # ====================== Plot: Asymmetry model ======================

  predict_phi = left_join(expand.grid(zIndex = seq(0, 1, length.out = 30), Season = params_asymm$Season), params_asymm, by = "Season") %>%
    mutate(phi = get_phi_z(z = zIndex, nu = nu))
  
  predict_mean_weights = left_join(expand.grid(zIndex = seq(0, 1, length.out = 30), Season = params_asymm$Season), params_asymm, by = "Season") %>%
    mutate(m = get_mean_z(z = zIndex, lambda = lambda))
  
  ggplot_asymm_ratio = ggplot()+
    geom_line(data = predict_phi, aes(x = zIndex, y = phi))+
    geom_point(data = params_MRC$phi_z, aes(x = zIndex, y = value))+
    scale_x_continuous(limits = c(0, 1))+
    scale_y_continuous(limits = c(0, 1))+
    labs(x = "z", 
         y = latex2exp::TeX(r'($\varphi$)'),
         title = "Probability asymmetry ratio")+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+ 
    theme_bw()
  
  ggplot_asymm_mean = ggplot()+
    geom_line(data = predict_mean_weights, aes(x = zIndex, y = m))+
    geom_point(data = params_MRC$emp_mean_z, aes(x = zIndex, y = value))+
    scale_x_continuous(limits = c(0, 1))+
    scale_y_continuous(limits = c(0, 1))+
    labs(x = "z", 
         y = "m",
         title = "Mean of the cascade weights")+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+ 
    theme_bw()
  
  list_plots[["Asymm_ratio"]] = ggplot_asymm_ratio
  list_plots[["Asymm_mean"]] = ggplot_asymm_mean
  
  params_all = as.data.frame(left_join(params_asymm, 
                                       left_join(fit_alpha_intensity_K, fit_intermittency_mu_sigma, by = "Season"),
                                       by = "Season"), row.names = NULL)

  params_all <- params_all %>%
    select(-c("I_min", "I_max"))
  
  params_all <- params_all %>%
    rename("mu" = "b_mu",
           "Sigma" = "b_sigma")
  
  return(list(params = params_all, fig_plots = list_plots))
}