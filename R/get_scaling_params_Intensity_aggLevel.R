#' get_scaling_params_Intensity_aggLevel
#' 
#' @description A function that estimates the scaling parameters of a MRC model, where the cascade weights are considered to depend on the temporal aggregation level and the intensity of the precipitation to be disaggregated.
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
get_scaling_params_Intensity_aggLevel = function(vecPrecip,
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
  lambda_start = 0.5
  nu_start = 0.5
  
  # Where to predict values 
  vecIntensity = exp(seq(from = log(0.003), to = log(150), length.out = 100))
  
  # ============================= Px model =============================
  
  # First step: estimate mu and sigma parameters for each temporal aggregation level
  fit_intermittency_mu_sigma = params_MRC$px_aggLevel_Intensity %>%
    group_by(Season, aggLevel) %>%
    dplyr::do(fit = robustbase::nlrob(value ~ get_Px_Intensity(Intensity = Intensity, mu = mu, Sigma = Sigma), 
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                      start = list(mu = mu_start, Sigma = sigma_start), maxit = 2000, data = .data)) %>%
    mutate(mu = fit$coefficients[1], Sigma = fit$coefficients[2]) %>% dplyr::select(-fit)
  
  # Second step: estimate mu and sigma scaling parameters over temporal aggregation levels, fit linear models
  fit_intermittency_mu = fit_intermittency_mu_sigma %>%
    group_by(Season) %>% dplyr::select(-Sigma) %>%
    dplyr::do(fit = robustbase::nlrob(mu ~ (a_mu*log(aggLevel/1440)+b_mu),
                                      start = list(a_mu = a_start, b_mu = b_start), maxit = 2000, data = .data)) %>%  
    mutate(a_mu = fit$coefficients[1], b_mu = fit$coefficients[2]) %>% dplyr::select(-fit)
  
  fit_intermittency_sigma = fit_intermittency_mu_sigma %>% 
    group_by(Season) %>% dplyr::select(-mu) %>%
    dplyr::do(fit = robustbase::nlrob(Sigma ~ (a_sigma*log(aggLevel/1440)+b_sigma), 
                                      start = list(a_sigma = a_start, b_sigma = b_start), maxit = 2000, data = .data)) %>%
    mutate(a_sigma = fit$coefficients[1], b_sigma = fit$coefficients[2]) %>%  dplyr::select(-fit)
  
  params_px = left_join(fit_intermittency_mu, fit_intermittency_sigma, by = "Season")
  
  # =========================== Alpha model ===========================
  
  # Fit scaling model over temporal aggregation levels
  fit_alpha_aggLevel = params_MRC$alpha_aggLevel %>% group_by(Season) %>%
    dplyr::do(fit = robustbase::nlrob(value ~ get_alpha_aggLevel(alpha0 = alpha0, H = H, res_aggLevel = aggLevel/1440),
                                      start = list(alpha0 = 1, H = -1), maxit = 2000, data = .data)) %>%
    mutate(alpha0 = fit$coefficients[1], H = fit$coefficients[2]) %>% dplyr::select(-fit)
  
  # Fit quadratic model to observed alpha star
  fit_quad_model = params_MRC$alpha_star_aggLevel_Intensity %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(log(value) ~ (c0+c1*log(Intensity)+c2*log(Intensity)*log(Intensity)),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                      start = list(c0 = 0.1, c1 = 2, c2 = 1), maxit = 1000, data = .data)) %>% 
    mutate(c0 = fit$coefficients[1], c1 = fit$coefficients[2], c2 = fit$coefficients[3]) %>% dplyr::select(-fit)
  
  params_alpha = left_join(fit_alpha_aggLevel, fit_quad_model, by = "Season")
  
  # ========================= Asymmetry model =========================
  
  # Fit linear model to estimated means for each z index class
  fit_mean_weights = params_MRC$emp_mean_z %>% dplyr::filter(!is.na(zIndex)) %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value ~ get_mean_z(z = zIndex, lambda = lambda),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                      start = list(lambda = lambda_start), maxit = 2000, data = .data)) %>%
    mutate(lambda = fit$coefficients[1]) %>% dplyr::select(-fit)

  # Fit erf model to estimated ratios for each z index class
  fit_intermittency_ratio_phi = params_MRC$phi_z %>% dplyr::filter(!is.na(zIndex)) %>% group_by(Season) %>% 
    dplyr::do(fit = robustbase::nlrob(value ~ get_phi_z(z = zIndex, nu = nu),
                                      control = list(maxiter = 1000, printEval = F, minFactor = 1/100000),
                                      start = list(nu = nu_start), maxit = 2000, data = .data)) %>%
    mutate(nu = fit$coefficients[1]) %>% dplyr::select(-fit)
  
  params_asymm = left_join(fit_mean_weights, fit_intermittency_ratio_phi, by = "Season")
  
  # ========================= Plot: Px model ==========================

  list_plots = list()

  # Predict values based on model parameters
  predict_Px = left_join(expand.grid(Intensity = vecIntensity, Season = params_px$Season, aggLevel = aggLevels), 
                         params_px, by = "Season") %>%
    mutate(Px = get_Px_Intensity_aggLevel(Intensity, aggLevel/1440, a_mu, b_mu, a_sigma, b_sigma))
  
  # Px obs
  ggplot_px_intensity_time = ggplot(data = params_MRC$px_aggLevel_Intensity)+
    geom_line(data = predict_Px, aes(x = Intensity, y = Px, color = factor(aggLevel)))+
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
         title = "Px as a function of temporal aggregation level and intensity")+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+
    theme_bw()
  
  list_plots[["Px_aggLevel_intensity"]] = ggplot_px_intensity_time
  
  int_pars_px = left_join(fit_intermittency_mu_sigma, params_px, by = "Season") %>% 
    mutate(pred_mu = a_mu*log(aggLevel/1440)+b_mu, pred_sigma = a_sigma*log(aggLevel/1440)+b_sigma)
  
  ggplot_px_params_aggLevel = ggplot(int_pars_px)+
    geom_line(aes(x = aggLevel, y = pred_mu))+
    geom_point(aes(x = aggLevel, y = mu, color = "mu"), fill = "black", shape = 15)+
    geom_line(aes(x = aggLevel, y = pred_sigma))+
    geom_point(aes(x = aggLevel, y = Sigma, color = "sigma"), shape = 17)+
    scale_x_continuous(trans = 'log10', breaks = aggLevels)+
    scale_color_manual(name = "", breaks = c("mu", "sigma"), values=c("black", "black"))+
    labs(x = "Temporal aggregation level", 
         y = "Value",
         title = "Px parameters as a function of temporal scale")+
    guides(colour = guide_legend(override.aes = list(shape = c(15, 17))))+
    facet_wrap(~Season)+
    theme_bw()
    
  list_plots[["Px_scaling_params"]] = ggplot_px_params_aggLevel
  
  # ======================== Plot: Alpha model ========================
  
  # Predict alpha star values by the model 
  predict_alpha = left_join(expand.grid(aggLevel = aggLevels, Season = params_alpha$Season), params_alpha, by = "Season") %>%
    mutate(alpha = get_alpha_aggLevel(alpha0, H, res_aggLevel = aggLevel/1440))
  
  # Alpha as a function of temporal scale
  ggplot_alpha_ts = ggplot(data = params_MRC$alpha_aggLevel)+
    geom_line(data = predict_alpha, aes(x = aggLevel, y = alpha))+
    geom_point(aes(x = aggLevel, y = value))+
    scale_x_continuous(trans = 'log10', breaks = aggLevels)+
    scale_y_continuous(trans = 'log10', limits = c(0.2, 12))+
    labs(x = "Temporal aggregation level", 
         y = latex2exp::TeX(r'($h(\tau)$)'))+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+
    theme_bw()
  
  list_plots[["Alpha_aggLevel"]] = ggplot_alpha_ts
  
  # Alpha temporal scale and intensity class 
  ggplot_alpha_ts_int = ggplot()+
    geom_point(data = params_MRC$alpha_aggLevel_Intensity, aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_x_continuous(trans = 'log10',
                       limits = range(vecIntensity),
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100"))+
    scale_y_continuous(trans = 'log10', limits = c(0.2, 12))+
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]")+ 
    labs(x = latex2exp::TeX(r'(Intensity [mm h$^{-1}$])'), 
         y = latex2exp::TeX(r'($\alpha$)'))+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+
    theme_bw()
  
  list_plots[["Alpha_aggLevel_intensity"]] = ggplot_alpha_ts_int
  
  # Predict alpha star values by the model 
  predict_alpha_star = left_join(expand.grid(Intensity = vecIntensity, Season = params_alpha$Season), params_alpha, by = "Season") %>%
    mutate(alpha_star = alpha_star_Intensity(Intensity, c0, c1, c2))
  
  # Alpha star
  ggplot_alpha_star = ggplot()+
    geom_line(data = predict_alpha_star, aes(x = Intensity, y = alpha_star))+
    geom_point(data = params_MRC$alpha_star_aggLevel_Intensity, aes(x = Intensity, y = value, color = factor(aggLevel)))+
    scale_x_continuous(trans = 'log10',
                       limits = range(vecIntensity),
                       breaks = c(0.01, 0.1, 1, 10, 100),
                       labels = c("0.01", "0.1", "1", "10", "100"))+
    scale_y_continuous(trans = 'log10', limits = c(0.2, 12))+ 
    scale_color_manual(values = get_col_aggLevel(aggLevels),
                       breaks = aggLevels,
                       name = "T.scale [min]")+ 
    labs(x = latex2exp::TeX(r'(Intensity [mm h$^{-1}$])'),
         y = latex2exp::TeX(r'(g(I))'))+
    facet_wrap(~Season, labeller = labeller(Season = c("1" = "DJF", "2" = "MAM", "3" = "JJA", "4" = "SON")))+
    theme_bw()
  
  list_plots[["Alpha_star_aggLevel_intensity"]] = ggplot_alpha_star
  
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
                                       left_join(params_alpha, params_px, by = "Season"),
                                       by = "Season"), row.names = NULL)
  
  return(list(params = params_all, fig_plots = list_plots))
}