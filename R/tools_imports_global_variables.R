#' @importFrom stats runif rbeta
#' @importFrom dplyr group_by mutate summarise_at select rename filter
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes labs guides guide_legend ggplot geom_line geom_point scale_x_continuous scale_y_continuous scale_color_manual facet_wrap theme_bw
#' @importFrom magrittr %>%
#' 
#' @export
NULL

globalVariables(c(".data", "H", "I_max", "I_min", "Intensity", "Intensity_min", "K", "P01", "P10", "Precip", "Px", "Season", "Sigma", "Weight", "a_mu", 
                  "a_sigma", "aggLevel", "alpha", "alpha0", "alpha_star", "alpha_t", "b_mu", "b_sigma", "bdc", "c0", "c1", "c2", "fit", "labeller", "lambda", "m", 
                  "mu", "n.y", "nu", "ones", "p01", "p10", "param", "params_px", "phi", "pred_mu", "pred_sigma", "recode", "rnd_unif", "value", "zIndex", "zeros"))