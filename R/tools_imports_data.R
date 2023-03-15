#' World Health Organization TB data
#'
#'A fictive 10 minute resolution time series of precipitation 
#'
#' @format ## `data_test_MRC`
#' A data frame of two columns and 1577952 rows:
#' \describe{
#'   \item{obs}{Precipitation amounts}
#'   \item{date}{Dates}
#' }
"data_test_MRC"


#' @importFrom stats runif rbeta
#' @importFrom dplyr group_by mutate summarise_at select rename filter
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes labs guides guide_legend ggplot geom_line geom_point scale_x_continuous scale_y_continuous scale_color_manual facet_wrap theme_bw
#' @importFrom magrittr %>% 
#' 
#' @export
NULL

globalVariables(c(".data", "H", "Intensity", "K", "P01", "P10", "Precip", "Px", "Season", "Sigma", "Weight", "a_mu",
                  "a_sigma", "aggLevel", "alpha", "alpha0", "alpha_star", "alpha_t", "b_mu", "b_sigma", "bdc", "c0",
                  "c1", "c2", "fit", "lambda", "m", "mu", "nu", "ones", "p01", "p10", "param", "phi", "pred_mu", "pred_sigma",
                  "rnd_unif", "small_int", "value", "zIndex", "zeros"))


