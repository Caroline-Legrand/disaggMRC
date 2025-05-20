#' compute_score
#' 
#' @description A function that calculates the score as assigned by the CASE framework
#' 
#' @param sim.mean A numerical value, the mean among simulated scenarios for a given statistic
#' @param sim.sd A numerical value, the standard deviation among simulated scenarios for a given statistic
#' @param q05 A numerical value, the 5th percentile among simulated scenarios for a given statistic
#' @param q95 A numerical value, the 95th percentile among simulated scenarios for a given statistic
#' @param stat.obs A numerical value, the statistic calculated from the observations
#' 
#' @return A numerical value, the score as assigned by the CASE framework
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#'
#' @export
compute.score = function(sim.mean, sim.sd, q05, q95, stat.obs){
  if (stat.obs >= q05 & stat.obs <= q95){
    score = 1
  } else if (stat.obs >= (sim.mean-3*sim.sd) & stat.obs <= (sim.mean+3*sim.sd)){
    score = 2
  } else if (abs(stat.obs-sim.mean)/stat.obs <= 0.05){
    score = 3
  } else{
    score = 4
  }
  return(score)
}