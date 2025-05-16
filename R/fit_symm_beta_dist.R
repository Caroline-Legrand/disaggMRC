#' fit_symm_beta_dist
#' 
#' @description A function that fits a symmetric Beta distribution, with the probability density function:
#' \deqn{f(x, \alpha) = \frac{1}{B(\alpha, \alpha)} {x}^{(\alpha-1)}(1-x)^{(\alpha-1)}}
#' 
#' @param x A vector of weights
#' 
#' @return The estimated \eqn{\alpha} parameter
#' 
#' @details If sample size is less than 10, return NA
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#' 
#' @export
fit_symm_beta_dist = function(x){
  sample_length = length(x)
  
  if (sample_length < 10){ # If sample size is less than 10, return NA (do not estimate alpha)
    alpha_par = NA
  } else { # Calculate the variance
    var_pw = stats::var(x, na.rm = T)
    
    # If all data have the same value, return NA
    var_pw[var_pw == 0] = NA
    
    # Estimate alpha
    alpha_par = 1/(8*var_pw)-0.5
    
    # If sample size is less than 20, perform a goodness of fit test
    if (sample_length < 20 & sample_length > 10 & !is.na(alpha_par)){
      probabilities = (1:sample_length)/(sample_length+1)
      pred = stats::qbeta(probabilities, shape1 = alpha_par, shape2 = alpha_par)
      k = stats::ks.test(x, pred)
      
      if (k$p.value <= 0.1){
        print(k)
        print(x)
        alpha_par = NA
      }
    }
  }
  
  return(alpha_par)
}