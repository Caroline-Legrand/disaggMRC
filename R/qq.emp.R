#' qq.emp
#' 
#' @description A function that calculates an empirical quantile based on a Gumbel transformation
#' 
#' @param x A vector of data
#' @param p A numerical value, the probability associated with the quantile to be calculated
#' @param a Parameter associated with the Gumbel distribution for probability fitting, based on a Gringorten approach
#' 
#' @return A numerical value, the empirical quantile
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#'
#' @export
qq.emp = function(x, p, a=0.44){
  # Remove missing values (NA)
  x = na.omit(x)
  
  # If the vector is too short or empty after removing NAs, return NA
  if(length(x) <= 2 | all(is.na(x))){
    x.out = NA
    
  } else {
    # Sort the data
    x.srt = sort(x)
    
    # Compute the empirical probabilities using L-moments
    p.srt = lmomco::pp(x.srt, a=a)
    
    # Apply the Gumbel transformation to the empirical probabilities
    u.srt = Q_Gumbel(p.srt)
    
    # Apply the Gumbel transformation to the provided probability
    u.out = Q_Gumbel(p)
    
    # Perform interpolation to get the empirical quantile corresponding to the given probability
    x.out = approx(x=u.srt, y=x.srt, xout=u.out, rule=2)$y
  }
  
  # Return the empirical quantile
  return(x.out)
}