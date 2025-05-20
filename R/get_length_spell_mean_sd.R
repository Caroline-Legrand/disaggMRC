#' get_length_spell_mean_sd
#' 
#' @description A function that calculates the mean and standard deviation of the length of dry and wet spells
#' 
#' @param x A vector of precipitation, in mm
#' 
#' @return A list containing the mean and standard deviation of the length of dry and wet spells
#' 
#' @author Kaltrina Maloku
#' 
#' @references Maloku, K., Hingray, B., and Evin, G. (2023).
#' Accounting for precipitation asymmetry in a multiplicative random cascade disaggregation model.
#' Hydrol. Earth Syst. Sci., 27, 3643–3661, https://doi.org/10.5194/hess-27-3643-2023
#'
#' @export
get_length_spell_mean_sd = function(x){
  # Transform x into a Boolean indicating dry state (= 0) and wet state
  x.states = x>0
  
  # Run length encoding
  rle.x = rle(x.states)
  
  # Dry spells
  is.dry = rle.x$values==F
  dry.len = rle.x$lengths[is.dry]
  dry.freq = table(dry.len)
  
  mean.dry = sum(as.numeric(names(dry.freq))*dry.freq)/sum(dry.freq)
  sd.dry = sqrt((1/(sum(dry.freq)-1))*sum(dry.freq*(as.numeric(names(dry.freq))-mean.dry)^2))
  
  # Wet spells
  is.wet = rle.x$values==T
  wet.len = rle.x$lengths[is.wet]
  wet.freq = table(wet.len)
  
  mean.wet = sum(as.numeric(names(wet.freq))*wet.freq)/sum(wet.freq)
  sd.wet = sqrt((1/(sum(wet.freq)-1))*sum(wet.freq*(as.numeric(names(wet.freq))-mean.wet)^2))
  
  return(list(mean.dry = mean.dry, mean.wet = mean.wet, sd.dry = sd.dry, sd.wet = sd.wet))
}