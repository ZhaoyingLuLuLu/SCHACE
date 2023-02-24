meansquared <- function(test.y, pred.y){
  
  diff <- (test.y - pred.y)^2
  return(diff)
  
}

absolute <- function(test.y, pred.y){
  
  diff <- abs(test.y - pred.y)
  return(diff)
  
}

trimmedmse <- function(test.y, pred.y, percent){
  
  diff <- tmspe(test.y, pred.y, trim = percent)
  return(diff)
  
}