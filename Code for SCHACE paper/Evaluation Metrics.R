
# informative set

Informative = function(changepoints, true_CP, range){
  
  normal_cp = changepoints/range
  normal_true_cp = true_CP/range
  
  tau = c(0, normal_cp, 1)
  k = length(tau)
  middle_CPs = (tau[1:(k-1)] + tau[2:k]) / 2
  
  infor_points = normal_cp[findInterval(normal_true_cp, middle_CPs)]
  #Given a vector of non-decreasing breakpoints in middle_CPs, find the interval containing each element of normal_true_cp
  Uninfor_points = normal_cp[!normal_cp %in% infor_points]
  
  in_points = changepoints[which(normal_cp %in% infor_points)]
  Un_points = changepoints[!changepoints %in% in_points]
  
  result = list("infor" = in_points, "uninfor" = Un_points)
  
  return(result)
}

# True Positives or not

catch = function(sample_x, informative, true_CP, h){
  
  number_info = length(informative)
  number_trueCP = length(true_CP)
  Estimated_CPs_location = rep(NA, length(true_CP))
  
  if(number_info == 0){
    
    number_catched = 0
    
  } else {
    
    number_catched = 0
    
    for(i in 1:number_info){
      
      #i = 2
      distance = abs(informative[i] - true_CP)
      catched = distance <= h 
      catched[-which.min(distance)] = FALSE
      
      number_catched = number_catched + sum(catched)
      
      location = which(sample_x %in% informative[i])
      Estimated_CPs_location[catched] = sample_x[location]
      
      true_CP[catched] = 1000000
      
    }
    
  }
  
  results = list("Number_of_Catched_CPs" = number_catched, "Estimated_CPs_location" = Estimated_CPs_location)
  return(results)
  
}
