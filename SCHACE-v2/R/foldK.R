foldK <- function(n,K){
  reminder <- n%%K
  if (reminder == 0) return(rep(c(1:K),n/K))
  return(c(rep(c(1:K),floor(n/K)),1:reminder))
}