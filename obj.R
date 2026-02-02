### calculate scaled objective function

obj = function(D, Wg = NULL, Ag = NULL, G = NULL, d = NULL, ...) {
  D = matrix(D, nrow = d, ncol = d)
  tD = t(D)
  z = rowSums(sapply(1:G, function(g) { 
    ( (Ag[ , g] * tD ) %*% Wg[ , , g] )
  }))
  val = sum( z * as.numeric( tD ) )
  return(val/(d*G)) ### edited Aug 20 to divide by dimension and number of groups 
}
