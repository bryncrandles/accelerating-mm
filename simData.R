### functions to simulate data used in simulation study

simAllAg = function(d = NULL, G = NULL, similarVec = NULL){
  
  # d : dimension 
  # G : number of groups
  # similarVec: logical vector of length G 
  
  # return 
  # matrix of dimension d by G, where each column has elements diag(Ag^{-1})
  
  Ag = matrix(0, nrow = d, ncol = G)
  for(g in 1:G){
    similar = similarVec[g]
    if(similar == TRUE){
      Ag[ , g] = 1 +  runif(d, -0.05, 0.05)
    } else{
      Ag[, g] = 1 +  c(runif(d/2, 0.05, 1) , -runif(d/2, 0.05, 1))
    }
  }
  return(Ag)
}


Wg_sim = function(d, similar = TRUE){
  # d : dimension 
  # similar: logical 
  
  # details
  # generates Wg with similar evalues if TRUE, otherwise generates Wg with dissimilar eigenvalues 
  
  # return 
  # matrix of dimension d by G, where each column has elements diag(Ag^{-1})
  
  temp = cov(matrix(rnorm((d + 1)*d), nrow = d + 1, ncol = d))
  P = eigen(temp)$vectors
  
  if(similar == TRUE){
    ev = 1 +  runif(d, -0.05, 0.05)
  } else{
    ev = 1 +  c(runif(d/2, 0.05, 1) , -runif(d/2, 0.05, 1))
  }
  out = P %*% diag(ev) %*% t(P)
  return(out)
}

simAllWg = function(d = NULL, G = NULL, similarVec = NULL){
  # d : dimension 
  # G : number of groups
  # similarVec: logical vector of length G 
  
  # return 
  # 3D array with dimensions d by d by G, where each "slice" is Wg 
  
  Wg = array(0, c(d, d, G))
  for (g in 1:G) {
    Wg[ , , g] = Wg_sim(d, similar = similarVec[g])
  }
  return(Wg)
}
