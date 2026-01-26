MMA = function(D, d = NULL, G = NULL, Wg = NULL, Ag = NULL, ...) {
  
  # D : Vector. Length d^2 of d x d matrix of current value of D.
  # p : Scalar. Dimension of D 
  # G : Scalar. Number of groups.
  # Wg : Array. d x d x G 
  # Ag: Matrix of dimension d x G. Each column contains diag(Ag^{-1}). 
  
  # details
  # MMA update from Browne 2013. Adapted from code written by Dr. Ryan Browne.
  
  # return 
  # Dnew : vector. Updated d x d matrix from MMA update. 
  
  
  D = matrix(D, nrow = d, ncol = d)

  Ftemp = rowSums(sapply(1:G, function(g) {
    omega.g = RSpectra::eigs_sym(Wg[ , , g], 1, which = "LM", opts = list(retvec = FALSE))$values 
    (Ag[ , g]) * ( t(D) %*% (Wg[ , , g]) )  - omega.g* (Ag[ , g]) * t(D)
  }))
  dim(Ftemp) = c(d, d)
  
  Ftemp.svd = svd(-Ftemp)
  Dnew  = tcrossprod(Ftemp.svd$v, Ftemp.svd$u)
  return(as.vector(Dnew))
  
}

MMW = function(D = NULL, d = NULL, G = NULL, Wg = NULL, Ag = NULL, ...) {
  
  # D : Vector of length d^2 of vec(D).
  # d : Scalar. Dimension of D.
  # G : Scalar. Number of groups.
  # Wg : Array. d x d x G 
  # Ag: Matrix of dimension d x G. Each column contains diag(Ag^{-1}). 
  # useGradU: Logical. TRUE uses unconstrained gradient, FALSE uses constrained gradient.
  
  # Details
  # MMW update from Browne 2013. Adapted from code written by Dr. Ryan Browne.
  
  # return 
  # Dnew : Vector of length d^2 of updated vec(D) from MMW update. 
  
  D = matrix(D, nrow = d, ncol = d)
  
  Ftemp = rowSums(sapply(1:G, function(g) {
    alpha.g = max(Ag[, g])
    (Ag[ , g]) * ( t(D) %*% (Wg[ , , g]) )  - alpha.g * ( t(D) %*% Wg[ , , g] )
  }))
  dim(Ftemp) = c(d, d)
  
  Ftemp.svd = svd(-Ftemp)
  Dnew = tcrossprod(Ftemp.svd$v, Ftemp.svd$u)
  return(as.vector(Dnew))
}


MMAW = function(D = NULL, d = NULL, G = NULL, Wg = NULL, Ag = NULL, iter = NULL, ...) {
  
  # D : Vector of length d^2 of vec(D).
  # p : Scalar. Dimension of D 
  # G : Scalar. Number of groups.
  # Wg : Array. d x d x G 
  # Ag: Matrix of dimension d x G. Each column contains diag(Ag^{-1}). 
  
  # details
  # MMAW update from Browne 2013. Adapted from code written by Dr. Ryan Browne.
  
  # return 
  # Dnew : Vector of length d^2 of updated vec(D) from MM3 update.
  
  
  if (iter %% 2) {
    Dnew = MMA(D = D, d = d, G = G, Wg = Wg, Ag = Ag)
  } else {
    Dnew = MMW(D = D, d = d, G = G, Wg = Wg, Ag = Ag)
  }
  return(Dnew)
}

MMlambda = function(D = NULL, d = NULL, G = NULL, Wg = NULL, Ag = NULL, fObj, ...) {
  
  # D : Vector of length d^2 of vec(D).
  # d : Scalar. Dimension of D 
  # G : Scalar. Number of groups.
  # Wg : Array. d x d x G 
  # Ag: Matrix of dimension d x G. Each column contains diag(Ag^{-1}). 
  
  # details
  # MMlambda update. Adapted from code written by Dr. Ryan Browne.
  
  # return 
  # Dnew : Vector of length d^2 of updated vec(D) from MMlambda update.
  
  
  Dnew1 = MMA(D = D, d = d, G = G, Wg = Wg, Ag = Ag)
  Dnew2 = MMW(D = D, d = d, G = G, Wg = Wg, Ag = Ag)
  if(fObj(D = Dnew1, Wg = Wg, Ag = Ag, G = G, d = d) <  fObj(D = Dnew2, Wg = Wg, Ag = Ag, G = G, d = d)){
    Dnew = Dnew1
  } else{
    Dnew = Dnew2
  }
  
  return(Dnew)
  
}

MMdelta = function(D = NULL, d = NULL, G = NULL, Wg = NULL, Ag = NULL, delta1 = NULL, delta2 = NULL, ...) {
  
  # D : Vector of length d^2 of vec(D).
  # p : Scalar. Dimension of D 
  # G : Scalar. Number of groups.
  # Wg : Array. d x d x G 
  # Ag: Matrix of dimension d x G. Each column contains diag(Ag^{-1}). 
  
  # details
  # MMdelta update. Adapted from code written by Dr. Ryan Browne.
  
  # return
  # Dnew : Vector of length d^2 of updated vec(D) from MMdelta update.
  
  if(delta1 > delta2){
    MM = "MMA"
    Dnew = MMA(D = D, d = d, G = G, Wg = Wg, Ag = Ag)
    
  } else{
    MM = "MMW"
    Dnew = MMW(D = D, d = d, G = G, Wg = Wg, Ag = Ag)
    }
  
  out = list(Dnew = Dnew, MM = MM)
  return(out)
  
}

