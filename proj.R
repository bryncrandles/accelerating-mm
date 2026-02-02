### projections considered in paper

projection.SVD = function(D){
  D = matrix(D, nrow = sqrt(length(D)), ncol = sqrt(length(D)))
  D.SVD = svd(D)
  D.proj = tcrossprod(D.SVD$u, D.SVD$v)
  return(as.vector(D.proj))
}

projection.QR = function(D){
  D = matrix(D, nrow = sqrt(length(D)), ncol = sqrt(length(D)))
  as.vector(QR_Rplus(D)$Q)
}