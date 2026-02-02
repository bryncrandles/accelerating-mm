### QR decomposition where R has positive diagonal entries 

QR_Rplus = function(X) {
  QR1 = qr(X)
  Q = qr.Q(QR1, TRUE)
  R = qr.R(QR1, TRUE)
  Q = Q %*% diag(sign(diag(R)))
  R = diag(sign(diag(R))) %*% R
  return(list(Q = Q, R = R)) 
}
