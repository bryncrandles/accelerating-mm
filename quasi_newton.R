#' Quasi Newton Method for Accelerating Slowly-Convergent Fixed-Point Iterations
#'
#' An implementation of quasi Newton method described in Zhou (2011). Including monotonicity control and projection.
#'
#' @param par Vector for initial parameters
#' @param fixptfn Fixed point updating function
#' @param objfn Objective function
#' @param control A list containing parameters controlling the algorithm
#' @param ... Other arguments required by \code{fixptfn} and \code{objfn}
#'
#' @details The task it to \strong{minimize} \code{objfn}. Default values of \code{control} are: \code{qn=3, mono.tol=1, projection=function(x) x, tol=1e-7, maxiter=2000, convtype="parameter", par.track=FALSE, conv.spec=NULL}.
#' \describe{
#'  \item{qn}{An integer variable indicating the order of Quasi-Newton algorithm used. Default is 3.}
#'  \item{mono.tol}{A non-negative scalar that dictates the degree of non-montonicity. Default is 1. Set objfn.inc = 0 to obtain monotone convergence. Setting objfn.inc = Inf gives a non-monotone scheme. In-between values result in partially-monotone convergence.}
#'  \item{projection}{A function projecting the parameter after each iteration. Default is identity function \eqn{f(x) = x}}
#'  \item{tol}{A small, positive scalar that determines when iterations should be terminated, see \code{convtype} for details. Default is \code{1e-7}}
#'  \item{maxiter}{An integer denoting the maximum limit on the number of evaluations of \code{fixptfn}. Default is 2000.}
#'  \item{convtype}{A string indicating the convergence criteria.
#'                 If it is "parameter", the algorithm will termenate when L2 norm of parameters difference \eqn{x_{new} - x_{old} < tol}.
#'                 If it is "objfn", the algorithm will terminate when the absolute difference of objective function \eqn{|L_{new} - L_{old}| < tol}.
#'                 If it is "user" or \code{conv.spec} is not \code{NULL}. Then the convergence is guided by the user defined function \code{conv.spec}.
#'                 Default is "parameter".}
#'  \item{par.track}{An bool value indicating whether to track parameters along the algorithm. \code{TRUE} for tracking and \code{FALSE} for not. Default is \code{FALSE}}
#'  \item{conv.spec}{A function for user specified convergence criteria. When using "parameter" or "objfn" option in \code{convtype}, this should be \code{NULL}.
#'                  The function should have the form \code{f(old_parameter, new_parameter, old_objective, new_objective, tolerance)} and return 1 if convergent, 0 if not.
#'                  Defalut is \code{NULL}.}
#' }
#'
#' @return A list of results
#'  \item{par}{Parameter values, x* that are the fixed-point of fixptfn F such that x*=F(x*) if convergence is successful.}
#'  \item{value.objfn}{The objective function value at termination.}
#'  \item{fpevals}{Number of times the fixed-point function \code{fixptfn} was evaluated.}
#'  \item{objfevals}{Number of times the objective function \code{objfn} was evaluated.}
#'  \item{iter}{Numbers of iteration used at termination. (for different algorithms, multiple fixed point iteration might be evaluated in one iteration)}
#'  \item{convergence}{An integer code indicating whether the algorithm converges. 1 for convergence and 0 denote failure.}
#'  \item{objfn.track}{An array tracking objective function values along the algorithm}
#'  \item{par.track}{A matrix tracking parameters along the algorithm, where each row is an array of parameters at some iteration. If not tracking paramters, this will be \code{NULL}}
#'
#' @references Zhou H, Alexander D, Lange K (2011). A quasi-Newton acceleration for high-dimensional optimization algorithms. Statistics and Computing, 21(2): 261–273.
#'
#' @examples
#' \dontrun{
#' set.seed(54321)
#' prob = lasso_task(lam=1)
#' quasi_newton(prob$initfn(), prob$fixptfn, prob$objfn, X=prob$X, y=prob$y)
#' }
#'
#' @export quasi_newton
quasi_newton = function(par, fixptfn_name, objfn, ..., control=list()){
  # fixptfn should return feasible par and no need for projection
  # minimizing task
  control.default <- list(
    convtype="parameter", tol=1.0e-07,
    maxiter=2000, par.track=FALSE,
    qn=3, projection=function(x) x,
    mono.tol=0.1, conv.spec=NULL
  )
  control.sub = control[names(control) %in% names(control.default)]
  ctrl = modifyList(control.default, control.sub)
  
  convergence = TRUE
  fpevals = 0
  objfevals = 0
  objfn.track = c()
  p0 = par
  
  par.track = c(par)
  
  QN = ctrl$qn
  convtype=ctrl$convtype
  convf <- ctrl$conv.spec
  if(!is.null(convf)) convtype="user"
  
  if(fixptfn_name == "MMdelta"){
    delta1 <- Inf
    delta2 <- Inf
    MM.track = numeric(ctrl$maxiter*3)
  } else{
    MM.track = NULL
  }
  
  # set up
  Y <- par
  U <- V <- matrix(0, length(par), ctrl$qn)
  
  iter = 1
  
  for(i in 1:QN){
    Y_old = Y
    # Y = fixptfn(Y_old, ...)
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        # fixptfn = "MMA"
        Y = MMA(D = Y_old, ...)
      } else{
        MM = "MMW"
        Y = MMW(D = Y_old, ...)
      }
      fpevals = fpevals + 1
      MM.track[fpevals] = MM
    } else if(fixptfn_name == "MMAW"){
      if (iter %% 2) {
        Y = MMA(D = Y_old, ...)
      } else {
        Y = MMW(D = Y_old, ...)
      }
      fpevals = fpevals + 1
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      Y = fixptfn(Y_old, ...)
      fpevals = fpevals + 1
    }
    # fpevals = fpevals + 1
    U[ , i] = Y - Y_old
  }
  
  
  V[, 1:(QN-1)] = U[, 2:QN]
  Y_old = Y
  # Y = fixptfn(Y_old, ...)
  if(fixptfn_name == "MMdelta"){
    if(delta1 > delta2){
      MM = "MMA"
      # fixptfn = "MMA"
      Y = MMA(D = Y_old, ...)
    } else{
      MM = "MMW"
      Y = MMW(D = Y_old, ...)
    }
    fpevals = fpevals + 1
    MM.track[fpevals] = MM
  } else if(fixptfn_name == "MMAW"){
    if (iter %% 2) {
      Y = MMA(D = Y_old, ...)
    } else {
      Y = MMW(D = Y_old, ...)
    }
    fpevals = fpevals + 1
  } else{
    if(fixptfn_name == "MMA"){
      fixptfn = MMA
    } else if(fixptfn_name == "MMW"){
      fixptfn = MMW
    }
    Y = fixptfn(Y_old, ...)
    fpevals = fpevals + 1
  }
  
  # fpevals = fpevals + 1
  V[, QN] = Y - Y_old
  
  if (ctrl$par.track) par.track = rbind(par.track, Y)
  Lold = try(objfn(p0, ...), silent=TRUE)
  Lnew = try(objfn(Y, ...), silent=TRUE)
  if(fixptfn_name == "MMdelta"){
    if(MM == "MMA"){
      delta1 = Lold - Lnew
    } else if(MM == "MMW"){
      delta2 = Lold - Lnew
    }
  }
  Lold = Lnew
  
  objfevals = objfevals + 1
  objfn.track = c(objfn.track, Lold)
  
  it = 0
  iter = 2
  
  
  
  while(fpevals <= ctrl$maxiter){
    Y_orig = Y
    Y_old = Y
    # Y = fixptfn(Y_old, ...)
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        # fixptfn = "MMA"
        Y = MMA(D = Y_old, ...)
      } else{
        MM = "MMW"
        Y = MMW(D = Y_old, ...)
      }
      fpevals = fpevals + 1
      MM.track[fpevals] = MM
    } else if(fixptfn_name == "MMAW"){
      if (iter %% 2) {
        Y = MMA(D = Y_old, ...)
      } else {
        Y = MMW(D = Y_old, ...)
      }
      fpevals = fpevals + 1
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      Y = fixptfn(Y_old, ...)
      fpevals = fpevals + 1
    }
    # fpevals = fpevals + 1
    # do not support QN = 0
    
    U[, 1:(QN-1)] = U[, 2:QN]
    U[, QN] = Y - Y_old
    
    Y_old = Y
    # Y = fixptfn(Y_old, ...) ## F(F(x))
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        # fixptfn = "MMA"
        Y = MMA(D = Y_old, ...)
      } else{
        MM = "MMW"
        Y = MMW(D = Y_old, ...)
      }
      fpevals = fpevals + 1
      MM.track[fpevals] = MM
    } else if(fixptfn_name == "MMAW"){
      if (iter %% 2) {
        Y = MMA(D = Y_old, ...)
      } else {
        Y = MMW(D = Y_old, ...)
      }
      fpevals = fpevals + 1
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      Y = fixptfn(Y_old, ...)
      fpevals = fpevals + 1
    }
    
    # fpevals = fpevals + 1
    
    V[, 1:(QN-1)] = V[, 2:QN]
    V[, QN] = Y - Y_old
    
    AA = crossprod(U, U - V)
    
    utdiff = crossprod(U, as.vector(U[, QN]))
    AinvU = try(solve(AA, utdiff), silent=TRUE)
    
    if(!inherits(AinvU, "try-error")){
      Y_QN = ctrl$projection( as.vector(Y_old + V %*% AinvU) )
    } else{
      Y_QN = Y
    }
    Lnew = try(objfn(Y_QN, ...), silent=TRUE)
    objfevals = objfevals + 1
    if(inherits(Lnew, "try-error") || !is.finite(Lnew) || Lnew > Lold + ctrl$mono.tol){
      Lnew = try(objfn(Y, ...), silent=TRUE)
      objfevals = objfevals + 1
    } else{
      Y = Y_QN
    }
    
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        delta1 = Lold - Lnew
      } else if(MM == "MMW"){
        delta2 = Lold - Lnew
      }
    }
    
    objfn.track = c(objfn.track, Lnew)
    
    if(convtype == "objfn" && is.finite(abs(Lnew - Lold)) && abs(Lnew - Lold) < ctrl$tol)
      break
    if(convtype == "parameter" && is.finite(c(crossprod(Y - Y_orig))) && sqrt(c(crossprod(Y - Y_orig))) < ctrl$tol)
      break
    if(convtype == "user")
      if (convf(Y_orig, Y, Lold, Lnew, ctrl$tol))
        break
    
    Lold = Lnew
    
    if (ctrl$par.track) par.track = rbind(par.track, Y)
    
    it = it + 1
    iter = iter + 1
  }
  
  if(fpevals > ctrl$maxiter){
    convergence = FALSE
    # warning("Algorithm did not converge")
  }
  
  rownames(par.track) = NULL
  
  # list(par = c(Y),
  #      value.objfn = Lold,
  #      iter = it,
  #      fpevals = fpevals,
  #      objfevals = objfevals,
  #      convergence = convergence,
  #      objfn.track = objfn.track,
  #      par.track = par.track)
  
  list(par = c(Y),
       value.objfn = Lold,
       iter = it,
       fpevals = fpevals,
       objfevals = objfevals,
       convergence = convergence,
       objfn.track = objfn.track,
       par.track = par.track, MM.track = MM.track[1:fpevals])
}
