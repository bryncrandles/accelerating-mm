#' Restart Nesterov Method for Accelerating Slowly-Convergent Fixed-Point Iterations
#'
#' Using restart Nesterov method O’donoghue (2015) to accelerate general fixed-point iteration problems.
#'
#' @param par Vector for initial parameters
#' @param fixptfn Fixed point updating function
#' @param objfn Objective function
#' @param control A list containing parameters controlling the algorithm
#' @param ... Other arguments required by \code{fixptfn} and \code{objfn}
#'
#' @details The task it to \strong{minimize} \code{objfn}. Default values of \code{control} are: \code{projection=function(x) x, tol=1e-7, maxiter=2000, convtype="parameter", par.track=FALSE, conv.spec=NULL}.
#' \describe{
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
#' @references O’donoghue B, Candes E (2015). Adaptive restart for accelerated gradient schemes. Foundations of Computational Mathematics, 15(3): 715–732.
#'
#' @examples
#' \dontrun{
#' set.seed(54321)
#' prob = lasso_task(lam=1)
#' reNesterov(prob$initfn(), prob$fixptfn, prob$objfn, X=prob$X, y=prob$y)
#' }
#'
#' @export reNesterov
reNesterov <- function (par, fixptfn_name, objfn, ..., control=list()){
  control.default <- list(tol=1e-7, maxiter=2000, par.track=FALSE,
                          convtype="parameter", projection = function(x) x,
                          conv.spec=NULL)
  namc <- names(control)
  ctrl <- modifyList(control.default, control[namc %in% names(control.default)])
  
  tol = ctrl$tol
  maxiter = ctrl$maxiter
  trace = ctrl$par.track
  convtype = ctrl$convtype
  proj = ctrl$projection
  
  convf <- control$conv.spec
  if(!is.null(convf)) convtype="user"
  
  par.track = c(par)
  
  if(fixptfn_name == "MMdelta"){
    delta1 <- Inf
    delta2 <- Inf
    MM.track = numeric(ctrl$maxiter*3)
  } else{
    MM.track = NULL
  }
  
  pp <- length(par)
  conv <- TRUE
  beta.vecy <- beta.old <- par
  objfn.val <- rep(NA, maxiter + 1)
  alpha <- 1
  objfn.val[1] <- objfn(beta.old, ...)
  obj.old = objfn.val[1]
  
  for (k in 1:maxiter) {
    # beta.new <- fixptfn(beta.vecy, ...)
    # oval <- objfn(beta.new, ...)
    
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        # fixptfn = "MMA"
        beta.new = MMA(D = beta.vecy, ...)
      } else{
        MM = "MMW"
        beta.new = MMW(D = beta.vecy, ...)
      }
      MM.track[k] = MM
    } else if(fixptfn_name == "MMAW"){
      if (k %% 2) {
        beta.new = MMA(D = beta.vecy, ...)
      } else {
        beta.new = MMW(D = beta.vecy, ...)
      }
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      beta.new = fixptfn(beta.vecy, ...)
    }
    oval <- objfn(beta.new, ...)
    obj.new = oval
    obj.old = objfn.val[k]
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        delta1 = obj.old - obj.new
      } else if(MM == "MMW"){
        delta2 = obj.old - obj.new
      }
    }
    if(convtype == "objfn" & abs(oval - objfn.val[k]) < tol){
      break
    }
    
    if (is.na(objfn.val[k]) | is.na(oval) | oval > objfn.val[k]) {
      alpha <- 1
      beta.vecy <- beta.new
    }
    else {
      alpha.new <- 1/2 + sqrt(1 + 4 * alpha * alpha)/2
      mom <- (alpha - 1)/alpha.new
      alpha <- alpha.new
      beta.vecy <- beta.new + mom * (beta.new - beta.old)
      ## projection
      beta.vecy <- proj(beta.vecy)
      ss.resids <- sqrt(crossprod(beta.new - beta.old))
      if (convtype == "parameter" & ss.resids < tol)
        break
    }
    
    if (convtype == "user")
      if (convf(beta.old, beta.new, objfn.val[k], oval, tol))
        break
    
    
    objfn.val[k + 1] <- oval
    beta.old <- beta.new
    
    if (trace) par.track = rbind(par.track, beta.new)
  }
  
  rownames(par.track) = NULL
  
  objfn.val <- objfn.val[!is.na(objfn.val)]
  if (k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  value.obj <- objfn.val[k]
  # return(list(par = beta.new, value.objfn = value.obj, objfn.track = objfn.val,
  #             fpevals = k, iter = k, convergence = conv, par.track=par.track))
  # 
  return(list(par = beta.new, value.objfn = value.obj, objfn.track = objfn.val,
              fpevals = k, iter = k, convergence = conv, par.track=par.track, MM.track = MM.track[1:k]))
}