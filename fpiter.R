### Adapted from fpiter function from AccelBenchmark. 


#' Fixed point iteration
#' 
#' Apply given fixed point function iteratively until convergence
#' 
#' @param par Vector for initial parameters
#' @param fixptfn Fixed point updating function
#' @param objfn Objective function
#' @param control A list containing parameters controlling the algorithm
#' @param ... Other arguments required by \code{fixptfn} and \code{objfn}
#' 
#' @details Default values of \code{control} are: \code{tol=1e-7, maxiter=2000, convtype="parameter", par.track=FALSE, conv.spec=NULL}
#' \describe{
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
#' @examples
#' \dontrun{
#' set.seed(54321)
#' prob = lasso_task(lam=1)
#' fpiter(prob$initfn(), prob$fixptfn, prob$objfn, X=prob$X, y=prob$y)
#' }
#' 
#' @export fpiter


fpiter <- function(par, fixptfn_name, objfn, control = list(), ...){
  
  if(fixptfn_name == "MMA"){
    fixptfn = MMA
  } else if(fixptfn_name == "MMW"){
    fixptfn = MMW
  } else if(fixptfn_name == "MMAW"){
    fixptfn = MMAW
  } else if(fixptfn_name == "MMdelta"){
    fixptfn = MMdelta
  } else if(fixptfn_name == "MMlambda"){
    fixptfn = MMlambda
  }
  
  control.default <- list(tol=1.e-07, maxiter=2000, convtype="parameter",
                          par.track=FALSE, conv.spec=NULL)
  namc <- names(control)
  ctrl <- modifyList(control.default, control[namc %in% names(control.default)])
  
  tol <- ctrl$tol
  maxiter <- ctrl$maxiter
  convtype <- ctrl$convtype
  missing_obj <- is.null(objfn)
  track <- ctrl$par.track
  if(missing(objfn)) {
    convtype <- "parameter"
  }
  convf <- control$conv.spec
  if(!is.null(convf)) convtype="user"
  
  iter <- 0 ### change to 0. 
  par.track = c(par)
  # obj <- obj.new <- NA
  obj = objfn(par, ...)
  leval <- 0
  objval.track <- c()
  
  if (!missing(objfn)) {
    obj = objfn(par, ...)
    leval = leval + 1
    objval.track = c(objval.track, obj)
  }
  conv <- FALSE
  
  if(fixptfn_name == "MMdelta"){
    delta1 = Inf
    delta2 = Inf
    MM.track = numeric(maxiter)
  } else{
    MM.track = NULL
  }
  
  
  while (iter < maxiter) {
    iter = iter + 1
    if(fixptfn_name == "MMdelta"){
      p.new_out <- fixptfn(par, ..., fObj = objfn, delta1 = delta1, delta2 = delta2, iter = iter)
      MM = p.new_out$MM
      MM.track[iter] = MM
      p.new <- p.new_out$Dnew
      if (!missing(objfn)) {
        obj.new <- objfn(p.new, ...)
        leval <- leval + 1
        res.obj <- abs(obj.new - obj)
      }
      if(MM == "MMA"){
        delta1 = obj - obj.new
      } else if(MM == "MMW"){
        delta2 = obj - obj.new
      }
    } else{
      p.new <- fixptfn(par, ..., fObj = objfn, iter = iter)
      if (!missing(objfn)) {
        obj.new <- objfn(p.new, ...)
        leval <- leval + 1
        res.obj <- abs(obj.new - obj)
      }
    }
    
    res <- c(sqrt(crossprod(p.new - par)))
    
    if (convtype == "parameter" & !is.na(res) & res < tol) {
      conv <- TRUE
      par = p.new
      break
    }
    if (convtype == "objfn" & !is.na(res.obj) & res.obj < tol) {
      conv <- TRUE
      par = p.new
      break
    }
    if (convtype == "user")
      if(convf(par, p.new, obj, obj.new, tol)) {
        conv <- TRUE
        par = p.new
        break
      }
    
    if (!missing(objfn)){
      obj = obj.new
      objval.track = c(objval.track, obj)
    }
    
    par <- p.new
    if (track) par.track = rbind(par.track, par)
  }
  
  rownames(par.track) = NULL
  
  
  return(list(par = par, value.objfn = obj, iter = (iter), fpevals = (iter), objfevals = leval,
              convergence = conv, objfn.track = objval.track, par.track = par.track, MM.track = MM.track[1:(iter)]))
}