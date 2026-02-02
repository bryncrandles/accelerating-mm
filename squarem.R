### Adapted from squarem function from AccelBenchmark. 


#' Squared Extrapolation Methods for Accelerating Slowly-Convergent Fixed-Point Iterations
#'
#' Globally-convergent, partially monotone, acceleration schemes for accelerating the convergence of *any* smooth, monotone, slowly-converging contraction mapping. It can be used to accelerate the convergence of a wide variety of iterations including the expectation-maximization (EM) algorithms and its variants, majorization-minimization (MM) algorithm, power method for dominant eigenvalue-eigenvector, Google's page-rank algorithm, and multi-dimensional scaling.
#' This is a modification of the original squarem in \code{SQUAREM} package, including choice of projection and user defined convergence function.
#'
#' @param par Vector for initial parameters
#' @param fixptfn Fixed point updating function
#' @param objfn Objective function
#' @param control A list containing parameters controlling the algorithm
#' @param ... Other arguments required by \code{fixptfn} and \code{objfn}
#'
#' @details The task it to \strong{minimize} \code{objfn}. Default values of \code{control} are: \code{method=3, step.min0=1, step.max0=1, mstep=4, objfn.inc=1, projection=function(x) x, tol=1e-7, maxiter=2000, convtype="parameter", par.track=FALSE, conv.spec=NULL}.
#' \describe{
#'  \item{method}{An integer variable that denotes the particular SQUAREM scheme to be used. Value should be either 1, 2, or 3. These correspond to the 3 schemes discussed in Varadhan and Roland (2008). Default is 3.}
#'  \item{step.min0}{A scalar denoting the minimum steplength taken by a SQUAREM algorithm. Default is 1. For contractive fixed-point iterations (e.g. EM and MM), this defualt works well. In problems where an eigenvalue of the Jacobian of $F$ is outside of the interval (0,1), step.min0 should be less than 1 or even negative in some cases.}
#'  \item{step.max0}{A positive-valued scalar denoting the initial value of the maximum steplength taken by a SQUAREM algorithm. Default is 1. When the steplength computed by SQUAREM exceeds step.max0, the steplength is set equal to step.max0, but then step.max0 is increased by a factor of mstep.}
#'  \item{mstep}{A scalar greater than 1. When the steplength computed by SQUAREM exceeds step.max0, the steplength is set equal to step.max0, but step.max0 is increased by a factor of mstep. Default is 4.}
#'  \item{objfn.inc}{A non-negative scalar that dictates the degree of non-montonicity. Default is 1. Set objfn.inc = 0 to obtain monotone convergence. Setting objfn.inc = Inf gives a non-monotone scheme. In-between values result in partially-monotone convergence.}
#'  \item{projection}{A function projecting the parameter after each iteration. Default is identity function \eqn{f(x) = x}}
#'  \item{tol}{A small, positive scalar that determines when iterations should be terminated, see \code{convtype} for details. Default is \code{1e-7}}
#'  \item{maxiter}{An integer denoting the maximum limit on the number of evaluations of \code{fixptfn}. Default is 2000.}
#'  \item{convtype}{A string indicating the convergence criteria.
#'                 If it is "parameter", the algorithm will terminate when L2 norm of parameters difference \eqn{x_{new} - x_{old} < tol}.
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
#' @references R Varadhan and C Roland (2008), Simple and globally convergent numerical schemes for accelerating the convergence of any EM algorithm, Scandinavian Journal of Statistics, 35:335-353.
#' @references C Roland, R Varadhan, and CE Frangakis (2007), Squared polynomial extrapolation methods with cycling: an application to the positron emission tomography problem, Numerical Algorithms, 44:159-172.
#' @references Y Du and R Varadhan (2020), SQUAREM: An R package for off-the-shelf acceleration of EM, MM, and other EM-like monotone algorithms, Journal of Statistical Software, 92(7): 1-41. <doi:10.18637/jss.v092.i07>
#'
#' @examples
#' \dontrun{
#' set.seed(54321)
#' prob = lasso_task(lam=1)
#' squarem(prob$initfn(), prob$fixptfn, prob$objfn, X=prob$X, y=prob$y)
#' }
#'
#' @importFrom utils modifyList
#'
#' @export squarem
squarem <- function(par, fixptfn_name, objfn, ... , control=list()) {
  control.default <- list(
    method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, moment=1,
    objfn.inc=1, tol=1.e-07, maxiter=2000, convtype="parameter",
    projection = function(x) x, par.track=FALSE, conv.spec=NULL
  )
  namc <- names(control)
  ctrl <- modifyList(control.default, control[namc %in% names(control.default)])
  
  if (!(ctrl$method %in% c(1,2,3,4))) ctrl$method <- 3
  
  if (missing(objfn)) {
    ctrl$missing_obj = TRUE
    ctrl$convtype = "parameter"
  } else
    ctrl$missing_obj = FALSE
  
  sqobj <- squarem_body(par, fixptfn_name, objfn, ... , control=ctrl)
  
  return(sqobj)
}

######################################################################
## main function for SQUAREM, requires objfn

squarem_body <- function(par, fixptfn_name, objfn, ... , control=list()) {
  # par = starting value of parameter vector
  ### CHANGED fixptfn_name = string of name of function of fixed-point iteration F(x)
  # for which the solution: F(x*) = x* is sought
  # objfn = underlying objective function which is minimized at x*
  
  
  #####
  # method = 1, 2, or 3, indicating the type of steplength to be used
  # A key parameter is `step.min0'. This can be either positive or negative if the eigenvalues of the Jacobian of fixed-point mapping are all positive;
  # We set it to +1 for contraction mappings such as EM and MM algorithms
  # Must be negative, e.g. -1, if the fixed-point mapping can have negative eigenvalues
  #####
  # parameter "objfn.inc" dictates the amount of non-monotonicity in the objective function
  # setting objfn.inc=0 would enforce monotonicity, whereas objfn.inc=Inf would be a non-monotonic scheme
  # The defalut objfn.inc=1 would enforce monotonicity far from solution, but allows for non-monotonicity closer to solution
  
  
  method <- control$method
  maxiter <- control$maxiter
  tol <- control$tol
  step.min <- control$step.min0
  step.max <- control$step.max0
  step.max0 <- control$step.max0
  mstep <- control$mstep
  objfn.inc <- control$objfn.inc
  moment <- control$moment
  convtype <- control$convtype
  proj <- control$projection
  missing_obj <- control$missing_obj
  track <- control$par.track
  convf <- control$conv.spec
  if(!is.null(convf)) convtype="user"
  
  #fixptfn_mm = function(par, ..., iter = iter) (1-moment) * par + moment * fixptfn(par, ..., iter = iter)
  #fixptfn_mm = fixptfn
  
  iter <- 1
  p <- par
  
  leval <- 0
  objval.track <- c()
  lold <- lnew <- NA
  
  if (!missing_obj) {
    lold = objfn(p, ...)
    leval = leval + 1
    objval.track = c(objval.track, lold)
  }
  
  feval <- 0
  par.track = c(par)
  conv <- TRUE
  
  alpha.track = c()
  
  if(fixptfn_name == "MMdelta"){
    delta1 <- Inf
    delta2 <- Inf
    MM.track = numeric(maxiter*3)
  } else{
    MM.track = NULL
  }
  
  while (feval < maxiter) {
    
    extrap <- TRUE
    
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        # fixptfn = "MMA"
        p1 = MMA(D = p, ...)
      } else{
        MM = "MMW"
        p1 = MMW(D = p, ...)
      }
      feval = feval + 1
      MM.track[feval] = MM
    } else if(fixptfn_name == "MMAW"){
      if (iter %% 2) {
        p1 = MMA(D = p, ...)
      } else {
        p1 = MMW(D = p, ...)
      }
      feval = feval + 1
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      p1 = fixptfn(p, ...)
      feval = feval + 1
    }

    q1 <- p1 - p
    sr2 <- c( crossprod(q1) )
    if (convtype=="parameter" & !is.nan(sr2) & sqrt(sr2) < tol) break
    
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        p2 = MMA(D = p1, ...)
      } else{
        p2 = MMW(D = p1, ...)
      }
      feval = feval + 1
      MM.track[feval] = MM
    } else if(fixptfn_name == "MMAW"){
      if (iter %% 2) {
        p2 = MMA(D = p1, ...)
      } else {
        p2 = MMW(D = p1, ...)
      }
      feval = feval + 1
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      p2 = fixptfn(p1, ...)
      feval = feval + 1
    }

    q2 <- p2 - p1
    sq2 <- c( sqrt(crossprod(q2)) )
    if (convtype=="parameter" & !is.na(sq2) & sq2 < tol) break
    
    sv2 <- c( crossprod(q2 - q1) )
    srv <- c( crossprod(q1, q2 - q1) )
    
    if(method == 4) {
      alphas = c(-srv/sv2, -sr2/srv, sqrt(sr2/sv2))
      steplen = c(
        norm(2*alphas[1]*q1 + alphas[1]^2*(q2-q1), "2"),
        norm(2*alphas[2]*q1 + alphas[2]^2*(q2-q1), "2"),
        norm(2*alphas[3]*q1 + alphas[3]^2*(q2-q1), "2")
      )
      alpha = alphas[which.min(steplen)]
    } else {
      alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
    }
    
    alpha <- max(step.min, min(step.max, alpha))
    alpha.track  = c(alpha.track, alpha)
    
    
    p.new <- p + 2*alpha*q1 + alpha^2*(q2-q1)
    
    ## projection
    p.new <- proj(p.new)
    
    ### compute obj function at squarem update
    # lnew = objfn(p.new, ...)
    # leval = leval + 1
    
    
    if (!is.finite(alpha) || abs(alpha - 1) > 0.01) {
      if(fixptfn_name == "MMdelta"){
        if(MM == "MMA"){
          p.new = MMA(D = p.new, ...)
        } else{
          p.new = MMW(D = p.new, ...)
        }
        feval = feval + 1
        MM.track[feval] = MM
      } else if(fixptfn_name == "MMAW"){
        if (iter %% 2) {
          p.new = MMA(D = p.new, ...)
        } else {
          p.new = MMW(D = p.new, ...)
        }
        feval = feval + 1
      } else{
        if(fixptfn_name == "MMA"){
          fixptfn = MMA
        } else if(fixptfn_name == "MMW"){
          fixptfn = MMW
        }
        p.new = fixptfn(p.new, ...)
        feval = feval + 1
      }
    }
    
    if (inherits(p.new, "try-error") | any(is.nan(p.new))) {
      p.new <- p2
      if (!missing_obj) {
        lnew <- try(objfn(p2, ...), silent=TRUE)
        leval <- leval + 1
      }
      if (is.finite(alpha) & alpha == step.max) step.max <- max(step.max0, step.max/mstep)
      alpha <- 1
      extrap <- FALSE
    } else if(!missing_obj) {
      if (is.finite(objfn.inc)) {
        lnew <- try(objfn(p.new, ...), silent=TRUE)
        leval <- leval + 1
      } else lnew <- lold
      # if (class(lnew) == "try-error" | is.nan(lnew) |
      if (inherits(lnew, "try-error") | is.nan(lnew) |
          (lnew > lold + objfn.inc)) {
        p.new <- p2
        lnew <- try(objfn(p2, ...), silent=TRUE)
        leval <- leval + 1
        if (is.finite(alpha) & alpha==step.max) step.max <- max(step.max0, step.max/mstep)
        alpha <- 1
        extrap <- FALSE
      }
    }
    
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        delta1 = lold - lnew
      } else if(MM == "MMW"){
        delta2 = lold - lnew
      }
    }
    
    if (is.finite(alpha) & alpha == step.max) step.max <- mstep*step.max
    if (step.min < 0 & is.finite(alpha) & alpha == step.min) step.min <- mstep*step.min
    
    if (convtype == "user")
      if (convf(p, p.new, lold, lnew, tol)) {
        lold = lnew
        p = p.new
        break
      }
    
    p <- p.new
    if(convtype == "objfn" & abs(lold - lnew) < tol) break
    if (!missing_obj) {
      lold <- lnew
      objval.track = c(objval.track, lold)
    }
    if (track) par.track = rbind(par.track, p)
    
    
    iter <- iter + 1
  }
  
  if (feval >= maxiter) conv <- FALSE
  if (!missing_obj & is.infinite(objfn.inc)) {
    lold <- objfn(p, ...)
    leval <- leval + 1
  }
  
  rownames(par.track) = NULL
  
  return(list(par = p, value.objfn = lold, iter = iter, fpevals = feval,
              objfevals = leval, convergence = conv, objfn.track = objval.track,
              par.track = par.track, alpha.track = alpha.track, MM.track = MM.track[1:(feval)]))
  
  # return(list(par = par, value.objfn = obj, iter = (iter), fpevals = (iter), objfevals = leval,
  #             convergence = conv, objfn.track = objval.track, par.track = par.track, MM.track = MM.track[1:(iter)]))
}