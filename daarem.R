### Adapted from daarem function from AccelBenchmark


#' Damped Anderson Acceleration with Restarts and Epsilon-Montonicity for Accelerating Slowly-Convergent, Monotone Fixed-Point Iterations
#'
#' An ‘off-the-shelf’ acceleration scheme for accelerating the convergence of any smooth, monotone, slowly-converging fixed-point iteration. It can be used to accelerate the convergence of a wide variety of montone iterations including, for example, expectation-maximization (EM) algorithms and majorization-minimization (MM) algorithms.
#' This is an modified version of the original \code{daarem} in \code{daarem} package, including projection and user defined convergence function.
#'
#' @param par Vector for initial parameters
#' @param fixptfn Fixed point updating function
#' @param objfn Objective function
#' @param control A list containing parameters controlling the algorithm
#' @param ... Other arguments required by \code{fixptfn} and \code{objfn}
#'
#' @details The task it to \strong{minimize} \code{objfn}. Default values of \code{control} are: \code{order=10, mon.tol=0.01, cycl.mon.tol=0.0, alpha=1.2, kappa=25, projection=function(x) x, tol=1e-7, maxiter=2000, convtype="parameter", par.track=FALSE, conv.spec=NULL}.
#' \describe{
#'  \item{order}{An integer >= 1 denoting the order of the DAAREM acceleration scheme. Default is 10.}
#'  \item{mon.tol}{A nonnegative scalar that determines whether the montonicity condition is violated. The monotonicity condition is violated whenver \eqn{L(x[k+1]) > L(x[k]) + mon.tol}. Such violations determine how much damping is to be applied on subsequent steps of the algorithm. Default value of mon.tol is 1.e-02.}
#'  \item{cycl.mon.tol}{A nonegative scalar that determines whether a montonicity condition is violated after the end of the cycle. This cycle-level monotonicity condition is violated whenver \eqn{L(x[end cycle]) > L(x[start cycle]) + cycl.mon.tol}. Here, x[start cycle] refers to the value of x at the beginning of the current cycle while x[end cycle] refers to the value of x at the end of the current cycle. Such violations also determine how much damping is to be applied on subsequent steps of the algorithm.}
#'  \item{kappa}{A nonnegative parameter which determines the “half-life” of relative damping and how quickly relative damping tends to one. In the absence of monotonicity violations, the relative damping factor is <= 1/2 for the first kappa iterations, and it is then greater than 1/2 for all subsequent iterations. The relative damping factor is the ratio between the norm of the unconstrained coefficients in Anderson acceleration and the norm of the damped coefficients. In the absence of any monotonicity violations, the relative damping factor in iteration k is \eqn{1/(1 + a^{(kappa - k)})}.}
#'  \item{alpha}{A parameter > 1 that determines the initial relative damping factor and how quickly the relative damping factor tends to one. The initial relative damping factor is \eqn{1/(1 + a^{kappa})}. In the absence of any monotonicity violations, the relative damping factor in iteration k is \eqn{1/(1 + a^{(kappa - k)})}.}
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
#' @references Henderson, N.C. and Varadhan, R. (2019) Damped Anderson acceleration with restarts and monotonicity control for accelerating EM and EM-like algorithms, Journal of Computational and Graphical Statistics, Vol. 28(4), 834-846.
#'
#' @examples
#' \dontrun{
#' set.seed(54321)
#' prob = lasso_task(lam=1)
#' daarem(prob$initfn(), prob$fixptfn, prob$objfn, X=prob$X, y=prob$y)
#' }
#'
#' @importFrom stats dnorm pnorm runif uniroot
#'
#' @export daarem
daarem <- function(par, fixptfn_name, objfn, ..., control=list()) {
  
  control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0,
                          kappa=25, alpha=1.2, resid.tol=0.95, convtype="parameter",
                          par.track=FALSE, projection=function(x) x, conv.spec=NULL)
  namc <- names(control)
  control <- modifyList(control.default, control[namc %in% names(control.default)])
  
  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa
  resid.tol <- control$resid.tol
  proj <- control$projection
  trace <- control$par.track
  convf <- control$conv.spec
  if(control$convtype=="parameter") {
    check.par.resid <- TRUE
  } else if(control$convtype=="objfn") {
    check.par.resid <- FALSE
  }
  
  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))
  
  if(!missing(objfn)) {
    ans <- daarem_base_objfn(
      par=par, fixptfn_name=fixptfn_name, objfn_neg=objfn, maxiter=maxiter,
      tol=tol, mon.tol=mon.tol, cycl.mon.tol=cycl.mon.tol,
      a1=a1, kappa=kappa, num.params=num.params, nlag=nlag,
      check.par.resid=check.par.resid, proj=proj, trace=trace,
      convf=convf, ...
    )
  } else {
    ans <- daarem_base_noobjfn(par, fixptfn_name, maxiter, tol, resid.tol,
                               a1, kappa, num.params, nlag, proj, trace, ...)
  }
  if(!ans$convergence) {
    warning("Algorithm did not converge")
  }
  return(ans)
}
