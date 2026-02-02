#' Parabolic-EM Method for Accelerating Slowly-Convergent Fixed-Point Iterations
#'
#' Using parabolic-EM method Berlinet (2009) to accelerate general fixed-point iteration problems.
#'
#' @param par Vector for initial parameters
#' @param fixptfn Fixed point updating function
#' @param objfn Objective function
#' @param control A list containing parameters controlling the algorithm
#' @param ... Other arguments required by \code{fixptfn} and \code{objfn}
#'
#' @details The task it to \strong{minimize} \code{objfn}. Default values of \code{control} are: \code{warmup=5, h=0.1, a=1.5, maxtry=Inf, version="geometric", objfn.inc=0, projection=function(x) x, tol=1e-7, maxiter=2000, convtype="parameter", par.track=FALSE, conv.spec=NULL}.
#' \describe{
#'  \item{warmup}{An integer variable indicating the number of \code{fixptfn} to be evaluated before starting this algorithm. Default is 5.}
#'  \item{a}{A positive real number in the line search of geometric method. Default is 1.5.}
#'  \item{h}{A positive real number indicating the step size in the line search step. Default is 0.1.}
#'  \item{maxtry}{An integer variable indicating maximum number of try when searching for optimal step length in ever iteration. Default is Inf.}
#'  \item{version}{A string indicating the method used in searching the step length \eqn{t}. \eqn{t = 1 + h \times a^i} for "geometric" and \eqn{t = 1 + i \times h} for "arithmetic". Default is "geometric".}
#'  \item{projection}{A function projecting the parameter after each iteration. Default is identity function \eqn{f(x) = x}.}
#'  \item{objfn.inc}{A non-negative scalar that dictates the degree of non-montonicity. Default is 0. Set objfn.inc = 0 to obtain monotone convergence. Setting objfn.inc = Inf gives a non-monotone scheme. In-between values result in partially-monotone convergence.}
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
#' @references Berlinet A, Roland C (2009). Parabolic acceleration of the EM algorithm. Statistics and Computing, 19(1): 35–47.
#'
#' @examples
#' \dontrun{
#' set.seed(54321)
#' prob = lasso_task(lam=1)
#' parabolic_em(prob$initfn(), prob$fixptfn, prob$objfn, X=prob$X, y=prob$y)
#' }
#'
#' @export parabolic_em
parabolic_em = function(par, fixptfn_name, objfn, ..., control=list()){
  # fixptfn should return feasible par and no need for projection
  # minimizing task
  control.default <- list(
    convtype="parameter", tol=1.0e-07,
    maxiter=2000, par.track=FALSE,
    projection=function(x) x,
    warmup=5, h=0.1, a=1.5, maxtry=Inf,
    version="geometric", objfn.inc=0,
    conv.spec=NULL
  )
  control.sub = control[names(control) %in% names(control.default)]
  ctrl = modifyList(control.default, control.sub)
  
  convergence = TRUE
  fpevals = 0
  objfevals = 0
  objfn.track = c()
  
  par.track = c(par)
  
  convtype=ctrl$convtype
  convf <- control$conv.spec
  if(!is.null(convf)) convtype="user"
  obj.old = objfn(par, ...)
  objfevals = objfevals + 1
  
  if(fixptfn_name == "MMdelta"){
    delta1 <- Inf
    delta2 <- Inf
    MM.track = numeric(ctrl$maxiter*3)
  } else{
    MM.track = NULL
  }
  
  # warm up
  for(i in 1:ctrl$warmup){
    # par = fixptfn(par, ...)
    # fpevals = fpevals + 1
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        par = MMA(D = par, ...)
      } else{
        MM = "MMW"
        par = MMW(D = par, ...)
      }
      fpevals = fpevals + 1
      MM.track[fpevals] = MM
    } else if(fixptfn_name == "MMAW"){
      if (i %% 2) {
        par = MMA(D = par, ...)
      } else {
        par = MMW(D = par, ...)
      }
      fpevals = fpevals + 1
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      par = fixptfn(par, ...)
      fpevals = fpevals + 1
    }
    obj.new = objfn(par, ...)
    objfevals = objfevals + 1
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        delta1 = obj.old - obj.new
      } else if(MM == "MMW"){
        delta2 = obj.old - obj.new
      }
    }
    obj.old = obj.new
  }
  
  if (ctrl$par.track) par.track = rbind(par.track, par)
  iter = i + 1
  # set up
  p0 = par
  # p1 = fixptfn(p0, ...)
  if(fixptfn_name == "MMdelta"){
    if(delta1 > delta2){
      MM = "MMA"
      p1 = MMA(D = p0, ...)
    } else{
      MM = "MMW"
      p1 = MMW(D = p0, ...)
    }
    fpevals = fpevals + 1
    MM.track[fpevals] = MM
  } else if(fixptfn_name == "MMAW"){
    if (iter %% 2) {
      p1 = MMA(D = p0, ...)
    } else {
      p1 = MMW(D = p0, ...)
    }
    fpevals = fpevals + 1
  } else{
    if(fixptfn_name == "MMA"){
      fixptfn = MMA
    } else if(fixptfn_name == "MMW"){
      fixptfn = MMW
    }
    p1 = fixptfn(p0, ...)
    fpevals = fpevals + 1
  }
  
  # p2 = fixptfn(p1, ...)
  if(fixptfn_name == "MMdelta"){
    if(delta1 > delta2){
      MM = "MMA"
      p2 = MMA(D = p1, ...)
    } else{
      MM = "MMW"
      p2 = MMW(D = p1, ...)
    }
    fpevals = fpevals + 1
    MM.track[fpevals] = MM
  } else if(fixptfn_name == "MMAW"){
    if (iter %% 2) {
      p2 = MMA(D = p1, ...)
    } else {
      p2 = MMW(D = p1, ...)
    }
    fpevals = fpevals + 1
  } else{
    if(fixptfn_name == "MMA"){
      fixptfn = MMA
    } else if(fixptfn_name == "MMW"){
      fixptfn = MMW
    }
    p2 = fixptfn(p1, ...)
    fpevals = fpevals + 1
  }
  # fpevals = fpevals + 2
  obj.new = objfn(p2, ...)
  L2 = -obj.new
  objfevals = objfevals + 1
  if(fixptfn_name == "MMdelta"){
    if(MM == "MMA"){
      delta1 = obj.old - obj.new
    } else if(MM == "MMW"){
      delta2 = obj.old - obj.new
    }
  }
  obj.old = obj.new
  
  objfn.track = c(objfn.track, -L2)
  
  i = 1; it = iter + 1
  
  while(fpevals <= ctrl$maxiter){
    p2old = p2; l2old = L2
    if(ctrl$version=="geometric") {
      t <- 1 + ctrl$a^(i-1) * ctrl$h
    } else if(ctrl$version=="arithmetic") {
      t <- 1 + i * ctrl$h
    } else {
      stop("for pem, version must be one of c('geometric','arithmetic')")
    }
    
    p_pem = ctrl$projection( (1-t)^2*p0 + 2*t*(1-t)*p1 + t^2*p2 )
    Lpem = try(-objfn(p_pem, ...), silent=TRUE)
    objfevals = objfevals + 1
    
    if(inherits(Lpem, "try-error") | is.nan(Lpem) | Lpem <= L2) {
      p0 <- p2
      if(fixptfn_name == "MMdelta"){
        if(delta1 > delta2){
          MM == "MMA"
          p1 = MMA(D = p0, ...)
        } else{
          MM = "MMW"
          p1 = MMW(D = p0, ...)
        }
        fpevals = fpevals + 1
        MM.track[fpevals] = MM
      } else if(fixptfn_name == "MMAW"){
        if (iter %% 2) {
          p1 = MMA(D = p0, ...)
        } else {
          p1 = MMW(D = p0, ...)
        }
        fpevals = fpevals + 1
      } else{
        if(fixptfn_name == "MMA"){
          fixptfn = MMA
        } else if(fixptfn_name == "MMW"){
          fixptfn = MMW
        }
        p1 = fixptfn(p0, ...)
        fpevals = fpevals + 1
      }
      
      if(fixptfn_name == "MMdelta"){
        if(delta1 > delta2){
          MM = "MMA"
          p2 = MMA(D = p1, ...)
        } else{
          MM = "MMW"
          p2 = MMW(D = p1, ...)
        }
        fpevals = fpevals + 1
        MM.track[fpevals] = MM
      } else if(fixptfn_name == "MMAW"){
        if (iter %% 2) {
          p2 = MMA(D = p1, ...)
        } else {
          p2 = MMW(D = p1, ...)
        }
        fpevals = fpevals + 1
      } else{
        if(fixptfn_name == "MMA"){
          fixptfn = MMA
        } else if(fixptfn_name == "MMW"){
          fixptfn = MMW
        }
        p2 = fixptfn(p1, ...)
        fpevals = fpevals + 1
      }
      # fpevals <- fpevals + 2
    } else {
      times = 0
      while(Lpem + ctrl$objfn.inc > L2 & times < ctrl$maxtry){
        pold = p_pem
        L2 = Lpem
        i = i + 1
        if(ctrl$version=="geometric") {
          t <- 1 + ctrl$a^(i - 1) * ctrl$h
        } else{
          t <- 1 + i * ctrl$h
        }
        p_pem = ctrl$projection( (1 - t)^2*p0 + 2*t*(1-t)*p1 + t^2*p2 )
        Lpem = try(-objfn(p_pem, ...), silent=TRUE)
        objfevals = objfevals + 1
        times = times + 1
      }
      p0 = p1
      p1 = p2
      
      # p2 = fixptfn(fixptfn(pold, ...), ...)
      if(fixptfn_name == "MMdelta"){
        if(delta1 > delta2){
          MM = "MMA"
          p2 = MMA(MMA(D = pold, ...), ...)
        } else{
          MM = "MMW"
          p2 = MMW(MMW(D = pold, ...), ...)
        }
        fpevals = fpevals + 2
        MM.track[fpevals] = MM
      } else if(fixptfn_name == "MMAW"){
        if (iter %% 2) {
          p2 = MMA(MMA(D = pold, ...), ...)
        } else {
          p2 = MMW(MMW(D = pold, ...), ...)
        }
        fpevals = fpevals + 2
      } else{
        if(fixptfn_name == "MMA"){
          fixptfn = MMA
        } else if(fixptfn_name == "MMW"){
          fixptfn = MMW
        }
        p2 = fixptfn(fixptfn(pold, ...), ...)
        fpevals = fpevals + 2
      }
      # fpevals <- fpevals + 2
    }
    L2 = try(-objfn(p2, ...), silent=TRUE)
    objfevals = objfevals + 1
    objfn.track = c(objfn.track, -L2)
    obj.new = -L2
    
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        delta1 = obj.old - obj.new
      } else if(MM == "MMW"){
        delta2 = obj.old - obj.new
      }
    }
    obj.old = obj.new
    
    if(convtype == "objfn" & abs(L2 - l2old) < ctrl$tol)
      break
    if(convtype == "parameter" & sqrt(crossprod(p2 - p2old)) < ctrl$tol)
      break
    if(convtype == "user")
      if (convf(p2old, p2, -l2old, -L2, ctrl$tol))
        break
    
    if (ctrl$par.track) par.track = rbind(par.track, p2)
    
    it = it + 1
    iter = iter + 1
  }
  
  rownames(par.track) = NULL
  
  if(fpevals > ctrl$maxiter){
    convergence = FALSE
    # warning("Algorithm did not converge")
  }
  
  # list(par = c(p2),
  #      value.objfn = -L2,
  #      iter = it,
  #      fpevals = fpevals,
  #      objfevals = objfevals,
  #      convergence = convergence,
  #      objfn.track = objfn.track,
  #      par.track = par.track)
  
  list(par = c(p2),
       value.objfn = -L2,
       iter = iter,
       fpevals = fpevals,
       objfevals = objfevals,
       convergence = convergence,
       objfn.track = objfn.track,
       par.track = par.track,
       MM.track = MM.track[1:(fpevals)])
}
