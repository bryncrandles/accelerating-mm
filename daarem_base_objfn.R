daarem_base_objfn <- function(par, fixptfn_name, objfn_neg, maxiter, tol, mon.tol,
                              cycl.mon.tol, a1, kappa, num.params, nlag,
                              check.par.resid, proj, trace, convf, ...) {
  
  ## objfn for minimizing
  objfn = function(par, ...) -objfn_neg(par, ...)
  
  par.track = c(par)
  
  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)
  
  if(fixptfn_name == "MMdelta"){
    delta1 <- Inf
    delta2 <- Inf
    MM.track = numeric(maxiter*3)
  } else{
    MM.track = NULL
  }
  
  xold <- par
  #xnew <- fixptfn(xold, ...)

  
  k = iter = 1
  
  if(fixptfn_name == "MMdelta"){
    if(delta1 > delta2){
      MM = "MMA"
      # fixptfn = "MMA"
      xnew = MMA(D = xold, ...)
    } else{
      MM = "MMW"
      xnew = MMW(D = xold, ...)
    }
    MM.track[iter] = MM
  } else if(fixptfn_name == "MMAW"){
    if (iter %% 2) {
      xnew = MMA(D = xold, ...)
    } else {
      xnew = MMW(D = xold, ...)
    }
  } else{
    if(fixptfn_name == "MMA"){
      fixptfn = MMA
    } else if(fixptfn_name == "MMW"){
      fixptfn = MMW
    }
    xnew = fixptfn(xold, ...)
  }
  
  lold = objfn(xold, ...)
  lnew = objfn(xnew, ...)
  
  # if(fixptfn_name == "MMdelta"){
  #   if(MM == "MMA"){
  #     delta1 = lold - lnew
  #   } else if(MM == "MMW"){
  #     delta2 = lold - lnew
  #   }
  # }
  
  obj_funvals[1] <- lold
  obj_funvals[2] <- lnew
  likchg <- obj_funvals[2] - obj_funvals[1]
  obj.evals <- 2
  if (trace) par.track = rbind(par.track, c(xnew))
  
  fold <- xnew - xold
  iter = k 
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- TRUE
  ell.star <- obj_funvals[2]
  
  lold = lnew
  
  while(k < maxiter) {
    count <- count + 1
    
    if(fixptfn_name == "MMdelta"){
      if(delta1 > delta2){
        MM = "MMA"
        # fixptfn = "MMA"
        xnew2 = MMA(D = xnew, ...)
      } else{
        MM = "MMW"
        xnew2 = MMW(D = xnew, ...)
      }
      MM.track[iter] = MM
    } else if(fixptfn_name == "MMAW"){
      if (iter %% 2) {
        xnew2 = MMA(D = xnew, ...)
      } else {
        xnew2 = MMW(D = xnew, ...)
      }
    } else{
      if(fixptfn_name == "MMA"){
        fixptfn = MMA
      } else if(fixptfn_name == "MMW"){
        fixptfn = MMW
      }
      xnew2 = fixptfn(xnew, ...)
    }
    fnew <- xnew2 - xnew
    ss.resids <- sqrt(crossprod(fnew))
    if(ss.resids < tol & check.par.resid & is.null(convf)) break
    #new.objective.val >= obj_funvals[k+1] - mon.tol
    
    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold
    
    np <- count
    if(np==1) {
      Ftmp <- matrix(Fdiff[,1], nrow=num.params, ncol=np)
      Xtmp <- matrix(Xdiff[,1], nrow=num.params, ncol=np)  ## is this matrix function needed?
      
    } else {
      Ftmp <- Fdiff[,1:np]
      Xtmp <- Xdiff[,1:np]
    }
    tmp <- La.svd(Ftmp)
    dvec <- tmp$d
    dvec.sq <- dvec*dvec
    uy <- crossprod(tmp$u, fnew)
    uy.sq <- uy*uy
    
    ### Still need to compute Ftf
    Ftf <- sqrt(sum(uy.sq*dvec.sq))
    tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr
    dd <- (dvec*uy)/(dvec.sq + lambda.ridge)
    gamma_vec <- crossprod(tmp$vt, dd)
    
    xbar <- xnew - drop(Xtmp%*%gamma_vec)
    fbar <- fnew - drop(Ftmp%*%gamma_vec)
    ## add projection
    x.propose <- proj( xbar + fbar )
    new.objective.val <- try(objfn(x.propose, ...), silent=TRUE)
    obj.evals <- obj.evals + 1
    
    
    
    if(class(new.objective.val)[1] != "try-error" & !is.na(obj_funvals[k+1]) &
       !is.nan(new.objective.val)) {
      if(new.objective.val >= obj_funvals[k+1] - mon.tol) {  ## just change this line in daarem_base_noobjfn
        ## Increase delta
        obj_funvals[k+2] <- new.objective.val
        fold <- fnew
        xold <- xnew
        
        xnew <- x.propose
        
        shrink.count <- shrink.count + 1
      } else {
        ## Keep delta the same
        fold <- fnew
        xold <- xnew
        
        xnew <- fold + xold
        obj_funvals[k+2] <- objfn(xnew, ...)
        obj.evals <- obj.evals + 1
      }
    } else {
      ## Keep delta the same
      fold <- fnew
      xold <- xnew
      
      xnew <- fold + xold
      obj_funvals[k+2] <- objfn(xnew, ...)
      obj.evals <- obj.evals + 1
      count <- 0
    }
    
    lnew = obj_funvals[k + 2]
    if(fixptfn_name == "MMdelta"){
      if(MM == "MMA"){
        delta1 = lnew - lold
      } else if(MM == "MMW"){
        delta2 = lnew - lold
      }
    }
    
    if(!is.null(convf))
      if (convf(xold, xnew, -obj_funvals[k+1], -obj_funvals[k+2], tol)) {
        break
      }
    if(count==nlag) {
      count <- 0
      ## restart count
      ## make comparison here l.star vs. obj_funvals[k+2]
      if(obj_funvals[k+2] < ell.star - cycl.mon.tol) {
        ## Decrease delta
        shrink.count <- max(shrink.count - nlag, -2*kappa)
      }
      ell.star <- obj_funvals[k+2]
    }
    if(abs(obj_funvals[k+2] - obj_funvals[k+1]) < tol & !check.par.resid & is.null(convf)) {
      break
    }
    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
    k <- k + 1
    iter = iter + 1
    if (trace) par.track = rbind(par.track, c(xnew))
    
  }
  rownames(par.track) = NULL
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]
  value.obj <- objfn(xnew, ...)
  if(k >= maxiter) {
    conv <- FALSE
  }
  
  return(list(par=c(xnew), value.objfn=-value.obj, iter = k, fpevals = k, objfevals=obj.evals,
              convergence=conv, objfn.track=-obj_funvals, par.track=par.track, MM.track = MM.track[1:k]))
  
#   return(list(par=c(xnew), fpevals = k, value.objfn=-value.obj, objfevals=obj.evals,
#               convergence=conv, objfn.track=-obj_funvals, par.track=par.track))
# 
}
