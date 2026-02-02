runAllMethods = function(Wg = NULL, Ag = NULL, d = NULL, G = NULL, D0 = NULL, MMfns = c("MMA" ,"MMW", "MMAW", "MMdelta", "MMlambda"), fObj, algorithms, control, control.spec = list()){
  

  ### adapted from benchmark function from AccelBenchmark
  ### run all methods on simulated data
  
  alg_dict = list(
    base = function(...) fpiter(...),
    squarem = function(...) squarem(...),
    daarem = function(...) daarem(...),
    qn = function(...) quasi_newton(...),
    pem = function(...) parabolic_em(...),
    nes = function(...) reNesterov(...)
  )
    
    result_table = data.frame(
      MM = NA, acc = NA, objval = NA, niter = NA,
      fpevals = NA, convergence = NA, user_time = NA, elapsed_time = NA
    )
    
    all_results = list()
    for(m in 1:length(MMfns)){
      fixptfn_name = MMfns[m]
        for(alg in algorithms){
          control.alg = control
          if(!is.null(control.spec[[alg]])){
            control.alg = modifyList(control.alg, control.spec[[alg]]) 
          }
          param_alg = list(par = as.vector(D0), fixptfn_name = fixptfn_name, objfn = fObj, control = control.alg, 
                           Wg = Wg, Ag = Ag, d = d, G = G)
          res = suppressWarnings(try({
            run.time = system.time({run.alg = do.call(alg_dict[[alg]], param_alg)}, gcFirst = TRUE)[c(1, 3)]
            run.alg$elapsed_time = run.time[2]
            run.alg$user_time = run.time[1]
            if(is.null(run.alg$iter)){
              run.alg$iter = run.alg$fpevals
            }
            if(MMfns[m] == "MMlambda"){
              run.alg$fpevals = 2*run.alg$fpevals
            }
            result_table = rbind(result_table,
                                c(MMfns[m], alg, run.alg$value.objfn, run.alg$iter,
                                  run.alg$fpevals, run.alg$convergence, run.time))
              
          }, TRUE))
          if(inherits(res, "try-error")){
            result_table = rbind(result_table, c(MMfns[m], alg, NA, NA, NA, FALSE, NA, NA))
            run.alg = list()
          }
          
          all_results = c(all_results, run.alg)
        }
      }
    #}
    result_table = result_table[-1, ] ## remove NA row
    return(list(result_table = result_table, all_results = all_results))
}
