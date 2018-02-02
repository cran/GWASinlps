GWASinlps = function( y, x, prior, tau, priorDelta = modelbbprior(1,1), k0, m, rxx, nskip = 3, niter = 2000, verbose = F, 
  skip.return = F, seed = NULL, tau.hs.method = "halfCauchy", sigma.hs.method = "Jeffreys" )
{
    if(!exists("time", mode="integer")) time = Sys.time()
    if(!is.null(seed)) set.seed(seed)

    # ### check if inputs are sufficient and consistent

    # if( !(prior %in% c("mom","imom","zellner","horseshoe") ) 
    #   stop("currently implemented priors are mom, imom, zellner and horseshoe")


    x = as.matrix(x)
    if( is.null( colnames(x) ) ) colnames(x) = paste0( "X", 1:ncol(x) )  # if x has no column names, give column names

    y1 = y
    vars = colnames(x)
    vars.remove = c()
    vars.final = c() 
    selected = list()
    
    max.nocollect = 0   # number of times no variables show up in the hppm 
    iter = 0
    while( length(vars.final) < m )
    {
      while( length( setdiff( vars, c( vars.final, vars.remove ) ) ) > 0  &&  max.nocollect < nskip )
      {  
        iter = iter + 1
        if( verbose == T) cat( "-------------", "\n", "Iteration ", iter, "\n", "-------------", "\n", sep = "") 
        run = nlpsLM( y = y1, x = x[ , colnames(x) %in% setdiff( vars, c( vars.final, vars.remove ) ) ], prior = prior, tau = tau, k0 = k0, rxx = rxx, verbose = verbose, tau.hs.method = tau.hs.method, sigma.hs.method = sigma.hs.method )

        collect =  run $ hppm 

        if( length(collect) > 0) break else 
        { 
          vars.remove = c( vars.remove, run $names.not.selected ) 
          #here you can choose to regress out not selected vars.
          #k0 = min( k0 + k0_inc, length( setdiff( vars, c( vars.final, vars.remove ) ) ) ) 
          max.nocollect = max.nocollect + 1
          if( verbose == T) {cat("***","nskip=",max.nocollect,"***","\n")} 
          if( length(vars.final) == 0 ) { ## which means in the very first iteration no variable is selected
            selected = c( selected, '') 
          } else selected = c( selected, list(vars.final[ 1 : min(length(vars.final),m) ]) )
        }
      }

      if( length(collect) == 0) break else
      {
        vars.final = c( vars.final, collect )  # include HPPM vars in selected set
        vars.remove = c( vars.remove, run $names.not.selected )

        #if( verbose == T) cat( "# selected variables up to now:", length( vars.final ), "\n" )  # print number of selected vars

        y1 = lm( y1 ~ as.matrix( x[ , colnames(x) %in% collect]) ) $ residuals  # regress out HPPM vars
      }

    }

    if( length(selected) == 0 ) { ## which means no iterations are skipped yet, but number of selected variables already reached m
      selected = c( selected, list(vars.final[ 1 : min(length(vars.final),m) ]) )
    }
    if(verbose == T) cat( "=================================", "\n","Number of selected variables: ", length(selected[[length(selected)]]), "\n", "Time taken: ", round(difftime(Sys.time(), time, units = "mins"),2), " min", "\n", "=================================", "\n", sep = "")  
    if(skip.return) return(selected) else return(selected[[length(selected)]])

}

