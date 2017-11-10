
GWASinlps = function( x, y, prior, tau, priorDelta = modelbbprior(1,1), k0, m, rxx, niter = 2000, nskip = 3, verbose = F )
{
    if(!exists("time", mode="integer")) time = Sys.time()

    x = as.matrix(x)
    if( is.null( colnames(x) ) ) colnames(x) = paste0( "X", 1:ncol(x) )  # if x has no column names, give column names

    y1 = y
    vars = colnames(x)
    vars.remove = c()
    vars.final = c() 
    
    max.nocollect = 0   # number of times no variables show up in the hppm 
    iter = 0
    while( length(vars.final) < m )
    {
      while( length( setdiff( vars, c( vars.final, vars.remove ) ) ) > 0  &&  max.nocollect < nskip )
      {  
        iter = iter + 1
        if( verbose == T) cat( "-------------", "\n", "Iteration ", iter, "\n", "-------------", "\n", sep = "") 
        run = nlpsLM( x = x[ , colnames(x) %in% setdiff( vars, c( vars.final, vars.remove ) ) ], y = y1, prior = prior, tau = tau, k0 = k0, rxx = rxx, verbose = verbose )

        collect =  run $ hppm 

        if( length(collect) > 0) break else 
        { 
          vars.remove = c( vars.remove, run $ names.not.selected ) 
          #here you can choose to regress out not selected vars.
          #k0 = min( k0 + k0_inc, length( setdiff( vars, c( vars.final, vars.remove ) ) ) ) 
          max.nocollect = max.nocollect + 1
          if( verbose == T) cat("***","nskip=",max.nocollect,"***","\n")
        }
      }

      if( length(collect) == 0) break else
      {
        vars.final = c( vars.final, collect )  # include HPPM vars in selected set

        #if( verbose == T) cat( "# selected variables up to now:", length( vars.final ), "\n" )  # print number of selected vars

        y1 = lm( y1 ~ as.matrix( x[ , colnames(x) %in% collect]) ) $ residuals  # regress out HPPM vars
      }

    }

    selected = vars.final[ 1 : min(length(vars.final),m) ]
    if(verbose == T) cat( "=================================", "\n","Number of selected variables: ", length(selected), "\n", "Time taken: ", round(difftime(Sys.time(), time, units = "mins"),2), " min", "\n", "=================================", "\n", sep = "")  
    return(selected)

}

