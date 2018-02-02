nlpsLM = function( y, x, prior, tau, priorDelta = modelbbprior(1,1), k0, rxx, niter = 2000, verbose = F, tau.hs.method = "halfCauchy", sigma.hs.method = "Jeffreys" )
{ 
  corr.xy = sapply( as.data.frame(x), function(x) cor(x, y, use = "pairwise.complete.obs") )  # find corr of all x's with y

  names(corr.xy) = colnames(x)

  names.sorted.corr.xy = names( sort( abs(corr.xy), decreasing = T ) [1:k0] )  # find x's with top k0 corrs

  hppm = list()

  names.xx.input.set = c()

  for(i in 1:k0)  # take the i'th of top k0 x's
  {
    corr.xx = as.vector( abs( cor( x[ , names.sorted.corr.xy[i] ], x, use = "pairwise.complete.obs" ) ) )   # find corrs of all x variables with xi, and allocate names

    names(corr.xx) = colnames(x) 

    names.xx.thres = na.omit( names( corr.xx [corr.xx >= rxx] ) ) # take x's with corr > rxx with xi

    names.xx.input =  setdiff( names.xx.thres, unlist(hppm) )  # exclude x's previously chosen in hppm and set the rest as the input x's 

    names.xx.input.set = c( names.xx.input.set, names.xx.input )  # store the input x's

    if( verbose == T) cat( "j =",i, "\ninput :", names.xx.input, "\n" )

    if( length(names.xx.input) != 0  )  # if there is some input x
    {
      if( prior == "mom" ) 
      {
        bb = modelSelection( y, x = x[ , colnames(x) %in% names.xx.input, drop=F ], priorCoef = momprior(tau=tau), priorDelta = priorDelta, niter = niter, center = T, scale = T, verbose=F )  # NLP-MCMC with those vars only
        hppm[[i]] = names.xx.input [ which( bb $ postMode == 1 ) ] # collect the HPPM vars
      } else
      if( prior == "imom" ) 
      {
        bb = modelSelection( y, x = x[ , colnames(x) %in% names.xx.input, drop=F ], priorCoef = imomprior(tau=tau), priorDelta = priorDelta, niter = niter, center = T, scale = T, verbose=F )  # NLP-MCMC with those vars only
        hppm[[i]] = names.xx.input [ which( bb $ postMode == 1 ) ] # collect the HPPM vars
      } else
      if( prior == "zellner" ) 
      {
        bb = modelSelection( y, x = x[ , colnames(x) %in% names.xx.input, drop=F ], priorCoef = zellnerprior(tau=tau), priorDelta = priorDelta, niter = niter, center = T, scale = T, verbose=F )  # NLP-MCMC with those vars only
        hppm[[i]] = names.xx.input [ which( bb $ postMode == 1 ) ] # collect the HPPM vars
      } else
      # emom not yet implemented in modelSelection. Also, emomLM won't run with only one x variable.
      if( prior == "horseshoe" ) 
      {
        fit = horseshoe(y, x[ , colnames(x) %in% names.xx.input, drop=F ], method.tau = tau.hs.method, method.sigma = sigma.hs.method)
        fitsel = HS.var.select(fit, y, "intervals")
        hppm[[i]] = names.xx.input [ which( fitsel == 1 ) ]
      }
      if( verbose == T) cat( "selected :", hppm[[i]], "\n")  # print the HPPM vars
    }
  }

  names.xx.not.selected = setdiff( names.xx.input.set, unlist(hppm) )

  return( list( hppm = unlist(hppm), names.not.selected = names.xx.not.selected ) )
}
