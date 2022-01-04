pars.estab <- function(hsamples,parameter){
  out <- do.call(rbind, lapply(hsamples, function(x) {as.vector(x$theta[, parameter, ])
  }))
  return(out)
}