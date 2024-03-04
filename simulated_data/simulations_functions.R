# -----------------------------
# Functions for simulations.R
# -----------------------------

## requires paleoTS version 0.5.3 


# simulate data

sim <- function(i, ns, nn, vs, vp) {
  
  i <- i
  ns <- ns
  nn <- nn
  vs <- vs
  vp <- vp
  
  data <- foreach(i=i) %do% {
    data <- sim.GRW(ns = ns, nn = nn, vs = vs, vp = vp)
    return(data)
  }
  return(data)
}


# calculate darwins *** DOESN'T WORK AS A FUNCTION, ONLY MANUALLY ***
sim_darwins <- function(data){
  
  # subset
  #df_darwins <- lapply(data, function(x) x[(names(x) %in% c("mm", "tt", "popID"))])
  df_darwins <- data
  
  # get interval length
  #for (i in 1:length(df_darwins)){
  #  df_darwins[[i]]$interval_MY <- tail(df_darwins[[i]][[4]], 1)
  #}
  
  df_darwins <- lapply(df_darwins, function(x) {
    x$interval_MY <- tail(x$tt, n = 1)
    return(x)
  })
  
  # calculate darwins for all datasets (***measurements are already log transformed)
  for (i in 1:length(df_darwins)){
    df_darwins[[i]]$darwins <- abs((tail(df_darwins[[i]][[1]], n = 1) - df_darwins[[i]]$mm[1]) / df_darwins[[i]]$interval_MY)
    df_darwins[[i]]$darwins <- log(df_darwins[[i]]$darwins)
  }
  
  #df_darwins <- lapply(df_darwins, function(x) {
  #  x$darwins <- abs((tail(x$mm, n = 1) - x$mm[1]) / x$interval_MY)
  #  x$darwins <- log(x$darwins)
  #  return(df_darwins)
  #})
  
  # transform time to log scale
  df_darwins <- lapply(df_darwins, function(x) {
    x$interval_MY <- log(x$interval_MY)
    return(x)
  })
  
  return(df_darwins)
}


