# -----------------------------
# Functions for URW simulations
# -----------------------------

## requires paleoTS 


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



# check if URW fits the simulated data best
check_fit_URW <- function(model){
  
  # extract AICc values from results
  aicc <- lapply(model, function(x) x[(names(x) %in% c("AICc"))])
  
  # check if the second AICc value (URW) is the lowest
  aicc <- lapply(aicc, function(x) {
    which.min(as.numeric(unlist(x)))
  })
  
  # filter for timeseries with second value as the lowest
  URW_ts <- Filter(function(x) x[[1]] == 2, aicc)
  
  # check if all data fitted URW best
  fit <- isTRUE(length(model) == length(URW_ts))
  
  # make results df
  fit_out <- matrix(nrow = 1, ncol = 2)
  colnames(fit_out) <- c("Fit", "N_URW")
  fit_out[1,1] <- fit
  fit_out[1,2] <- as.character(length(URW_ts))
  return(fit_out)
}

# check how much the AICc discrepancy is between GRW and URW

check_GRW_URW <- function(model, tolerance){
  
  # extract AICc values from results
  aicc <- lapply(model, function(x) x[(names(x) %in% c("AICc"))])
  
  # check if the second or first AICc value (URW) is the lowest
  aicc_char <- lapply(aicc, function(x) {
    which.min(as.numeric(unlist(x)))
  })
  
  # filter for timeseries with second value as the lowest
  aicc_GRW <- Filter(function(x) x[[1]] == 1, aicc_char)
  aicc_URW <- Filter(function(x) x[[1]] == 2, aicc_char)
  
  # keep only data sets where the second AICc value is the lowest
  aicc_add <- mapply(c, aicc, aicc_char, SIMPLIFY = FALSE)
    #adds index of lowest AICc as column in ln_data
  
  GRW_best <- Filter(function(x) x[[2]] == 1, aicc_add)
  URW_best <- Filter(function(x) x[[2]] == 2, aicc_add)
    # keeps only the data frames with index 1 or 2
  
  # check if GRW and URW amounts to total N time series
  if (length(GRW_best) + length(URW_best) == length(model)){
    print("Timeseries are either GRW or URW, not stasis.")
  } else{
    print("Warning: Some timeseries are best fitted to stasis.")
  }
  
  # check how much GRW deviates from URW
  count <- 0
  
  for (i in 1:length((GRW_best))){
    if (GRW_best[[i]]$AICc[2] - tolerance > GRW_best[[i]]$AICc[1]){
      count <- count + 1
    } else {
      count <- count
    }
  }

  output <- NULL
  if(isTRUE(count == 0) == TRUE){
    output <- "GRWs AICc are not unreasonably smaller than URWs"
  } else{
    output <- "GRWs AICc are unreasonably smaller than URWs"
  }
  
  return(output)
}

# check if the range of variances is ok
var_ok <- function(data, vs, rng){
  
  # get vstep
  variance <- lapply(data, function(x) {
    x$vstep = as.numeric(toString(x$parameters[2]))
    return(x)
  })
  
  variance <- lapply(variance, function(x) x[(names(x) %in% c("vstep"))])

  # set range values
  rng_min <- vs - rng
  rng_max <- vs + rng
  
  # check range
  for (i in 1:length(variance)) {
    variance[[i]]$OK <- isTRUE(variance[[i]] > rng_min & variance[[i]] < rng_max)
  }

  # get ok
  ok <- list()
  ok <- lapply(variance, function(x) {
    y <- append(ok, x$OK)
    return(y)
  })
  ok <- unlist(ok)
  
  # check if all are TRUE
  true <- all(ok)
  false <- list()
  
  for (i in 1:length(variance)){
    if (variance[[i]]$OK == FALSE){
      false <- append(false, paste0(i, " is not within tolerance range"))
    } else {
      false <- false
    }
  }
  
  if(length(false) == 0){
    return(true)
  } else {
    return(false)
  }
  
}


# remove timeseries with too much deviation in variance
var_remove <- function(model, data, vs, rng){
  
  data <- data
  
  # get vstep
  variance <- lapply(model, function(x) {
    x$vstep = as.numeric(toString(x$parameters[2]))
    return(x)
  })
  
  variance <- lapply(variance, function(x) x[(names(x) %in% c("vstep"))])
  
  # set range values
  rng_min <- vs - rng
  rng_max <- vs + rng
  
  # check range
  for (i in 1:length(variance)) {
    variance[[i]]$OK <- isTRUE(variance[[i]] > rng_min & variance[[i]] < rng_max)
  }
  
  # get ok
  ok <- list()
  ok <- lapply(variance, function(x) {
    y <- append(ok, x$OK)
    return(y)
  })
  ok <- unlist(ok)
  
  # check if all are TRUE
  new_data <- list()
  false <- list()
  
  for (i in 1:length(variance)){
    if (variance[[i]]$OK == TRUE){
      new_data[[length(new_data)+1]] <- data[[i]]
    } else {
      false <- false
    }
  }
  
  return(new_data)
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


