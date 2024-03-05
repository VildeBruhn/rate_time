# -----------------------------
# Functions for simulations.R
# -----------------------------

## requires paleoTS version 0.5.3 


# ----------------------------------------- #
# Simulate unbiased random walk time series #
# ----------------------------------------- #

sim <- function(i, ns, nn, vs, vp) {
  
  i <- i
  ns <- ns
  nn <- nn
  vs <- vs
  vp <- vp
  
  data <- foreach(i=i) %do% {
    data <- sim.GRW(ms = 0, ns = ns, nn = nn, vs = vs, vp = vp)
    return(data)
  }
  return(data)
}

#----------------- #
# Bind data (bind) #
#----------------- #

bind <- function(data, variance_term, unit_list){
  if (variance_term == "darwins"){
    
    binded <- lapply(data, function(x) x[(names(x) %in% unit_list)])
    binded <- bind_rows(binded, .id="data_frame")
    
    # filter out infinite values
    binded <- binded %>% filter_all(all_vars(!is.infinite(.)))
    
    # remove duplicated data
    binded <- binded[!duplicated(binded),]
    
  } else if (variance_term == "vstep"){
    
    # replace vstep with parameters
    unit_ls <- unit_list[unit_list != "vstep"]
    unit_ls <- append(unit_ls, "parameters")
    new_data <- lapply(data, function(x) x[(names(x) %in% unit_ls)])
    
    # get only chosen units in list
    new_data <- lapply(new_data, function(x) {
      x$vstep = as.numeric(toString(x$parameters[2]))
      return(x)
    })
    
    new_data <- lapply(new_data, function(x) x[(names(x) %in% unit_list)])
    
    # bind to data frame
    binded <- bind_rows(new_data, .id="data_frame")
    
    # filter out infinite values
    binded <- binded %>% filter_all(all_vars(!is.infinite(.))) 
    
    # remove duplicated data
    binded <- binded[!duplicated(binded),]
    
  }else {
    
    print("The chosen unit has do be darwins or vstep")
    
  }
  
  return(binded)
  
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




