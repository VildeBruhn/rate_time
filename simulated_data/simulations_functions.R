# --------------------------- #
# FUNCTIONS FOR simulations.R #
# --------------------------- #

# ----------------------------------------------- #
# Simulate unbiased random walk time series (sim) #
# ----------------------------------------------- #

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


# --------------------------- #
# Calculate darwins (darwins) #
# --------------------------- #

darwins <- function(data) {
  data <- data
  for (i in 1:length(data)){
    data[[i]]$darwins <- abs((tail(data[[i]]$mm, n = 1) - 
                                data[[i]]$mm[1]) / data[[i]]$interval_MY)
    # log transform darwins
    data[[i]]$darwins <- log(data[[i]]$darwins)
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