#------------------------
# Data manipulation (dt)
#------------------------

dt <- function(df, tsID){
 
  # split into two lists; one with log transformed and one without
  df_log <- list()
  df_not_log <- list()
  
  df_log <- lapply(df, function(x){
    if(isTRUE(x$log_transformed[1] == "yes") == TRUE){
      df_log$x = x
    }
  })
  
  df_log = df_log[-which(sapply(df_log, is.null))]
  
  df_not_log <- lapply(df, function(x){
    if(isTRUE(x$log_transformed[1] == "no") == TRUE){
      df_not_log$x = x
    }
  })
  
  df_not_log = df_not_log[-which(sapply(df_not_log, is.null))]
  
  # check if oldest is first, change those with youngest first to oldest first
  df_log <- lapply(df_log, function(x) {
    if (isTRUE(x$oldest_first[1] == "no") == TRUE){
      arrange(x, desc(age_MY))
    } else if (isTRUE(x$oldest_first[1] == "yes") == TRUE){
      x <- x
    }
  })
  
  df_not_log <- lapply(df_not_log, function(x) {
    if (isTRUE(x$oldest_first[1] == "no") == TRUE){
      arrange(x, desc(age_MY))
    } else if (isTRUE(x$oldest_first[1] == "yes") == TRUE){
      x <- x
    }
  })
  
  ## CREATE paleoTS OBJECTS AND TRANSFORM DATA ##
  
  # apply the as.paleoTS function to all data frames in the list
  data_log <- lapply(df_log, function(x) {
    as.paleoTS(mm = x$trait_mean, vv = x$trait_var, nn = x$N, tt = x$age_MY, oldest = "first")
  })
  
  data_not_log <- lapply(df_not_log, function(x) {
    as.paleoTS(mm = x$trait_mean, vv = x$trait_var, nn = x$N, tt = x$age_MY, oldest = "first")
  })
  
  # approximate log-transformation of the data sets that do not have log transformed values
  data_not_log_log <- lapply(data_not_log, ln.paleoTS)
  #nan <- as.data.frame(which(rapply(data_not_log_log, is.nan)))
  
  # join log and not log data sets again now that all are on log scale
  data_log <- mapply(c, data_log, df_log, SIMPLIFY = FALSE)
  data_not_log_log <- mapply(c, data_not_log_log, df_not_log, SIMPLIFY = FALSE)
  ln_data <- c(data_log, data_not_log_log)
  
  return(ln_data)
}

#-----------------------------
# Calculate darwins (darwins)
#-----------------------------

darwins <- function(data){
 
  # subset
  #df_darwins <- lapply(data, function(x) x[(names(x) %in% c("mm", "tt", "popID"))])
  
  # get interval length
  df_darwins <- lapply(data, function(x) {
    x$interval_MY <- x$interval_MY[1]
    return(x)
  })

  # calculate darwins for all datasets (***measurements are already log transformed)
  df_darwins <- lapply(df_darwins, function(x) {
    x$darwins <- abs((tail(x$mm, n = 1) - x$mm[1]) / x$interval_MY)
    x$darwins <- log(x$darwins)
    return(x)
  })

  # transform time to log scale
  df_darwins <- lapply(df_darwins, function(x) {
    x$interval_MY <- log(x$interval_MY)
    return(x)
  })
  
  return(df_darwins)
}

#-----------------------------
# Bind data (bind)
#-----------------------------

bind <- function(data, variance_term, unit_list){
  if (variance_term == "darwins"){
    binded <- lapply(data, function(x) x[(names(x) %in% unit_list)])
    # change nn to mean nn
    binded <- lapply(binded, function(x) {
      x$nn = sum(x$nn)/length(x$nn)
      return(x)
    })
    binded <- bind_rows(binded, .id="data_frame")
    # filter out infinite values
    binded <- binded %>% filter_all(all_vars(!is.infinite(.)))
    # remove duplicated data
    binded <- binded[!duplicated(binded),]
  } else if (variance_term == "vstep"){
    # repalce vstep with parameters
    unit_ls <- unit_list[unit_list != "vstep"]
    unit_ls <- append(unit_ls, "parameters")
    new_data <- lapply(data, function(x) x[(names(x) %in% unit_ls)])
    
    # get only chosen units in list
    new_data <- lapply(new_data, function(x) {
      x$vstep = as.numeric(toString(x$parameters[2]))
      return(x)
    })
    new_data <- lapply(new_data, function(x) x[(names(x) %in% unit_list)])
    
    # change nn to mean nn
    new_data <- lapply(new_data, function(x) {
      x$nn = sum(x$nn)/length(x$nn)
      return(x)
    })
    
    # bind to dataframe
    binded <- bind_rows(new_data, .id="data_frame")
    # filter out infinite values
    binded <- binded %>% filter_all(all_vars(!is.infinite(.))) 
    # remove duplicated data
    binded <- binded[!duplicated(binded),]
  } else if (variance_term == "omega"){
    # keep only parameters, time interval and popID
    unit_ls <- unit_list[unit_list != "omega"]
    unit_ls <- append(unit_ls, "parameters")
    new_data <- lapply(data, function(x) x[(names(x) %in% unit_ls)])
    
    # get only chosen units in list
    new_data <- lapply(new_data, function(x) {
      x$omega = as.numeric(toString(x$parameters[2]))
      return(x)
    })
    new_data <- lapply(new_data, function(x) x[(names(x) %in% unit_list)])
    
    # change nn to mean nn
    new_data <- lapply(new_data, function(x) {
      x$nn = sum(x$nn)/length(x$nn)
      return(x)
    })
    
    # bind to dataframe
    binded <- bind_rows(new_data, .id="data_frame")
    # filter out infinite values
    binded <- binded %>% filter_all(all_vars(!is.infinite(.)))
    # remove duplicated data
    binded <- binded[!duplicated(binded),]
  } else {
    print("The chosen unit has do be darwins, vstep or omega")
  }
  
  return(binded)

}

#----------------------------
# Unbiased random walk (URW)
#----------------------------

URW <- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", 
                 hess = FALSE) 
{
  if (pool) 
    y <- pool.var(y, ret.paleoTS = TRUE)
  if (y$tt[1] != 0) 
    stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  p0 <- array(dim = 2)
  p0[1] <- y$mm[1]
  p0[2] <- min(c(mle.URW(y), 1e-07))
  names(p0) <- c("anc", "vstep")
  if (is.null(cl$ndeps)) 
    cl$ndeps <- abs(p0/10000)
  cl$ndeps[cl$ndeps == 0] <- 1e-08
  if (meth == "L-BFGS-B") 
    w <- optim(p0, fn = paleoTS:::logL.joint.URW, control = cl, method = meth, 
               lower = c(NA, 0), hessian = hess, y = y)
  else w <- optim(p0, fn = paleoTS:::logL.joint.URW, control = cl, method = meth, 
                  hessian = hess, y = y)
  if (hess) 
    w$se <- sqrt(diag(-1 * solve(w$hessian)))
  else w$se <- NULL
  wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = "URW", 
                      method = "Joint", K = 2, n = length(y$mm), se = w$se)
  return(wc)
}


#----------------------------------
# Get timeseries that fit URW best
#----------------------------------

get_URW_best <- function(data){

  # extract AICc values from URW on all results
  aicc <- lapply(data, function(x) x[(names(x) %in% c("AICc"))])
  
  # check if the second AICc value (URW) is the lowest
  aicc <- lapply(aicc, function(x) {
    which.min(as.numeric(unlist(x)))
  })
  
  # filter for timeseries with second value as the lowest
  aicc_URW <- Filter(function(x) x[[1]] == 2, aicc)
  
  # keep only data sets where the second AICc value is the lowest
  data_AICc <- mapply(c, data, aicc, SIMPLIFY = FALSE)
    #adds index of lowest AICc as column in ln_data
  
  data_AICc_URW <- Filter(function(x) x[[52]] == 2, data_AICc)
  # keeps only the data frames with index 2 (URW)
  
  data_URW <- lapply(data_AICc_URW, function(x) { x[[52]] <- NULL; x })
  # final list of data sets with lowest AICc for URW
  
  return(data_URW)
}

#-------------------------------------
# Get timeseries that fit stasis best
#-------------------------------------

get_stasis_best <- function(model_test){
  
  data <- model_test_meta
  
  # extract AICc values from URW on all results
  aicc <- lapply(data, function(x) x[(names(x) %in% c("AICc"))])
  
  # check if the third AICc value (URW) is the lowest
  aicc <- lapply(aicc, function(x) {
    which.min(as.numeric(unlist(x)))
  })
  
  # filter for timeseries with third value as the lowest
  aicc_URW <- Filter(function(x) x[[1]] == 3, aicc)
  
  # keep only data sets where the third AICc value is the lowest
  data_AICc <- mapply(c, data, aicc, SIMPLIFY = FALSE)
  # adds index of lowest AICc as column in ln_data
  
  data_AICc_URW <- Filter(function(x) x[[52]] == 3, data_AICc)
  # keeps only the data frames with index 2 (URW)
  
  data_URW <- lapply(data_AICc_URW, function(x) { x[[52]] <- NULL; x })
  # final list of data sets with lowest AICc for URW
  
  return(data_URW)
}


#---------------------------------------------------
# Get timeseries that fit URW best and are adequate
#---------------------------------------------------

adequate <- function(data){
  
  # get only passed or failed adequacy test into list
  data_URW_adeq <- lapply(data, function(x) {
    x$result = toString(x$summary[5])
    return(x)
  })
  
  # filter out those that did not pass the adequacy test
  data_URW_adeq_passed <- Filter(function(x) x$result ==  "c(\"PASSED\", \"PASSED\", \"PASSED\")", data_URW_adeq)
  
  return(data_URW_adeq_passed)
  
}

adequate_stasis <- function(data){
  
  # get only passed or failed adequacy test into list
  data_stasis_adeq <- lapply(data, function(x) {
    x$results = toString(x$summary[5])
    return(x)
  })
  
  # filter out those that did not pass the adequacy test
  data_stasis_adeq_passed <- Filter(function(x) x$results ==  "c(\"PASSED\", \"PASSED\", \"PASSED\", \"PASSED\")", data_stasis_adeq)
  
  return(data_stasis_adeq_passed)
  
}

#---------------
# Extras
#---------------

# feilsøking
#for (i in 1:length(ln_data)){
#  print(i)
#  fit.all.univariate(ln_data[[i]], pool = TRUE)
#}
#ln_data[253]


