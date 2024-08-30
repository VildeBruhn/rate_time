# --------------------------- #
# FUNCTIONS FOR model_ident.R #
# --------------------------- #


# ---------------------------------- #
# Likelihood surface (like_surf_URW) #
# ---------------------------------- #

# define function with lower and upper output (modified from evoTS::loglik.surface.URW)
like_surf_URW <- function (y, vstep.vec, header_plot, save_plot, pool = TRUE) 
{
  x <- NULL
  if (min(vstep.vec) < 0) 
    stop("the vstep parameter cannot take negative values. The smallest value is 0")
  y <- paleoTS::pool.var(y, ret.paleoTS = TRUE)
  y$tt <- y$tt - min(y$tt)
  n <- length(y$mm)
  anc <- paleoTS::opt.joint.URW(y)$parameters[1]
  loglik <- rep(NA, length(vstep.vec))
  vstep.limit <- rep(NA, length(vstep.vec))
  for (i in 1:length(vstep.vec)) {
    vstep <- vstep.vec[i]
    VV <- vstep * outer(y$tt, y$tt, FUN = pmin)
    diag(VV) <- diag(VV) + y$vv/y$nn
    M <- rep(anc, n)
    loglik[i] <- mvtnorm::dmvnorm(y$mm, mean = M, sigma = VV, 
                                  log = TRUE)
  }
  loglik <- loglik - max(loglik)
  loglik <- loglik + 2
  loglik[loglik < 0] <- 0
  for (i in 1:length(vstep.vec)) {
    if (all(loglik[i] == 0)) 
      vstep.limit[i] <- NA
    else vstep.limit[i] <- vstep.vec[i]
  }
  out <- matrix(c(min(na.exclude(vstep.limit)), max(na.exclude(vstep.limit))), 
                ncol = 2, byrow = TRUE)
  colnames(out) <- c("lower", "upper")
  rownames(out) <- "vstep"
  #print(out)
  #pdf(file = save_plot)
  like_surf_plot <- plot(loglik ~ vstep.vec, type = "l", col = "black", lwd = 3, 
       xlab = "vstep", ylab = "log-likelihood", cex.lab = 1.2,
       main = paste("tsID = ", header_plot))
  like_surf_plot
  #dev.off()
  like_surf_plot
  return(out)
}


# --------------------------------------------------------------- #
# Search likelihood space and plot landscapes (search_likelihood) #
# --------------------------------------------------------------- #

# search likelihood space with plot output
search_likelihood <- function(data, data_list, header_plot, save_plot){
  
  # set up dataframe
  likelihood <- as.data.frame(matrix(nrow = length(data$data_frame), ncol = 3))
  colnames(likelihood) <- c("ts", "lower", "upper")
  likelihood$ts <- data$data_frame
  
  # search
  for (i in 1:length(data$vstep)){
    if (data$vstep[i] < 0.01){
      u <- 0.01
      l <- 0
        vstep.vec <- seq(l, u, 0.00001)
      lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
      upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
      while(upper == u){
        u <- u + 0.01
        vstep.vec <- seq(l, u, 0.00001)
        upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf") )[1,2]
        like_surf_plot <- 
        if (upper < u){
          break
        }
      }
      likelihood$lower[i] <- lower
      likelihood$upper[i] <- upper
    } else if (between(data$vstep[i], 0.01, 0.1) == TRUE){
      u <- 0.1
      l <- 0
      vstep.vec <- seq(l, u, 0.00001)
      lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
      upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
      while(upper == u){
        u <- u + 0.3
        vstep.vec <- seq(l, u, 0.00001)
        upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
        if (upper < u){
          break
        }
      }
      likelihood$lower[i] <- lower
      likelihood$upper[i] <- upper
    }
      else if (between(data$vstep[i], 0.1, 1) == TRUE){
      u <- 1
      l <- 0
      vstep.vec <- seq(l, u, 0.0001)
      lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
      upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
      while(upper == u){
        u <- u + 1
        vstep.vec <- seq(l, u, 0.0001)
        upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
        if (upper < u){
          break
        }
      }
      likelihood$lower[i] <- lower
      likelihood$upper[i] <- upper
    } else if (between(data$vstep[i], 1, 10) == TRUE){
      u <- 2
      l <- 1
      vstep.vec <- seq(l, u, 0.0001)
      lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
      upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
      while(lower == l){
        if (l > 0){
          l <- l - 0.1
          vstep.vec <- seq(l, u, 0.0001)
          lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
        }
        if (lower == 0 | lower > l){
          break
        }
      }
      while(upper == u){
        u <- u + 2
        vstep.vec <- seq(l, u, 0.0001)
        upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
        if (upper < u){
          break
        }
      }
      likelihood$lower[i] <- lower
      likelihood$upper[i] <- upper
    } else if (between(data$vstep[i], 10, 100) == TRUE){
      u <- 100
      l <- 10
      vstep.vec <- seq(l, u, 0.001)
      lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
      upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
      while(lower == l){
        if (l > 0){
          l <- l - 1
          vstep.vec <- seq(l, u, 0.001)
          lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
        }
        if (lower == 0 | lower > l){
          break
        }
      }
      while(upper == u){
        u <- u + 100
        vstep.vec <- seq(l, u, 0.001)
        upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
        if (upper < u){
          break
        }
      }
      likelihood$lower[i] <- lower
      likelihood$upper[i] <- upper
    } else if (data$vstep[i] >= 100){
      u <- 200
      l <- 100
      vstep.vec <- seq(l, u, 0.01)
      lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
      upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
      while(lower == l){
        if (l > 0){
          l <- l - 20
          vstep.vec <- seq(l, u, 0.01)
          lower <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,1]
        }
        if (lower == 0 | lower > l){
          break
        }
      }
      while(upper == u){
        u <- u + 200
        vstep.vec <- seq(l, u, 0.01)
        upper <- like_surf_URW(y = data_list[[i]], vstep.vec = vstep.vec, header_plot = header_plot[i,1], save_plot = paste0(save_plot, "/", header_plot[i,1], ".pdf"))[1,2]
        if (upper < u){
          break
        }
      }
      likelihood$lower[i] <- lower
      likelihood$upper[i] <- upper
    } else {
      likelihood <- likelihood
    }
  }
  return(likelihood)
}


# ----------------------------------------------------------------------- #
# Estimate vsteps and perform mixed effect linear regressions (vstep_reg) #
# ----------------------------------------------------------------------- #

# estimate new vsteps and regressions
vstep_reg <- function(list_data, df, lw_up_data){
  
  for (i in 1:length(list_data)){
    
    list_data[[i]] <- df
    x <- i
    
    # generate random vstep between lower and upper
    for (i in 1:length(lw_up_data$ts)){
      list_data[[x]]$vstep[i] <- runif(1, min=lw_up_data$lower[i], max=lw_up_data$upper[i])
    }
    # log transform tt and vstep
    list_data[[x]]$tt <- log(list_data[[x]]$tt)
    list_data[[x]]$vstep <- log(list_data[[x]]$vstep)
    
    # remove infinite values
    list_data[[x]] <- list_data[[x]] %>% filter_all(all_vars(!is.infinite(.)))
    
    # mixed effect linear regression
    lmer <- lmer(vstep ~ tt + (1|popID), list_data[[x]], weights = 1/nn)
    list_data[[x]]$estimate_1[1] <- tidy(lmer, conf.int = TRUE)$estimate[1]
    list_data[[x]]$estimate_2[1] <- tidy(lmer, conf.int = TRUE)$estimate[2]
    list_data[[x]]$SE[1] <- tidy(lmer, conf.int = TRUE)$std.error[2]
  }
  
  return(list_data)
  
}


#----------------- #
# Bind data (bind) #
#----------------- #

bind <- function(data, variance_term, variables){
  
    # replace vstep with parameters
    variables_ls <- variables[variables != "vstep"]
    variables_ls <- append(variables_ls, "parameters")
    new_data <- lapply(data, function(x) x[(names(x) %in% variables_ls)])
    
    # get only chosen variables in list
    new_data <- lapply(new_data, function(x) {
      x$vstep = as.numeric(toString(x$parameters[2]))
      return(x)
    })
    
    new_data <- lapply(new_data, function(x) x[(names(x) %in% variables)])
    
    # change nn to mean nn (for regression)
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
    
  
  return(binded)
  
}
