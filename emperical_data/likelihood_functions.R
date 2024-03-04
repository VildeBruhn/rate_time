# define function with lower and upper output (from evoTS::loglik.surface.URW)
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