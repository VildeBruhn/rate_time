adeq_stasis <- function (y, nrep = 1000, conf = 0.95, plot = FALSE, omega = NULL) 
{
  x <- y$mm
  v <- y$vv
  n <- y$nn
  tt <- y$tt
  theta <- opt.joint.Stasis(y)$parameters[1]
  if (is.null(omega)) 
    omega <- opt.joint.Stasis(y)$parameters[2]
  lower <- (1 - conf)/2
  upper <- (1 + conf)/2
  obs.auto.corr <- auto.corr(x, model = "stasis")
  obs.runs.test <- runs.test(x, model = "stasis", theta = theta)
  obs.slope.test <- slope.test(x, tt, model = "stasis", theta = theta)
  obs.net.change.test <- net.change.test(x, model = "stasis")
  out.auto <- auto_corr_test_stasis(y, nrep, conf, plot = FALSE, 
                                    theta, omega)
  out.runs <- runs.test.stasis(y, nrep, conf, plot = FALSE, 
                               theta, omega)
  out.slope <- slope.test.stasis(y, nrep, conf, plot = FALSE, 
                                 theta, omega)
  out.net <- net.change.test.stasis(y, nrep, conf, plot = FALSE, 
                                    theta, omega)
  output <- c(as.vector(matrix(unlist(out.auto[[3]]), ncol = 5, 
                               byrow = FALSE)), as.vector(matrix(unlist(out.runs[[3]]), 
                                                                 ncol = 5, byrow = FALSE)), as.vector(matrix(unlist(out.slope[[3]]), 
                                                                                                             ncol = 5, byrow = FALSE)), as.vector(matrix(unlist(out.net[[3]]), 
                                                                                                                                                         ncol = 5, byrow = FALSE)))
  output <- as.data.frame(cbind(c(output[c(1, 6, 11, 16)]), 
                                c(output[c(2, 7, 12, 17)]), c(output[c(3, 8, 13, 18)]), 
                                c(output[c(4, 9, 14, 19)]), c(output[c(5, 10, 15, 20)])), 
                          ncol = 5)
  rownames(output) <- c("auto.corr", "runs.test", "slope.test", 
                        "net.change.test")
  colnames(output) <- c("estimate", "min.sim", "max.sim", "p-value", 
                        "result")
  if (plot == TRUE) {
    par(mfrow = c(2, 2))
    model.names <- c("auto.corr", "runs.test", "slope.test", 
                     "net.change.test")
    plotting.distributions(out.auto$replicates, obs.auto.corr, 
                           model.names[1], xlab = "Simulated data", main = "Autocorrelation")
    plotting.distributions(out.runs$replicates, obs.runs.test, 
                           model.names[2], xlab = "Simulated data", main = "Runs")
    plotting.distributions(out.slope$replicates, obs.slope.test, 
                           model.names[3], xlab = "Simulated data", main = "Fixed variance")
    plotting.distributions(out.net$replicates, obs.net.change.test, 
                           model.names[4], xlab = "Simulated data", main = "Net evolution")
  }
  summary.out <- as.data.frame(c(nrep, conf))
  rownames(summary.out) <- c("replications", "confidence level")
  colnames(summary.out) <- ("Value")
  out <- list(info = summary.out, summary = output)
  return(out)
}
