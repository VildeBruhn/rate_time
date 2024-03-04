
# ----------------- #
#  SIMULATED DATA   #
# ----------------- #

# packages needed
#install.packages(evoTS)        # version 1.0.2
#install.packages(paleoTS)      # version 0.5.3
#install.packages(wesanderson)  # colors for figures
#install.packages(tidyverse)
#install.packages(ggpmisc)

library(evoTS)
library(paleoTS)
library(wesanderson)
library(tidyverse)
library(ggpmisc)


# set working directory and import function script
setwd("[working directory]")
source(file = "[working directory]/simulations_functions.R")

#*source(file = "/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/rates_time_functions.R")


#### SIMULATE DATA ###

# example of how to simulate data (sim function from simulations_functions.R)
data <- sim(i = 1:100, ns = 100, nn = rep(50,100), vs = 1, vp = 0.1)

## load simulated data used in the article
load("[working directory]/simulations.Rdata")


### CUT SIMULATIONS RANDOMLY IN TWO ###

## example of how to cut data (load article data below)

# duplicate data list
new_data <- c(data, data)

# generate random sequences
start <- rep(0, length(new_data))
stop <- rep(0, length(new_data))
range_ts <- length(new_data[[1]]$tt)
len_data <- length(new_data)/2
for (i in 1:len_data){
  cut <- sample(10:(range_ts-10), 1, replace=T)
  start[i] <- 0
  stop[i] <- cut
  start[i+len_data] <- cut + 1
  stop[i+len_data] <- range_ts
}

# cut data
cut_data <- list()
for (i in 1:length(new_data)){
  cut_data[[i]] <- sub.paleoTS(new_data[[i]], ok = start[i]:stop[i], reset.time = FALSE)
}

## load cut data used in the article
load("[working directory]/cut_data.Rdata")


### MAKE INCOMPLETE DATA BY REMOVING POPULATIONS RANDOMLY ###


# generate random k for sub.paleoTS() between 0.1-1
k <- runif(length(cut_data), min=0.1, max=1)

# make empty list
incompl_data <- list()

# remove populations
for (i in 1:length(cut_data)){
  incompl_data[[i]] <- sub.paleoTS(cut_data[[i]], k = k[i], reset.time = FALSE)
}

## load cut and incomplete data used in the article
load("[working directory]/incompl_data.Rdata")


### MAKE BIASED DATA BY REMOVING TIME CHUNKS OF POPULATIONS RANDOMLY ###

# biased sub sampling
biased_data <- vector(mode = "list", length = length(incompl_data))
for (i in 1:length(incompl_data)){
  
  if (length(incompl_data[[i]]$tt) < 10){
    # < 10 to assure enough populations to estimate vstep
    # keep time series with little data as is
    biased_data[[i]]$tt <- incompl_data[[i]]$tt
    biased_data[[i]]$mm <- incompl_data[[i]]$mm
    biased_data[[i]]$nn <- incompl_data[[i]]$nn
    biased_data[[i]]$vv <- incompl_data[[i]]$vv
  } else {
    # generate random start
    start <- floor(runif(1, 1, length(incompl_data[[i]]$tt)))
    while(start > (length(incompl_data[[i]]$tt) - 5)){
      start <- floor(runif(1, 1, length(incompl_data[[i]]$tt)))
      if (start <= (length(incompl_data[[i]]$tt) - 5)){
        break
      }
    }
  
    # generate random stop
    stop <- floor(runif(1, (start + 1), length(sub_data2[[i]]$tt)))
  
    # cut out data points
    sub_data3[[i]]$tt <- sub_data2[[i]]$tt[-(start:stop)]
    sub_data3[[i]]$mm <- sub_data2[[i]]$mm[-(start:stop)]
    sub_data3[[i]]$nn <- sub_data2[[i]]$nn[-(start:stop)]
    sub_data3[[i]]$vv <- sub_data2[[i]]$vv[-(start:stop)]
  }
}


#save(sub_data3, file = "./data_temp/sub_data3.Rdata")

# Load subset data 
load("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/data_temp/sub_data.Rdata")
load("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/data_temp/sub_data2.Rdata")
load("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/data_temp/sub_data3.Rdata")


# --------------------------
# Run darwins on subset data
# --------------------------
darwins_sub <-sub_data
darwins_sub <- lapply(darwins_sub, function(x) {
  x$interval_MY <- (tail(x$tt, n = 1)) - x$tt[1]
  return(x)
})
# calculate darwins for all datasets (***measurements are already log transformed)
for (i in 1:length(darwins_sub)){
  darwins_sub[[i]]$darwins <- abs((tail(darwins_sub[[i]]$mm, n = 1) - darwins_sub[[i]]$mm[1]) / darwins_sub[[i]]$interval_MY)
  darwins_sub[[i]]$darwins <- log(darwins_sub[[i]]$darwins)
}

# transform time to log scale
darwins_sub <- lapply(darwins_sub, function(x) {
  x$interval_MY <- log(x$interval_MY)
  return(x)
})

darwins_sub2 <-sub_data2
darwins_sub2 <- lapply(darwins_sub2, function(x) {
  if (length(x$tt) == 1){
      x$interval_MY <- x$tt
    } else{
      x$interval_MY <- (tail(x$tt, n = 1)) - x$tt[1]
    }
  return(x)
})

# calculate darwins for all datasets (***measurements are already log transformed)
for (i in 1:length(darwins_sub2)){
  darwins_sub2[[i]]$darwins <- abs((tail(darwins_sub2[[i]]$mm, n = 1) - darwins_sub2[[i]]$mm[1]) / darwins_sub2[[i]]$interval_MY)
  darwins_sub2[[i]]$darwins <- log(darwins_sub2[[i]]$darwins)
}
# transform time to log scale
darwins_sub2 <- lapply(darwins_sub2, function(x) {
  x$interval_MY <- log(x$interval_MY)
  return(x)
})

darwins_sub3 <-sub_data3
darwins_sub3 <- lapply(darwins_sub3, function(x) {
  x$interval_MY <- (tail(x$tt, n = 1)) - x$tt[1]
  return(x)
})
# calculate darwins for all datasets (***measurements are already log transformed)
for (i in 1:length(darwins_sub3)){
  darwins_sub3[[i]]$darwins <- abs((tail(darwins_sub3[[i]]$mm, n = 1) - darwins_sub3[[i]]$mm[1]) / darwins_sub3[[i]]$interval_MY)
  darwins_sub3[[i]]$darwins <- log(darwins_sub3[[i]]$darwins)
}
# transform time to log scale
darwins_sub3 <- lapply(darwins_sub3, function(x) {
  x$interval_MY <- log(x$interval_MY)
  return(x)
})


# --------------------------
# Run URW on subset data
# --------------------------

#doParallel::registerDoParallel(cl = my_cluster)

# make paleoTS objects
sub_data_paleo <- lapply(sub_data, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt, oldest = "first")
})

#save(file = "./data_temp/sub_data_paleo.Rdata", sub_data_paleo)

sub_data2_paleo <- lapply(sub_data2, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt, oldest = "first")
})

#save(file = "./data_temp/sub_data2_paleo.Rdata", sub_data2_paleo)

sub_data3_paleo <- lapply(sub_data3, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt, oldest = "first")
})

#save(file = "./data_temp/sub_data3_paleo.Rdata", sub_data3_paleo)


#URW_sub <- mclapply(sub_data, opt.joint.URW) # best to do on HPC with high number of simulations
#URW_sub2 <- lapply(sub_data2, opt.joint.URW) # best to do on HPC with high number of simulations
#URW_sub3 <- lapply(sub_data3_paleo, opt.joint.URW) # best to do on HPC with high number of simulations


# Load URW on subset data from HPC
load("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/data_temp/URW_sub.Rdata")
load("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/data_temp/URW_sub2.Rdata")
load("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/data_temp/URW_sub3.Rdata")


# --------------------------
# Prepare plotting and plot
# --------------------------

# append time to data
#URW <- mapply(c, URW_joint, data, SIMPLIFY = FALSE)
URW_sub <- mapply(c, URW_sub, sub_data, SIMPLIFY = FALSE)
URW_sub2 <- mapply(c, URW_sub2, sub_data2, SIMPLIFY = FALSE)
URW_sub3 <- mapply(c, URW_sub3, sub_data3, SIMPLIFY = FALSE)

# change time to length of time interval
URW_sub <- lapply(URW_sub, function(x) {
  x$tt = tail(x$tt, n = 1) - x$tt[1]
  return(x)
})

URW_sub2 <- lapply(URW_sub2, function(x) {
  if (length(x$tt) == 1){
    x$tt <- x$tt
  } else {
    x$tt <- tail(x$tt, n = 1) - x$tt[1]
  }
  return(x)
})

URW_sub3 <- lapply(URW_sub3, function(x) {
  x$tt = tail(x$tt, n = 1) - x$tt[1]
  return(x)
})

# bind to dataframe
#URW_bind <- bind(URW, variance_term = "vstep", unit_list = c("tt", "vstep"))
URW_sub_bind <- bind(URW_sub, variance_term = "vstep", unit_list = c("tt", "vstep"))
URW_sub2_bind <- bind(URW_sub2, variance_term = "vstep", unit_list = c("tt", "vstep"))
URW_sub3_bind <- bind(URW_sub3, variance_term = "vstep", unit_list = c("tt", "vstep"))
darwins_sub_bind <- bind(darwins_sub, variance_term = "darwins", unit_list = c("interval_MY", "darwins"))
darwins_sub2_bind <- bind(darwins_sub2, variance_term = "darwins", unit_list = c("interval_MY", "darwins"))
darwins_sub3_bind <- bind(darwins_sub3, variance_term = "darwins", unit_list = c("interval_MY", "darwins"))

# Only to plot increase in trait variance across time for supplementary figure
#data <- sim(i = 1:100, ns = 1000, nn = rep(50,1000), vs = 1, vp = 0.1)
#URW_bind <- bind(data, variance_term = "vstep", unit_list = c("mm", "tt"))

#pdf(file = "./results/URW_sim_increased_var.pdf")
#ggplot() +
#  geom_line(data = URW_bind, aes(x = tt, y = mm, group = data_frame), color = c(wes_palette("Rushmore1")[3])) +
#  theme_classic() + ylab("TRAIT MEAN") + xlab("TIME") +
#  theme(axis.title = element_text(size = 15)) +
#  theme(axis.text = element_text(size = 13))
#dev.off()


# log-transform tt and vstep
#URW_bind$tt <- log(URW_bind$tt)
#URW_bind$vstep <- log(URW_bind$vstep)
URW_sub_bind$tt <- log(URW_sub_bind$tt)
URW_sub_bind$vstep <- log(URW_sub_bind$vstep)
URW_sub2_bind$tt <- log(URW_sub2_bind$tt)
URW_sub2_bind$vstep <- log(URW_sub2_bind$vstep)
URW_sub3_bind$tt <- log(URW_sub3_bind$tt)
URW_sub3_bind$vstep <- log(URW_sub3_bind$vstep)


# linerar regression
summary(lm(vstep ~ tt, URW_sub_bind))
summary(lm(vstep ~ tt, URW_sub2_bind))
summary(lm(vstep ~ tt, URW_sub3_bind))
summary(lm(darwins ~ interval_MY, darwins_sub_bind))
summary(lm(darwins ~ interval_MY, darwins_sub2_bind))
summary(lm(darwins ~ interval_MY, darwins_sub3_bind))

# plot and write to pdf 

pdf(file = "./results/URW_simulations_cut.pdf")
plot_sub <- ggplot(URW_sub_bind, aes(tt, vstep)) + ylim(c(-1.4,1)) + xlim(c(1.1,7)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.010, slope = -0.004, linewidth = 0.7) + 
  ggtitle(expression(bold("Cut time series"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = 2.5, y = -0.6, parse = TRUE, label="italic(y)==-0.010-0.002~italic(x)", size = 5) +
  annotate("text", x = 2.5, y = -0.8, parse = TRUE, label="italic(SE)=='' %+-% '0.003'", size = 5) +
  annotate("text", x = 2.5, y = -1, parse = TRUE, label = "italic(R)^2== 0.000", size = 5) +
  annotate("text", x = 2.5, y = -1.2, parse = TRUE, label = "italic(n)== 2000", size = 5) +
  geom_rect(aes(xmin = 1.1, xmax = 3.88, ymin = -1.4, ymax = -0.4), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))
 
print(plot_sub)
dev.off()



pdf(file = "./results/URW_simulations_mis.pdf")
plot_sub2 <- ggplot(URW_sub2_bind, aes(tt, vstep)) + ylim(c(-3.4,1.7)) + xlim(c(0.2,7.3)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.052, slope = 0.006, linewidth = 0.7) + 
  ggtitle(expression(bold("Incomplete time series"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = 2, y = -1.5, parse = TRUE, label="italic(y)==-0.052+0.006~italic(x)", size = 5) +
  annotate("text", x = 2, y = -2, parse = TRUE, label="italic(SE)=='' %+-% '0.005'", size = 5) +
  annotate("text", x = 2, y = -2.5, parse = TRUE, label = "italic(R)^2== 0.000", size = 5) +
  annotate("text", x = 2, y = -3, parse = TRUE, label = "italic(n)== 2000", size = 5) +
  geom_rect(aes(xmin = 0.2, xmax = 3.81, ymin = -1.05, ymax = -3.4), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_sub2)
dev.off()



pdf(file = "./results/URW_simulations_biased.pdf")
plot_sub3 <- ggplot(URW_sub3_bind, aes(tt, vstep)) + ylim(c(-6.5,1.7)) + xlim(c(1.089,7)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.185, slope = 0.028, linewidth = 0.7) + 
  ggtitle(expression(bold("Biased time series"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = 5.5, y = -3.7, parse = TRUE, label="italic(y)==-0.185+0.028~italic(x)", size = 5) +
  annotate("text", x = 5.5, y = -4.4, parse = TRUE, label="italic(SE)=='' %+-% '0.005'", size = 5) +
  annotate("text", x = 5.5, y = -5.1, parse = TRUE, label = "italic(R)^2== 0.013", size = 5) +
  annotate("text", x = 5.5, y = -5.85, parse = TRUE, label = "italic(n)== 1984", size = 5) +
  geom_rect(aes(xmin = 4, xmax = 7, ymin = -3, ymax = -6.5), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_sub3)
dev.off()


# print together
library(gridExtra)
pdf(width = 11.50, height = 5.50, file = "./results/URW_simulations_cut_mis.pdf")
grid.arrange(plot_sub, plot_sub2, nrow = 1)
dev.off()


pdf(file = "./results/darwins_simulations_cut.pdf")
plot_sub4 <- ggplot(darwins_sub_bind, aes(interval_MY, darwins)) + ylim(c(-13.9,0)) + xlim(c(1.2,7.1)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.609, slope = -0.499, linewidth = 0.7) + 
  ggtitle(expression(bold("Cut time series"))) +
  ylab(expression(paste("Log ", italic(darwins)))) + xlab(expression("Log time")) +
  annotate("text", x = 2.7, y = -9, parse = TRUE, label="italic(y)==-0.609-0.499~italic(x)", size = 5) +
  annotate("text", x = 2.7, y = -10.2, parse = TRUE, label="italic(SE)=='' %+-% '0.026'", size = 5) +
  annotate("text", x = 2.7, y = -11.5, parse = TRUE, label = "italic(R)^2== 0.155", size = 5) +
  annotate("text", x = 2.7, y = -12.9, parse = TRUE, label = "italic(n)== 2000", size = 5) +
  geom_rect(aes(xmin = 1.2, xmax = 4.2, ymin = -13.9, ymax = -8.1), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_sub4)
dev.off()


pdf(file = "./results/darwins_simulations_mis.pdf")
plot_sub5 <- ggplot(darwins_sub2_bind, aes(interval_MY, darwins)) + ylim(c(-13.9,0)) + xlim(c(0.2,7.1)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.754, slope = -0.479, linewidth = 0.7) + 
  ggtitle(expression(bold("Incomplete time series"))) +
  ylab(expression(paste("Log ", italic(darwins)))) + xlab(expression("Log time")) +
  annotate("text", x = 2, y = -9, parse = TRUE, label="italic(y)==-0.754-0.479~italic(x)", size = 5) +
  annotate("text", x = 2, y = -10.2, parse = TRUE, label="italic(SE)=='' %+-% '0.026'", size = 5) +
  annotate("text", x = 2, y = -11.5, parse = TRUE, label = "italic(R)^2== 0.141", size = 5) +
  annotate("text", x = 2, y = -12.9, parse = TRUE, label = "italic(n)== 1999", size = 5) +
  geom_rect(aes(xmin = 0.2, xmax = 3.8, ymin = -13.9, ymax = -8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_sub5)
dev.off()

pdf(file = "./results/darwins_simulations_biased.pdf")
plot_sub6 <- ggplot(darwins_sub3_bind, aes(interval_MY, darwins)) + ylim(c(-16,0.24)) + xlim(c(2,7.1)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.781, slope = -0.475, linewidth = 0.7) + 
  ggtitle(expression(bold("Biased time series"))) +
  ylab(expression(paste("Log ", italic(darwins)))) + xlab(expression("Log time")) +
  annotate("text", x = 3.5, y = -10.5, parse = TRUE, label="italic(y)==-0.781-0.475~italic(x)", size = 5) +
  annotate("text", x = 3.5, y = -12, parse = TRUE, label="italic(SE)=='' %+-% '0.027'", size = 5) +
  annotate("text", x = 3.5, y = -13.5, parse = TRUE, label = "italic(R)^2== 0.130", size = 5) +
  annotate("text", x = 3.5, y = -15.15, parse = TRUE, label = "italic(n)== 1984", size = 5) +
  geom_rect(aes(xmin = 2.1, xmax = 4.9, ymin = -16, ymax = -9.5), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_sub6)
dev.off()

# print together
pdf(width = 11.50, height = 5.50, file = "./results/darwins_simulations_cut_mis.pdf")
grid.arrange(plot_sub4, plot_sub5, nrow = 1)
dev.off()

# print together
pdf(width = 11, height = 8.50, file = "./results/URW_darwins_simulations_cut_mis.pdf")
grid.arrange( plot_sub4, plot_sub5, plot_sub, plot_sub2, nrow = 2)
dev.off()

# print together
pdf(width = 11.50, height = 5.50, file = "./results/URW_darwins_simulations_biased.pdf")
grid.arrange(plot_sub3, plot_sub6, nrow = 1)
dev.off()


# all plots
pdf(width = 12, height = 12, file = "./results/darwins_URW_sim.pdf")
grid.arrange(plot_sub4, plot_sub, plot_sub5, plot_sub2,
             plot_sub6, plot_sub3, nrow = 3)
dev.off()


