
# ----------------- #
#  SIMULATED DATA   #
# ----------------- #

# packages needed
#install.packages("evoTS")        # version 1.0.3
#install.packages("paleoTS")      # version 1.0.6
#install.packages("wesanderson")  # colors for figures
#install.packages("tidyverse")
#install.packages("ggpmisc")
#install.packages("foreach")
#install.packages("gridExtra")

library(evoTS)
library(paleoTS)
library(wesanderson)
library(tidyverse)
library(ggpmisc)
library(foreach)
library(gridExtra)


######################################
## REMEMBER TO CHANGE PATH TO FILES ## 
######################################

PATH = "./GitHub/rate_time/simulated_data/"

# import functions
source(paste0(PATH, "/simulations_functions.R"))

# ------------- #
# Simulated data #
# ------------- #

# example of how to simulate unbiased random walk time series 
# (sim function from simulations_functions.R)
data_sim <- sim(i = 1:1000, ns = 1000, nn = rep(50,1000), vs = 1, vp = 0.1)

## load simulated data used in the article
load(paste0(PATH, "/data_sim.Rdata"))


# --------------------------------------- #
# Cut simulations to mimic time averaging #
# --------------------------------------- #

# duplicate data list
new_data <- c(data_sim)

# generate time averaged intervals

start_list <- rep(0,1000)
stop_list <- rep(0,1000)

periods <- c(1,201,401,601,801,1001) 
intervals <- c(10,20,50,80,100)

for (i in 1:length(intervals)){
  start <- periods[i] + intervals[i]
  stop <- periods[i+1]-1
  cut <- sample(start:stop, 200, replace=T)
  stop_list[(start-intervals[i]):stop] <- cut
  start_list[(start-intervals[i]):stop] <- cut-intervals[i]
}

# cut data
cut_data <- list()
for (i in 1:length(new_data)){
  cut_data[[i]] <- sub.paleoTS(new_data[[i]], ok = start_list[i]:stop_list[i], 
                               reset.time = FALSE)
}

# ------------------------------------- #
#  Calculate darwins for simulated data #
# ------------------------------------- #

# copy cut dataset
darwins_cut <- cut_data

# calculate length of time interval
darwins_cut <- lapply(darwins_cut, function(x) {
  x$interval_MY <- (tail(x$tt, n = 1)) - x$tt[1]
  return(x)
})

# calculate darwins (function darwins() from simulations_functions.R)
darwins_cut <- darwins(darwins_cut)

# log transform time
darwins_cut <- lapply(darwins_cut, function(x) {
  x$interval_MY <- log(x$interval_MY)
  return(x)
})


# -------------------------------------------------------- #
# Fit an unbiased random walk to the simulated time series #
# -------------------------------------------------------- #

## This is faster to do on a high performance computer. 
## It is also better to parallelize

## load data used in the article below the example code

# make paleoTS objects
cut_data_paleo <- lapply(cut_data, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt, oldest = "first")
})

# fit unbiased random walk (best to do on HPC with high number of simulations)
URW_cut <- lapply(cut_data_paleo, opt.joint.URW)


# ------------------------- #
# Prepare plotting and plot #
# ------------------------- #

# append time to unbiased random walk data
URW_cut <- mapply(c, URW_cut, cut_data, SIMPLIFY = FALSE)

# calculate time interval
URW_cut <- lapply(URW_cut, function(x) {
  x$tt = tail(x$tt, n = 1) - x$tt[1]
  return(x)
})


# bind to data frame (the bind function is defined in simulations_functions.R)
URW_cut_bind <- bind(URW_cut, variance_term = "vstep", unit_list = c("tt", "vstep"))
darwins_cut_bind <- bind(darwins_cut, variance_term = "darwins", unit_list = c("interval_MY", "darwins"))

# log transform tt and vstep for unbiased random walk fits
URW_cut_bind$tt <- log(URW_cut_bind$tt)
URW_cut_bind$vstep <- log(URW_cut_bind$vstep)

# linerar regressions (manually written in to plots below)
summary(lm(vstep ~ tt, URW_cut_bind))
summary(lm(darwins ~ interval_MY, darwins_cut_bind))

# plot and write to pdf 
plot_URW_cut <- ggplot(URW_cut_bind, aes(tt, vstep)) + ylim(c(-4.7,1)) + xlim(c(1.1,7)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.010, slope = -0.004, linewidth = 0.7) + 
  ggtitle(expression(bold("Cut time series"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = 2.5, y = -1.9, parse = TRUE, label="italic(y)==-0.031+0.005~italic(x)", size = 5) +
  annotate("text", x = 2.5, y = -2.6, parse = TRUE, label="italic(SE)=='' %+-% '0.002'", size = 5) +
  annotate("text", x = 2.5, y = -3.3, parse = TRUE, label = "italic(R)^2== 0.002", size = 5) +
  annotate("text", x = 2.5, y = -4.05, parse = TRUE, label = "italic(n)== 2000", size = 5) +
  geom_rect(aes(xmin = 1.1, xmax = 3.88, ymin = -1.2, ymax = -4.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_URW_cut)



plot_URW_incompl <- ggplot(URW_incompl_bind, aes(tt, vstep)) + ylim(c(-6.7,1.7)) + xlim(c(0.2,7.3)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.052, slope = 0.006, linewidth = 0.7) + 
  ggtitle(expression(bold("Incomplete time series"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = 2, y = -3.9, parse = TRUE, label="italic(y)==-0.146+0.023~italic(x)", size = 5) +
  annotate("text", x = 2, y = -4.6, parse = TRUE, label="italic(SE)=='' %+-% '0.004'", size = 5) +
  annotate("text", x = 2, y = -5.3, parse = TRUE, label = "italic(R)^2== 0.013", size = 5) +
  annotate("text", x = 2, y = -6.05, parse = TRUE, label = "italic(n)== 1997", size = 5) +
  geom_rect(aes(xmin = 0.2, xmax = 3.81, ymin = -3.2, ymax = -6.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_URW_incompl)



plot_URW_biased <- ggplot(URW_biased_bind, aes(tt, vstep)) + ylim(c(-6.5,1.7)) + xlim(c(1.089,7)) +
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = -0.185, slope = 0.028, linewidth = 0.7) + 
  ggtitle(expression(bold("Biased time series"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = 2.55, y = -3.7, parse = TRUE, label="italic(y)==-0.186+0.028~italic(x)", size = 5) +
  annotate("text", x = 2.55, y = -4.4, parse = TRUE, label="italic(SE)=='' %+-% '0.005'", size = 5) +
  annotate("text", x = 2.55, y = -5.1, parse = TRUE, label = "italic(R)^2== 0.013", size = 5) +
  annotate("text", x = 2.55, y = -5.85, parse = TRUE, label = "italic(n)== 1984", size = 5) +
  geom_rect(aes(xmin = 1.1, xmax = 4, ymin = -3, ymax = -6.5), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 15))

print(plot_URW_biased)


plot_darwins_cut <- ggplot(darwins_cut_bind, aes(interval_MY, darwins)) + ylim(c(-13.9,0)) + xlim(c(1.2,7.1)) +
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

print(plot_darwins_cut)


plot_darwins_incompl <- ggplot(darwins_incompl_bind, aes(interval_MY, darwins)) + ylim(c(-13.9,0)) + xlim(c(0.2,7.1)) +
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

print(plot_darwins_incompl)


plot_darwins_biased <- ggplot(darwins_biased_bind, aes(interval_MY, darwins)) + ylim(c(-16,0.24)) + xlim(c(2,7.1)) +
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

print(plot_darwins_biased)


# print all plots together to file
pdf(width = 12, height = 12, file = "[PATH_TO_FOLDER]/simulations.pdf")
grid.arrange(plot_darwins_cut, plot_URW_cut, plot_darwins_incompl, plot_URW_incompl,
             plot_darwins_biased, plot_URW_biased, nrow = 3)
dev.off()


