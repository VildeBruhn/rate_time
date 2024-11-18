# -------------- #
# EMPIRICAL DATA #
# -------------- #

# packages needed
#install.packages("evoTS")        # version 1.0.3
#install.packages("paleoTS")      # version 0.6.1
#install.packages("adePEM")
#install.packages("tidyverse")
#install.packages("wesanderson")  # colors for figures
#install.packages("Matrix")
#install.packages("lme4")
#install.packages("MuMIn")
#install.packages("gridExtra")
#install.packages("broom.mixed")

library(evoTS)
library(paleoTS)
library(adePEM)
library(tidyverse)
library(wesanderson)
library(Matrix)
library(lme4)
library(MuMIn)
library(gridExtra)
library(broom.mixed)


######################################
## REMEMBER TO CHANGE PATH TO FILES ## 
######################################

PATH = "[PATH_TO_DATA_FOLDER]/"

# import functions
source(paste0(PATH, "empirical_functions.R"))



#--------------------------------------- #
# Import and process files from database #  
#--------------------------------------- #

# import time series and metadata
timeseries <- read_delim(paste0(PATH, "timeseries.txt"), col_names = TRUE, delim = "\t")
metadata <- read_delim(paste0(PATH, "metadata.txt"), col_names = TRUE, delim = "\t")

# join data frames
df <- left_join(timeseries, metadata, by = c("tsID"))

# remove modern time series
df <- subset(df, sediment!="none")
df <- subset(df, trait_type!="complex")

# remove time series with less than 5 steps
df <- subset(df, steps >= 5)

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))


# --------------------- #
# Make complete dataset #
# --------------------- #

# process data into right form (function dt from empirical_functions.R)
# with metadata
complete_meta <- dt(df, "tsID")

# without metadata
complete <- lapply(complete_meta, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})


# ------------------------- #
# Make absolute fit dataset #
# ------------------------- #

# Example of how to run adequacy test for an unbiased random walk on complete dataset
# (the warnings are ok)
# The output of the test can vary with ca. +/- one time series (i.e., you will
# not always get the exact same number of time series that pass the adequacy test)
# set.seed to get the same output every time
# (takes time, load model test used in article below)
set.seed(1)
adequacy <- lapply(complete, fit3adequacy.RW, plot = FALSE)

## this will give some error messages with paleoTS v0.6.1,
## circumvent the errors with this approach:
adequacy <- list()
for(i in 1:length(complete)){
  try(adequacy[[i]] <- fit3adequacy.RW(complete[[i]], plot = FALSE))
}

# add metadata
adequacy_meta <- mapply(c, adequacy, complete_meta, SIMPLIFY = FALSE)

# add time series IDs
names_list <- names(complete)
names(adequacy) <- names_list
names(adequacy_meta) <- names_list

# remove time series that cannot be processed by the loglikelihood function
adequacy_meta <-  adequacy_meta[-which(sapply(adequacy, is.null))]
adequacy <- adequacy[-which(sapply(adequacy, is.null))]

## load adequacy test used in article
load(file = paste0(PATH, "adequacy.Rdata"))
load(file = paste0(PATH, "adequacy_meta.Rdata"))

# get only adequate unbiased random walk time series 
# (adequate function from empirical_functions.R)
absolute <- adequate(adequacy_meta)

#--------------------------------------- #
# Calculate darwins for complete dataset #
#--------------------------------------- #

darwins_compl <- complete_meta

# calculate darwins (puts time and darwins on log scale)
darwins_compl <- darwins(darwins_compl)

# bind data and choose variables from metadata (nn is changed to mean population nn for regression)
# (bind function from empirical_functions.R)
bind_darwins_compl <- bind(darwins_compl, variance_term = "darwins", 
                           variables = c("popID","darwins", "interval_MY", "nn"))

# mixed effect linear regression
darwins_compl_lmer <- lmer(darwins ~ interval_MY + (1|popID), bind_darwins_compl, weights = 1/nn)

# summary statistics (written manually into plot)
summary(darwins_compl_lmer)
r.squaredGLMM(darwins_compl_lmer)
print(paste("Cor =", cor(bind_darwins_compl$interval_MY, bind_darwins_compl$darwins)))

# tidy regression results
darwins_compl_tidy <- tidy(darwins_compl_lmer)

# plot regression
darwins_compl_lmer_plot <- ggplot(bind_darwins_compl, aes(interval_MY, darwins)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = darwins_compl_tidy$estimate[1], slope = darwins_compl_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("A")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  annotate("text", x = -5, y = -3.5, parse = TRUE, label="italic(y)==-2.389-0.856~italic(x)", size = 5) +
  annotate("text", x = -5, y = -5, parse = TRUE, label="italic(SE)=='' %+-% '0.044'", size = 5) +
  annotate("text", x = -5, y = -6.5, parse = TRUE, label = "italic(R)^2== 0.795", size = 5) +
  annotate("text", x = -5, y = -8, parse = TRUE, label = "italic(n)== 634", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -2.8, ymin = -9, ymax = -2.5), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19))
print(darwins_compl_lmer_plot)


#------------------------------------------- #
# Calculate darwins for absolute fit dataset # 
#------------------------------------------- #

# calculate darwins
darwins_abs <- darwins(absolute)

# bind data and choose variables
bind_darwins_abs <- bind(darwins_abs, variance_term = "darwins",
                         variables = c("popID","darwins", "interval_MY", "nn"))

# mixed effect linear regression
darwins_abs_lmer <- lmer(darwins ~ interval_MY + (1|popID), bind_darwins_abs, weights = 1/nn)

# summary statistics (written manually into plot)
summary(darwins_abs_lmer)
r.squaredGLMM(darwins_abs_lmer)
print(paste("Cor =", cor(bind_darwins_abs$interval_MY, bind_darwins_abs$darwins)))

# tidy regression result
darwins_abs_tidy <- tidy(darwins_abs_lmer)

# plot
darwins_abs_lmer_plot <- ggplot(bind_darwins_abs, aes(interval_MY, darwins)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = darwins_abs_tidy$estimate[1], slope = darwins_abs_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("C")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  annotate("text", x = -5, y = -4, parse = TRUE, label="italic(y)==-2.439-0.831~italic(x)", size = 5) +
  annotate("text", x = -5, y = -5.5, parse = TRUE, label="italic(SE)=='' %+-% '0.046'", size = 5) +
  annotate("text", x = -5, y = -7, parse = TRUE, label = "italic(R)^2== 0.786", size = 5) +
  annotate("text", x = -5, y = -8.7, parse = TRUE, label = "italic(n)== 492", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -2.8, ymin = -9.8, ymax = -2.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19))
print(darwins_abs_lmer_plot)


# ----------------------------------------------------------------------- #
# Estimate vstep from complete dataset, not accounting for sampling error #
# ----------------------------------------------------------------------- #

# set variance in empirical data to near zero
complete_no_vv <- lapply(complete, function(x){
  len <- length(x$vv)
  x$vv <- rep(0.00000001, len)
  return(x)
})

# fit an unbiased random walk
URW_compl_no_error <- lapply(complete_no_vv, opt.joint.URW, pool = TRUE)

# append parameters from test to complete_meta
URW_compl_no_error <- mapply(c, complete_meta, URW_compl_no_error, SIMPLIFY = FALSE)

# bind data and choose variables
bind_URW_compl_no_error <- bind(data = URW_compl_no_error, variance_term = "vstep",
                                variables = c("popID","vstep", "interval_MY", "nn"))

# put tt and vstep on log scale
bind_URW_compl_no_error$interval_MY <- log(bind_URW_compl_no_error$interval_MY)
bind_URW_compl_no_error$vstep <- log(bind_URW_compl_no_error$vstep)

# mixed effect linear regression
URW_compl_no_error_lmer <- lmer(vstep ~ interval_MY + (1|popID),
                                bind_URW_compl_no_error, weights = 1/nn)

# summary statistics (written manually into plot)
summary(URW_compl_no_error_lmer)
r.squaredGLMM(URW_compl_no_error_lmer)
print(paste("Cor =", cor(bind_URW_compl_no_error$interval_MY, bind_URW_compl_no_error$vstep)))

# tidy regression results
URW_compl_no_error_tidy <- tidy(URW_compl_no_error_lmer)

# plot
URW_compl_no_error_lmer_plot <- ggplot(bind_URW_compl_no_error, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_compl_no_error_tidy$estimate[1], slope = URW_compl_no_error_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("D")))) +
  ylab(expression(paste("Log ", italic("v")[bold("step")]))) + xlab("Log time") +
  annotate("text", x = -5, y = -6, parse = TRUE, label="italic(y)==-2.323-0.707~italic(x)", size = 5) +
  annotate("text", x = -5, y = -8, parse = TRUE, label="italic(SE)=='' %+-% '0.062'", size = 5) +
  annotate("text", x = -5, y = -10, parse = TRUE, label = "italic(R)^2== 0.533", size = 5) +
  annotate("text", x = -5, y = -12.2, parse = TRUE, label = "italic(n)== 643", size = 5) +
  geom_rect(aes(xmin = -7.3, xmax = -2.7, ymin = -13.7, ymax = -4.5), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))
print(URW_compl_no_error_lmer_plot)


# --------------------------------------------------------------------------- #
# Estimate vstep from absolute fit dataset, not accounting for sampling error #
# --------------------------------------------------------------------------- #

# set variance to near zero in the absolute fit dataset
absolute_no_vv <- lapply(absolute, function(x){
  len <- length(x$vv)
  x$vv <- rep(0.00000001, len)
  return(x)
})

# make paleoTS objects for the model fit to work
absolute_no_vv_paleo <- lapply(absolute_no_vv, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

# fit an unbiased random walk 
URW_abs_no_error <- lapply(absolute_no_vv_paleo, opt.joint.URW)

# append metadata 
URW_abs_no_error <- mapply(c, URW_abs_no_error, absolute_no_vv, SIMPLIFY = FALSE)

# bind and choose variables 
bind_URW_abs_no_error <- bind(URW_abs_no_error, variance_term = "vstep",
                              variables = c("popID","vstep", "interval_MY", "nn"))

# log transform tt and vstep
bind_URW_abs_no_error$interval_MY <- log(bind_URW_abs_no_error$interval_MY)
bind_URW_abs_no_error$vstep <- log(bind_URW_abs_no_error$vstep)

# mixed effect linear regression
URW_abs_no_error_lmer <- lmer(vstep ~ interval_MY + (1|popID), bind_URW_abs_no_error, weights = 1/nn)

# summary statistics (written manually into plot)
summary(URW_abs_no_error_lmer)
r.squaredGLMM(URW_abs_no_error_lmer)
print(paste("Cor =", cor(bind_URW_abs_no_error$interval_MY, bind_URW_abs_no_error$vstep)))

# tidy regression results
URW_abs_no_error_tidy <- tidy(URW_abs_no_error_lmer, conf.int = TRUE)

# plot
URW_abs_no_error_lmer_plot <- ggplot(bind_URW_abs_no_error, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_abs_no_error_tidy$estimate[1], slope = URW_abs_no_error_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(bold("F"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = -5, y = -4.5, parse = TRUE, label="italic(y)==-2.570-0.658~italic(x)", size = 5) +
  annotate("text", x = -5, y = -6, parse = TRUE, label="italic(SE)=='' %+-% '0.064'", size = 5) +
  annotate("text", x = -5, y = -7.5, parse = TRUE, label = "italic(R)^2== 0.495", size = 5) +
  annotate("text", x = -5, y = -9.2, parse = TRUE, label = "italic(n)== 494", size = 5) +
  geom_rect(aes(xmin = -7.4, xmax = -2.7, ymin = -10.5, ymax = -3), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_abs_no_error_lmer_plot)


# ------------------------------------------------------------- #
# Estimate vstep from the complete dataset, with sampling error #
# ------------------------------------------------------------- #

# example of how to fit an unbiased random walk
URW_compl <- lapply(complete, opt.joint.URW, pool = TRUE)
## this will give some error messages with paleoTS v0.6.1,
## circumvent the errors with this approach:

URW_compl_fit <- list()
for(i in 1:length(complete)){
  try(URW_compl_fit[[i]] <- opt.joint.URW(complete[[i]], pool = TRUE))
}

# append metadata
URW_compl <- mapply(c, complete_meta, URW_compl_fit, SIMPLIFY = FALSE)

# remove time series that cannot be processed by the loglikelihood function
URW_compl = URW_compl[-which(sapply(URW_compl_fit, is.null))]

# bind data and choose variables
bind_URW_compl <- bind(data = URW_compl, variance_term = "vstep",
                       variables = c("popID","vstep", "interval_MY", "nn"))

# remove time series with estimated rates very close to 0
bind_URW_compl <- bind_URW_compl[!bind_URW_compl$vstep <= 1.000000e-06, ]

# put tt and vstep on log scale
bind_URW_compl$interval_MY <- log(bind_URW_compl$interval_MY)
bind_URW_compl$vstep <- log(bind_URW_compl$vstep)

# mixed effect linear regression
URW_compl_lmer <- lmer(vstep ~ interval_MY + (1|popID), bind_URW_compl, weights = 1/nn)

# summary statistict (written manually into plot)
summary(URW_compl_lmer)
r.squaredGLMM(URW_compl_lmer)
print(paste("Cor =", cor(bind_URW_compl$interval_MY, bind_URW_compl$vstep)))

# tidy regression results
URW_compl_tidy <- tidy(URW_compl_lmer)

# plot
URW_compl_lmer_plot <- ggplot(bind_URW_compl, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_compl_tidy$estimate[1], slope = URW_compl_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("G")))) +
  ylab(expression(paste("Log ", italic("v")["step"]))) + xlab("Log time") +
  annotate("text", x = -5, y = -6.2, parse = TRUE, label="italic(y)==-3.467-0.809~italic(x)", size = 5) +
  annotate("text", x = -5, y = -8.2, parse = TRUE, label="italic(SE)=='' %+-% '0.066'", size = 5) +
  annotate("text", x = -5, y = -10.2, parse = TRUE, label = "italic(R)^2== 0.611", size = 5) +
  annotate("text", x = -5, y = -12.4, parse = TRUE, label = "italic(n)== 537", size = 5) +
  geom_rect(aes(xmin = -7.3, xmax = -2.7, ymin = -13.9, ymax = -4.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_compl_lmer_plot)


# ----------------------------------------------------------------- #
# Estimate vstep from the absolute fit dataset, with sampling error #
# ----------------------------------------------------------------- #

# make paloTS objects for the model fit
absolute_paleo <- lapply(absolute, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

# fit unbiased random walk 
URW_abs <- lapply(absolute_paleo, opt.joint.URW)

# append metadata  
URW_abs <- mapply(c, URW_abs, absolute, SIMPLIFY = FALSE)

# bind and choose variables 
bind_URW_abs <- bind(URW_abs, variance_term = "vstep",
                     variables = c("popID","vstep", "interval_MY", "nn"))

# remove time series with estimated rates very close to 0
bind_URW_abs <- bind_URW_abs[!bind_URW_abs$vstep <= 1.000000e-06, ]

# log transform tt and vstep
bind_URW_abs$interval_MY <- log(bind_URW_abs$interval_MY)
bind_URW_abs$vstep <- log(bind_URW_abs$vstep)

# mixed effect linear regression
URW_abs_lmer <- lmer(vstep ~ interval_MY + (1|popID), bind_URW_abs, weights = 1/nn)

# summary statistics (written manually into plot)
summary(URW_abs_lmer)
r.squaredGLMM(URW_abs_lmer)
print(paste("Cor =", cor(bind_URW_abs$interval_MY, bind_URW_abs$vstep)))

# tidy regression results
URW_abs_tidy <- tidy(URW_abs_lmer, conf.int = TRUE)

# plot 
URW_abs_lmer_plot <- ggplot(bind_URW_abs, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_abs_tidy$estimate[1], slope = URW_abs_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(bold("I"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = -5, y = -6, parse = TRUE, label="italic(y)==-3.690-0.752~italic(x)", size = 5) +
  annotate("text", x = -5, y = -7.5, parse = TRUE, label="italic(SE)=='' %+-% '0.070'", size = 5) +
  annotate("text", x = -5, y = -9.2, parse = TRUE, label = "italic(R)^2== 0.591", size = 5) +
  annotate("text", x = -5, y = -11, parse = TRUE, label = "italic(n)== 415", size = 5) +
  geom_rect(aes(xmin = -7.4, xmax = -2.7, ymin = -12.5, ymax = -4.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_abs_lmer_plot)


# ------------------------- #
# Plot regressions together #
# ------------------------- #

# plot and write to file
pdf(width = 20.5, height = 11.5, file = "[PATH_TO_RESULTS_FOLDER]/empirical.pdf")
grid.arrange(darwins_compl_lmer_plot, darwins_abs_lmer_plot,
             URW_compl_no_error_lmer_plot, URW_abs_no_error_lmer_plot,
             URW_compl_lmer_plot, URW_abs_lmer_plot, nrow = 3)
dev.off()
