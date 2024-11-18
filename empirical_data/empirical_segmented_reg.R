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
#install.packages("segmented")

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
library(segmented)


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
# Make relative fit dataset #
# ------------------------- #

# load model test used in article
load(paste0(PATH, "model_test.Rdata"))
load(paste0(PATH, "model_test_meta.Rdata"))

# get only time series that fit an unbiased random walk best according to AICc
# (relative_fit function from empirical_functions.R)
relative <- relative_fit(model_test_meta)


# ------------------------- #
# Make absolute fit dataset #
# ------------------------- #

# make paleoTS objects for the adequacy test
relative_paleo <- lapply(relative, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

# run adequacy test for an unbiased random walk on relative fit dataset
# (the warnings are ok)
# The output of the test can vary with ca. +/- one time series (i.e., you will
# not always get the exact same number of time series that pass the adequacy test)
# set.seed to get the same output every time
set.seed(1)
adequacy <- lapply(relative_paleo, fit3adequacy.RW, plot = FALSE)

# append results from adequacy test to relative 
absolute <- mapply(c, relative, adequacy, SIMPLIFY = FALSE)

# get only adequate unbiased random walk time series 
# (adequate function from empirical_functions.R)
absolute <- adequate(absolute)


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

# segmented glm regression
darwins_compl_glm <- glm(darwins ~ interval_MY, family = gaussian, bind_darwins_compl, weights = 1/nn)
darwins_compl_seg <- segmented(darwins_compl_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(darwins_compl_seg)
slope(darwins_compl_seg)

# plot
data1 <- data.frame(x = bind_darwins_compl$interval_MY, y = bind_darwins_compl$darwins)
data2 <- data.frame(x = bind_darwins_compl$interval_MY, y = broken.line(darwins_compl_seg)$fit)
darwins_compl_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("A")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -5.5, parse = TRUE, label="beta[1]==-1.677 %+-% '0.734' ", size = 5) +
  annotate("text", x = -4.5, y = -7, parse = TRUE, label="beta[2]==-0.875 %+-% '0.044' ", size = 5) +
  geom_rect(aes(xmin = -7, xmax = -2, ymin = -7.8, ymax = -4.6), 
            fill = "white", alpha = 0, color = "black")
plot(darwins_compl_seg_plot)


#------------------------------------------- #
# Calculate darwins for relative fit dataset #  
#------------------------------------------- #

# calculate darwins for relative fit
darwins_rel <- darwins(relative)

# bind data and choose variables
bind_darwins_rel <- bind(darwins_rel, variance_term = "darwins", 
                         variables = c("popID","darwins", "interval_MY", "nn"))

# segmented glm regression
darwins_rel_glm <- glm(darwins ~ interval_MY, family = gaussian, bind_darwins_rel, weights = 1/nn)
darwins_rel_seg <- segmented(darwins_rel_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(darwins_rel_seg)
slope(darwins_rel_seg)

# plot
data1 <- data.frame(x = bind_darwins_rel$interval_MY, y = bind_darwins_rel$darwins)
data2 <- data.frame(x = bind_darwins_rel$interval_MY, y = broken.line(darwins_rel_seg)$fit)
darwins_rel_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("B")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -5.5, parse = TRUE, label="beta[1] ==-0.800 %+-% '0.068' ", size = 5) +
  annotate("text", x = -4.5, y = -7, parse = TRUE, label="beta[2] ==-1.369 %+-% '0.756' ", size = 5) +
  geom_rect(aes(xmin = -7, xmax = -2, ymin = -7.8, ymax = -4.6), 
            fill = "white", alpha = 0, color = "black")
plot(darwins_rel_seg_plot)


#------------------------------------------- #
# Calculate darwins for absolute fit dataset # 
#------------------------------------------- #

# calculate darwins
darwins_abs <- darwins(absolute)

# bind data and choose variables
bind_darwins_abs <- bind(darwins_abs, variance_term = "darwins",
                         variables = c("popID","darwins", "interval_MY", "nn"))

# segmented glm regression
darwins_abs_glm <- glm(darwins ~ interval_MY, family = gaussian, bind_darwins_abs, weights = 1/nn)
darwins_abs_seg <- segmented(darwins_abs_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(darwins_abs_seg)
slope(darwins_abs_seg)

# plot
data1 <- data.frame(x = bind_darwins_abs$interval_MY, y = bind_darwins_abs$darwins)
data2 <- data.frame(x = bind_darwins_abs$interval_MY, y = broken.line(darwins_abs_seg)$fit)
darwins_abs_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("C")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -5.5, parse = TRUE, label="beta[1] ==-0.816 %+-% '0.072' ", size = 5) +
  annotate("text", x = -4.5, y = -7, parse = TRUE, label="beta[2] ==-1.282 %+-% '1.245' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -7.8, ymax = -4.6), 
            fill = "white", alpha = 0, color = "black")
plot(darwins_abs_seg_plot)


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

# segmented glm regression
URW_compl_no_error_glm <- glm(vstep ~ interval_MY, family = gaussian, bind_URW_compl_no_error, weights = 1/nn)
URW_compl_no_error_seg <- segmented(URW_compl_no_error_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(URW_compl_no_error_seg)
slope(URW_compl_no_error_seg)

# plot
data1 <- data.frame(x = bind_URW_compl_no_error$interval_MY, y = bind_URW_compl_no_error$vstep)
data2 <- data.frame(x = bind_URW_compl_no_error$interval_MY, y = broken.line(URW_compl_no_error_seg)$fit)
URW_compl_no_error_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("D")))) +
  ylab(expression(paste("Log ", italic("vstep")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -5.5, parse = TRUE, label="beta[1] ==-0.373 %+-% '0.049' ", size = 5) +
  annotate("text", x = -4.5, y = -7, parse = TRUE, label="beta[2] ==-2.707 %+-% '0.424' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -7.8, ymax = -4.6), 
            fill = "white", alpha = 0, color = "black")
plot(URW_compl_no_error_seg_plot)


# --------------------------------------------------------------------------- #
# Estimate vstep from relative fit dataset, not accounting for sampling error #
# --------------------------------------------------------------------------- #

# set variance to near zero in the relative fit dataset
relative_no_vv <- lapply(relative, function(x){
  len <- length(x$vv)
  x$vv <- rep(0.00000001, len)
  return(x)
})

# make paleoTS objects for the model fit to work
relative_no_vv_paleo <- lapply(relative_no_vv, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

# fit an unbiased random walk 
URW_rel_no_error <- lapply(relative_no_vv_paleo, opt.joint.URW)

# append metadata 
URW_rel_no_error <- mapply(c, URW_rel_no_error, relative_no_vv, SIMPLIFY = FALSE)

# bind and choose variables
bind_URW_rel_no_error <- bind(URW_rel_no_error, "vstep",
                              variables = c("popID","vstep", "interval_MY", "nn"))

# put tt and vstep on log scale
bind_URW_rel_no_error$interval_MY <- log(bind_URW_rel_no_error$interval_MY)
bind_URW_rel_no_error$vstep <- log(bind_URW_rel_no_error$vstep)

# segmented glm regression
URW_rel_no_error_glm <- glm(vstep ~ interval_MY, family = gaussian, bind_URW_rel_no_error, weights = 1/nn)
URW_rel_no_error_seg <- segmented(URW_rel_no_error_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(URW_rel_no_error_seg)
slope(URW_rel_no_error_seg)

# plot
data1 <- data.frame(x = bind_URW_rel_no_error$interval_MY, y = bind_URW_rel_no_error$vstep)
data2 <- data.frame(x = bind_URW_rel_no_error$interval_MY, y = broken.line(URW_rel_no_error_seg)$fit)
URW_rel_no_error_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("E")))) +
  ylab(expression(paste("Log ", italic("vstep")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -5.5, parse = TRUE, label="beta[1] ==-0.630 %+-% '0.086' ", size = 5) +
  annotate("text", x = -4.5, y = -7, parse = TRUE, label="beta[2] ==-2.708 %+-% '0.863' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -7.8, ymax = -4.6), 
            fill = "white", alpha = 0, color = "black")
plot(URW_rel_no_error_seg_plot)


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

# segmented glm regression
URW_abs_no_error_glm <- glm(vstep ~ interval_MY, family = gaussian, bind_URW_abs_no_error, weights = 1/nn)
URW_abs_no_error_seg <- segmented(URW_abs_no_error_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(URW_abs_no_error_seg)
slope(URW_abs_no_error_seg)

# plot
data1 <- data.frame(x = bind_URW_abs_no_error$interval_MY, y = bind_URW_abs_no_error$vstep)
data2 <- data.frame(x = bind_URW_abs_no_error$interval_MY, y = broken.line(URW_abs_no_error_seg)$fit)
URW_abs_no_error_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("F")))) +
  ylab(expression(paste("Log ", italic("vstep")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -5.5, parse = TRUE, label="beta[1] ==-0.600 %+-% '0.093' ", size = 5) +
  annotate("text", x = -4.5, y = -7, parse = TRUE, label="beta[2] ==-2.820 %+-% '1.018' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -7.8, ymax = -4.6), 
            fill = "white", alpha = 0, color = "black")
plot(URW_abs_no_error_seg_plot)


# ------------------------------------------------------------- #
# Estimate vstep from the complete dataset, with sampling error #
# ------------------------------------------------------------- #

# example of how to fit an unbiased random walk
#URW_compl <- lapply(complete, opt.joint.URW, pool = TRUE)
## this will give some error messages with paleoTS v0.6.1,
## circumvent the errors with this approach:

URW_compl_fit <- list()
for(i in 1:length(complete)){
  print(i)
  try(URW_compl_fit[[i]] <- opt.joint.URW(complete[[i]], pool = TRUE))
}

# append metadata
URW_compl <- mapply(c, complete_meta, URW_compl_fit, SIMPLIFY = FALSE)

# remove time series that didn't work with paleoTS v0.6.1
URW_compl = URW_compl[-which(sapply(URW_compl_fit, is.null))]

# bind data and choose variables
bind_URW_compl <- bind(data = URW_compl, variance_term = "vstep",
                       variables = c("popID","vstep", "interval_MY", "nn"))

# remove time series with estimated rates very close to 0
bind_URW_compl <- bind_URW_compl[!bind_URW_compl$vstep <= 1.000000e-06, ]

# put tt and vstep on log scale (warning ok)
bind_URW_compl$interval_MY <- log(bind_URW_compl$interval_MY)
bind_URW_compl$vstep <- log(bind_URW_compl$vstep)

# segmented glm regression
URW_compl_glm <- glm(vstep ~ interval_MY, family = gaussian, bind_URW_compl, weights = 1/nn)
URW_compl_seg <- segmented(URW_compl_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(URW_compl_seg)
slope(URW_compl_seg)

# plot
data1 <- data.frame(x = bind_URW_compl$interval_MY, y = bind_URW_compl$vstep)
data2 <- data.frame(x = bind_URW_compl$interval_MY, y = broken.line(URW_compl_seg)$fit)
URW_compl_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("G")))) +
  ylab(expression(paste("Log ", italic("vstep")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -9.7, parse = TRUE, label="beta[1] ==-0.846 %+-% '0.098' ", size = 5) +
  annotate("text", x = -4.5, y = -11.5, parse = TRUE, label="beta[2] ==-0.365 %+-% '0.173' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -12.4, ymax = -8.7), 
            fill = "white", alpha = 0, color = "black")
plot(URW_compl_seg_plot)


# ----------------------------------------------------------------- #
# Estimate vstep from the relative fit dataset, with sampling error #
# ----------------------------------------------------------------- #

# fit unbiased random walk 
URW_rel <- lapply(relative_paleo, opt.joint.URW)

# append metadata  
URW_rel <- mapply(c, URW_rel, relative, SIMPLIFY = FALSE)

# bind and choose variables
bind_URW_rel <- bind(URW_rel, "vstep",
                     variables = c("popID","vstep", "interval_MY", "nn"))

# log transform tt and vstep
bind_URW_rel$interval_MY <- log(bind_URW_rel$interval_MY)
bind_URW_rel$vstep <- log(bind_URW_rel$vstep)

# segmented glm regression
URW_rel_glm <- glm(vstep ~ interval_MY, family = gaussian, bind_URW_rel, weights = 1/nn)
URW_rel_seg <- segmented(URW_rel_glm, seg.Z = ~interval_MY)

# summary stats (written manually into plot)
summary(URW_rel_seg)
slope(URW_rel_seg)

# plot
data1 <- data.frame(x = bind_URW_rel$interval_MY, y = bind_URW_rel$vstep)
data2 <- data.frame(x = bind_URW_rel$interval_MY, y = broken.line(URW_rel_seg)$fit)
URW_rel_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("H")))) +
  ylab(expression(paste("Log ", italic("vstep")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -7.5, parse = TRUE, label="beta[1] ==-0.767 %+-% '0.087' ", size = 5) +
  annotate("text", x = -4.5, y = -9, parse = TRUE, label="beta[2] ==-1.895 %+-% '0.652' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -9.9, ymax = -6.6), 
            fill = "white", alpha = 0, color = "black")
plot(URW_rel_seg_plot)


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

# log transform tt and vstep
bind_URW_abs$interval_MY <- log(bind_URW_abs$interval_MY)
bind_URW_abs$vstep <- log(bind_URW_abs$vstep)

# segmented glm regression
URW_abs_glm <- glm(vstep ~ interval_MY, family = gaussian, bind_URW_abs, weights = 1/nn)
URW_abs_seg <- segmented(URW_abs_glm, seg.Z = ~interval_MY)

# summary stats(written manually into plot)
summary(URW_abs_seg)
slope(URW_abs_seg)

# plot
data1 <- data.frame(x = bind_URW_abs$interval_MY, y = bind_URW_abs$vstep)
data2 <- data.frame(x = bind_URW_abs$interval_MY, y = broken.line(URW_abs_seg)$fit)
URW_abs_seg_plot <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) +
  geom_line(data = data2, linewidth = 0.7) +
  theme_classic() +
  ggtitle(expression(paste(bold("I")))) +
  ylab(expression(paste("Log ", italic("vstep")))) + xlab(("Log time")) +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19)) +
  annotate("text", x = -4.5, y = -7.5, parse = TRUE, label="beta[1] ==-0.762 %+-% '0.091' ", size = 5) +
  annotate("text", x = -4.5, y = -9, parse = TRUE, label="beta[2] ==-2.180 %+-% '1.578' ", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -1.8, ymin = -9.9, ymax = -6.6), 
            fill = "white", alpha = 0, color = "black")
plot(URW_abs_seg_plot)


# ------------------------- #
# Plot regressions together #
# ------------------------- #

# plot and write to file
pdf(width = 20.5, height = 11.5, file = "[PATH_TO_RESULTS_FOLDER]/empirical_segmented.pdf")
grid.arrange (darwins_compl_seg_plot, darwins_rel_seg_plot, darwins_abs_seg_plot,
             URW_compl_no_error_seg_plot, URW_rel_no_error_seg_plot, URW_abs_no_error_seg_plot,
             URW_compl_seg_plot, URW_rel_seg_plot, URW_abs_seg_plot, nrow = 3)
dev.off()
