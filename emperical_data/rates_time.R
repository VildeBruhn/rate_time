#########################################
## Evolutionary rates and time scaling ##
#########################################

#paleoTS.v.0.5.3
#evoTS GitHub version

rm(list = ls())

library(parallel)
library(doParallel)
library(broom.mixed)
library(Matrix)
library(ggpmisc)
library(evoTS)
library(paleoTS)
library(adePEM)
library(tidyverse)
library(data.table)
library(wesanderson)
library(lme4)
library(MuMIn)

source("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/rates_time_functions.R")

# set working directory for database
setwd("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/timeseries/")

# -------------------------
# Set up for parallel runs
# -------------------------

n_cores <- parallel::detectCores() - 1

# create the cluster
my_cluster <- parallel::makeCluster(
  n_cores, 
  type = "FORK"
)

# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my_cluster)

#------------------------------------------
# IMPORT AND EDIT FILES FROM DATABASE 
#------------------------------------------

# import
timeseries <- read_delim("timeseries.txt", col_names = TRUE, delim = "\t")
metadata <- read_delim("metadata.txt", col_names = TRUE, delim = "\t")

# join dataframes
df <- left_join(timeseries, metadata, by = c("tsID"))

# remove modern timeseries
df <- subset(df, sediment!="none")
df <- subset(df, trait_type!="complex")

# remove time series with less than 5 steps
df <- subset(df, steps >= 5)

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))

# process data
ln_data_meta <- dt(df, "tsID")
ln_data <- lapply(ln_data_meta, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})

#-------------------------------------
# Calculate darwins for all datasets 
#-------------------------------------

# make dataframe with ln(mm) and interval length MY
df_darwins <- ln_data_meta

# calculate darwins (puts time and darwins on log scale)
df_darwins <- darwins(df_darwins)

# bind data and choose units (will but interval_MY and darwins on log scale)
binded_darwins <- bind(df_darwins, variance_term = "darwins", unit_list = c("popID","darwins", "trait_type", 
                                                                            "steps", "microfossil", "lat", "lon", "interval_MY", "nn"))
  
# linear regression
darwins_lmer <- lmer(darwins ~ interval_MY + (1|popID), binded_darwins, weights = 1/nn)

# summary
summary(darwins_lmer)
r.squaredGLMM(darwins_lmer)
print(paste("Cor =", cor(binded_darwins$interval_MY, binded_darwins$darwins)))

# save results
sink("../results/darwins_summaries.txt")
summary(darwins_lmer)
r.squaredGLMM(darwins_lmer)
print(paste("Cor =", cor(binded_darwins$interval_MY, binded_darwins$darwins)))
sink()

# tidy lmer result with lowest BIC
darwins_tidy <- tidy(darwins_lmer)

# plot and write lmer to pdf
pdf(file = "../results/darwins_lmer.pdf")
darwins_plot_lmer <- ggplot(binded_darwins, aes(interval_MY, darwins)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = darwins_tidy$estimate[1], slope = darwins_tidy$estimate[2], linewidth = 0.7) + 
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
print(darwins_plot_lmer)
dev.off() 


########################################
## RUN URW ON ALL DATASETS AND PLOT   ##
##        vstep VS. TIME,             ##
## not accunting for measurement error ##
########################################

# fit URW to the data
ln_data_no_vv <- lapply(ln_data, function(x){
  len <- length(x$vv)
  x$vv <- rep(0.00000001, len)
  return(x)
})

URW_no_error <- mclapply(ln_data_no_vv, URW, pool = TRUE)

# append parameters from test to ln_data_meta
URW_no_error <- mapply(c, ln_data_meta, URW_no_error, SIMPLIFY = FALSE)

# bind data and choose units
binded_URW_no_error <- bind(data = URW_no_error, variance_term = "vstep", unit_list = c("popID","vstep", "trait_type", 
                                                                              "steps", "microfossil", "lat", "lon", 
                                                                              "interval_MY", "nn"))

# put tt and vstep on log scale
binded_URW_no_error$interval_MY <- log(binded_URW_no_error$interval_MY)
binded_URW_no_error$vstep <- log(binded_URW_no_error$vstep)

# remove infinite values
binded_URW_no_error <- binded_URW_no_error %>% filter_all(all_vars(!is.infinite(.)))
binded_URW_no_error <- binded_URW_no_error %>% filter_all(all_vars(!is.na(vstep)))

# linear regression
URW_no_error_lmer <- lmer(vstep ~ interval_MY + (1|popID), binded_URW_no_error, weights = 1/nn)

# summary
summary(URW_no_error_lmer)
r.squaredGLMM(URW_no_error_lmer)
print(paste("Cor =", cor(binded_URW_no_error$interval_MY, binded_URW_no_error$vstep)))

# save results
sink("../results/URW_no_error_summaries.txt")
summary(URW_no_error_lmer)
r.squaredGLMM(URW_no_error_lmer)
print(paste("Cor =", cor(binded_URW_no_error$interval_MY, binded_URW_no_error$vstep)))
sink()

# tidy lmer results
URW_no_error_tidy <- tidy(URW_no_error_lmer)

# plot and write to file
pdf(file = "../results/log_URW_all_lmer.pdf")
URW_no_error_plot_lmer <- ggplot(binded_URW_no_error, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_no_error_tidy$estimate[1], slope = URW_no_error_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("D")))) +
  ylab(expression(paste("Log ", italic("v")[bold("step")]))) + xlab("Log time") +
  annotate("text", x = -5, y = -6, parse = TRUE, label="italic(y)==-2.385-0.681~italic(x)", size = 5) +
  annotate("text", x = -5, y = -8, parse = TRUE, label="italic(SE)=='' %+-% '0.061'", size = 5) +
  annotate("text", x = -5, y = -10, parse = TRUE, label = "italic(R)^2== 0.515", size = 5) +
  annotate("text", x = -5, y = -12.2, parse = TRUE, label = "italic(n)== 643", size = 5) +
  geom_rect(aes(xmin = -7.3, xmax = -2.7, ymin = -13.7, ymax = -4.5), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_no_error_plot_lmer)
dev.off() 



######################################
## RUN URW ON ALL DATASETS AND PLOT ##
##        vstep VS. TIME            ##
######################################

# fit URW to the data
URW_all <- mclapply(ln_data, URW, pool = TRUE)

# append parameters from test to ln_data_meta
URW_all <- mapply(c, ln_data_meta, URW_all, SIMPLIFY = FALSE)

# bind data and choose units
binded_URW_all <- bind(data = URW_all, variance_term = "vstep", unit_list = c("popID","vstep", "trait_type", 
                                                                       "steps", "microfossil", "lat", "lon", 
                                                                       "interval_MY", "nn"))

# put tt and vstep on log scale
binded_URW_all$interval_MY <- log(binded_URW_all$interval_MY)
binded_URW_all$vstep <- log(binded_URW_all$vstep)

# remove infinite values
binded_URW_all <- binded_URW_all %>% filter_all(all_vars(!is.infinite(.)))
binded_URW_all <- binded_URW_all %>% filter_all(all_vars(!is.na(vstep)))

# linear regression
URW_all_lmer <- lmer(vstep ~ interval_MY + (1|popID), binded_URW_all, weights = 1/nn)

# summary
summary(URW_all_lmer)
r.squaredGLMM(URW_all_lmer)
print(paste("Cor =", cor(binded_URW_all$interval_MY, binded_URW_all$vstep)))

# save results
sink("../results/URW_all_summaries.txt")
summary(URW_all_lmer)
r.squaredGLMM(URW_all_lmer)
print(paste("Cor =", cor(binded_URW_all$interval_MY, binded_URW_all$vstep)))
sink()

# tidy lmer results
URW_all_tidy <- tidy(URW_all_lmer)

# plot and write to file
pdf(file = "../results/log_URW_all_lmer.pdf")
URW_all_plot_lmer <- ggplot(binded_URW_all, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_all_tidy$estimate[1], slope = URW_all_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("G")))) +
  ylab(expression(paste("Log ", italic("v")["step"]))) + xlab("Log time") +
  annotate("text", x = -5, y = -6.2, parse = TRUE, label="italic(y)==-3.481-0.780~italic(x)", size = 5) +
  annotate("text", x = -5, y = -8.2, parse = TRUE, label="italic(SE)=='' %+-% '0.065'", size = 5) +
  annotate("text", x = -5, y = -10.2, parse = TRUE, label = "italic(R)^2== 0.584", size = 5) +
  annotate("text", x = -5, y = -12.4, parse = TRUE, label = "italic(n)== 529", size = 5) +
  geom_rect(aes(xmin = -7.3, xmax = -2.7, ymin = -13.9, ymax = -4.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_all_plot_lmer)
dev.off() 

#################################################################### 
## CHECK IF URW IS THE BEST MODEL FOR THE DATA SETS,              ##
## RUN darwins and URW ON DATA SETS THAT HAD LOWEST BICc FOR URW, ## 
## EXTRACT RATE (vstep) AND TIME INTERVAL AND PLOT                ##
####################################################################

#-----------------------------------------
# Calculate darwins for URW best datasets 
#-----------------------------------------

model_test <- mclapply(ln_data, fit.all.univariate)
save(file = "../model_test_lim5.Rdata", model_test)

load("../model_test_lim5.Rdata")

# add metadata
model_test_meta <- mapply(c, model_test, ln_data_meta, SIMPLIFY = FALSE)

# get only timeseries that fit URW best
URW_best <- get_URW_best(model_test_meta)

# calculate darwins (puts time and darwins on log scale)
df_darwins2 <- darwins(URW_best)

# bind data and choose units (will but interval_MY and darwins on log scale)
binded_darwins2 <- bind(df_darwins2, variance_term = "darwins", unit_list = c("popID","darwins", "trait_type", 
                                                                            "steps", "microfossil", "lat", "lon", "interval_MY", "nn"))
# linear regression
darwins2_lmer <- lmer(darwins ~ interval_MY + (1|popID), binded_darwins2, weights = 1/nn)

# summary
summary(darwins2_lmer)
r.squaredGLMM(darwins2_lmer)
print(paste("Cor =", cor(binded_darwins2$interval_MY, binded_darwins2$darwins)))

# save results
sink("../results/darwins_URWbest_summaries.txt")
summary(darwins2_lmer)
r.squaredGLMM(darwins2_lmer)
print(paste("Cor =", cor(binded_darwins2$interval_MY, binded_darwins2$darwins)))
sink()

# tidy lmer result with lowest BIC
darwins2_tidy <- tidy(darwins2_lmer)

# plot and write lmer to pdf
pdf(file = "../results/darwins_URWbest_lmer.pdf")
darwins2_plot_lmer <- ggplot(binded_darwins2, aes(interval_MY, darwins)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = darwins2_tidy$estimate[1], slope = darwins2_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("B")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  annotate("text", x = -5, y = -3, parse = TRUE, label="italic(y)==-2.063-0.872~italic(x)", size = 5) +
  annotate("text", x = -5, y = -4.5, parse = TRUE, label="italic(SE)=='' %+-% '0.065'", size = 5) +
  annotate("text", x = -5, y = -6, parse = TRUE, label = "italic(R)^2== 0.736", size = 5) +
  annotate("text", x = -5, y = -7.7, parse = TRUE, label = "italic(n)== 164", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -2.8, ymin = -8.7, ymax = -1.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19))
print(darwins2_plot_lmer)
dev.off() 

### URW without measurement error ###

# test all possible univariate models from evoTS on timeseries
model_test_meta <- model_test_meta

# get only timeseries that fit URW best
URW_best <- get_URW_best(model_test_meta)

# remove vv
URW_best_no_error <- lapply(URW_best, function(x){
  len <- length(x$vv)
  x$vv <- rep(0.00000001, len)
  return(x)
})

# need to make new paleoTS objects for the model test to work
data_URW_paleo_no_error <- lapply(URW_best_no_error, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

# fit URW to the data sets that had URW as the best model
URW_best_no_error_model <- mclapply(data_URW_paleo_no_error, URW)

# append parameters from model test to URW_best  
URW_best_no_error <- mapply(c, URW_best_no_error, URW_best_no_error_model, SIMPLIFY = FALSE)

# bind and choose units
binded_URW_best_no_error <- bind(URW_best_no_error, "vstep", unit_list = c("popID","vstep", "trait_type", 
                                                         "steps", "microfossil", "lat", "lon", 
                                                         "interval_MY", "nn"))

# put tt and vstep on log scale
binded_URW_best_no_error$interval_MY <- log(binded_URW_best_no_error$interval_MY)
binded_URW_best_no_error$vstep <- log(binded_URW_best_no_error$vstep)

# remove infinite values
binded_URW_best_no_error <- binded_URW_best_no_error %>% filter_all(all_vars(!is.infinite(.)))

# linear regression
URW_best_no_error_lmer <- lmer(vstep ~ interval_MY + (1|popID), binded_URW_best_no_error, weights = 1/nn)

# summary
summary(URW_best_no_error_lmer)
r.squaredGLMM(URW_best_no_error_lmer)
print(paste("Cor =", cor(binded_URW_best_no_error$interval_MY, binded_URW_best_no_error$vstep)))

# save results
sink("../results/URW_best_no_error_summaries.txt")
summary(URW_best_no_error_lmer)
r.squaredGLMM(URW_best_no_error_lmer)
print(paste("Cor =", cor(binded_URW_best_no_error$interval_MY, binded_URW_best_no_error$vstep)))
sink()

# tidy lmer results
URW_best_no_error_tidy <- tidy(URW_best_no_error_lmer)

# plot and write lm to pdf 
pdf(file = "../results/log_URW_best_no_error_lmer.pdf")
URW_best_no_error_plot_lmer <- ggplot(binded_URW_best_no_error, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_best_no_error_tidy$estimate[1], slope = URW_best_no_error_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(bold("E"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = -5, y = -4.5, parse = TRUE, label="italic(y)==-2.640-0.815~italic(x)", size = 5) +
  annotate("text", x = -5, y = -6, parse = TRUE, label="italic(SE)=='' %+-% '0.091'", size = 5) +
  annotate("text", x = -5, y = -7.5, parse = TRUE, label = "italic(R)^2== 0.603", size = 5) +
  annotate("text", x = -5, y = -9.2, parse = TRUE, label = "italic(n)== 165", size = 5) +
  geom_rect(aes(xmin = -7.3, xmax = -2.7, ymin = -10.7, ymax = -3), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_best_no_error_plot_lmer)
dev.off() 

### URW with measurement error ###

# test all possible univariate models from evoTS on timeseries
model_test_meta <- model_test_meta

# get only timeseries that fit URW best
URW_best <- get_URW_best(model_test_meta)

# need to make new paleoTS objects for the model test to work
data_URW_paleo <- lapply(URW_best, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

# fit URW to the data sets that had URW as the best model
URW_best_model <- mclapply(data_URW_paleo, URW)

# append parameters from model test to URW_best  
URW_best <- mapply(c, URW_best, URW_best_model, SIMPLIFY = FALSE)

# bind and choose units
binded_URW_best <- bind(URW_best, "vstep", unit_list = c("popID","vstep", "trait_type", 
                                                                           "steps", "microfossil", "lat", "lon", 
                                                                           "interval_MY", "nn"))

# put tt and vstep on log scale
binded_URW_best$interval_MY <- log(binded_URW_best$interval_MY)
binded_URW_best$vstep <- log(binded_URW_best$vstep)

# remove infinite values
binded_URW_best <- binded_URW_best %>% filter_all(all_vars(!is.infinite(.)))

# linear regression
URW_best_lmer <- lmer(vstep ~ interval_MY + (1|popID), binded_URW_best, weights = 1/nn)

# summary
summary(URW_best_lmer)
r.squaredGLMM(URW_best_lmer)
print(paste("Cor =", cor(binded_URW_best$interval_MY, binded_URW_best$vstep)))

# save results
sink("../results/URW_best_summaries.txt")
summary(URW_best_lmer)
r.squaredGLMM(URW_best_lmer)
print(paste("Cor =", cor(binded_URW_best$interval_MY, binded_URW_best$vstep)))
sink()

# tidy lmer results
URW_best_tidy <- tidy(URW_best_lmer)

# plot and write lm to pdf 
pdf(file = "../results/log_URW_best_lmer.pdf")
URW_best_plot_lmer <- ggplot(binded_URW_best, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_best_tidy$estimate[1], slope = URW_best_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(bold("H"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = -5, y = -5, parse = TRUE, label="italic(y)==-3.461-0.855~italic(x)", size = 5) +
  annotate("text", x = -5, y = -6.5, parse = TRUE, label="italic(SE)=='' %+-% '0.087'", size = 5) +
  annotate("text", x = -5, y = -8, parse = TRUE, label = "italic(R)^2== 0.630", size = 5) +
  annotate("text", x = -5, y = -9.7, parse = TRUE, label = "italic(n)== 165", size = 5) +
  geom_rect(aes(xmin = -7.3, xmax = -2.7, ymin = -10.8, ymax = -3.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(URW_best_plot_lmer)
dev.off() 


############################################
## TEST ADEQUACY of URW model WITH adePEM ##
############################################

#--------------------------------------------
# Calculate darwins for URW adequate datasets 
#--------------------------------------------

# run adequacy test on datasets that fitted URW best
adeq_URW <- mclapply(data_URW_paleo, fit3adequacy.RW, plot = FALSE)

# append results from adequacy test to URW_best 
URW_adeq <- mapply(c, URW_best, adeq_URW, SIMPLIFY = FALSE)

# get only adequate timeseries
URW_adeq_passed <- adequate(URW_adeq) 

 # calculate darwins (puts time and darwins on log scale)
df_darwins3 <- darwins(URW_adeq_passed)

# bind data and choose units (will but interval_MY and darwins on log scale)
binded_darwins3 <- bind(df_darwins3, variance_term = "darwins", unit_list = c("popID","darwins", "trait_type", 
                                                                              "steps", "microfossil", "lat", "lon", "interval_MY", "nn"))
# linear regression
darwins3_lmer <- lmer(darwins ~ interval_MY + (1|popID), binded_darwins3, weights = 1/nn)

# summary
summary(darwins3_lmer)
r.squaredGLMM(darwins3_lmer)
print(paste("Cor =", cor(binded_darwins3$interval_MY, binded_darwins3$darwins)))

# save results
sink("../results/darwins_URWadeq_summaries.txt")
summary(darwins3_lmer)
r.squaredGLMM(darwins3_lmer)
print(paste("Cor =", cor(binded_darwins3$interval_MY, binded_darwins3$darwins)))
sink()

# tidy lmer result with lowest BIC
darwins3_tidy <- tidy(darwins3_lmer)

# plot and write lmer to pdf
pdf(file = "../results/darwins_URWadeq_lmer.pdf")
darwins3_plot_lmer <- ggplot(binded_darwins3, aes(interval_MY, darwins)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = darwins3_tidy$estimate[1], slope = darwins3_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(paste(bold("C")))) +
  ylab(expression(paste("Log ", italic("darwins")))) + xlab(("Log time")) +
  annotate("text", x = -5, y = -3, parse = TRUE, label="italic(y)==-2.074-0.889~italic(x)", size = 5) +
  annotate("text", x = -5, y = -4.5, parse = TRUE, label="italic(SE)=='' %+-% '0.067'", size = 5) +
  annotate("text", x = -5, y = -6, parse = TRUE, label = "italic(R)^2== 0.735", size = 5) +
  annotate("text", x = -5, y = -7.7, parse = TRUE, label = "italic(n)== 148", size = 5) +
  geom_rect(aes(xmin = -7.2, xmax = -2.8, ymin = -8.8, ymax = -1.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15), ) +
  theme(title = element_text(size = 19))
print(darwins3_plot_lmer)
dev.off() 



### URW no measurement error ###

# run adequacy test on datasets that fitted URW best
adeq_URW_no_error <- mclapply(data_URW_paleo_no_error, fit3adequacy.RW, plot = FALSE)

# append results from adequacy test to URW_best 
URW_adeq_no_error <- mapply(c, URW_best_no_error, adeq_URW_no_error, SIMPLIFY = FALSE)

# get only adequate timeseries
URW_adeq_no_error_passed <- adequate(URW_adeq_no_error)

# append results from adequacy test to URW_best 
URW_adeq_no_error <- mapply(c, URW_best_no_error, adeq_URW_no_error, SIMPLIFY = FALSE)

# get only adequate timeseries
URW_adeq_no_error_passed <- adequate(URW_adeq_no_error)

# bind and choose unit 
binded_adeq_no_error <- bind(URW_adeq_no_error_passed, "vstep", c("popID","vstep", "trait_type", 
                                                "steps", "microfossil", "lat", "lon", 
                                                "interval_MY", "nn"))

# log-transform tt and vstep
binded_adeq_no_error$interval_MY <- log(binded_adeq_no_error$interval_MY)
binded_adeq_no_error$vstep <- log(binded_adeq_no_error$vstep)

# remove infinite values
binded_adeq_no_error <- binded_adeq_no_error %>% filter_all(all_vars(!is.infinite(.)))

# linear regression
URW_adeq_no_error_lmer <- lmer(vstep ~ interval_MY + (1|popID), binded_adeq_no_error, weights = 1/nn)

# summary
summary(URW_adeq_no_error_lmer)
r.squaredGLMM(URW_adeq_no_error_lmer)
print(paste("Cor =", cor(binded_adeq_no_error$interval_MY, binded_adeq_no_error$vstep)))

# save results
sink("../results/URW_adeq_no_error_summaries.txt")
summary(URW_adeq_no_error_lmer)
r.squaredGLMM(URW_adeq_no_error_lmer)
print(paste("Cor =", cor(binded_adeq_no_error$interval_MY, binded_adeq_no_error$vstep)))
sink()

# tidy lmer results
URW_adeq_no_error_tidy <- tidy(URW_adeq_no_error_lmer, conf.int = TRUE)


# plot and write lm to pdf
pdf(file = "../results/log_URW_adeq_no_error_lmer.pdf")
adeq_URW_no_error_lmer <- ggplot(binded_adeq_no_error, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_adeq_no_error_tidy$estimate[1], slope = URW_adeq_no_error_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(bold("F"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = -5, y = -4.5, parse = TRUE, label="italic(y)==-2.791-0.832~italic(x)", size = 5) +
  annotate("text", x = -5, y = -6, parse = TRUE, label="italic(SE)=='' %+-% '0.096'", size = 5) +
  annotate("text", x = -5, y = -7.5, parse = TRUE, label = "italic(R)^2== 0.614", size = 5) +
  annotate("text", x = -5, y = -9.2, parse = TRUE, label = "italic(n)== 131", size = 5) +
  geom_rect(aes(xmin = -7.4, xmax = -2.7, ymin = -10.5, ymax = -3), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(adeq_URW_no_error_lmer)
dev.off()


### Adequate URW with measurement error ###

# bind and choose unit 
binded_adeq <- bind(URW_adeq_passed, "vstep", c("popID","vstep", "trait_type", 
                                                "steps", "microfossil", "lat", "lon", 
                                                "interval_MY", "nn"))

# log-transform tt and vstep
binded_adeq$interval_MY <- log(binded_adeq$interval_MY)
binded_adeq$vstep <- log(binded_adeq$vstep)

# remove infinite values
binded_adeq <- binded_adeq %>% filter_all(all_vars(!is.infinite(.)))

# linear regression
URW_adeq_lmer <- lmer(vstep ~ interval_MY + (1|popID), binded_adeq, weights = 1/nn)

# summary
summary(URW_adeq_lmer)
r.squaredGLMM(URW_adeq_lmer)
print(paste("Cor =", cor(binded_adeq$interval_MY, binded_adeq$vstep)))
 
# save results
sink("../results/URW_adeq_summaries.txt")
summary(URW_adeq_lmer)
r.squaredGLMM(URW_adeq_lmer)
print(paste("Cor =", cor(binded_adeq$interval_MY, binded_adeq$vstep)))
sink()

# tidy lmer results
URW_adeq_tidy <- tidy(URW_adeq_lmer, conf.int = TRUE)


# plot and write lm to pdf
pdf(file = "../results/log_URW_adeq_lmer.pdf")
adeq_URW_lmer <- ggplot(binded_adeq, aes(interval_MY, vstep)) + 
  geom_point(color = c(wes_palette("Rushmore1")[3])) + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = URW_adeq_tidy$estimate[1], slope = URW_adeq_tidy$estimate[2], linewidth = 0.7) + 
  ggtitle(expression(bold("I"))) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  annotate("text", x = -5, y = -5, parse = TRUE, label="italic(y)==-3.514-0.855~italic(x)", size = 5) +
  annotate("text", x = -5, y = -6.5, parse = TRUE, label="italic(SE)=='' %+-% '0.091'", size = 5) +
  annotate("text", x = -5, y = -8.2, parse = TRUE, label = "italic(R)^2== 0.621", size = 5) +
  annotate("text", x = -5, y = -10, parse = TRUE, label = "italic(n)== 147", size = 5) +
  geom_rect(aes(xmin = -7.4, xmax = -2.7, ymin = -11.5, ymax = -3.7), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(title = element_text(size = 19))

print(adeq_URW_lmer)
dev.off() 


########################################
## PLOT ALL REGRESSION LINES TOGETHER ##
########################################
library(gridExtra)
pdf(width = 20.5, height = 11.5, file = "../results/darwins_URW_best_adeq_error.pdf")
grid.arrange(darwins_plot_lmer, darwins2_plot_lmer, darwins3_plot_lmer,
             URW_no_error_plot_lmer, URW_best_no_error_plot_lmer, adeq_URW_no_error_lmer,
             URW_all_plot_lmer, URW_best_plot_lmer, adeq_URW_lmer, nrow = 3)
dev.off()


### EXTRA ###

test_data <- lapply(URW_adeq_passed, function(x){
  x$mm <- exp(x$mm)
  x$mm = x$mm - x$mm[1]
  return(x)
})
test_data <- lapply(test_data, function(x){
  x$var_z <- c()
  x$var_z[1] <- 0
  for (i in 2:length(x$mm)){
    x$var_z[i] <- var(x$mm[1:i])
  }
  return(x)
})

#which((test$mm > 100) == TRUE)
test_data <- test_data[names(test_data) != 452]
test_data <- test_data[names(test_data) != 608]
test_data <- test_data[names(test_data) != 318]


test <- bind(test_data, "vstep", c("popID","vstep", "trait_type", 
                                   "steps", "microfossil", "lat", "lon", 
                                   "interval_MY", "nn", "mm", "tt", "var_z"))


pdf(file = "../results/trait_mean_time.pdf")
ggplot() +
  geom_line(data = test, aes(x = tt, y = mm, group = data_frame), color = c(wes_palette("Rushmore1")[3])) +
  theme_classic() + ylab("TRAIT MEAN") + xlab("TIME") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13))
dev.off()

pdf(file = "../results/var_trait_time.pdf")
ggplot() +
  geom_line(data = test, aes(x = tt, y = var_z, group = data_frame), color = c(wes_palette("Rushmore1")[3])) +
  theme_classic() + ylab("VARIANCE OF TRAIT MEAN") + xlab("TIME") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13))
dev.off()



