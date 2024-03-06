#########################################
##    Explore likelihood landscape     ##
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
source("/Users/vildeki/Dropbox (UiO)/PhD/Evo. rates and time scaling/likelihood_functions.R")

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
ln_data_all <- lapply(ln_data_meta, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$N, tt = x$tt, oldest = "first")
})

#save(file = "../data_temp/ln_data_all.Rdata", ln_data_all)

######################################
##    RUN URW ON ALL DATASETS       ##
######################################

# fit URW to the data
URW_all <- mclapply(ln_data_all, URW, pool = TRUE)

# append parameters from test to ln_data_meta
URW_all <- mapply(c, ln_data_meta, URW_all, SIMPLIFY = FALSE)

# bind data and choose units
binded_URW_all <- bind(data = URW_all, variance_term = "vstep", unit_list = c("vstep", "interval_MY", "popID", "nn"))

# set time series with neg vstep to 0
binded_URW_all$vstep[444] <- 0
#save(file = "../data_temp/binded_URW_all.Rdata", binded_URW_all)

#-----------------------------------------
# URW AICc
#-----------------------------------------

load("../model_test.Rdata")

# add metadata
model_test_meta <- mapply(c, model_test, ln_data_meta, SIMPLIFY = FALSE)

# get only timeseries that fit URW best
URW_best <- get_URW_best(model_test_meta)

# need to make new paleoTS objects for the model test to work
ln_data_aicc <- lapply(URW_best, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

#save(file = "../data_temp/ln_data_aicc.Rdata", ln_data_aicc)

# fit URW to the data sets that had URW as the best model
URW_aicc <- mclapply(ln_data_aicc, URW)

# append parameters from test to ln_data_meta
URW_aicc <- mapply(c, URW_best, URW_aicc, SIMPLIFY = FALSE)

# bind and choose units
binded_URW_aicc <- bind(URW_aicc, "vstep", unit_list = c("vstep", "interval_MY", "popID", "nn"))

#save(file = "../data_temp/binded_URW_aicc.Rdata", binded_URW_aicc)


#-----------------------------------------
# URW adequate
#-----------------------------------------

# run adequacy test on datasets that fitted URW best
adeq_URW <- mclapply(ln_data_aicc, fit3adequacy.RW, plot = FALSE)

# append results from adequacy test to URW_best 
URW_adeq <- mapply(c, URW_best, adeq_URW, SIMPLIFY = FALSE)

# get only adequate timeseries
URW_adeq_passed <- adequate(URW_adeq)

# paleoTS object adequate
ln_data_adeq <- lapply(URW_adeq_passed, function(x) {
  as.paleoTS(mm = x$mm, vv = x$vv, nn = x$nn, tt = x$tt)
})

#save(file = "../data_temp/ln_data_adeq.Rdata", ln_data_adeq)

# fit URW to the data sets that had URW as the best model
URW_adeq_passed_fit <- mclapply(ln_data_adeq, URW)

# append parameters from test to ln_data_meta
URW_adeq_passed <- mapply(c, URW_adeq_passed_fit, URW_adeq_passed, SIMPLIFY = FALSE)

# bind and choose unit 
binded_URW_adeq <- bind(URW_adeq_passed, "vstep", c("vstep", "interval_MY", "popID", "nn"))

#save(file = "../data_temp/binded_URW_adeq.Rdata", binded_URW_adeq)


load("../data_temp/ln_data_all.Rdata")
load("../data_temp/binded_URW_all.Rdata")

load("../data_temp/ln_data_aicc.Rdata")
load("../data_temp/binded_URW_aicc.Rdata")

load("../data_temp/ln_data_adeq.Rdata")
load("../data_temp/binded_URW_adeq.Rdata")

######################################
##    Explore goodness of fit       ##
######################################

test_data <- binded_URW_aicc[c(149,164),]
test_ln <- ln_data_aicc[c(149,164)]

test <- search_likelihood(binded_URW_all, ln_data_all, binded_URW_all, save_plot = "../data_temp/likelihood_plots/")

# find upper and lower vstep
test <- search_likelihood(binded_URW_aicc, ln_data_aicc, binded_URW_aicc)

load("../data_temp/l_u_all.Rdata")
load("../data_temp/l_u_aicc.Rdata")
load("../data_temp/l_u_adeq.Rdata")

# check difference between data sets
#x <- setdiff(binded_URW_all, binded_URW_aicc)
#mean(x$interval_MY)
#mean(binded_URW_aicc$interval_MY)
#mean(binded_URW_adeq$interval_MY)
#mean(x$vstep)
#mean(binded_URW_aicc$vstep)
#mean(binded_URW_adeq$vstep)

####################################################
##    Select random vstep and run regression      ##
####################################################

# set up dataframe
df_new <- as.data.frame(matrix(nrow = length(l_u_adeq$ts), ncol = 7))
colnames(df_new) <- c("ts", "vstep", "tt", "nn", "popID", "estimate_1", "estimate_2")
df_new$ts <- l_u_adeq$ts
df_new$tt <- binded_URW_adeq$interval_MY
df_new$nn <- binded_URW_adeq$nn
df_new$popID <- binded_URW_adeq$popID

# get new vsteps and lmer
list_adeq <- vector(mode = "list", length = 1000)
for (i in 1:length(list_adeq)){
  list_adeq[[i]] <- df_new
  x <- i
  # generate random vstep between lower and upper
  for (i in 1:length(l_u_adeq$ts)){
    list_adeq[[x]]$vstep[i] <- runif(1, min=l_u_adeq$lower[i], max=l_u_adeq$upper[i])
  }
  # log-transform tt and vstep
  list_adeq[[x]]$tt <- log(list_adeq[[x]]$tt)
  list_adeq[[x]]$vstep <- log(list_adeq[[x]]$vstep)
  # remove infinite values
  list_adeq[[x]] <- list_adeq[[x]] %>% filter_all(all_vars(!is.infinite(.)))
  
  # linear regression
  lmer <- lmer(vstep ~ tt + (1|popID), list_adeq[[x]], weights = 1/nn)
  list_adeq[[x]]$estimate_1[1] <- tidy(lmer, conf.int = TRUE)$estimate[1]
  list_adeq[[x]]$estimate_2[1] <- tidy(lmer, conf.int = TRUE)$estimate[2]
  list_adeq[[x]]$SE[1] <- tidy(lmer, conf.int = TRUE)$std.error[2]
}

# Plot
plot_adeq <- bind(list_adeq, variance_term = "vstep", unit_list = c("estimate_1", "estimate_2", "SE"))
plot_adeq <- plot_adeq[,-5]
plot_adeq <- na.omit(plot_adeq)

binded_URW_adeq_log <- binded_URW_adeq
binded_URW_adeq_log$vstep <- log(binded_URW_adeq_log$vstep)
binded_URW_adeq_log$interval_MY <- log(binded_URW_adeq_log$interval_MY)

pdf(file = "../results/adeq_nonident.pdf")
adeq <- ggplot(binded_URW_adeq_log, aes(interval_MY, vstep)) +
  geom_point(color = "white") + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = plot_adeq$estimate_1, slope = plot_adeq$estimate_2, linewidth = 0.1, color = c(wes_palette("Rushmore1")[3]), alpha = 0.1) +
  geom_abline(intercept = -3.552, slope = -0.837, linewidth = 0.9, color = c(wes_palette("IsleofDogs1")[4])) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  ggtitle(expression(bold("Absolute fit"))) +
  annotate("text", x = -4, y = -8, parse = TRUE, label="italic(y)[orig]==-3.514-0.855~italic(x)", size = 5) +
  annotate("text", x = -4, y = -9.5, parse = TRUE, label="italic(SE)[orig]=='' %+-% '0.091'", size = 5) +
  annotate("text", x = -4, y = -11, parse = TRUE, label= paste0("mean~italic(y)[new]==", round(mean(plot_adeq$estimate_1),3), 
                                                                round(mean(plot_adeq$estimate_2), 3), "~italic(x)"), size = 5) +
  annotate("text", x = -4, y = -12.61, parse = TRUE, label= paste0("mean~italic(SE)[new]=='' %+-%", round(mean(plot_adeq$SE), 3)), size = 5) +
  geom_rect(aes(xmin = -7.5, xmax = -0.5, ymin = -13.7, ymax = -6.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13)) +
  theme(title = element_text(size = 17))

adeq
dev.off()


### AICc ###


# set up dataframe
df_new <- as.data.frame(matrix(nrow = length(l_u_aicc$ts), ncol = 7))
colnames(df_new) <- c("ts", "vstep", "tt", "nn", "popID", "estimate_1", "estimate_2")
df_new$ts <- l_u_aicc$ts
df_new$tt <- binded_URW_aicc$interval_MY
df_new$nn <- binded_URW_aicc$nn
df_new$popID <- binded_URW_aicc$popID

# get new vsteps and lmer
list_aicc <- vector(mode = "list", length = 1000)
for (i in 1:length(list_aicc)){
  list_aicc[[i]] <- df_new
  x <- i
  # generate random vstep between lower and upper
  for (i in 1:length(l_u_aicc$ts)){
    list_aicc[[x]]$vstep[i] <- runif(1, min=l_u_aicc$lower[i], max=l_u_aicc$upper[i])
  }
  # log-transform tt and vstep
  list_aicc[[x]]$tt <- log(list_aicc[[x]]$tt)
  list_aicc[[x]]$vstep <- log(list_aicc[[x]]$vstep)
  # remove infinite values
  list_aicc[[x]] <- list_aicc[[x]] %>% filter_all(all_vars(!is.infinite(.)))
  
  # linear regression
  lmer <- lmer(vstep ~ tt + (1|popID), list_aicc[[x]], weights = 1/nn)
  list_aicc[[x]]$estimate_1[1] <- tidy(lmer, conf.int = TRUE)$estimate[1]
  list_aicc[[x]]$estimate_2[1] <- tidy(lmer, conf.int = TRUE)$estimate[2]
  list_aicc[[x]]$SE[1] <- tidy(lmer, conf.int = TRUE)$std.error[2]
}


# Plot
plot_aicc <- bind(list_aicc, variance_term = "vstep", unit_list = c("estimate_1", "estimate_2", "SE"))
plot_aicc <- plot_aicc[,-5]
plot_aicc <- na.omit(plot_aicc)

binded_URW_aicc_log <- binded_URW_aicc
binded_URW_aicc_log$vstep <- log(binded_URW_aicc_log$vstep)
binded_URW_aicc_log$interval_MY <- log(binded_URW_aicc_log$interval_MY)

pdf(file = "../results/aicc_nonident.pdf")
aicc <- ggplot(binded_URW_aicc_log, aes(interval_MY, vstep)) +
  geom_point(color = "white")  + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = plot_aicc$estimate_1, slope = plot_aicc$estimate_2, linewidth = 0.1, color = c(wes_palette("Rushmore1")[3]), alpha = 0.1) +
  geom_abline(intercept = -3.461, slope = -0.855, linewidth = 0.9, color = c(wes_palette("IsleofDogs1")[4])) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  ggtitle(expression(bold("Relative fit"))) +
  annotate("text", x = -4, y = -8, parse = TRUE, label="italic(y)[orig]==-3.461-0.855~italic(x)", size = 5) +
  annotate("text", x = -4, y = -9.5, parse = TRUE, label="italic(SE)[orig]=='' %+-% '0.087'", size = 5) +
  annotate("text", x = -4, y = -11, parse = TRUE, label= paste0("mean~italic(y)[new]==", round(mean(plot_aicc$estimate_1),3), 
                                                                round(mean(plot_aicc$estimate_2), 3), "~italic(x)"), size = 5) +
  annotate("text", x = -4, y = -12.61, parse = TRUE, label= paste0("mean~italic(SE)[new]=='' %+-%", round(mean(plot_aicc$SE), 3)), size = 5) +
  geom_rect(aes(xmin = -7.5, xmax = -0.5, ymin = -13.7, ymax = -6.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13)) +
  theme(title = element_text(size = 17))

aicc
dev.off()



### All time series ###


# set up dataframe
df_new <- as.data.frame(matrix(nrow = length(l_u_all$ts), ncol = 7))
colnames(df_new) <- c("ts", "vstep", "tt", "nn", "popID", "estimate_1", "estimate_2")
df_new$ts <- l_u_all$ts
df_new$tt <- binded_URW_all$interval_MY
df_new$nn <- binded_URW_all$nn
df_new$popID <- binded_URW_all$popID

# get new vsteps and lmer
list_all <- vector(mode = "list", length = 1000)
for (i in 1:length(list_all)){
  list_all[[i]] <- df_new
  x <- i
  # generate random vstep between lower and upper
  for (i in 1:length(l_u_all$ts)){
    list_all[[x]]$vstep[i] <- runif(1, min=l_u_all$lower[i], max=l_u_all$upper[i])
  }
  # log-transform tt and vstep
  list_all[[x]]$tt <- log(list_all[[x]]$tt)
  list_all[[x]]$vstep <- log(list_all[[x]]$vstep)
  # remove infinite values
  list_all[[x]] <- list_all[[x]] %>% filter_all(all_vars(!is.infinite(.)))
  
  # linear regression
  lmer <- lmer(vstep ~ tt + (1|popID), list_all[[x]], weights = 1/nn)
  list_all[[x]]$estimate_1[1] <- tidy(lmer, conf.int = TRUE)$estimate[1]
  list_all[[x]]$estimate_2[1] <- tidy(lmer, conf.int = TRUE)$estimate[2]
  list_all[[x]]$SE[1] <- tidy(lmer, conf.int = TRUE)$std.error[2]
}


# Plot
plot_all <- bind(list_all, variance_term = "vstep", unit_list = c("estimate_1", "estimate_2", "SE"))
plot_all <- plot_all[,-5]
plot_all <- na.omit(plot_all)

binded_URW_all_log <- binded_URW_all
binded_URW_all_log$vstep <- log(binded_URW_all_log$vstep)
binded_URW_all_log$interval_MY <- log(binded_URW_all_log$interval_MY)

pdf(file = "../results/all_nonident.pdf")
all <- ggplot(binded_URW_all_log, aes(interval_MY, vstep)) +
  geom_point(color = "white")  + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = plot_all$estimate_1, slope = plot_all$estimate_2, linewidth = 0.1, color = c(wes_palette("Rushmore1")[3]), alpha = 0.1) +
  geom_abline(intercept = -3.481, slope = -0.780, linewidth = 0.9, color = c(wes_palette("IsleofDogs1")[4])) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  ggtitle(expression(bold("Complete"))) +
  annotate("text", x = -4, y = -8, parse = TRUE, label="italic(y)[orig]==-3.481-0.780~italic(x)", size = 5) +
  annotate("text", x = -4, y = -9.5, parse = TRUE, label="italic(SE)[orig]=='' %+-% '0.065'", size = 5) +
  annotate("text", x = -4, y = -11, parse = TRUE, label= paste0("mean~italic(y)[new]==", round(mean(plot_all$estimate_1),3), 
                                                                round(mean(plot_all$estimate_2), 3), "~italic(x)"), size = 5) +
  annotate("text", x = -4, y = -12.61, parse = TRUE, label= paste0("mean~italic(SE)[new]=='' %+-%", round(mean(plot_all$SE), 3)), size = 5) +
  geom_rect(aes(xmin = -7.5, xmax = -0.5, ymin = -13.7, ymax = -6.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13)) +
  theme(title = element_text(size = 17))

all
dev.off()


#------------------
# Find vstep range
#------------------

# find range aicc
aicc_range <- l_u_aicc
aicc_range$range <- rep(NA, length(aicc_range$ts))
aicc_range$range <- aicc_range$upper - aicc_range$lower

# find length of range above certain limit
limit <- 1
above_lim_aicc <- length(which(aicc_range$range > limit))

# find range adeq
adeq_range <- l_u_adeq
adeq_range$range <- rep(NA, length(adeq_range$ts))
adeq_range$range <- adeq_range$upper - adeq_range$lower

# find length of range above certain limit
limit <- 1
above_lim_adeq <- length(which(adeq_range$range > limit))

# find range all
all_range <- l_u_all
all_range$range <- rep(NA, length(all_range$ts))
all_range$range <- all_range$upper - all_range$lower

# find length of range above certain limit
limit <- 1
above_lim_all <- length(which(all_range$range > limit))

# make output
above_lim <- as.data.frame(matrix(nrow = 3, ncol = 4))
colnames(above_lim) <- c("data_set", "limit", "count", "total")
above_lim$limit <- c(1,1,1)
above_lim$data_set <- c("all", "aicc", "adequate")
above_lim$count <- c(above_lim_all, above_lim_aicc, above_lim_adeq)
above_lim$total <- c(length(all_range$ts), length(aicc_range$ts), length(adeq_range$ts))


# save tables
write.table(aicc_range, file = "../results/l_u_aicc.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(adeq_range, file = "../results/l_u_adeq.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(all_range, file = "../results/l_u_all.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)
write.table(above_lim, file = "../results/above_lim.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)


### Box plot vstep ###

ggplot(binded_URW_adeq, aes(x,vstep)) + geom_boxplot() + theme_classic() +
  ylim(0,0.5)
 