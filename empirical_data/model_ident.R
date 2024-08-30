# --------------------- #
# MODEL IDENTIFIABILITY #
# --------------------- #


#install.packages("evoTS")        version 0.1.3
#install.packages("paleoTS")      version 0.6.1
#install.packages("adePEM")
#install.packages("tidyverse")
#install.packages("Matrix")
#install.packages("lme4")
#install.packages("broom.mixed")
#install.packages("wesanderson")  colors for plots 

library(evoTS)
library(paleoTS)
library(adePEM)
library(tidyverse)
library(Matrix)
library(lme4)
library(broom.mixed)
library(wesanderson)


######################################
## REMEMBER TO CHANGE PATH TO FILES ## 
######################################

PATH = "[PATH_TO_DATA_FOLDER]"

# import functions
source(paste0(PATH, "/model_ident_functions.R"))


# --------------- #
# Import datasets #  
# --------------- #

## The data is processed in the script empirical.R

# import complete, relative fit and absolute fit datasets
load(paste0(PATH, "complete.Rdata"))
load(paste0(PATH, "relative.Rdata"))
load(paste0(PATH, "absolute.Rdata"))

# import dataframes with estimated vstep and some metadata
load(paste0(PATH, "bind_URW_compl.Rdata"))
load(paste0(PATH, "bind_URW_rel.Rdata"))
load(paste0(PATH, "bind_URW_abs.Rdata"))


# ----------------------- #
# Explore goodness of fit #
# ----------------------- #

## Plots from the articles are saved in the supplementary_material_1_(SM1) folder.
# They are stored as zip-files and are called likelihood_compl_plots.zip,
## likelihood_rel_plots.zip and likelihood_abs_plots.zip.

# examples of how to explore likelihood space and find lower and upper vsteps
# (function from model_ident_functions.R)
explore_compl <- search_likelihood(bind_URW_compl, complete, bind_URW_compl, save_plot = "./GitHub/rate_time/empirical_data/test1")
explore_rel <- search_likelihood(bind_URW_rel, relative, bind_URW_rel, save_plot = "./GitHub/rate_time/empirical_data/test2")
explore_abs <- search_likelihood(bind_URW_abs, absolute, bind_URW_abs, save_plot = "./GitHub/rate_time/empirical_data/test3")

# load data with lower and upper vsteps from article
load("[PATH_TO_DATA]/l_u_compl.Rdata")
load("[PATH_TO_DATA]/l_u_rel.Rdata")
load("[PATH_TO_DATA]/l_u_abs.Rdata")


# ----------------------------------------------------------- #
# Select random vstep and run regression for complete dataset #
# ----------------------------------------------------------- #

# set up dataframe
df_compl <- as.data.frame(matrix(nrow = length(l_u_compl$ts), ncol = 7))
colnames(df_compl) <- c("ts", "vstep", "tt", "nn", "popID", "estimate_1", "estimate_2")

# add data
df_compl$ts <- l_u_compl$ts
df_compl$tt <- bind_URW_compl$interval_MY
df_compl$nn <- bind_URW_compl$nn
df_compl$popID <- bind_URW_compl$popID

# define list for function
list_compl <- vector(mode = "list", length = 1000)

# find new vsteps and estimate regressions (function from model_ident_functions.R)
list_compl <- vstep_reg(list_data = list_compl, df = df_compl, lw_up_data = l_u_compl)

# bind data into correct form for plotting (function from model_ident_functions.R)
compl <- bind(list_compl, variance_term = "vstep", variables = c("estimate_1", "estimate_2", "SE"))
compl <- compl[,-5]
compl <- na.omit(compl)

# log transform complete data
bind_URW_compl_log <- bind_URW_compl
bind_URW_compl_log$vstep <- log(bind_URW_compl_log$vstep)
bind_URW_compl_log$interval_MY <- log(bind_URW_compl_log$interval_MY)

# plot and write to file

## The estimates will differ slightly from the article because vstep is chosen 
## randomly between the lower and upper value.

pdf(file = "[PATH_TO_RESULTS]/compl_ident.pdf")
compl_plot <- ggplot(bind_URW_compl_log, aes(interval_MY, vstep)) +
  geom_point(color = "white")  + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = compl$estimate_1, slope = compl$estimate_2, linewidth = 0.1, color = c(wes_palette("Rushmore1")[3]), alpha = 0.1) +
  geom_abline(intercept = -3.481, slope = -0.780, linewidth = 0.9, color = c(wes_palette("IsleofDogs1")[4])) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  ggtitle(expression(bold("Complete"))) +
  annotate("text", x = -4, y = -8, parse = TRUE, label="italic(y)[orig]==-3.481-0.780~italic(x)", size = 5) +
  annotate("text", x = -4, y = -9.5, parse = TRUE, label="italic(SE)[orig]=='' %+-% '0.065'", size = 5) +
  annotate("text", x = -4, y = -11, parse = TRUE, label= paste0("mean~italic(y)[new]==", round(mean(compl$estimate_1),3), 
                                                                round(mean(compl$estimate_2), 3), "~italic(x)"), size = 5) +
  annotate("text", x = -4, y = -12.61, parse = TRUE, label= paste0("mean~italic(SE)[new]=='' %+-%", round(mean(compl$SE), 3)), size = 5) +
  geom_rect(aes(xmin = -7.5, xmax = -0.5, ymin = -13.7, ymax = -6.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13)) +
  theme(title = element_text(size = 17))

compl_plot
dev.off()

# --------------------------------------------------------------- #
# Select random vstep and run regression for relative fit dataset #
# --------------------------------------------------------------- #

# set up dataframe
df_rel <- as.data.frame(matrix(nrow = length(l_u_rel$ts), ncol = 7))
colnames(df_rel) <- c("ts", "vstep", "tt", "nn", "popID", "estimate_1", "estimate_2")

# add data
df_rel$ts <- l_u_rel$ts
df_rel$tt <- bind_URW_rel$interval_MY
df_rel$nn <- bind_URW_rel$nn
df_rel$popID <- bind_URW_rel$popID

# define list for function
list_rel <- vector(mode = "list", length = 1000)

# find new vsteps and estimate regressions (function from model_ident_functions.R)
list_rel <- vstep_reg(list_data = list_rel, df = df_rel, lw_up_data = l_u_rel)

# bind data into correct form for plotting (function from model_ident_functions.R)
rel <- bind(list_rel, variance_term = "vstep", variables = c("estimate_1", "estimate_2", "SE"))
rel <- rel[,-5]
rel <- na.omit(rel)

# log transform relative fit data
bind_URW_rel_log <- bind_URW_rel
bind_URW_rel_log$vstep <- log(bind_URW_rel_log$vstep)
bind_URW_rel_log$interval_MY <- log(bind_URW_rel_log$interval_MY)

# plot and write to file

## The estimates will differ slightly from the article because vstep is chosen 
## randomly between the lower and upper value.

pdf(file = "[PATH_TO_RESULTS]/rel_ident.pdf")
rel_plot <- ggplot(bind_URW_rel_log, aes(interval_MY, vstep)) +
  geom_point(color = "white")  + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = rel$estimate_1, slope = rel$estimate_2, linewidth = 0.1, color = c(wes_palette("Rushmore1")[3]), alpha = 0.1) +
  geom_abline(intercept = -3.481, slope = -0.780, linewidth = 0.9, color = c(wes_palette("IsleofDogs1")[4])) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  ggtitle(expression(bold("Relative fit"))) +
  annotate("text", x = -4, y = -8, parse = TRUE, label="italic(y)[orig]==-3.481-0.780~italic(x)", size = 5) +
  annotate("text", x = -4, y = -9.5, parse = TRUE, label="italic(SE)[orig]=='' %+-% '0.065'", size = 5) +
  annotate("text", x = -4, y = -11, parse = TRUE, label= paste0("mean~italic(y)[new]==", round(mean(rel$estimate_1),3), 
                                                                round(mean(rel$estimate_2), 3), "~italic(x)"), size = 5) +
  annotate("text", x = -4, y = -12.61, parse = TRUE, label= paste0("mean~italic(SE)[new]=='' %+-%", round(mean(rel$SE), 3)), size = 5) +
  geom_rect(aes(xmin = -7.5, xmax = -0.5, ymin = -13.7, ymax = -6.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13)) +
  theme(title = element_text(size = 17))

rel_plot
dev.off()


# --------------------------------------------------------------- #
# Select random vstep and run regression for absolute fit dataset #
# --------------------------------------------------------------- #

# set up dataframe
df_abs <- as.data.frame(matrix(nrow = length(l_u_abs$ts), ncol = 7))
colnames(df_abs) <- c("ts", "vstep", "tt", "nn", "popID", "estimate_1", "estimate_2")

# add data
df_abs$ts <- l_u_abs$ts
df_abs$tt <- bind_URW_abs$interval_MY
df_abs$nn <- bind_URW_abs$nn
df_abs$popID <- bind_URW_abs$popID

# define list for function
list_abs <- vector(mode = "list", length = 1000)

# find new vsteps and estimate regressions (function from model_ident_functions.R)
list_abs <- vstep_reg(list_data = list_abs, df = df_abs, lw_up_data = l_u_abs)

# bind data into correct form for plotting (function from model_ident_functions.R)
abs <- bind(list_abs, variance_term = "vstep", variables = c("estimate_1", "estimate_2", "SE"))
abs <- abs[,-5]
abs <- na.omit(abs)

# log transform absolute fit data
bind_URW_abs_log <- bind_URW_abs
bind_URW_abs_log$vstep <- log(bind_URW_abs_log$vstep)
bind_URW_abs_log$interval_MY <- log(bind_URW_abs_log$interval_MY)

# plot and write to file

## The estimates will differ slightly from the article because vstep is chosen 
## randomly between the lower and upper value.

pdf(file = "[PATH_TO_RESULTS]/abs_ident.pdf")
abs_plot <- ggplot(bind_URW_abs_log, aes(interval_MY, vstep)) +
  geom_point(color = "white")  + theme_classic() + theme(legend.position="none") +
  geom_abline(intercept = abs$estimate_1, slope = abs$estimate_2, linewidth = 0.1, color = c(wes_palette("Rushmore1")[3]), alpha = 0.1) +
  geom_abline(intercept = -3.481, slope = -0.780, linewidth = 0.9, color = c(wes_palette("IsleofDogs1")[4])) +
  ylab(expression(paste("Log ", italic(v)["step"]))) + xlab(expression("Log time")) +
  ggtitle(expression(bold("Absolute fit"))) +
  annotate("text", x = -4, y = -8, parse = TRUE, label="italic(y)[orig]==-3.481-0.780~italic(x)", size = 5) +
  annotate("text", x = -4, y = -9.5, parse = TRUE, label="italic(SE)[orig]=='' %+-% '0.065'", size = 5) +
  annotate("text", x = -4, y = -11, parse = TRUE, label= paste0("mean~italic(y)[new]==", round(mean(abs$estimate_1),3), 
                                                                round(mean(abs$estimate_2), 3), "~italic(x)"), size = 5) +
  annotate("text", x = -4, y = -12.61, parse = TRUE, label= paste0("mean~italic(SE)[new]=='' %+-%", round(mean(abs$SE), 3)), size = 5) +
  geom_rect(aes(xmin = -7.5, xmax = -0.5, ymin = -13.7, ymax = -6.8), 
            fill = "white", alpha = 0, color = "black") +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text = element_text(size = 13)) +
  theme(title = element_text(size = 17))

abs_plot
dev.off()
