library(deSolve)
library(dplyr)
library(rstan)
library(gridExtra)
library(outbreaks)
library(bayesplot)
library(data.table)
library(knitr)
library(kableExtra)
library(tidyverse)

source("load_data.R")
source("init_SIR.R")
source("mcmc_stain_summary.R")
source("model_fit.R")
source("divergence_plots.R")
source("mcmc_post_analysis.R")

####################

df_PSCL4 <- read_csv("/home/gabrielsalcedo/Documentos/TYLCVD-fit/model3/data/interpolated_Ramsh_data2.csv")
df_PSCL4_real <- read_csv("/home/gabrielsalcedo/Documentos/TYLCVD-fit/model3/data/datos_squares.csv")

onset_real <- df_PSCL4_real$t
cum_cases_real <- df_PSCL4_real$infected

onset <- df_PSCL4$time
cum_cases <- as.integer(df_PSCL4$Interpolate_Ip)
cum_cases <- unlist(cum_cases, use.names = FALSE)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# data_tomato <- load_data()
# onset <- data_tomato[[1]]
# cum_cases <- data_tomato[[3]]
#cum_cases <- unlist(cum_cases, use.names = FALSE)
#
N = length(onset) # Number of days observed throughout the outbreak
pop = 1000         # Plant Population 
sample_time=1:N

#### SIR rstan model ####
model <- stan_model("SIR_Tomato3.stan")
#

#### mcmc parameters ####
n_chains <- 4
n_warmups <- 500
n_iter <- 25500
n_thin <- 50
set.seed(972198)

Tomato_data <- list(n_obs = N,
                     n_theta = 9,
                     n_difeq = 9,
                     n_pop = pop,
                     y = cum_cases,
                     t0 = 0,
                     ts = sample_time)
#parameters = c("y_init", "theta", "R_0")
parameters <- c("theta", "y_hat", "y_init", "R_0")

time.start_nuts <- Sys.time()
nuts_fit <-
  sampling(model,
           data = Tomato_data,
           pars = parameters,
           init = init_sir,
           chains = n_chains,
           warmup = n_warmups,
           iter = n_iter,
           thin = n_thin#,
           #control = list(adapt_delta = 0.86)
           )

#parameters = c("y_init", "theta", "R_0")
parameters = c("theta[1]","y_hat")#,"y_init", "R_0")
nuts_fit_summary <- summary(nuts_fit, pars = parameters)$summary
print(nuts_fit_summary,
      scientific = FALSE,
      digits = 4) 
time.end_nuts <- Sys.time()
duration_nuts <- time.end_nuts - time.start_nuts
parameters = c("theta[1]", "y_hat", "y_init", "R_0")
#parameters = c("y_hat", "y_init", "R_0")
nuts_fit_summary <- summary(nuts_fit, pars = parameters)$summary
print(nuts_fit_summary,
      scientific = FALSE,
      digits = 4)

sub_path_1 <- "/home/gabrielsalcedo"
sub_path_2 <- "Documentos/TYLCVD-fit"
sub_path_3 <- "model3/runs"

prefix_time <- Sys.time()
prefix_time <- paste(as.Date(prefix_time),
                     hour(prefix_time),
                     minute(prefix_time), sep="_")
file_name <- paste("run-", prefix_time,".RData", sep="")
runs_path <- 
  paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/") 
save.image(file=runs_path)
#
#### Post analysis ####
#
mcmcm_post_analysis(nuts_sample = nuts_fit)
#
### Model Fit ####
# Model fitted values across the observed time period
#
model_fit()
#### Divergence analysis ####
divergence_analysis()
###
##
#
