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

df_data_real <- read_csv("datos_squares.csv")
df_data <- read_csv("interpolated_Ramsh_data2.csv")

t_real <- df_data_real$t
cases_real <- df_data_real$y

cases_real <- 1000*cases_real
theme_set(theme_bw())
ggplot(data = df_data) + 
  geom_point(mapping = aes(x = t_real, y = cases_real ), color ='red') + 
  labs(y = "Number of infected plants")

rstan_options (auto_write = TRUE)
options (mc.cores = parallel::detectCores ())



# time series of cases
cases <- round(df_data$Interpolate_Ip, digits = 0)
dates <- df_data$time
# total count
N <- 1000;

# times
n_days <- length(cases) 
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions
lp0 <- 0.003
ip0 <- 0.001
sp0 <- 1 - lp0 - ip0 
sv0 <- 0.4
iv0 <- 0.3

y0 = c(Sp = sp0, Lp = lp0, Ip = ip0, Sv = sv0, Iv = iv0, Ipc = 0.001)

# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases)
model <- stan_model("SIR_tomato.stan")

fit_sir_negbin <- sampling(model,
                           data = data_sir,
                           iter = 2500,
                           chains = 4,
                           control = list(adapt_delta = 0.9),
                           seed = 0)


check_hmc_diagnostics(fit_sir_negbin)

posts <-  rstan::extract(fit_sir_negbin)
mod_diagnostics  <- rstan::get_sampler_params(fit_sir_negbin)
color_scheme_set("viridisE")
fit_CIS <- posts$y[, ,6]*1000
fit_SIR <- fit_CIS
median_I = apply(fit_SIR, 2, median)
low_I = apply(fit_SIR, 2, quantile, probs = c(0.025))
high_I = apply(fit_SIR, 2, quantile, probs = c(0.975))
df_sample_N = data.frame(cases, dates)
df_sample_N_real = data.frame(cases_real, t_real)
df_fit_CIS = data.frame(median_I, low_I, high_I, dates)

plt <- ggplot(df_sample_N,aes(x = dates, y = cases))+ 
  geom_ribbon(aes(x = dates, ymin = low_I, ymax = high_I), fill = "orange", 
              alpha = 0.6) +  
  geom_line(data = df_fit_CIS, aes(x = dates, y = median_I, 
                                   color = "Median"), size = 1.3) + 
  geom_point(data = df_sample_N_real, 
             aes(x = t_real, y = cases_real), color = "red") +
  geom_point(shape = 1, size = 2, (aes(color = "Data"))) + 
  scale_colour_manual(name = '', values = c('Data' = 'black', 
                                            'Median' = 'darkorange3')) + 
  guides(colour = guide_legend(override.aes = list(shape = c(16, NA), 
                                                   linetype = c(0, 1)))) +
  labs(x = "Time (days)", y = "Cumulative Infected Cases") #+ xlim(0, 33)

show(plt)
nuts_fit_summary <- summary(fit_sir_negbin, pars = c('beta_p', 'beta_v'))$summary

print(nuts_fit_summary,
      scientific = FALSE,
      digits = 4) 