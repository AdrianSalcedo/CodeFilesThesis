library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)
library(tidyverse)
library(gridExtra)

rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

df_PSCL4 <- read_csv("datos_circle.csv")

cases <- df_PSCL4$infected
days <- df_PSCL4$t

theme_set(theme_bw())
ggplot(data = df_PSCL4) + 
  geom_point(mapping = aes(x = days, y = cases)) + 
  labs(y = "Number of infected plants")
###############################################
N = length(days) # Number of days observed throughout the outbreak
pop = 1000         # Population 
sample_time=1:N

# Modify data into a form suitable for Stan
tomato_data = list(n_obs = N,
                n_theta = 9,
                n_difeq = 6,
                n_pop = pop,
                y = cases,
                t0 = 0,
                ts = sample_time)

# Specify parameters to monitor

parameters = c("betap", "b", "R0", "y_hat")

m1 <- stan_model("Tomato_sys2.stan")

### Fit and sample from the posterior using Hamiltonian Monte Carlo-NUTS:

n_chains=4
n_warmups=500
n_iter=2000
n_thin=50
set.seed(1234)
# Set initial values:
ini_1 = function(){
  list(theta=c(0.1, 0.01, 0.01, 0.075, 0.01, 0.06, 0.06, 0.4, 1.0))  
}

time.start_nuts1 <- Sys.time()
nuts_fit_1 = sampling(m1, 
                      data = tomato_data, 
                      pars = parameters, 
                      init = ini_1, 
                      chains = n_chains, 
                      warmup = n_warmups,
                      iter = n_iter, 
                      thin=n_thin, 
                      seed=0)
time.end_nuts1 <- Sys.time()
duration_nuts1<- time.end_nuts1 - time.start_nuts1

nuts_fit_1_summary <- summary(nuts_fit_1, pars = c("lp__", "betap", "b", "R0"))$summary
print(nuts_fit_1_summary,scientific=FALSE,digits=2)
posts_1 <-  rstan::extract(nuts_fit_1)

### Check HMC diagnostics:

mod1_diagnostics <-rstan::get_sampler_params(nuts_fit_1)

# Check for divergent transitions
rstan::check_divergences(nuts_fit_1)

posterior_1 <- as.array(nuts_fit_1)
color_scheme_set("viridis")
# Markov chain traceplots
mcmc_trace(posterior_1, pars="lp__")
mcmc_trace(posterior_1, pars="b")
mcmc_trace(posterior_1, pars="R_0")

# Univariate and bivariate marginal posterior distributions
pairs(nuts_fit_1, pars = "b", labels = "b", 
      cex.labels=1.5, font.labels=9, condition = "accept_stat__")  

# Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior_1, pars="b")

#Central posterior uncertainty intervals
mcmc_intervals(posterior_1,pars = "b")














