library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)
library(tidyverse)
library(gridExtra)

rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

df_PSCL4 <- read_csv("/home/gabrielsalcedo/Documentos/TYLCVD-fit/model3/datos_circle.csv")

onset <- df_PSCL4$t   
cases <- df_PSCL4$infected

N = length(onset) # Number of days observed throughout the outbreak
pop = 1000         # Plant Population 
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

parameters = c("y_hat", "y_init", "theta",  "R_0")  #deterministic models (Model 1, Model 2)

m1 <- stan_model("/home/gabrielsalcedo/Documentos/TYLCVD-fit/model3/Stan_tomato.stan")

n_chains=4
n_warmups=500
n_iter=5000
n_thin=50
set.seed(1234)
# Set initial values:
ini_1 = function(){
 list(theta=c(runif(1,0,1), runif(1,0.008,0.1), runif(1,0.008,0.1),
      runif(1,0.05,0.1), runif(1,0,1), runif(1,0.03,1.0),
      runif(1,0,0.5), runif(1,0,1), runif(1,0,1)),
      S0=runif(1,(pop-3)/pop,(pop-1)/pop))
}
# ini_1 = function(){
#     list(theta=c(0.1, 0.01, 0.01, 0.075,0.01,0.06,0.06,0.4,1.0),
#          S0=runif(1,(pop-3)/pop,(pop-1)/pop),L0 = 1- S0)  
# }


time.start_nuts1 <- Sys.time()
nuts_fit_1 = sampling(m1, data = tomato_data, pars = parameters, init = ini_1, 
                      chains = n_chains, warmup = n_warmups, iter = n_iter, 
                      control=list(adapt_delta=0.9),
                      thin = n_thin, seed=13219)
time.end_nuts1 <- Sys.time()
duration_nuts1<- time.end_nuts1 - time.start_nuts1
nuts_fit_1_summary <- summary(nuts_fit_1, pars = c("lp__", "theta[3]", 
                                                   "theta[4]", 
                                                   "y_init[1]", "R_0"))$summary
print(nuts_fit_1_summary,scientific=FALSE,digits=2)
posts_1 <-  rstan::extract(nuts_fit_1)

#Check HMC diagnostics:

mod1_diagnostics <-rstan::get_sampler_params(nuts_fit_1)

# Check for divergent transitions
rstan::check_divergences(nuts_fit_1)

posterior_1 <- as.array(nuts_fit_1)
color_scheme_set("viridis")
# Markov chain traceplots
mcmc_trace(posterior_1, pars="lp__")
mcmc_trace(posterior_1, pars=c("theta[3]", "theta[4]", "y_init[1]"))
mcmc_trace(posterior_1, pars="R_0")

# Univariate and bivariate marginal posterior distributions
pairs(nuts_fit_1, pars = c("theta[3]", "theta[4]", "y_init[1]"), 
      labels = c("r_2", "b", "Sp(0)"), 
      cex.labels=1.5, font.labels=9, condition = "accept_stat__")  

# Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior_1, pars=c("theta[3]", "theta[4]", "y_init[1]"))

#Central posterior uncertainty intervals
mcmc_intervals(posterior_1,pars = c("theta[3]", "theta[4]", "y_init[1]"))


##Plot model fit, median and 95% credible interval:

# Model fitted values across the observed time period
fit_I_1 <- posts_1$y_hat[,,6]    # Fitted fraction of infected 
fit_SIR_1 <- fit_I_1*pop         # Fitted number of infected
median_I_1 = apply(fit_SIR_1, 2, median)
low_I_1 = apply(fit_SIR_1, 2, quantile, probs=c(0.025))
high_I_1 = apply(fit_SIR_1, 2, quantile, probs=c(0.975))
df_sample_N = data.frame(cases, sample_time)
df_fit_I_1 = data.frame(median_I_1, low_I_1, high_I_1, sample_time)

save(df_sample_N,file="data.Rda")
save(df_fit_I_1,file="df_I_det_Poiss.Rda")

ggplot(df_sample_N, aes(x=sample_time, y=cases)) +
  geom_ribbon(aes(x=sample_time, ymin = low_I_1, ymax = high_I_1), fill = "orange", alpha = 0.6) +
  geom_line(data = df_fit_I_1, aes(x=sample_time, y=median_I_1, color = "Median"), size = 1.3) +
  geom_point(shape = 19, size = 3, (aes(color="Data"))) +
  scale_colour_manual(name='', values=c('Data'='black', 'Median'='darkorange3'))+
  guides(colour = guide_legend(override.aes = list(shape=c(16,NA),  linetype=c(0,1))))+
  labs(x = "Time (days)", y = "Number of Infected students") + 
  scale_x_continuous(limits=c(0, 6), breaks=c(0,2,4,6)) +
  scale_y_continuous(limits=c(0,1000), breaks=c(0,200,400,600,800,1000)) +
  theme_bw()+ theme(text = element_text(size=20))





















