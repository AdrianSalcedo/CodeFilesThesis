functions {
  real[] SIR(real t,  // time
  real[] y,           // system state {susceptible,infected,recovered}
  real[] theta,       // parameters
  real[] x_r,
  int[] x_i) {
  
  real dy_dt[6];
      real betap = theta[1];
      real r_1 = theta[2];
      real r_2 = theta[3];
      real b = theta[4];
      real betav = theta[5];
      real gamma = theta[6];
      real gammaf = theta[7];
      real theta1 = theta[8];
      real mu = theta[9];
      real Nv = mu/(gamma + gammaf);
      
      //real init[5] = {1, 0, 0, 0, 0}; // initial values
      real Sp = y[1];
      real Lp = y[2];
      real Ip = y[3];
      real Sv = y[4];
      real Iv = y[5];
      real Ip_aux = y[6];
 
      real dSp_dt = -betap * Sp * Iv + r_1 * Lp + r_2 * Ip;
      real dLp_dt =  betap * Sp * Iv - (b + r_1) * Lp;
      real dIp_dt = b * Lp - r_2 * Ip;
      real dSv_dt =  -betav * Sv * Ip  - (gamma + gammaf) * Sv + (1 - theta1) * mu / Nv;
      real dIv_dt =  betav * Sv * Ip - (gamma + gammaf) * Iv + theta1 * mu / Nv;
      real dIp_aux_dt = b * Lp;
      dy_dt = {dSp_dt, dLp_dt, dIp_dt, dSv_dt, dIv_dt, dIp_aux_dt};
  return dy_dt;
  }
  
  }
  data {
  int<lower = 1> n_obs;       // number of days observed
  int<lower = 1> n_theta;     // number of model parameters
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_pop;       // population 
  int y[n_obs];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_obs];         // time points observed
  }
  
  transformed data {
  real x_r[0];
  int x_i[0];
  //real<lower=0, upper = 1> betap;
  real<lower=0, upper = 1> r_1;
  real<lower=0, upper = 1> r_2;
  real<lower=0, upper = 1> betav;
  real<lower=0, upper = 1> gamma;
  real<lower=0, upper = 1> theta1;
  real<lower=0, upper = 1> mu;
  real<lower=0, upper = 1> gammaf;
  real<lower=0, upper = 100000> Nv;
  //betap = 0.1;
  r_1 = 0.01;
  r_2 = 0.01;
  betav = 0.01;
  gamma = 0.06;
  gammaf = 0.06;
  theta1 = 0.4;
  mu = 1.0;
  Nv = mu/(gamma + gammaf);
  }
  
  parameters {
  real<lower = 0> theta[n_theta]; // model parameters {beta,gamma}
  real<lower=0, upper = 1> b;
  real<lower=0, upper = 1> betap;
  //real<lower = 0, upper = 1> S0;  // initial fraction of susceptible individuals
  }
  
  transformed parameters{
  real y_hat[n_obs, n_difeq]; // solution from the ODE solver
  real y_init[n_difeq];     // initial conditions for both fractions of S and I
  
  y_init[1] = 1;
  y_init[2] = 0;
  y_init[3] = 0;
  y_init[4] = 0;
  y_init[5] = 0;
  y_init[6] = 0;
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  
  }
 model {
  real lambda[n_obs];      //poisson parameter
  
  //priors
  betap ~ uniform(0,1);
  b ~ uniform(0,1);
  //likelihood
  //for (i in 1:n_obs){
  //lambda[i] = n_pop * y_hat[i,3];
  //}
  lambda = y_hat[,6];
  y ~ poisson(lambda);
  }
  
generated quantities {
  real R0 = (betap * betav * b) / ((b + r_1) * (gamma + gammaf) * r_2) ;
  //real replanting_time = 1 / r_2;
  //real incubation_time = 1 / b;
  //real pred_cases[n_days-1];//-1
  //pred_cases = poisson_rng(incidence, phi);
}

