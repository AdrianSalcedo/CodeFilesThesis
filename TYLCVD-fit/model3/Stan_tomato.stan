functions {
  real[] SIR(real t,  // time
  real[] y,           // system state {susceptible,infected,recovered}
  real[] theta,       // parameters
  real[] x_r,
  int[] x_i) {
  
  real dy_dt[6];
  
  dy_dt[1] = - theta[1] * y[1] * y[5] + theta[2] * y[2] + theta[3] * y[3] ;
  dy_dt[2] = theta[1] * y[1] * y[5] - theta[2] * y[2] - theta[4] * y[2];
  dy_dt[3] = theta[4] * y[2] - theta[3] * y[3];
  dy_dt[4] = - theta[5] * y[4] * y[3] - (theta[6] + theta[7]) * y[4] + (1 - theta[8]) * theta[9] / (theta[9] / (theta[6] + theta[7]));
  dy_dt[5] = theta[5] * y[4] * y[3] - (theta[6] + theta[7]) * y[4] + theta[8] * theta[9]/ (theta[9] / (theta[6] + theta[7]));
  dy_dt[6] = theta[4] * y[2];
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
  }
  
  parameters {
    // model parameters {betap, r_1, r_2, b, betav, gamma, gammaf, theta, mu}
    real<lower = 0,upper =1> theta[n_theta];
  // real<lower = 0,upper =1> theta[1];
  // real<lower = 0.008,upper = 0.1> theta[2];
  // real<lower = 0.008,upper = 0.1> theta[3];
  // real<lower = 0.05,upper =0.1> theta[4];
  // real<lower = 0,upper =1> theta[5];
  // real<lower = 0.03,upper =1> theta[6];
  // real<lower = 0,upper = 0.5> theta[7];
  // real<lower = 0,upper =1> theta[8];
  // real<lower = 0,upper =1> theta[9];
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
  theta[1] ~ lognormal(0,1); //betap
  theta[2] ~ gamma(0.008,0.1); //r_1
  theta[3] ~ gamma(0.008,0.1); //r_2
  theta[4] ~ gamma(0.05,0.1); //b
  theta[5] ~ lognormal(0,1);//betav
  theta[6] ~ uniform(0.03,1.0);//gamma
  theta[7] ~ uniform(0,0.5);//gammaf
  theta[8] ~ uniform(0,1);//theta
  theta[9] ~ uniform(0,1);//mu
  //S0 ~ uniform(0, 1);
  
  //likelihood
  for (i in 1:n_obs){
  lambda[i] = y_hat[i,6];
  }
  y ~ poisson(lambda);
  //y ~ exponential_lcdf(lambda);
}
  
  generated quantities {
  real R_0;      // Basic reproduction number
  R_0 = (theta[1] * theta[4] * theta[5])/ ((theta[2] + theta[4]) * (theta[6] + theta[7]) * theta[2]);
  }
  
  