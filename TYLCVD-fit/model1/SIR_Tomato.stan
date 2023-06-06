functions {
  real[] SIR(real t,  // time
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {

    real dy_dt[6];
    real Sp = y[1];
    real Lp = y[2];
    real Ip = y[3];
    real Sv = y[4];
    real Iv = y[5];
    real N_v;
//
    real beta_p = theta[1];
    real r_1 = theta[2];
    real r_2 = theta[3];
    real b = theta[4];
    real beta_v = theta[5];
    real gamma = theta[6];
    real gamma_f = theta[7];
    real theta_1 = theta[8];
    real mu = theta[9];
//
    N_v = mu /(gamma + gamma_f);

    dy_dt[1] = -beta_p * Sp * Iv + r_1 * Lp + r_2 * Ip;
    dy_dt[2] = beta_p * Sp * Iv - (b + r_1) * Lp;
    dy_dt[3] = b * Lp - r_2 * Ip;
    dy_dt[4] = -beta_v * Sv * Ip - (gamma + gamma_f) * Sv + (1 - theta_1) * mu;
    dy_dt[5] = beta_v * Sv * Ip - (gamma + gamma_f) * Iv + theta_1 * mu;
    dy_dt[6] = b * Lp;
    return dy_dt;
  }
  
}
data {
    int<lower = 1> n_obs;       // number of days observed
    int<lower = 1> n_theta;     // number of model parameters
    int<lower = 1> n_difeq;     // number of differential equations
    int<lower = 1> n_pop;       // population
    int y[n_obs];               // data, total number of infected
    real t0;                    // initial time point (zero)
    real ts[n_obs];             // time points observed
}

transformed data {
  real x_r[0];
  int x_i[0];
  // real<lower = 0.0001, upper = 0.1> beta_p;
  // real<lower = 0, upper = 1> r_1;
  // real<lower = 0, upper = 1> r_2;
  // real<lower = 0, upper = 1> b;
  // real<lower = 0, upper = 1> beta_v;
  // real<lower = 0, upper = 1> gamma;
  // real<lower = 0, upper = 1> gamma_f;
  // real<lower = 0, upper = 1> theta_1;
  // real<lower = 0, upper = 1> mu;
  // 
   // beta_p = 0.1;
   // r_1 = 0.05; 
   // r_2 = 0.05;
   // b = 0.05;
   // beta_v = 0.01;
   // gamma = 0.005;
   // gamma_f = 0.001;
   // theta_1 = 0.4;
   // mu =1.0;
  
}
parameters {
    //real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real<lower = 0, upper = 1> beta_p;
    real<lower = 0.0083, upper = 0.013> r_1;
    real<lower = 0.0083, upper = 0.013> r_2;
    real<lower = 0, upper = 1> b;
    real<lower = 0.0001, upper = 1> beta_v;
    real<lower = 0.03, upper = 1> gamma;
    real<lower = 0.0001, upper = 0.5> gamma_f;
    real<lower = 0, upper = 1> theta_1;
    real<lower = 0, upper = 0.75> mu;
    //real <lower = 0, upper =1000> Sp0;
    //real <lower = 0, upper= 0.0001> Lp0;
    //real <lower = 0, upper= 0.0001> Ip0;
    //real <lower = 0, upper= 0.0001> Sv0;
    //real <lower = 0, upper= 0.0001> Iv0;
    //real <lower = 0> N_v;
    real<lower=0> phi_inv;
}

transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];

    //initial conditions
    real Lp0 = 0.04;
    real Ip0 = 0.04;
    real Sv0 = 0.2;
    real Iv0 = 0.1;
    real theta[n_theta];
    
    real phi = 1. / phi_inv;
    
    theta[1] = beta_p;
    theta[2] = r_1;
    theta[3] = r_2;
    theta[4] = b;
    theta[5] = beta_v;
    theta[6] = gamma;
    theta[7] = gamma_f;
    theta[8] = theta_1;
    theta[9] = mu;
   
    y_init[1] = 1.0 - (Lp0 + Ip0);
    y_init[2] = Lp0;
    y_init[3] = Ip0;
    y_init[4] = Sv0;
    y_init[5] = Iv0;
    y_init[6] = 0.04;
    
    y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
    
}

model {
    real lambda[n_obs];      //poisson parameter
    //priors
    beta_p ~ normal(0.07,0.005);
    r_1 ~ normal(0.01,0.01);
    r_2 ~ normal(0.01,0.001);
    b ~ inv_gamma(16,1.2);
    beta_v ~ normal(0.1,0.01);
    gamma ~ inv_gamma(12,1.0);
    gamma_f ~ inv_gamma(24,3);
    theta_1 ~ normal(0.0028, 0.015);
    mu ~ normal(0.39, 0.005);
    phi_inv ~ exponential(5);
    //    likelihood
     //print("betap = ", beta_p);
    for (i in 1:n_obs){
       lambda[i] = (y_hat[i, 3]) * n_pop;
    }
    y ~ poisson_log(lambda);//neg_binomial_2(col(to_matrix(y_hat),3),phi);//poisson(lambda);

}

generated quantities {
    real R_0;      // Basic reproduction number
    R_0 = beta_p * beta_v * b /((b + r_1) * (gamma + gamma_f)* r_2);
   
}

