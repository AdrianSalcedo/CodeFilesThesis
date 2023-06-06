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
  //real<lower = 0.0083, upper = 0.013> r_1;
  real<lower = 0.0083, upper = 0.013> r_2;
  real<lower = 0.03, upper = 1> gamma;
  real<lower = 0.0001, upper = 0.5> gamma_f;
  
  //r_1 = 0.01; 
  r_2 = 0.01;
  gamma = 0.06;
  gamma_f = 0.06;
  
}
parameters {
    //real <lower = 0> theta[n_theta]; // model parameters {beta,gamma}
    real<lower = 0.0001, upper = 0.1> beta_p;
    real<lower = 0.0083, upper = 0.013> r_1;
    //real<lower = 0.0083, upper = 0.013> r_2;
    real<lower = 0.05, upper = 0.1> b;
    real<lower = 0.0001, upper = 0.2> beta_v;
    //real<lower = 0.03, upper = 1> gamma;
    //real<lower = 0.0001, upper = 0.5> gamma_f;
    real<lower = 0.0001, upper = 1> theta_1;
    real<lower = 0.0001, upper = 1> mu;

}
transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];

    // initial conditions
    real Lp0 = 150.0;
    real Ip0 = 10.0;
    real Sv0 = 6.0;
    real Iv0 = 1.0;
    real theta[n_theta];
    
    theta[1] = beta_p;
    theta[2] = r_1;
    theta[3] = r_2;
    theta[4] = b;
    theta[5] = beta_v;
    theta[6] = gamma;
    theta[7] = gamma_f;
    theta[8] = theta_1;
    theta[9] = mu;
   
    y_init[1] = n_pop - (Lp0 + Ip0);
    y_init[2] = Lp0;
    y_init[3] = Ip0;
    y_init[4] = Sv0;
    y_init[5] = Iv0;
    y_init[6] = Ip0;
    
    y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
}
model {
    real lambda[n_obs];      //poisson parameter
    //priors
    beta_p ~ lognormal(0.01,0.03);//(0.0001, 0.1); 0.01 0.1
    r_1 ~ uniform(0.0083, 0.013);
    //r_2 ~ uniform(0.0083, 0.013);
    b ~ gamma(0.05, 1);
    beta_v ~ lognormal(0,1);//normal(0.0001, 0.2);
    //gamma ~ uniform(0.03, 1.0);
    //gamma_f ~ uniform(0.0001, 0.5);
    theta_1 ~ uniform(0.0001, 1);
    mu ~ uniform(0.0001, 1);
    //E0 ~ uniform(1.038561e-06, 1.038561e-05);
    //IA0 ~ uniform(1.038561e-06, 1.038561e-05);
//    likelihood
    for (i in 1:n_obs){
       lambda[i] = y_hat[i, 6];
    }
    
    y ~ poisson(lambda);
    
}
generated quantities {
    real R_0;      // Basic reproduction number
    // real beta_p = theta[1];
    // real r_1 = theta[2];
    // real r_2 = theta[3];
    // real b = theta[4];
    // real beta_v = theta[5];
    // real gamma = theta[6];
    // real gamma_f = theta[7];

    R_0 = beta_p * beta_v * b /((b + r_1) * (gamma + gamma_f)* r_2);
}

