functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

    //real dy_dt[6];
    real Sp = y[1];
    real Lp = y[2];
    real Ip = y[3];
    real Sv = y[4];
    real Iv = y[5];
    real Ipc = y[6];
    real N = x_i[1];
      
    real beta_p = theta[1];
    real r_1 = theta[2];
    real r_2 = theta[3];
    real b = theta[4];
    real beta_v = theta[5];
    real gamma = theta[6];
    real gamma_f = theta[7];
    real theta_1 = theta[8];
    real mu = theta[9];
    
    real dSp_dt = -beta_p * Sp * Iv + r_1 * Lp + r_2 * Ip;
    real dLp_dt = beta_p * Sp * Iv - (b + r_1) * Lp;
    real dIp_dt = b * Lp - r_2 * Ip;
    real dSv_dt = -beta_v * Sv * Ip - (gamma + gamma_f) * Sv + (1 - theta_1) * mu;
    real dIv_dt = beta_v * Sv * Ip - (gamma + gamma_f) * Iv + theta_1 * mu;
    real dIpc_dt = b * Lp;
    
    return {dSp_dt, dLp_dt, dIp_dt, dSv_dt, dIv_dt, dIpc_dt};
  }
}
data {
  int<lower=1> n_days;
  real y0[6];
  real t0;
  real ts[n_days];
  int N;
  int cases[n_days];
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower = 0, upper = 0.1> beta_p;
  real<lower = 0.0083, upper = 0.013> r_1;
  real<lower = 0.0083, upper = 0.013> r_2;
  real<lower = 0.05, upper = 0.1> b;
  real<lower = 0, upper = 0.2> beta_v;
  real<lower = 0.03, upper = 1> gamma;
  real<lower = 0, upper = 1> gamma_f;
  real<lower = 0, upper = 1> theta_1;
  real<lower = 0, upper = 1> mu;
  //real<lower = 0, upper = 1> Lp0;
  real<lower=0> phi_inv;
}
transformed parameters{
  real<lower = 0, upper = 1.0> y[n_days, 6];
  real incidence[n_days-1];
  real phi = 1. / phi_inv;
  {
    real theta[9];
    theta[1] = beta_p;
    theta[2] = r_1;
    theta[3] = r_2;
    theta[4] = b;
    theta[5] = beta_v;
    theta[6] = gamma;
    theta[7] = gamma_f;
    theta[8] = theta_1;
    theta[9] = mu;

    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
    //print(y);
    for (i in 1:(n_days-1)){
    incidence[i] =-(y[i+1, 2] - y[i, 2] + y[i+1, 1] - y[i, 1]) ; //-(E(t+1) - E(t) + S(t+1) - S(t))
    }
  }
}
model {
  real lambda[n_days];
  //priors
  beta_p ~ normal(0.01,0.001);//uniform(0,1);
  r_1 ~ exponential(0.01);//
  r_2 ~ exponential(0.01);//
  b ~ exponential(0.075);//
  beta_v ~ normal(0.003,0.001);
  gamma ~ uniform(0.03, 1);//
  gamma_f ~ uniform(0, 1);
  theta_1 ~ uniform(0, 1);
  mu ~ uniform(0, 1);
  //Lp0 ~ uniform(0.001, 0.003);
  // Sv0 ~ uniform(0.0, 0.001);
  // Iv0 ~ uniform(0.0, 0.001);
  phi_inv ~ exponential(5);
  
  //sampling distribution
  for (i in 1:n_days){
       lambda[i] = y[i, 6] * N;
       // lambda[i] = y_hat[i, 6];
  }
    cases ~ poisson(lambda);
    print(cases);
  //cases[1:(n_days-1)] ~ poisson(incidence);
  //cases[1:(n_days)] ~ neg_binomial_2(incidence, phi);
  //cases[1:(n_days)]~poisson(col(to_matrix(y), 6)*N, phi);
  
}
generated quantities {
  // real R0 = beta / gamma;
  // real recovery_time = 1 / gamma;
  real pred_cases[n_days-1];
  // 
  // //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  pred_cases = poisson_rng(incidence);
  //pred_cases = neg_binomial_2_rng(incidence, phi);
  //pred_cases = neg_binomial_2_rng(col(to_matrix(y), 6)*N, phi);
}

