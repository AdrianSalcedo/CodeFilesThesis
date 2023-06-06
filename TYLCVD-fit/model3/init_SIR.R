init_sir <- function(){
  list(
    theta = c(
      beta_p = runif(1, 0.0001, 0.1),
      r_1 = runif(1, 0.0083, 0.013), 
      r_2 = runif(1, 0.0083, 0.013),
      b = runif(1, 0.05, 0.1),
      beta_v = runif(1, 0.0001, 0.2),
      gamma = runif(1, 0.03, 1.0),
      gamma_f = runif(1, 0.001, 0.5),
      theta_1 = runif(1, 0.0001, 0.75),
      mu = runif(1, 0.0001, 0.75)
    )
  )
}