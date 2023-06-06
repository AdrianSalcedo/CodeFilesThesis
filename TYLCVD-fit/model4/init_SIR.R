init_sir <- function(){
  list(theta = c(
    beta_p = runif(1, 0, 0.1),
    r_1 = runif(1, 0, 0.2),
    r_2 = runif(1, 0, 0.2),
    b = runif(1, 0, 0.1),
    beta_v = runif(1, 0, 0.1),
    gamma = runif(1, 0, 0.5),
    gamma_f = runif(1, 0, 0.5),
    theta_1 = runif(1, 0, 1),
    mu =runif(1, 0, 1)
  )#,
  #E0 = runif(1,  74, 2100),
  #IA0 = runif(1, 74, 2100)
  )
}