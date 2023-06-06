library(sde)
library(plotly)
library(ggplot2)
library(MASS)
library(deSolve)
library(stats4)
library(fda)
library(splines)
library(fds)
library(Matrix)
library(rainbow)
library(pcaPP)
library(RCurl)

df<-as.matrix(read.csv("simulatedTrajectory.csv",header = T))
T = 70
N = length(df[,1])-1
dt = T/N

y_init = df[1,2:6]

beta_p = 0.1
r_1 = 0.01
r_2 = 0.01
b = 0.075
beta_v = 0.01
gamma = 0.03
gamma_f = 0.03
theta = 0.4
mu = 1.0
sigma_L= 0.118417
sigma_I = 0.990762
sigma_V = 0.652238
N_p = y_init[1]+y_init[2]+y_init[3]


#R_0 = sqrt(beta_v * mu * b * beta_p / (r_l * r_l * ( r_i + b) * gamma))
#R_0_s
#print(R_0)
#print(R_0_s)
################################################################################
Drift<-function(y,t){
  s_p = y[1]
  l_p = y[2]
  i_p = y[3]
  s_v = y[4]
  i_v = y[5]

  s_p_prime_d = -beta_p * s_p * i_v/(s_v+i_v) + r_1 * l_p + r_2 * i_p
  l_p_prime_d = beta_p * s_p * i_v/(s_v+i_v) - b * l_p - r_1 * l_p
  i_p_prime_d = b * l_p - r_2 * i_p
  s_v_prime_d = - beta_v * s_v * i_p/(s_p+l_p+i_p) - (gamma+gamma_f) * s_v 
    + (1 - theta) * mu
  i_v_prime_d = beta_v * s_v * i_p/(s_p+l_p+i_p) - (gamma+gamma_f) * i_v 
    + theta * mu
  rhs_d_np_array = c(s_p_prime_d, l_p_prime_d, i_p_prime_d, s_v_prime_d, 
                             i_v_prime_d)
  return(rhs_d_np_array)
  }

################################################################################
Diffusion <- function(y,t){
  s_p = y[1]
l_p = y[2]
i_p = y[3]
s_v = y[4]
i_v = y[5]
G_matrix = matrix(c((s_p/N_p) *(sigma_L*l_p+sigma_I*i_p),-sigma_L *(s_p/N_p) * l_p,-sigma_I * (s_p/N_p) *i_p, 0,
                            0,0,0,0, -sigma_V * s_v, -sigma_V * i_v),5,2)
  
  #matrix(c((s_p/N_p) *(sigma_L*l_p+sigma_I*i_p), 0,
            #          -sigma_L *(s_p/N_p) * l_p, 0,
             #        -sigma_I * (s_p/N_p) *i_p, 0,
              #        0, -sigma_V * s_v,
               #       0, -sigma_V * i_v),5,2)
# s_p_prime_s = sigma * (N_p - s_p)
# l_p_prime_s = -sigma_l * l_p
#i_p_prime_s = -sigma_i * i_p
#s_v_prime_s = -sigma_v * s_v
# i_v_prime_s = -sigma_v * i_v
rhs_s_np_array = G_matrix #np.array([s_p_prime_s, l_p_prime_s, i_p_prime_s, s_v_prime_s, i_v_prime_s])
return(rhs_s_np_array)
}


Bt_1 = df[,7]#brownian_path_sampler(dt,N) #generating the path
Bt_2 = df[,8]#brownian_path_sampler(dt,N)
  dB_1 = diff(Bt_1)
  dB_2 = diff(Bt_2)

dB = matrix(0,2,N+1)
dB[1,2:length(df[,1])] = dB_1
dB[2,2:length(df[,1])] = dB_2

ts = seq(0, T, length.out = N+1)
ys = matrix(0,5,N+1)
TL <- length(ts)
ys[,1] = y_init

#for i in range(num_sims):
  
  for (i in 2:TL){
  t = (i-1) * dt
  y = ys[,i-1]
  ys[,i] = y + Drift(y,t) * dt + Diffusion(y,t) %*% dB[,i-1]
}
ggplot() + 
  geom_line(aes(ts, ys[1,], colour=g1), ys) + 
  geom_line(aes(ts, df[,1], colour=g2), df)
plot(ts, ys[3,], type = "l", lty = 1) + plot(df[,1],df[,4],type = "l", lty = 1)

matplot(ys[,], type = "l", lwd=1)
data_ys <- t(ys)
write.csv(data_ys, "No_controlSDE1.csv")






