function mup(y,u,t)
    dy = zeros(5,1)
    Nv = y[4] + y[5]
    dy[1] = - beta_p * y[1] * y[5]/Nv + (r_1 + u[1]) * y[2] + (r_2 + u[2]) * y[3]
    dy[2] = beta_p * y[1] * y[5]/Nv - (b + r_1 + u[1]) * y[2]
    dy[3] = b * y[2] - (r_2 + u[2]) * y[3]
    dy[4] = - beta_v * y[4] * y[3]/Np - (gamma + gamma_f + u[3]) * y[4] + 
            (1 - theta) * mu
    dy[5] = beta_v * y[4] * y[3]/Np - (gamma + gamma_f + u[3]) * y[5] + 
            theta * mu
    return dy
end # function

function sigma(y,u,t)
    Mdy=zeros(5,2)
    Mdy[1,1] = y[1] * (sigma_L * y[2] + sigma_I * y[3])/Np
    Mdy[1,2] = 0
    Mdy[2,1] = - sigma_L * y[1] * y[2]/Np
    Mdy[2,2] = 0
    Mdy[3,1] = - sigma_I * y[1] * y[3]/Np
    Mdy[3,2] = 0
    Mdy[4,1] = 0
    Mdy[4,2] = -sigma_v * y[4]
    Mdy[5,1] = 0
    Mdy[5,2] = -sigma_v * y[5]
    return Mdy
end # function

function sigma_prime(y, t)
    return 0
end # function

function dW(step_size,number_max_of_steps,scale)
    normal_sampler = sqrt(scale * step_size) * randn(number_max_of_steps-1)
    wt = zeros(number_max_of_steps)
    wt[2:end] = cumsum(normal_sampler)
    wt = wt / sqrt(scale)
    return wt
end # function