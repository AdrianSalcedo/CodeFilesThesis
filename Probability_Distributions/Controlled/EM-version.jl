function EulerMaruyama(iterations)
    Datos=DataFrame()
    j = 0
    i = 1
    trajectories = 1
#    for i in 1:num_sims
#         Btp = dW(dt,N+1,scale)
#         Btv = dW(dt,N+1,scale)
#         for j in 2:size(ts)[1]
#             t = t_init + (j - 1) * dt
#             y = ys[j - 1,:]
#             u = us[j-1,:]
#             dBtp = Btp[j] - Btp[j-1]
#             dBtv = Btv[j] - Btv[j-1]
#             dBt = zeros(2,1)
#             dBt[1,1] =  dBtp
#             dBt[2,1] =  dBtp
#             ys[j,:] = y + mup(y,u,t) * dt + sigma(y,u,t) * dBt 
#             #+(1 / 2) * sigma(y,u,t) * sigma_prime(y,u,t) * (dBt ^ 2 - dt)
#             #plot!(ts, ys)
#         end # for
#     #plot!(ts, ys[:,3])
#     end # for 
#     return ys
    while j <= iterations
        Btp = dW(dt,N+1,scale)
        Btv = dW(dt,N+1,scale)
        for j in 2:size(ts)[1]
            t = t_init + (j - 1) * dt
            y = ys[j - 1,:]
            u = us[j-1,:]
            dBtp = Btp[j] - Btp[j-1]
            dBtv = Btv[j] - Btv[j-1]
            dBt = zeros(2,1)
            dBt[1,1] =  dBtp
            dBt[2,1] =  dBtp
            ys[j,:] = y + mup(y,u,t) * dt + sigma(y,u,t) * dBt 
         end # for
        cond1 = all(ys[:,1].>=0)
        cond2 = all(ys[:,2].>=0)
        cond3 = all(ys[:,3].>=0)
        cond4 = all(ys[:,4].>=0)
        cond5 = all(ys[:,5].>=0)
        condition6 = cond1 && cond2 && cond3 && cond4 && cond5
        condition1 =  ys[:,1] + ys[:,2] + ys[:,3]
        condition2 = condition1-Np*ones(size(condition1))
        condition3 = maximum(abs.(condition2))
        condition4 = condition3 <= tol
        if ((condition4 == true) && (condition6 == true))
            Datos_aux = DataFrame(ys,:auto)
            Datos = append!(Datos, Datos_aux) #append the data in the loop
            j+=1
            println("acepted =",j)
        end
        println(i)
        i+=1
    end
    return Datos
end
