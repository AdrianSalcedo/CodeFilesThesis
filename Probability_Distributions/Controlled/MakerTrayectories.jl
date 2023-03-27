function MakerTrayectories(T,iterations)

    Datos=DataFrame()
    j = 0
    i = 1
    trajectories = 1
    t_s = range(0.0,T, step=1.0)
    while j <= iterations
        Ensemble_prob = EnsembleProblem(prob_sde_tomato_sys)
        sim = solve(Ensemble_prob, SROCKC2(),dt= dt,EnsembleThreads(),trajectories=trajectories)
        component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
        component = transpose(component) #transpose to obtain any*5 data matrix
        component = vcat(component...) #to obtain shape for dataframe
        component = vcat(component...) # again do a reshape
        variables = DataFrame(component,:auto) # define first data frame
        time_aux = t_s
        Datos_aux = DataFrame(t = time_aux, S_p = variables[:,1], L_p = variables[:,2],
            I_p = variables[:,3], S_v = variables[:,4], I_v = variables[:,5]) #only some variables
        cond1 = all(variables[:,1].>=0)
        cond2 = all(variables[:,2].>=0)
        cond3 = all(variables[:,3].>=0)
        cond4 = all(variables[:,4].>=0)
        cond5 = all(variables[:,5].>=0)
        condition6 = cond1 && cond2 && cond3 && cond4 && cond5
        condition1 =  variables[:,1] + variables[:,2] + variables[:,3]
        condition2 = condition1-N_p*ones(size(condition1))
        condition3 = maximum(abs.(condition2))
        condition4 = condition3 <= tol
        if ((condition4 == true) && (condition6 == true))
            Datos = append!(Datos, Datos_aux) #append the data in the loop
            j+=1
            println("acepted =",j)
        end
        println(i)
        i+=1
    end
    return Datos
end