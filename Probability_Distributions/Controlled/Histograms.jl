using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV, StatsPlots, Random
using DataFrames, Images, PGFPlotsX
using StochasticDiffEq
PLOTS_DEFAULTS = Dict(:dpi => 1200,:width => 2226, :height => 1492)
path = "C:/Users/Usuario1/Desktop/Probability_Distributions/No_control/"

include(path * "Dynamics.jl")
include(path * "Laws-Paths.jl")
include(path * "MakerTrayectories.jl")
include(path * "PlotSolution.jl")
#control = CSV.read(path * "Controls.csv", DataFrame)

beta_p = 0.1
r_1 = 0.01
r_2 = 0.01
b = 0.075
beta_v = 0.01
theta = 0.4
mu = 1.0
gamma = 0.03
gamma_f = 0.03
sigma_L = 0.1118417
sigma_I = 0.990762
sigma_v = 0.652238

u_0 = [650.0,250.0,100.0,300.0,200.0]
T = 20.0
time = (0.0,T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt=0.1

tol = 1e-08

################################################################################

R0 = (beta_p * beta_v * b) / ( (gamma + gamma_f) * ( b + r_1) * r_2)
println("R0=",R0);
Rs0 = R0 -
    (1 / 2) * (
        (sigma_L + sigma_I) ^ 2 -
            sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f))
        )
println("Rs0=",Rs0);
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift, G_Diffusion, u_0, time,
                    noise_rate_prototype = zeros(5,2))
# ####################################################################################
# ###########################   Parallel   Ensamble   ################################
# ####################################################################################
Datos = MakerTrayectories(T,1000)
########################## Here we compute Probability Laws ####################
t_s = range(0.0,T, step=0.1)
time = DataFrame(t = t_s)
Nrow = size(Datos)[2]
t_steps = size(time)[1]
Ncol = Int(size(Datos)[1]/t_steps)
i=4
    Datos_State = reshape(Datos[:,i],t_steps,Ncol)
    Datos_State = DataFrame(Datos_State, :auto)
    Datos_State = hcat(time,Datos_State)
    D1_last = DataFrame(Datos_State[201,2:end])
    CSV.write(path * "DataState$(i-1)NC.csv",Datos_State)
    CSV.write(path * "DataState$(i-1)NC_last.csv",D1_last)
    ProbabilityLaw(201,Matrix(Datos_State[:,2:end]),i,"Ip")
    TrajectoriesPlot(Matrix(Datos_State),i,"Time (Days)","Ip")

################ Controlled ###########################

path2 = "C:/Users/Usuario1/Desktop/BocopHJB-Tomato-adimensional/Trayectorias/"
Tiempo = CSV.read(path2 * "simulatedTrajectorytimes.csv",DataFrame)
#States = CSV.read(path1 * "simulatedTrajectory.csv",DataFrame)
Tiempo = Matrix(Tiempo)
#States = Matrix(States)
Datos2 = zeros(size(Tiempo)[1],1006)
Datos2[:,1] .= Tiempo
for i in 1:1005
    path1 = "C:/Users/Usuario1/Desktop/BocopHJB-Tomato-adimensional/Trayectorias/T-$i/"
    States = CSV.read(path1 * "simulatedTrajectory.csv",DataFrame)
    Datos2[:,i+1] .= States[:,3]
end

########################## Here we compute Probability Laws ####################
t_s = range(0.0,T, step=dt)
time = DataFrame(t = t_s)
Nrow = size(Datos2)[1]
t_steps = size(time)[1]
Ncol = Int(size(Datos2)[1]/t_steps)
##
plotlyjs()
Datos_State2 = DataFrame(Datos2,:auto)
D_aux2 = 1000*Datos2[2:end,201]
i=4
#for i in 2:6
    #Datos_State = reshape(Datos[:,i-1],t_steps,Ncol)
    #Datos_State = DataFrame(Datos_State, :auto)
    #Datos_State = hcat(time,Datos_State)
    D2_last = permutedims(D_aux2)
    CSV.write(path * "DataState$(i-1)C.csv",Datos_State2)
    CSV.write(path * "DataState$(i-1)C_last.csv", D2_last)
    ProbabilityLaw2(201,Matrix(Datos_State[:,2:end]),1000*Matrix(Datos_State2),i,"Ip")
    TrajectoriesPlot(Matrix(Datos_State),i,"Time (Days)","Ip")
#end