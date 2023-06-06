using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV, StatsPlots, Random
using DataFrames, Images
using StochasticDiffEq
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

path = "D:/CodeFilesThesis/Probability_Distributions/Controlled/"

include(path * "Dynamics_Controlleds.jl")
include(path * "Laws-Paths.jl")
include(path * "EM-version.jl")
include(path * "PlotSolution.jl")
control = CSV.read(path * "Controls.csv", DataFrame)

################################################################################
################################# Parameters ###################################
################################################################################
# beta_p = 0.1
# r_1 = 0.01
# r_2 = 0.01
# b = 0.075
# beta_v = 0.01
# theta = 0.4
# mu = 1.0
# gamma = 0.03
# gamma_f = 0.03
# sigma_L = 0.1118417
# sigma_I = 0.990762
# sigma_v = 0.652238

t_init = 0
t_end  = 20
N      = 200 # number of time steps
dt     = float(t_end - t_init) / N # time step
# y_init = [650.0,250.0,100.0,300.0,200.0]
# Np = y_init[1] + y_init[2]+ y_init[3]
# tol = 1e-08

scale = 1.0
step_size = dt
number_max_of_steps = N
T = scale*step_size*number_max_of_steps
time = LinRange(0,T,number_max_of_steps+1)
ts = LinRange(t_init, t_end,N+1)
ys = zeros(N+1,5)
ys[1,:] = y_init
us = Matrix(control)
################################################################################

R0 = (beta_p * beta_v * b) / ( (gamma + gamma_f) * ( b + r_1) * r_2)
println("R0=",R0);
Rs0 = R0 -
    (1 / 2) * (
        (sigma_L + sigma_I) ^ 2 -
            sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f))
        )
println("Rs0=",Rs0);
######################### Solution computation #################################
########################## Stochastis Solution #################################


################################################################################
############################ PLot variables ####################################
#PlotSolution(det_sol,stc_sol)
################################################################################
######################    data  Det Solution    ################################
################################################################################
#Data_Det = rename!(DataFrame(det_sol),[:t, :S_p,:Lp,:Ip,:Sv,:Iv])
#CSV.write(path * "Det_Solution.csv",Data_Det)
# ################################################################################
# ######################    data  Sto Solution    ################################
# ################################################################################
#Data_Stc = rename!(DataFrame(stc_sol),[:t, :S_p,:Lp,:Ip,:Sv,:Iv])
#CSV.write(path * "Stc_Solution.csv",Data_Stc)

# ####################################################################################
# ###########################   Parallel   Ensamble   ################################
# ####################################################################################
#Datos = EulerMaruyama(10000)
########################## Here we compute Probability Laws ####################
t_s = range(0.0,T, step=dt)
time = DataFrame(t = t_s)
Nrow = size(Datos2)[1]
t_steps = size(time)[1]
Ncol = Int(size(Datos2)[1]/t_steps)
##
Datos_State2 = DataFrame(Datos2,:auto)
i=4
#for i in 2:6
    #Datos_State = reshape(Datos[:,i-1],t_steps,Ncol)
    #Datos_State = DataFrame(Datos_State, :auto)
    #Datos_State = hcat(time,Datos_State)
    CSV.write(path * "DataState$(i-1).csv",Datos_State2)
    ProbabilityLaw(201,1000*Matrix(Datos_State2),i,"Ip")
    TrajectoriesPlot(Matrix(Datos_State),i,"Time (Days)","Ip")
#end