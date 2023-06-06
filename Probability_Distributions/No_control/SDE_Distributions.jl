using DifferentialEquations
using Plots; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV, StatsPlots
using DataFrames
using StochasticDiffEq

path = "C:/Users/Usuario1/Desktop/Probability_Distributions/No_control/"
path_wrt = "D:/CodeFilesThesis/Probability_Distributions/No_control/"
include(path_wrt * "Dynamics.jl")
include(path_wrt * "Laws-Paths.jl")
include(path_wrt * "MakerTrayectories.jl")
include(path_wrt * "PlotSolution.jl")
#Parameters = CSV.read(path * "Parameter_extinction_noise.csv", DataFrame)
################################################################################
################################# Parameters ###################################
################################################################################
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

u_0 = [650.0,250.0,100.0,5000.0,5000.0]

T = 20.0
time = (0.0,T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt=0.1
u_0_adim = [650.0/N_p,250.0/N_p,100.0/N_p, 0.5,0.5]
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
######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Holt, u_0, time)
det_sol = solve(prob_det, Tsit5(), dt = dt)

########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Drift_adim, G_Diffusion_adim, u_0_adim, time,
                    noise_rate_prototype = zeros(5,2))
stc_sol = solve(prob_sde_tomato_sys, EM(),abstol =tol, dt = dt)
plot(stc_sol)
################################################################################
############################ PLot variables ####################################
PlotSolution(det_sol,stc_sol)
################################################################################
######################    data  Det Solution    ################################
################################################################################
Data_Det = rename!(DataFrame(det_sol),[:t, :S_p,:Lp,:Ip,:Sv,:Iv])
CSV.write(path * "Det_Solution.csv",Data_Det)
# ################################################################################
# ######################    data  Sto Solution    ################################
# ################################################################################
Data_Stc = rename!(DataFrame(stc_sol),[:t, :S_p,:Lp,:Ip,:Sv,:Iv])
CSV.write(path * "Stc_Solution.csv",Data_Stc)

# ####################################################################################
# ###########################   Parallel   Ensamble   ################################
# ####################################################################################
Datos = MakerTrayectories(T,1000)
########################## Here we compute Probability Laws ####################
t_s = range(0.0,T, step=0.1)
time = DataFrame(t = t_s)
Nrow = size(Datos)[2]
t_steps = size(time)[1]
Ncol = Int(round(size(Datos)[1]/t_steps))

for i in 2:6
    Datos_State = reshape(Datos[:,i],t_steps,Ncol)
    Datos_State = DataFrame(Datos_State, :auto)
    Datos_State = hcat(time,Datos_State)
    CSV.write(path_wrt * "DataState$(i-1).csv",Datos_State)
    ProbabilityLaw(t_steps,Matrix(Datos_State[:,2:end]),i,"State$(i-1)")
    TrajectoriesPlot(Matrix(Datos_State),i,"Time (Days)","State $(i-1)")
end