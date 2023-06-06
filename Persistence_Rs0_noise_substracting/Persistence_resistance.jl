using DifferentialEquations
using Plots, PlotlyBase; plotly()
using DifferentialEquations.EnsembleAnalysis
using CSV
using IterableTables, DataFrames, DataTables
using StochasticDiffEq
using Distributions

path1 = "/home/gabrielsalcedo/Documentos/Julia_code_for_tomate_SDE_paper-main-master/"
path2 = "Persistence_Rs0_noise_substracting/"
path = path1 * path2

include(path * "Dynamics.jl")
circles = CSV.read(path * "resistence_data.csv", DataFrame)
#squarts = CSV.read(path * "datos_squares.csv", DataFrame)
#Parameters = CSV.read(path * "Parameter_Persistence_1.csv", DataFrame)
#=
beta_p = Parameters.beta_p[1]
r_1 =  Parameters.r_1[1]
r_2 =  Parameters.r_2[1]
b =  Parameters.b[1]
beta_v =  Parameters.beta_v[1]
theta =  Parameters.theta[1]
mu =  Parameters.mu[1]
gamma =  Parameters.gamma[1]
gamma_f =  Parameters.gamma_f[1]
N_v =  Parameters.N_v[1]
sigma_L =  Parameters.sigma_L[1]
sigma_I =  Parameters.sigma_I[1]
sigma_v = Parameters.sigma_v[1]
=#

beta_p = 5.0018440/100#0.0966796875
r_1 = 0.01#6.945409e-02#0.00018310546875
r_2 = 8.740866/10000#0.0068359375
b = 5.000/100#0.0546875
beta_v = 7.058968/1000#0.0113525
theta = 0.4
mu = 1.0
gamma = 0.16#0.03125
gamma_f = 0.17#0.00048828125
N_v =  mu / (gamma + gamma_f)
sigma_L = 0.350
sigma_I = 0.3
sigma_v = 0.3


#u_0 = [ 1.0, 0.0, 0.0, 0.0 , 0.0]
#u_0 = [907.9833883/1000,87.01661167/1000,5/1000,1000/50000,40000/50000,5.576870231/1000]
u_0 = [1000/1000,0/1000,0/1000,29892/51000,1000/51000,0/1000]

T = 70.0
time = (21.0, T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt = 0.01
t_s = range( 21.0, T, step = 1.0)

R0 = (beta_p * beta_v * b) / ((gamma + gamma_f) * (b + r_1) * r_2)
Rs0 = R0 -
    (1 / 2) * (
        (sigma_L + sigma_I) ^ 2 -
            sigma_v ^ 2 / (beta_v + sigma_v ^ 2 + theta * (gamma + gamma_f))
        )

################################################################################
######################### Solution computation #################################
########################## Deterministic Solution ##############################

prob_det = ODEProblem(F_Holt,u_0,time)
det_sol = solve(prob_det, AutoTsit5(Rosenbrock23()), dt = dt)
########################## Stochastis Solution #################################
prob_sde_tomato_sys = SDEProblem(F_Holt,G_Diffusion,u_0,time,
noise_rate_prototype=zeros(6,2))
stc_sol = solve(prob_sde_tomato_sys,SROCKC2(), dt = dt)
################################################################################
############################ PLot variables ####################################
################################################################################

title = plot(title = "R_s =$Rs0", grid = false, showaxis = false,
    bottom_margin = -50Plots.px)
p1 = plot(
det_sol, vars=(6), color="red", title = "Holt Sol", legend = false
)
p1 = plot!(
stc_sol, vars=(6), color="blue", title = "stoc Sol", legend = false
)
p1 = plot!(
    circles.Time_LA_1582, circles.LA_1582,
    seriestype = :scatter, color="red", legend = false
)
# p1 = plot!(
# squarts.t, squarts.y, seriestype = :scatter, color="black", legend = false
# )
################################################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################

Datos=DataFrame()
j = 0
trajectories = 1
while j <= 10000
    monte_prob = EnsembleProblem(prob_sde_tomato_sys)
    sim = solve(monte_prob, SROCKC2(),dt = dt,EnsembleThreads(),
        trajectories = trajectories)
     component = componentwise_vectors_timepoint(sim,t_s) #gives all solution in time vector t_s
     component = transpose(component) #transpose to obtain any*5 data matrix
     component = vcat(component...) #to obtain shape for dataframe
     component = vcat(component...) # again do a reshape
     variables = DataFrame(component, :auto)# define first data frame
     Y_p = variables[:,6]
     Y_p = min.(1.0,Y_p)
     Datos_aux = DataFrame(t = t_s, Sp = variables[:,1], Ip = variables[:,3],
      Iv = variables[:,5], Yp = Y_p) #only some variables
     Datos = append!(Datos, Datos_aux) #append the data in the loop
     j+=1
     println("acepted =",j)
end
# summ = EnsembleSummary(sim)
# plot!(summ;idxs=6)

# ################################################################################
# ######################    data  Det Solution    ################################
# ################################################################################
# det_Time = det_sol.t
# det_xu = det_sol.u
# det_xu_glued = hcat(det_xu...)
# Xu1 = det_xu_glued[1:5:end]
# Xu2 = det_xu_glued[2:5:end]
# Xu3 = det_xu_glued[3:5:end]
# Xu4 = det_xu_glued[4:5:end]
# Xu5 = det_xu_glued[5:5:end]
#
#
# det_DF1 = DataFrame(t = det_Time,S_p = Xu1, L_p =Xu2, I_p = Xu3, S_v = Xu4,
#     I_v = Xu5)
# det_DF1_red = det_DF1#[1:10:end,1:end]
# CSV.write(path * "Det_solution_Persistence_1.csv",det_DF1_red)
#
# ################################################################################
# ######################    data  Sto Solution    ################################
# ################################################################################
# stc_Time = stc_sol.t
# det_yu = stc_sol.u
# det_yu_glued = hcat(det_yu...)
# Yu1 = det_yu_glued[1:5:end]
# Yu2 = det_yu_glued[2:5:end]
# Yu3 = det_yu_glued[3:5:end]
# Yu4 = det_yu_glued[4:5:end]
# Yu5 = det_yu_glued[5:5:end]
#
#
# stc_DF1 = DataFrame(t = stc_Time, S_p = Yu1, L_p = Yu2, I_p = Yu3, S_v = Yu4,
#     I_v = Yu5)
# stc_DF1_red = stc_DF1[1:20:end,1:end]
# CSV.write(path * "Stc_Solution_Persistence_1.csv",stc_DF1_red)

################################################################################
########################## Monte  Carlo Ensamble ###############################
################################################################################


CSV.write(path * "Data_persistence_incidence_resistance_LA_1582_substracting.csv",Datos)
