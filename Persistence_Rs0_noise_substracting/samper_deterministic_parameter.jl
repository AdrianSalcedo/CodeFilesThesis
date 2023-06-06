using CSV, IterableTables, DataFrames, DataTables, Distributions
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis
using  StochasticDiffEq, Plots; plotly()

path1 = "/home/gabrielsalcedo/Documentos/Julia_code_for_tomate_SDE_paper-main-master/"
path2 = "Persistence_Rs0_noise_substracting/"
path = path1 * path2
include(path * "Compute_deterministic_R0.jl")
include(path * "Dynamics.jl")
circles = CSV.read(path * "datos_circle.csv", DataFrame)
squarts = CSV.read(path * "datos_squares.csv", DataFrame)


function Sampler_persistence_parameters(N_p)
    # Input:
    #   N_p: plant size
    # Output:
    #   parameters: vector with parameters
Test1 = false
while Test1 == false
    beta_p = rand(Uniform(0.023860126,0.044936614))
    r_1 = rand(Uniform(0.0083965,0.012871061))
    b = rand(Uniform(0.050616588,0.078723884))
    r_2 = rand(Uniform(0.008457268,0.012892753))
    beta_v = rand(Uniform(0.008038328,0.011953985))
    theta = 0.4
    mu = 1.0
    gamma = 0.06
    gamma_f = rand(Uniform(0,0.5))

    par = DataFrame(beta_p = beta_p, r_1 = r_1, b = b, r_2 = r_2,
     beta_v = beta_v, theta = theta, mu = mu, gamma = gamma,
     gamma_f = gamma_f)
     cond1 = Compute_deterministic_R0(par) > 1

    Test1 = cond1

    if Test1 == true
        return par
    end
end
end
par = Sampler_persistence_parameters(1)

u_0 = [907.9833883/1000,87.01661167/1000,5/1000,1000/50000,40000/50000,5.576870231/1000]
T = 70.0
time = (7.0,T)
N_p = u_0[1] + u_0[2] + u_0[3]
dt = 0.01
t_s = range(0.0,T, step = 1.0)
prob_det = ODEProblem(F_Holt,u_0,time)
det_sol = solve(prob_det, Tsit5(), dt = dt)
prob_det2 = ODEProblem(F_Drift,u_0,time)
det_sol2 = solve(prob_det2, Tsit5(), dt = dt)
p1 = plot(
det_sol, vars=(6), color="blue", title = "Holt System", legend = false
)
p1 = plot!(
det_sol2, vars=(6), color="blue", title = "Holt System", legend = false
)
p1 = plot!(
circles.t, circles.y, seriestype = :scatter, color="red", legend = false
)
p1 = plot!(
squarts.t, squarts.y, seriestype = :scatter, color="black", legend = false
)
