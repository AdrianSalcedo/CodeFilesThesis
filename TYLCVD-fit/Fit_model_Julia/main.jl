using Downloads, DataFrames, CSV, Chain, Dates
using Plots, StatsPlots, LaTeXStrings
using DifferentialEquations
using Turing, AbstractMCMC
using LazyArrays
using Random:seed!
using StatsPlots
using Plots
using Distributed
addprocs(4)

path_reader = "/home/gabrielsalcedo/Documentos/TYLCVD-fit/Fit_model_Julia/"

df = CSV.File(path_reader * "interpolated_data.csv") |> DataFrame

# Here is a plot of the data:
dates = df.Time
data = df.Interpolate_Ip


@df df plot(:Time,
            :Interpolate_Ip,
            xlab=L"t", ylab="infected obs",
            label=false)
savefig(path_reader * "infected_obs.svg"); # hide
###############################################################################
include(path_reader*"det_tomato_tyclv.jl")

################################################
function NegativeBinomial2(μ₁, ϕ)
    p = 1 / ϕ
    r = μ₁

    return NegativeBinomial(r, p)
end
function Poisson2(λ)

    return Poisson(λ)
end

seed!(123)
setprogress!(true) # hide

@model function bayes_seiv(infected_obs, I₀, Nₚ, Nᵥ)
    # Prior distributions.

    βₚ ~ truncated(Normal(0.07, 0.005), 1e-4, 0.1)
    r₁ ~ truncated(Normal(0.01, 0.01), 1e-4, 0.1)
    r₂ ~ truncated(Normal(0.01, 0.001), 1e-4, 0.03)
    b  ~ InverseGamma(16.0, 1.2)
    βᵥ ~ truncated(Normal(0.1, 0.01), 1e-5, 0.2)
    γ  ~ InverseGamma(12, 1.0)
    γᵥ ~ InverseGamma(24, 3)
    θ  ~ truncated(Normal(0.0028, 0.015), 1e-4, 1.0)
    μ ~  truncated(Normal(0.39, 0.005), 1e-4, 0.5)

    # ϕ⁻ ~ truncated(Exponential(5), 0, 1e5)
    # ϕ = 1.0 / ϕ⁻
    # initial Conditions
    E₀ ~ truncated(Normal(9, 1), 1e-4, 10)
    Iᵥ ~ truncated(Normal(160, 2), 1e-4, 170)
    Dᵥ ~ truncated(Normal(1183, 1), 1e-4, 1200)
    ϕ⁻ ~ Exponential(5)
    ϕ = 1 / ϕ⁻

    par = [βₚ, r₁, r₂, b, βᵥ, γ, γᵥ, θ, μ]
    I = I₀
    E = E₀
    Yp = I
    Sv = Nᵥ - (Iᵥ + Dᵥ)
    Iv = Iᵥ
    Dv = Nᵥ - (Sv + Iv)
    l = length(infected_obs)
    u0 = [Nₚ - (E + I), E, I, Yp, Sv, Iv, Dv]
    tspan = (1.0, float(l))
    prob_ = ODEProblem(ode_rhs!, u0, tspan, par)
    sol = solve(prob_, Vern9(), saveat=1.0)
    #sol = solve(prob_, Vern7(lazy=false), saveat=1.0)
    solᵢ = sol[4, :]
    solᵢ = max.(1e-4, solᵢ)
    infected_obs ~ arraydist(LazyArray(@~ NegativeBinomial2.(solᵢ, ϕ)))
    #infected_obs ~ arraydist(LazyArray(@~ Poisson2.(solᵢ)))
end;

########################################################################

I₀ = first(df[:, :Interpolate_Ip])
Nₚ = 800
Nᵥ = Nₚ * 5
infected_obs = df[:, :Interpolate_Ip]
model_seiv =  bayes_seiv(infected_obs, I₀, Nₚ, Nᵥ)
chain_seiv = sample(model_seiv, NUTS(0.65), 1000)

sample(
    model_seiv,
    NUTS(0.65),
    MCMCDistributed(),
    100, 2; progress=true
)
summarystats(chain_seiv[[:βₚ, :r₁, :r₂, :b, :βᵥ, :γ, :γᵥ, :θ, :μ]])

plot(chain_seiv)
posterior_samples =
    chain_seiv[[:βₚ, :r₁,:r₂, :b, :βᵥ, :γ, :γᵥ, :θ, :μ]][rand(1:1000,300),:,:]

plot(; legend=false)
for p in eachrow(Array(posterior_samples))
    sol_p = solve(prob, Tsit5(); p=p, saveat=1)
    plot!(sol_p.t, sol_p[4,:]; alpha=0.1, color="#BBBBBB")
end
plot!(sol_p.t, infected_obs; alpha=1, color="#000000")
savefig("fit.svg"); # hide
