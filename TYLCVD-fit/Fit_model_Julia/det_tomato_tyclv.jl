function ode_rhs!(du,u,p,t)
    βₚ = p[1]
    r₁ = p[2]
    r₂ = p[3]
    b = p[4]
    βᵥ = p[5]
    γ = p[6]
    γᵥ = p[7]
    θ = p[8]
    μ = p[9]
    Nₚ = u[1] + u[2] + u[3]
    Nᵥ = u[5] + u[6]
 @inbounds begin
        du[1] = - βₚ * u[1] * u[6] / Nᵥ + r₁ * u[2] + r₂ * u[3]
        du[2] = βₚ * u[1] * u[6] / Nᵥ - (b + r₁) * u[2]
        du[3] = b * u[2] - r₂ * u[3]
        du[4] = b * u[2]
        du[5] = - βᵥ * u[5] * u[3] / Nₚ - (γ + γᵥ) * u[5] + (1 - θ) * μ
        du[6] = βᵥ * u[5] * u[3] / Nₚ - (γ + γᵥ) * u[6] + θ * μ
        du[7] = (γ + γᵥ) * (u[5] + u[6]) - μ
    end
    nothing
end
