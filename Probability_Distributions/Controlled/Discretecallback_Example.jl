using DifferentialEquations
using Plots

function odm2prod(dx, x, params, t)
    k_1, f_1, V_liq, X_in, Y_in, q_in = params

    rho_1 = k_1*x[1]
    q_prod = 0.52*f_1*x[1]
    # Differential Equations
    dx[1] = q_in/V_liq*(X_in - x[1]) - rho_1
    dx[2] = q_in/V_liq*(Y_in - x[2])
end

x0      = [3.15, 1.5]
tspan   = (0.0, 7.0)
params  = [0.22, 43, 155, 249, 58, 0]
prob    = ODEProblem(odm2prod, x0, tspan, params)

input   = [1.0 60; 1.1 0; 2.0 60; 2.3 0; 4.0 430; 4.05 0]
dosetimes = input[:,1]
function affect!(integrator)
    ind_t = findall(t -> t==integrator.t, dosetimes)
    integrator.p[6] = input[ind_t[1], 2]
end
# function affect!(integrator)
#     ind_t = findall(integrator.t == dosetimes)
#     integrator.p[6] = input[ind_t, 2]
# end
cb = PresetTimeCallback(dosetimes, affect!)
sol = solve(prob, Tsit5(), callback=cb, saveat=1/12)

plot(sol, vars=[1, 2])