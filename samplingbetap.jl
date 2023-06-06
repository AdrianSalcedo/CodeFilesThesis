Betap = LinRange(0.01,0.12,1000)
len = length(Betap)
Len_series = 1:len
p1 = plot!(
    circles.t, circles.y, seriestype = :scatter, color="red", legend = false
        )
p1 = plot!(
    squarts.t, squarts.y, seriestype = :scatter, color="black", legend = false
        )
ylims!(-0.01, 1.05)
for i in Len_series
    beta_p = Betap[i]
    prob_Holt = ODEProblem(F_Holt,u_0,time)
    det_Holt = solve(prob_Holt, Tsit5(), dt = dt)
    p1 = plot!(
        det_Holt,vars=(6), color="blue", title = "Holt System", legend = false
            )
end
