function ProbabilityLaw(time_idx,df,i::Int,statename::String)
    ProbLaw = histogram(df[time_idx,:],title = "Time = 20", xlabel = statename, 
    ylabel = "Density of Probability", normalize=:probability, legend = false, grid = false, 
    fg_legend = :transparent, size = (718, 446),dpi=1600)
    histogram!(size=(3*718, 3*446),dpi=600)
    figsize=(3*718, 3*446)
    savefig(ProbLaw,path * "ProbLaw$(i-1).png")
end
function ProbabilityLaw2(time_idx,df1,df2,i::Int,statename::String)
    PyPlot.plt[:legend](bbox_to_anchor = (1.05,1), loc =1)
    PLOTS_DEFAULTS = Dict(:dpi => 1200,:width => 2226, :height => 1492)
    Time_Value = df1[time_idx,1]
    ProbLaw = histogram(df1[time_idx,2:end],title = "Time = 20", xlabel = statename, 
    ylabel = "Density of Probability", normalize=:probability, legend=:outertopleft, 
    label = "Counterfactual", grid = false, size = (2126,1392),dpi=1200)
    ProbLaw = histogram!(df2[time_idx,2:end],title = "Time = 20", xlabel = statename, 
    fg_legend = :transparent, ylabel = "Density of Probability", normalize=:probability, 
    legend=:outertopleft, label = "Optimal Control", grid = false, size = (2126,1392), dpi=1200)
    figsize=(2126,1392)
    savefig(ProbLaw, path * "ProbLaws$(i-1).png")
end
###
function TrajectoriesPlot(df,i,xlabel::String,ylabel::String)
    Trayectory = plot(df[:,1], df[:,2:end], xlabel = xlabel, 
    ylabel = ylabel, normalize=:pdf, legend = false, grid = false, 
    size = (718, 446))
    Trayectory
    savefig(path * "Trayectory$(i-1).png")
end