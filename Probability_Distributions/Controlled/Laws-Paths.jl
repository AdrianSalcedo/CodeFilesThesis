function ProbabilityLaw(time_idx,df1,df2,i::Int,statename::String)
    Time_Value = df1[time_idx,1]
    ProbLaw = histogram(df1[time_idx,2:end],title = "Time = 20", xlabel = statename, 
    ylabel = "Density of Probability", normalize=:probability, legend = true, 
    fg_legend = :transparent, label = "Contrafactual", grid = false, 
    size = (718, 446),dpi=300)
    ProbLaw = histogram!(df2[time_idx,2:end],title = "Time =20", xlabel = statename, 
    fg_legend = :transparent, ylabel = "Density of Probability", 
    normalize=:probability, legend = true, label = "Controlled", grid = false, 
    size = (718, 446),dpi=300)
    figsize=(718,446)

    savefig(ProbLaw, path * "ProbLaws$(i-1).png")
end
###
function TrajectoriesPlot(df,i,xlabel::String,ylabel::String)
    Trayectory = plot(df[:,1], df[:,2:end], xlabel = xlabel, 
    ylabel = ylabel, legend = false, grid = false, size = (718, 446))
    Trayectory
    savefig(path * "Trayectory$(i-1).png")
end