function ProbabilityLaw(time_idx,df,i::Int,statename::String)
    ProbLaw = histogram(df[time_idx,:],title = "Time =20", xlabel = statename, 
    ylabel = "Probability Law", normalize=:probability, legend = false, grid = false, 
    size = (718, 446))
    #xlims!(0, 1)
    ProbLaw
    savefig(path_wrt * "ProbLaw$(i-1).png")
end
###
function TrajectoriesPlot(df,i,xlabel::String,ylabel::String)
    Trayectory = plot(df[:,1], df[:,2:end], xlabel = xlabel, 
    ylabel = ylabel, legend = false, grid = false, 
    size = (718, 446))
    Trayectory
    savefig(path_wrt * "Trayectory$(i-1).png")
end