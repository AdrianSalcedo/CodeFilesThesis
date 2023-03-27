function ProbabilityLaw(time_idx,df,i::Int,statename::String)
    ProbLaw = histogram(df[time_idx,:],title = "Time =70", xlabel = statename, 
    ylabel = "Probability Law", normalize=:pdf, legend = false, grid = false, 
    size = (718, 446))
    ProbLaw
    savefig(path * "ProbLaw$(i-1).png")
end
###
function TrajectoriesPlot(df,i,xlabel::String,ylabel::String)
    Trayectory = plot(df[:,1], df[:,2:end], xlabel = xlabel, 
    ylabel = ylabel, normalize=:pdf, legend = false, grid = false, 
    size = (718, 446))
    Trayectory
    savefig(path * "Trayectory$(i-1).png")
end