function ProbabilityLaw(time_idx,df,i::Int,statename::String)
    Time_Value = df[time_idx,1]
    ProbLaw = histogram(df[time_idx,2:end],title = "Time =$(Time_Value)", xlabel = statename, 
    ylabel = "Probability Law", normalize=:probability, legend = false, grid = false, 
    size = (718, 446))
    ProbLaw
    savefig(path * "ProbLaw$(i-1).png")
end
###
function TrajectoriesPlot(df,i,xlabel::String,ylabel::String)
    Trayectory = plot(df[:,1], df[:,2:end], xlabel = xlabel, 
    ylabel = ylabel, legend = false, grid = false, size = (718, 446))
    Trayectory
    savefig(path * "Trayectory$(i-1).png")
end