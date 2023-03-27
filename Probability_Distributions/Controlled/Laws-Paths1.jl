function ProbabilityLaw(time_idx,df,statename::String)
    ProbLaw = histogram(df[time_idx,:], xlabel = "statename", 
    ylabel = "Probability Law", normalize=:pdf, legend = false, grid = false)
    ProbLaw
   # savefig()
end
function TrajectoriesPlot(time,df,xlabel::String,ylabel::String)
    Trayectory = plot(time, df[time_idx,:], xlabel = "xlabel", 
    ylabel = "ylabel", normalize=:pdf, legend = false, grid = false)
    Trayectory
   # savefig()
end