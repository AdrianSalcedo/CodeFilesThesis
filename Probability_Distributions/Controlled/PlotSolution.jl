function PlotSolution(det_sol,stc_sol)
    title = plot(grid = false, showaxis = false, 
            bottom_margin = -50Plots.px)
    p1=plot(det_sol,vars=(1),color="blue")
    p1=plot!(stc_sol,vars=(1),color="darkgreen",title="Susc. p.")
    p2=plot(det_sol,vars=(2),color="blue")
    p2=plot!(stc_sol,vars=(2),color="darkorange", title ="Lat. p.")
    p3=plot(det_sol,vars=(3),color="blue")
    p3=plot!(stc_sol,vars=(3),color="darkred",title = "Infec. p.")
    p4=plot(det_sol,vars=(4),color="blue")
    p4=plot!(stc_sol,vars=(4),color="green", title = "Susc. v.")
    p5=plot(det_sol,vars=(5),color="blue")
    p5=plot!(stc_sol,vars=(5),color="red",title ="Infec. v.")
    plot(p1,p2,p3,p4,p5,layout = @layout([[A B C]; [D E]]),label = "")
end

