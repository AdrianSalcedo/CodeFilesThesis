using RecipesBase   # a tiny package defining the @recipe macro
sims = hcat((x->cumsum(randn(200))).(1:100)...) # produce 100 random walks with 200 steps
    @userplot SimPlot    # defines a plotting function called "simplot"
    
    @recipe function f(h::SimPlot; xlabel1 = "", xlabel2 = "") # define extra keywords to use in the plotting
        mat = h.args[1]      # the x, y, z data to be plotted are stored in the args array
    
        legend := false       # specify the plot attributes
        link := :y
        grid := false
        layout := grid(1, 2, widths = [0.7, 0.3])
    
        @series begin         # send the different data to the different subplots
            subplot := 2
            seriestype := :histogram
            orientation := :h
            xlabel := xlabel2
            title := ""
            ylabel := ""
            mat[end,:]
        end
    
        linealpha --> 0.4    # this (specifying the opacity of the line) can be overridden by the user
        seriestype := :path
        subplot := 1
        xlabel := xlabel1
        mat                         # the recipe returns the data to be plotted
    end

    using Plots; gr()
simplot(sims)