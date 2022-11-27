export plot_secondary_variables
function plot_secondary_variables(model::SimulationModel)
    plot_secondary_variables(MultiModel((model = model, )))
end

function plot_secondary_variables(model::MultiModel; resolution = default_jutul_resolution(), kwarg...)
    data = Dict{String, Any}()
    for (k, m) in pairs(model.models)
        for (vname, var) in Jutul.get_secondary_variables(m)
            d = JutulLinePlotData(var)
            if !isnothing(d)
                if d isa JutulLinePlotData
                    d = [d]
                end
                data["$k.$vname"] = d
            end
        end
    end
    fig = Figure(resolution = resolution)
    ax = []
    dkeys = keys(data)
    default = first(dkeys)
    m = Menu(fig[1, 1], options = dkeys, default = default)
    # slider = Slider(fig[1:2, 2], range = 0.5:0.1:10, startvalue = 1, horizontal = false)

    function plot_var!(datakey)
        d = data[datakey]
        n = length(d)
        plot_grid = fig[2, 1] = GridLayout(1, n)
        for (i, di) in enumerate(d)
            ax = Axis(plot_grid[1, i], xlabel = di.xlabel, ylabel = di.ylabel, title = di.title)
            for (x, y, lbl) in zip(di.xdata, di.ydata, di.datalabels)
                lines!(ax, x, y; label = lbl, kwarg...)
            end
            axislegend()
        end
    end
    on(m.selection) do datakey
        plot_var!(datakey)
    end
    plot_var!(default)
    display(fig)
    return fig
end
