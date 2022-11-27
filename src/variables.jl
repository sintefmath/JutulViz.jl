export plot_secondary_variables
function plot_secondary_variables(model::SimulationModel; kwarg...)
    plot_secondary_variables(MultiModel((model = model, )); kwarg...)
end

function plot_secondary_variables(model::MultiModel; resolution = default_jutul_resolution(), linewidth = 4, kwarg...)
    data = Dict{String, Any}()
    for (k, m) in pairs(model.models)
        for (vname, var) in Jutul.get_secondary_variables(m)
            d = line_plot_data(m, var)
            if !isnothing(d)
                if d isa JutulLinePlotData
                    d = [d]
                end
                data["$k.$vname"] = d
            end
        end
    end
    fig = Figure(resolution = resolution)
    dkeys = keys(data)
    default = first(dkeys)
    m = Menu(fig[1, 1], options = dkeys, default = default)
    ax = Axis(fig[2, 1])
    colors = Makie.wong_colors()
    # slider = Slider(fig[1:2, 2], range = 0.5:0.1:10, startvalue = 1, horizontal = false)
    function plot_var!(datakey)
        di = first(data[datakey])
        empty!(ax)
        ax.xlabel = di.xlabel
        ax.ylabel = di.ylabel
        ax.title = di.title
            # ax = Axis(plot_grid[1, i], xlabel = di.xlabel, ylabel = di.ylabel, title = di.title)
        labels_found = false
        ix = 1
        for (x, y, lbl) in zip(di.xdata, di.ydata, di.datalabels)
            c = colors[mod(ix, 7) + 1]
            lines!(ax, x, y; color = c, linewidth = linewidth, label = lbl, kwarg...)
            # text!(ax, x[end], y[end], text = lbl)
            ix += 1
        end
        # axislegend()
    end
    on(m.selection) do datakey
        plot_var!(datakey)
    end
    plot_var!(default)
    display(fig)
    return fig
end
