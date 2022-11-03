export plot_solve_breakdown, plot_cumulative_solve
function plot_solve_breakdown(allreports, names; per_it = false, include_local_solves = nothing)
    if per_it
        plot_title = "Time per iteration"
        to_plot = x -> [x.assembly/x.its, x.subdomains/x.its, x.solve/x.its, x.total/x.its]
    else
        plot_title = "Total time"
        to_plot = x -> [x.assembly, x.subdomains, x.solve, x.total]
    end
    labels = ["Assembly", "Local solves", "Linear solve", "Total"]

    D = map(x -> to_plot(timing_breakdown(x)), allreports)
    ndata = length(D)
    nel = length(D[1])

    if isnothing(include_local_solves)
        include_local_solves = sum(x -> x[2], D) > 0.0
    end
    colors = Makie.wong_colors(1.0)
    if include_local_solves
        colors = colors[[1, 2, 3, 6]]
    else
        subs = [1, 3, 4]
        colors = colors[[1, 2, 3]]
        labels = labels[subs]
        nel = length(subs)
        D = map(x -> x[subs], D)
    end

    h = vcat(D...)
    x = ones(nel)
    for i = 2:ndata
        x = vcat(x, i*ones(nel))
    end
    grp = repeat(1:nel, ndata)

    fig = Figure()
    ax = Axis(fig[1,1], xticks = (1:ndata, names), ylabel = "Time [s]", title = plot_title)
    barplot!(ax, x, h,
            dodge = grp,
            color = colors[grp])

    elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    title = nothing# "Legend"
    Legend(fig[2,1], elements, labels, title, orientation = :horizontal)
    display(fig)
    return (fig, D)
end

function plot_cumulative_solve(allreports, dt, names = nothing; use_time = false)
    r_rep = map(x -> timing_breakdown(x, reduce = false), allreports)
    if use_time
        F = D -> map(x -> x.total, D)
        yl = "Wall time [s]"
    else
        F = D -> map(x -> x.its, D)
        yl = "Iterations"
    end
    t = cumsum(vcat(0, dt))/(3600*24*365)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time [years]", ylabel = yl)
    get_data = x -> cumsum(vcat(0, F(x)))

    if isnothing(names)
        names = map(x -> "dataset $x", eachindex(allreports))
    end
    alldata = []
    for i in eachindex(allreports)
        data_i = get_data(r_rep[i])
        push!(alldata, data_i)
        lines!(ax, t, data_i, label = names[i], linewidth = 3.5)
        scatter!(ax, t, data_i)
    end
    axislegend(ax, position = :lt)
    display(fig)
    return (fig, alldata, t)
end