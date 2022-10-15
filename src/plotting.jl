export plot_well!, plot_interactive, plot_reservoir, plot_well_results
using .GLMakie
using ColorSchemes

function plot_reservoir(g, states; wells = nothing, kwarg...)
    if haskey(first(states), :Reservoir)
        states = map(x -> x[:Reservoir], states)
    end
    fig, ax = plot_interactive(g, states; kwarg...)
    if !isnothing(wells)
        # Precompute geometry
        geo = tpfv_geometry(g)
        for w in wells
            if w["sign"] > 0
                c = :midnightblue
            else
                c = :firebrick
            end
            plot_well!(ax, g, w, color = c, geometry = geo)
        end
    end
    return (fig, ax)
end

function plot_interactive(model::MultiModel, states, model_key = nothing; kwarg...)
    if states isa AbstractDict
        states = [states]
    end
    if isnothing(model_key)
        model_key = first(keys(model.models))
    end
    if haskey(states[1], model_key)
        model_states = map(x -> x[model_key], states)
    else
        # Hope that the user knew what they sent in
        model_states = states
    end
    return plot_interactive(model[model_key], model_states; kwarg...)
end


function plot_interactive(model::SimulationModel, states; kwarg...)
    if states isa AbstractDict
        states = [states]
    end
    mesh = model.plot_mesh
    if isnothing(mesh)
        @warn "No plotting possible. SimulationModel has .plot_mesh = nothing." 
    else
        return plot_interactive(mesh, states; kwarg...)
    end
end

default_jutul_resolution() = (1600, 900)

function plot_interactive(grid, states; plot_type = nothing, wells = nothing, resolution = default_jutul_resolution(), alpha = 1.0, colormap = :viridis, alphamap = :no_alpha_map, kwarg...)
    pts, tri, mapper = triangulate_mesh(grid)

    fig = Figure(resolution = resolution)
    if states isa AbstractDict
        states = [states]
    end
    data = states[1]
    labels = Vector{String}()
    limits = Dict()
    for k in keys(data)
        d = data[k]
        if eltype(d)<:Real
            push!(labels, "$k")
            mv = Inf
            Mv = -Inf
            for s in states
                di = s[k]
                mv = min(minimum(di), mv)
                Mv = max(maximum(di), Mv)
            end
            limits["$k"] = (mv, Mv)
        else
            @debug "Skipping $k: Non-numeric type" eltype(d)
        end
    end
    function get_valid_rows(s)
        sample = data[Symbol(s)]
        if sample isa AbstractVector
            n = 1
        else
            n = size(sample, 1)
        end
        return ["$x" for x in 1:n]
    end
    datakeys = labels
    initial_prop = datakeys[1]
    state_index = Observable{Int64}(1)
    row_index = Observable{Int64}(1)
    prop_name = Observable{Any}(initial_prop)
    lims = Observable(limits[initial_prop])
    menu = Menu(fig, options = datakeys, prompt = initial_prop)
    menu_2 = Menu(fig, options = get_valid_rows("$initial_prop"), prompt = "1", width = 60)

    nstates = length(states)

    function change_index(ix)
        tmp = max(min(ix, nstates), 1)
        sl_x.selected_index = tmp
        state_index[] = tmp
        notify(state_index)
        return tmp
    end

    function increment_index(inc = 1)
        change_index(state_index.val + inc)
    end

    fig[5, 3] = hgrid!(
        menu,
        menu_2,
        ; tellheight = false, width = 300)

    sl_x = Slider(fig[5, 2], range = 1:nstates, value = state_index, snap = true)

    low = Observable{Float64}(0.0)
    hi = Observable{Float64}(1.0)

    rs_v = IntervalSlider(fig[4, :], range = LinRange(0, 1, 1000))

    on(rs_v.interval) do x
        low[] = x[1]
        hi[] = x[2]
    end
    # point = sl_x.value
    on(sl_x.selected_index) do n
        val = sl_x.selected_index.val
        state_index[] = val
    end
    if size(pts, 2) == 3
        ax = Axis3(fig[2, 1:3])
    else
        ax = Axis(fig[2, 1:3])
    end

    is_3d = size(pts, 2) == 3

    # Selection of data
    ys = @lift(
                mapper.Cells(
                    select_data(states[$state_index], Symbol($prop_name), $row_index, $low, $hi, limits[$prop_name])
                )
            )
    # Selection of colormap
    colormap_name = Observable(colormap)
    alphamap_name = Observable(alphamap)
    cmap = @lift(generate_colormap($colormap_name, $alphamap_name, alpha, $low, $hi))
    # Actual plotting call
    scat = Makie.mesh!(ax, pts, tri, color = ys, colorrange = lims, size = 60; shading = is_3d, colormap = cmap, kwarg...)
    Colorbar(fig[3, 1:3], scat, vertical = false)
    # Menu for field to plot
    on(menu.selection) do s
        rows = get_valid_rows(s)
        msel =  menu_2.selection[]
        if isnothing(msel)
            old = 1
        else
            old = parse(Int64, msel)
        end
        nextn = min(old, length(rows))
        prop_name[] = s
        row_index[] = nextn
        notify(prop_name)
        notify(menu_2.selection)
        menu_2.options = rows
        menu_2.selection[] = "$nextn"
        lims[] = limits[s]
    end
    # Row of dataset selector
    on(menu_2.selection) do s
        if isnothing(s)
            s = "1"
        end
        row_index[] = parse(Int64, s)
    end

    # Colormap selector
    colormaps = ["viridis", "jet", "balance", "autumn1", "hot", "winter", "terrain", "turbo", "gnuplot", "ocean", "vik"]
    cmap_str = "$colormap"
    if !(cmap_str in colormaps)
        push!(colormaps, cmap_str)
    end
    menu_cmap = Menu(fig[1, 3], options = colormaps, prompt = cmap_str)
    on(menu_cmap.selection) do s
        colormap_name[] = Symbol(s)
    end
    # Alpha map selector
    alphamaps = ["no_alpha_map", "linear", "linear_scaled", "inv_linear", "inv_linear_scaled"]
    amap_str = "$alphamap"
    if !(amap_str in alphamaps)
        push!(alphamaps, cmap_str)
    end
    menu_amap = Menu(fig[1, 1], options = alphamaps, prompt = amap_str)
    on(menu_amap.selection) do s
        alphamap_name[] = Symbol(s)
    end

    function loopy()
        start = state_index.val
        if start == nstates
            increment_index(-nstates)
            start = 1
        end
        previndex = start
        for i = start:nstates
            newindex = increment_index()
            if newindex > nstates || previndex != newindex-1
                break
            end
            notify(state_index)
            previndex = newindex
            sleep(1/30)
        end
    end

    fig[5, 1] = buttongrid = GridLayout()
    rewind = Button(fig, label = "⏪")
    on(rewind.clicks) do n
        increment_index(-nstates)
    end
    prev = Button(fig, label = "◀️")
    on(prev.clicks) do n
        increment_index(-1)
    end

    play = Button(fig, label = "⏯️")
    on(play.clicks) do n
        @async loopy()
    end
    next =   Button(fig, label = "▶️")
    on(next.clicks) do n
        increment_index()
    end
    ffwd = Button(fig, label = "⏩")
    on(ffwd.clicks) do n
        increment_index(nstates)
    end
    buttongrid[1, 1:5] = [rewind, prev, play, next, ffwd]
    display(GLMakie.Screen(), fig)
    return fig, ax
end

function select_data(state, fld, ix, low, high, limits)
    d = unpack(state[fld], ix)
    if low > 0.0 || high < 1.0
        d = copy(d)
        L, U = limits
        for i in eachindex(d)
            val = (d[i] - L)/(U - L)
            if val < low || val > high
                d[i] = NaN
            end
        end
    end
    return d
end

unpack(x, ix) = x[ix, :]
unpack(x::AbstractVector, ix) = x


function generate_colormap(colormap_name, alphamap_name, base_alpha, low, high)
    cmap = to_colormap(colormap_name)
    n = length(cmap)
    if alphamap_name != :no_alpha_map
        if alphamap_name == :linear
            F = x -> x
        elseif alphamap_name == :inv_linear
            F = x -> 1.0 - x
        elseif alphamap_name == :linear_scaled
            F = x -> clamp((x - low)./(high-low), 0.0, 1.0)
        elseif alphamap_name == :inv_linear_scaled
            F = x -> clamp(((1.0 - x) - high)./(low - high), 0.0, 1.0)
        else
            error()
        end
        u = range(0, 1, length = n)
        for (i, c) in enumerate(cmap)
            cmap[i] = GLMakie.RGBA{Float64}(c.r, c.g, c.b, base_alpha*F(u[i]))
        end
    end
    return cmap
end

function plot_well!(ax, g, w; color = :darkred, textcolor = nothing, name = nothing, linewidth = 5, top_factor = 0.2, textsize = 18, geometry = tpfv_geometry(g), kwarg...)
    if isnothing(textcolor)
        textcolor = color
    end
    centers = geometry.cell_centroids
    coord_range(i) = maximum(centers[i, :]) - minimum(centers[i, :])

    if size(centers, 1) == 3
        z = centers[3, :]
    else
        z = [0.0, 1.0]
    end
    bottom = maximum(z)
    top = minimum(z)

    # xrng = coord_range(1)
    # yrng = coord_range(2)

    rng = top - bottom
    s = top + top_factor*rng

    c = well_cells_for_plot(w)
    pts = centers[:, [c[1], c...]]
    if size(pts, 1) == 2
        # 2D grid, add some zeros to make logic work
        pts = vcat(pts, zeros(1, size(pts, 2)))
    end
    pts[3, 1] = s

    l = pts[:, 1]
    text!(well_name_for_plot(w, name), position = Tuple([l[1], l[2], -l[3]]), space = :data, color = textcolor, align = (:center, :baseline), textsize = textsize)
    lines!(ax, vec(pts[1, :]), vec(pts[2, :]), -vec(pts[3, :]), linewidth = linewidth, color = color, kwarg...)
end

well_name_for_plot(w::Dict, ::Nothing) = w["name"]
well_name_for_plot(w, s::String) = s
well_name_for_plot(w, s::Nothing) = String(w.name)

function well_cells_for_plot(w::Dict)
    wc = w["cells"]
    if !isa(wc, AbstractArray)
        wc = [wc]
    end
    return vec(Int64.(wc))
end

function well_cells_for_plot(w)
    return w.perforations.reservoir
end

export plot_well_results
function plot_well_results(well_data::Dict, arg...; name = "Data", kwarg...)
    plot_well_results([well_data], arg...; names = [name], kwarg...)
end

function plot_well_results(well_data::Vector, time = nothing; start_date = nothing,
                                                              names =["Dataset $i" for i in 1:length(well_data)], 
                                                              linewidth = 3,
                                                              cmap = nothing, 
                                                              dashwidth = 1,
                                                              styles = [:solid, :scatter, :dash, :dashdot, :dot, :dashdotdot],
                                                              resolution = default_jutul_resolution(),
                                                              kwarg...)
    # Figure part
    names = Vector{String}(names)
    ndata = length(well_data)
    is_inj = is_injectors(first(well_data))
    @assert ndata <= length(styles) "Can't plot more datasets than styles provided"
    fig = Figure(resolution = resolution)
    if isnothing(time)
        t_l = "Time-step"
        @assert isnothing(start_date) "start_date does not make sense in the absence of time-steps"
    elseif isnothing(start_date)
        t_l = "Time [days]"
    else
        t_l = "Date"
    end
    ax = Axis(fig[1, 1], xlabel = t_l)

    wd = first(well_data)
    # Selected well
    wells = sort!(collect(keys(wd)))
    nw = length(wells)
    if isnothing(cmap)
        if nw > 20
            c_key = :turbo
        elseif nw > 10
            c_key = :tab20
        else
            c_key = :Paired_10
        end
        cmap = cgrad(c_key, nw, categorical=true)
    end
    wellstr = [String(x) for x in wells]

    # Type of plot (bhp, rate...)
    responses = collect(keys(wd[first(wells)]))
    respstr = [String(x) for x in responses]
    response_ix = Observable(1)
    is_accum = Observable(false)
    is_abs = Observable(false)
    type_menu = Menu(fig, options = respstr, prompt = respstr[1])

    on(type_menu.selection) do s
        val = findfirst(isequal(s), respstr)
        response_ix[] = val
        autolimits!(ax)
    end
    use_two_cols = nw > 5
    if use_two_cols
        right_block = 2:3
    else
        right_block = 2:2
    end
    fig[2, right_block] = hgrid!(
        type_menu)

    b_xlim = Button(fig, label = "Reset x")
    on(b_xlim.clicks) do n
        reset_limits!(ax; xauto = true, yauto = false)
    end
    b_ylim = Button(fig, label = "Reset y")
    on(b_ylim.clicks) do n
        reset_limits!(ax; xauto = false, yauto = true)
    end
    toggle_abs = Toggle(fig, active = false)
    connect!(is_abs, toggle_abs.active)
    toggle_accum = Toggle(fig, active = false)
    connect!(is_accum, toggle_accum.active)

    buttongrid = GridLayout(tellwidth = false)
    buttongrid[1, 1] = toggle_abs
    buttongrid[1, 2] = Label(fig, "Absolute")
    buttongrid[1, 3] = toggle_accum
    buttongrid[1, 4] = Label(fig, "Cumulative")
    buttongrid[1, 5] = b_xlim
    buttongrid[1, 6] = b_ylim

    # Lay out and do plotting
    fig[2, 1] = buttongrid
    function get_data(time, wix, rix, dataix, use_accum, use_abs)
        tmp = well_data[dataix][wells[wix]][responses[rix]]
        if use_accum && respstr[dataix] != "Bottom hole pressure"
            T = [0.0, time[dataix]...]
            tmp = cumsum(tmp.*diff(T))
        end
        if use_abs
            tmp = abs.(tmp)
        end
        return tmp
    end

    sample = map(x -> get_data([], 1, 1, x, false, false), 1:ndata)
    nsample = map(length, sample)
    if isnothing(time)
        time = map(s -> 1:length(s), sample)
    else
        if eltype(time)<:AbstractFloat || eltype(time)<:Date
            time = repeat([time], ndata)
        end
    end

    newtime = []
    for i = 1:ndata
        T = time[i]
        nt = length(T)
        ns = nsample[i]
        @assert nt == ns "Series $i: Recieved $nt steps, but wells had $ns results."
        if eltype(T)<:AbstractFloat
            # Scale to days
            if isnothing(start_date)
                @. T /= (3600*24)
            else
                T = @. Microsecond(ceil(T*1e6)) + start_date
            end
            push!(newtime, T)
        end
    end
    time = newtime

    lighten(x) = GLMakie.ARGB(x.r, x.g, x.b, 0.2)
    toggles = Vector{Any}([Toggle(fig, active = true, buttoncolor = cmap[i], framecolor_active = lighten(cmap[i])) for i in eachindex(wells)])

    b_inj_on = Button(fig, label = "✔️ I")
    b_inj_off = Button(fig, label = "❌ I")

    labels = Vector{Any}([Label(fig, w) for w in wellstr])

    b_prod_on = Button(fig, label = "✔️ P")
    b_prod_off = Button(fig, label = "❌ P")

    buttongrid[1, 7] = b_inj_on
    buttongrid[1, 8] = b_inj_off
    buttongrid[1, 9] = b_prod_on
    buttongrid[1, 10] = b_prod_off


    function toggle_wells(do_injectors, status)
        for (i, w) in enumerate(wellstr)
            if is_inj[Symbol(w)] == do_injectors
                toggles[i].active[] = status
            end
        end
    end
    on(b_inj_on.clicks) do n
        toggle_wells(true, true)
    end
    on(b_inj_off.clicks) do n
        toggle_wells(true, false)
    end
    on(b_prod_on.clicks) do n
        toggle_wells(false, true)
    end
    on(b_prod_off.clicks) do n
        toggle_wells(false, false)
    end
    tmp = hcat(toggles, labels)
    bgrid = tmp
    N = size(bgrid, 1)

    if use_two_cols
        M = div(N, 2, RoundUp)
        fig[1, 2] = grid!(bgrid[1:M, :], tellheight = false)
        fig[1, 3] = grid!(bgrid[(M+1):N, :], tellheight = false)
    else
        fig[1, 2] = grid!(bgrid, tellheight = false)
    end

    lineh = []
    for dix = 1:ndata
        T = time[dix]
        for i in 1:nw
            d = @lift(get_data(time, i, $response_ix, dix, $is_accum, $is_abs))
            style = styles[dix]
            if style == :scatter
                h = scatter!(ax, T, d, color = cmap[i], linewidth = linewidth, marker = :circle)
            else
                if style == :dash || style == :dashdot
                    lw = dashwidth
                else
                    lw = linewidth
                end
                h = lines!(ax, T, d, linewidth = linewidth, linestyle = style, color = cmap[i])
            end
            t = toggles[i]
            connect!(h.visible, t.active)
            push!(lineh, h)
        end
    end
    x = (tmp, tmp2) -> autolimits!(ax)
    @lift(x($is_abs, $is_accum))
    if ndata > 1
        elems = []
        for i = 1:ndata
            style = styles[i]
            if style == :scatter
                el = MarkerElement(color = :black, linewidth = linewidth, marker = :circle)
            else
                if style == :dash || style == :dashdot
                    lw = dashwidth
                else
                    lw = linewidth
                end
                el = LineElement(color = :black, linestyle = style, linewidth = lw)
            end
            push!(elems, el)
        end
        fig[1, 1] = Legend(fig, elems, names,
                        tellheight = false,
                        tellwidth = false,
                        margin = (10, 10, 10, 10),
                        halign = :left, valign = :top, orientation = :horizontal
        )
    end
    display(GLMakie.Screen(), fig)
    return fig
end

function is_injectors(well_data)
    D = Dict{Symbol, Bool}()
    for (k, v) in well_data
        D[k] = sum(v[Symbol("Surface total rate")]) > 0
    end
    return D
end

function basic_3d_figure()
    fig = Figure()
    ax = Axis3(fig[1, 1])
    return (fig, ax)
end

export plot_mesh, plot_mesh!, plot_cell_data

function plot_mesh(m; kwarg...)
    fig, ax = basic_3d_figure()
    p = plot_mesh!(ax, m; kwarg...)
    display(fig)
    return (fig, ax, p)
end

function plot_mesh!(ax, m; color = :lightblue, kwarg...)
    pts, tri, mapper = triangulate_mesh(m)
    f = mesh!(ax, pts, tri; color = color, kwarg...)
    return f
end

function plot_cell_data(m, data; colorbar = :vertical, kwarg...)
    fig, ax = basic_3d_figure()
    p = plot_cell_data!(ax, m, data; kwarg...)
    if !isnothing(colorbar) && maximum(data) != minimum(data)
        if colorbar == :vertical
            Colorbar(fig[2, 1], p, vertical = false)
        else
            Colorbar(fig[1, 2], p, vertical = true)
        end
    end
    display(fig)
    return (fig, ax, p)
end

function plot_cell_data!(ax, m, data, kwarg...)
    pts, tri, mapper = triangulate_mesh(m)
    return mesh!(ax, pts, tri; color = mapper.Cells(data), kwarg...)
end
