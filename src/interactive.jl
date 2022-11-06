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

function plot_interactive(grid, states; plot_type = nothing,
                                        primitives = nothing,
                                        transparency = false,
                                        resolution = default_jutul_resolution(),
                                        alpha = 1.0,
                                        title = "",
                                        colormap = :viridis,
                                        alphamap = :no_alpha_map,
                                        kwarg...)
    has_primitives = !isnothing(primitives)
    if grid isa Integer
        # Assume that someone figured out primitives already...
        nc = grid
        @assert has_primitives
        if isnothing(plot_type)
            plot_type = :mesh
        end
    else
        nc = number_of_cells(grid)
    end
    if !has_primitives
        if isnothing(plot_type)
            plot_candidates = [:mesh, :meshscatter, :lines]
            for p in plot_candidates
                primitives = plot_primitives(grid, p)
                if !isnothing(primitives)
                    plot_type = p
                    break
                end
            end
            if isnothing(primitives)
                @warn "No suitable plot found for mesh of type $(typeof(grid)). I tried $plot_candidates"
                return
            end
        else
            primitives = plot_primitives(grid, plot_type)
            if isnothing(primitives)
                @warn "Mesh of type $(typeof(grid)) does not support plot_type :$plot_type"
                return
            end
        end
    end
    pts = primitives.points
    mapper = primitives.mapper

    fig = Figure(resolution = resolution)
    if states isa AbstractDict
        states = [states]
    end
    if eltype(states)<:Number && (length(states) == nc || size(states, 1) == nc)
        states = [Dict(:Data => states)]
    end
    data = states[1]
    labels = Vector{String}()
    limits = Dict()
    for k in keys(data)
        d = data[k]
        is_valid_vec = d isa AbstractVector && length(d) == nc
        is_valid_mat = d isa AbstractMatrix && size(d, 2) == nc
        if eltype(d)<:Real && (is_valid_vec || is_valid_mat) 
            push!(labels, "$k")
            mv = Inf
            Mv = -Inf
            for s in states
                di = s[k]
                mv = min(minimum(x -> isnan(x) ? Inf : x, di), mv)
                Mv = max(maximum(x -> isnan(x) ? -Inf : x, di), Mv)
            end
            if mv == Mv
                Mv = 1.01*mv + 1e-12
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
    is_3d = size(pts, 2) == 3
    if is_3d
        make_axis = Axis3
    else
        make_axis = Axis
    end
    ax = make_axis(fig[2, 1:3], title = title)

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
    fig[1, 1] = lmap = GridLayout()
    menu_amap = Menu(lmap[1, 1], options = alphamaps, prompt = amap_str)
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

    # Actual plotting call
    if plot_type == :mesh
        tri = primitives.triangulation
        scat = Makie.mesh!(ax, pts, tri; color = ys,
                                        colorrange = lims,
                                        size = 60,
                                        shading = is_3d,
                                        colormap = cmap,
                                        transparency = transparency,
                                        kwarg...)
    elseif plot_type == :meshscatter
        sz = 0.8.*primitives.sizes
        npts, d = size(pts)
        if d < 3
            pts = hcat(pts, zeros(npts, 3 - d))
            sz = hcat(sz, ones(npts, 3 - d))
        end 
        sizes = zeros(GLMakie.Vec3f, size(sz, 1))
        for i in eachindex(sizes)
            sizes[i] = GLMakie.Vec3f(sz[i, 1], sz[i, 2], sz[i, 3])
        end
        scat = Makie.meshscatter!(ax, pts; color = ys,
                                        colorrange = lims,
                                        markersize = sizes,
                                        shading = is_3d,
                                        colormap = cmap,
                                        transparency = transparency,
                                        kwarg...)
    elseif plot_type == :lines
        x = pts[:, 1]
        y = pts[:, 2]
        z = pts[:, 3]
        scat = Makie.lines!(ax, x, y, z, color = ys,
                                                    linewidth = 15,
                                                    transparency = transparency,
                                                    colormap = cmap,
                                                    colorrange = lims)
        txt = primitives.top_text
        if !isnothing(txt)
            top = vec(pts[1, :])
            text!(txt,
                    position = Tuple([top[1], top[2], top[3] + 2.0]),
                    space = :data,
                    align = (:center, :baseline)
                    )
        end
        if primitives.marker_size > 0
            Makie.scatter!(ax, x, y, z, marker_size = primitives.marker_size, color = :black, alpha = 0.5, overdraw = true)
        end
    else
        error("Unsupported plot_type $plot_type")
    end

    Colorbar(fig[3, 1:3], scat, vertical = false)

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

function basic_3d_figure()
    fig = Figure()
    ax = Axis3(fig[1, 1])
    return (fig, ax)
end

export plot_multimodel_interactive
function plot_multimodel_interactive(model, states, model_keys = keys(model.models); plot_type = :mesh, shift = Dict(), kwarg...)
    n = length(model_keys)
    primitives = Vector{Any}(undef, n)
    ncells = zeros(Int64, n)
    active = BitArray(undef, n)
    active .= false
    for (i, k) in enumerate(model_keys)
        p = model[k].plot_mesh
        if isnothing(p)
            keep = false
        else
            primitive = plot_primitives(p, plot_type)
            keep = !isnothing(primitive)
            if keep
                nc = maximum(primitive.mapper.indices.Cells)
                ncells[i] = nc
                primitives[i] = primitive
                keep = keep && nc > 0
            end
        end
        active[i] = keep
    end
    model_keys = model_keys[active]
    primitives = primitives[active]
    ncells = ncells[active]
    # Remap states so that we have NaN padded versions
    offsets = cumsum(vcat([1], ncells))
    states_mapped = Vector{Dict{Symbol, Any}}()
    # Find all possible state fields
    all_state_fields = []
    state = states[1]
    for (i, model_key) in enumerate(model_keys)
        nc = ncells[i]
        for (k, v) in state[model_key]
            valid_vector = v isa AbstractVector && length(v) == nc
            valid_matrix = v isa AbstractMatrix && size(v, 2) == nc

            if valid_vector 
                push!(all_state_fields, (k, 1))
            elseif valid_matrix
                push!(all_state_fields, (k, size(v, 1)))
            end
        end
    end
    # Create flattened states with NaN for missing data
    total_number_of_cells = sum(ncells)
    new_states = Vector{Dict{Symbol, Any}}()
    for state in states
        new_state = Dict{Symbol, Any}()
        for (state_field, d) in all_state_fields
            data = zeros(d, total_number_of_cells)
            data .= NaN
            for (i, model_key) in enumerate(model_keys)
                state_m = state[model_key]
                if haskey(state_m, state_field)
                    old_data = state_m[state_field]
                    data[:, offsets[i]:(offsets[i+1]-1)] = old_data
                end
            end
            new_state[state_field] = data
        end
        push!(new_states, new_state)
    end
    # Merge triangulations
    face_index = Vector{Int64}()
    points = map(x -> x.points, primitives)
    tri = map(x -> x.triangulation, primitives)
    cell_index = map(x -> x.mapper.indices.Cells, primitives)
    offset = 0
    for (pts, T, k) in zip(points, tri, model_keys)
        @. T += offset
        if haskey(shift, k)
            for (i, dx) in enumerate(shift[k])
                @. pts[:, i] += dx
            end
        end
        offset += size(pts, 1)
    end
    offset = 0
    for (nc, cix) in zip(ncells, cell_index)
        @. cix += offset
        offset += nc
    end
    tri = vcat(tri...)
    points = vcat(points...)
    cell_index = vcat(cell_index...)

    mapper = (
                Cells = (cell_data) -> cell_data[cell_index],
                Faces = (face_data) -> face_data[face_index],
                indices = (Cells = cell_index, Faces = face_index)
              )
    acc_primitives = (points = points, triangulation = tri, mapper = mapper)
    plot_interactive(total_number_of_cells, new_states, primitives = acc_primitives; kwarg...)
end
