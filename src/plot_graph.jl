function plot_variable_graph(model)
    graph, nodes, = build_variable_graph(model, to_graph = true)
    colors = Vector{Symbol}()
    for n in nodes
        if n in keys(model.primary_variables)
            c = :red
        elseif n in keys(model.secondary_variables)
            c = :blue
        else
            c = :black
        end
        push!(colors, c)
    end
    plot_jutul_graph(graph, nodes, colors)
end

function plot_jutul_graph(graph, nodes, colors = [:black for _ in nodes]; kwarg...)
    xs, ys, paths = solve_positions(Zarate(), graph)
    lay = _ -> Point.(zip(xs,ys))
    # create a vector of Point2f per edge
    wp = [Point2f.(zip(paths[e]...)) for e in Graphs.edges(graph)]
    alignments = []
    for i in xs
        if i == 1
            al = (:right, :center)
        else
            al = (:center, :bottom)
        end
        push!(alignments, al)
    end
    N = length(nodes)
    graphplot(graph, layout=lay,
                     waypoints=wp,
                     nlabels_distance=10,
                     nlabels_textsize=20,
                     node_size = [20 for i in 1:N],
                     edge_width= [3 for i in 1:ne(graph)],
                     edge_color = :grey80,
                     node_color = colors,
                     nlabels_align= alignments,
                     nlabels = map(String, nodes))
end

export plot_model_graph
function plot_model_graph(model; kwarg...)
    plot_variable_graph(model; kwarg...)
end

function plot_model_graph(model::MultiModel)
    equation_to_label(k, eq_name) = "$k: $eq_name"
    edges = Vector{String}()
    nodes = Vector{String}()

    for (k, m) in pairs(model.models)
        for (e, eq) in m.equations
            push!(nodes, equation_to_label(k, e))
        end
    end

    n = length(nodes)
    directed = true
    if directed
        graph = SimpleDiGraph(n)
    else
        graph = SimpleGraph(n)
    end

    to_index(k, eq) = findfirst(isequal(equation_to_label(k, eq)), nodes)
    function add_cross_term_edge!(ct, target, source, equation)
        T = to_index(target, equation)
        F = to_index(source, equation)
        ct_bare_type = Base.typename(typeof(ct)).name
        if !has_edge(graph, T, F)
            add_edge!(graph, T, F)
            push!(edges, "$ct_bare_type"[1:end-2])
        end
    end

    for ctp in model.cross_terms
        (; cross_term, target, source, equation) = ctp
        add_cross_term_edge!(cross_term, target, source, equation)
        if Jutul.has_symmetry(cross_term) && directed
            add_cross_term_edge!(cross_term, source, target, equation)
        end
    end
    layout = SFDP(Ptype=Float32, tol=0.01, C=0.01, K=1)
    layout = Shell()
    # layout = SquareGrid()
    return graphplot(graph, nlabels = nodes,
                            elabels = edges,
                            nlabels_distance=15,
                            nlabels_textsize=20,
                            node_size = [30 for i in 1:n],
                            node_color = collect([Float64(i) for i in 1:n]),
                            layout=layout
                            )
end
