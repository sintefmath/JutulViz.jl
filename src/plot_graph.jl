function plot_variable_graph(model)
    graph, nodes, = build_variable_graph(model, to_graph = true)
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


