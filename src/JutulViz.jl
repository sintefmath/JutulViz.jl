module JutulViz
    using Jutul
    using GLMakie
    using ColorSchemes
    using GraphMakie, Graphs, LayeredLayouts, NetworkLayout
    using .GLMakie
    using ColorSchemes
    export plot_well!, plot_interactive, plot_reservoir, plot_well_results
    export plot_variable_graph
    include("3d.jl")
    include("interactive.jl")
    include("wells.jl")
    # Plotting of graphs (experimental)
    include("plot_graph.jl")
    # Plots of iteration counts, timing etc
    include("performance.jl")
end # module
