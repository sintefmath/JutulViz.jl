module JutulViz
    using Jutul
    using GLMakie
    using ColorSchemes
    # Conditional code loading
    using Requires

    include("plotting.jl")
    # Plotting of graphs (experimental)
    function __init__()
        @require GraphRecipes="bd48cda9-67a9-57be-86fa-5b3c104eda73" begin
            @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plot_graph.jl")
        end
    end
end # module
