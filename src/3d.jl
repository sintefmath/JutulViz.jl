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
