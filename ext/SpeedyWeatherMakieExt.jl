module SpeedyWeatherMakieExt

using SpeedyWeather
using Makie

function Makie.heatmap(
    grid::RingGrids.AbstractGrid;
    title::String = "$(RingGrids.get_nlat(grid))-ring $(typeof(grid))",
)
    full_grid = RingGrids.interpolate(RingGrids.full_grid(typeof(grid)), grid.nlat_half, grid)
    heatmap(full_grid; title, kwargs...)
end

function Makie.heatmap(
    grid::RingGrids.AbstractFullGrid;
    title::String = "$(RingGrids.get_nlat(grid))-ring $(typeof(grid))",
)
    mat = Matrix(grid)
    lond = RingGrids.get_lond(grid)
    latd = RingGrids.get_latd(grid)

    fig = Figure(size = (600, 300), figure_padding = 10)
    ax = Axis(fig[1, 1], aspect = 2,
        title = title,
        titlesize = 10,
        xticks = 0:60:360,
        yticks = -60:30:60,
        xticklabelsize = 10,
        yticklabelsize = 10,
        xtickformat = values -> ["$(round(Int,value))˚E" for value in values],
        ytickformat = values -> ["$(round(Int,value))˚N" for value in values],
    )

    hm = heatmap!(ax, lond, latd, mat)
    Colorbar(fig[1, 2], hm, ticklabelsize = 10)
    colsize!(fig.layout, 1, Aspect(1, 2.0))

    resize_to_layout!(fig)
    return fig
end