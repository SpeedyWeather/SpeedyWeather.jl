module SpeedyWeatherMakieExt

using SpeedyWeather, Makie
using DocStringExtensions

"""
$(TYPEDSIGNATURES)
Defines Makie's `heatmap` function for a`grid::AbstractGrid` via interpolation
to `::AbstractFullGrid` (which can be reshaped into a matrix.)"""
function Makie.heatmap(
    grid::RingGrids.AbstractGrid;
    title::String = "$(RingGrids.get_nlat(grid))-ring $(typeof(grid))",
    kwargs...   # pass on to Makie.heatmap
)
    full_grid = RingGrids.interpolate(RingGrids.full_grid(typeof(grid)), grid.nlat_half, grid)
    heatmap(full_grid; title, kwargs...)
end

"""
$(TYPEDSIGNATURES)
Defines Makie's `heatmap` function for a`grid::AbstractGrid` via interpolation
to `::AbstractFullGrid` (which can be reshaped into a matrix.)"""
function Makie.heatmap(
    grid::RingGrids.AbstractFullGrid;
    title::String = "$(RingGrids.get_nlat(grid))-ring $(typeof(grid))",
    size = (600,300),
    kwargs...
)

    mat = Matrix(grid)              # reshapes a full grid into a matrix
    lond = RingGrids.get_lond(grid) # get lon, lat axes in degrees
    latd = RingGrids.get_latd(grid)

    fig = Figure(size = size, figure_padding = 10)
    ax = Axis(fig[1, 1],
        aspect = 2,             # 0-360˚E -90-90˚N maps have an aspect of 2:1
        title = title,
        titlesize = 10,
        xticks = 0:60:360,      # label 0˚E, 60˚E, 120˚E, ...
        yticks = -60:30:60,     # label -60˚N, -30˚N, 0˚N, ... 
        xticklabelsize = 10,
        yticklabelsize = 10,
        xtickformat = values -> ["$(round(Int,value))˚E" for value in values],
        ytickformat = values -> ["$(round(Int,value))˚N" for value in values],
    )

    hm = heatmap!(ax, lond, latd, mat; kwargs...)
    Colorbar(fig[1, 2], hm, ticklabelsize = 10)
    colsize!(fig.layout, 1, Aspect(1, 2.0))

    resize_to_layout!(fig)
    return fig
end

end # module