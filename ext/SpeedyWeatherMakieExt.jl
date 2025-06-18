module SpeedyWeatherMakieExt

using SpeedyWeather, Makie
using DocStringExtensions

"""
$(TYPEDSIGNATURES)
Defines Makie's `heatmap` function for a`field::AbstractField2D` via interpolation
to `::AbstractFullField2D` (which can be reshaped into a matrix.)"""
function Makie.heatmap(
    field::RingGrids.AbstractField2D;
    title::String = "$(RingGrids.get_nlat(field))-ring $(typeof(field))",
    kwargs...   # pass on to Makie.heatmap
)
    full_field = RingGrids.interpolate(RingGrids.full_grid_type(field.grid), field.grid.nlat_half, field)
    heatmap(full_field; title, kwargs...)
end

"""
$(TYPEDSIGNATURES)
Defines Makie's `heatmap` function for a `field::AbstractFullField2D` which can be reshaped into a matrix."""
function Makie.heatmap(
    field::RingGrids.AbstractFullField2D;
    title::String = "$(RingGrids.get_nlat(field))-ring $(typeof(field))",
    size = (600,300),
    kwargs...
)

    mat = Matrix(field)                 # reshapes a full field into a matrix
    lond = RingGrids.get_lond(field)    # get lon, lat axes in degrees
    latd = RingGrids.get_latd(field)

    fig = Figure(size = size, figure_padding = 10)
    ax = Axis(fig[1, 1],
        aspect = 2,             # 0-360˚E -90-90˚N maps have an aspect of 2:1
        title = title,
        titlesize = 10,
        xticks = 0:60:360,      # label 0˚E, 60˚E, 120˚E, ...
        yticks = -60:30:60,     # label -60˚N, -30˚N, 0˚N, ... 
        xticklabelsize = 10,
        yticklabelsize = 10,
        xtickformat = values -> ["$(round(Int, value))˚E" for value in values],
        ytickformat = values -> ["$(round(Int, value))˚N" for value in values],
    )

    hm = heatmap!(ax, lond, latd, mat; kwargs...)
    Colorbar(fig[1, 2], hm, ticklabelsize = 10)
    colsize!(fig.layout, 1, Aspect(1, 2.0))

    resize_to_layout!(fig)
    return fig
end

end # module