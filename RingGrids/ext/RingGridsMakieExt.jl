module RingGridsMakieExt

using RingGrids
using Makie
using DocStringExtensions

function default_title(field::RingGrids.AbstractField)
    Grid = RingGrids.nonparametric_type(field.grid)
    NF = eltype(field)
    return "$(RingGrids.get_nlat(field))-ring Field{$NF} on $Grid"
end

# ============================================================================
# Non-mutating heatmap variants
# ============================================================================

function Makie.heatmap(
        field::RingGrids.AbstractField;
        title::String = default_title(field),
        kwargs...   # pass on to Makie.heatmap
    )
    @warn "Field of size $(size(field)) provided, 2D horizontal Field{T, 1} expected, selecting first indices of additional dimensions."
    inds = (Colon(), ntuple(_ -> 1, ndims(field) - 1)...)
    return heatmap(RingGrids.field_view(field, inds...))
end

"""
$(TYPEDSIGNATURES)
Defines Makie's `heatmap` function for a`field::AbstractField2D` via interpolation
to `::AbstractFullField2D` (which can be reshaped into a matrix.)
"""
function Makie.heatmap(
        field::RingGrids.AbstractField2D;
        title::String = default_title(field),
        kwargs...   # pass on to Makie.heatmap
    )
    full_field = RingGrids.interpolate(RingGrids.full_grid_type(field.grid), field.grid.nlat_half, field)
    return heatmap(full_field; title, kwargs...)
end


"""
$(TYPEDSIGNATURES)
Defines Makie's `heatmap` function for a `field::AbstractFullField2D` which can be reshaped into a matrix.
"""
function Makie.heatmap(
        field::RingGrids.AbstractFullField2D;
        title::String = default_title(field),
        size = (600, 300),
        kwargs...
    )
    fig = Figure(size = size, figure_padding = 10)
    _, hm = heatmap!(fig[1, 1], field; title, kwargs...)
    Colorbar(fig[1, 2], hm, ticklabelsize = 10)
    colsize!(fig.layout, 1, Aspect(1, 2.0))
    resize_to_layout!(fig)
    return fig
end

# ============================================================================
# Mutating heatmap variants
# ============================================================================

function Makie.heatmap!(pos::Makie.GridPosition, field::RingGrids.AbstractField2D; kwargs...)
    full_field = RingGrids.interpolate(RingGrids.full_grid_type(field.grid), field.grid.nlat_half, field)
    return heatmap!(pos, full_field; kwargs...)
end

function Makie.heatmap!(pos::Makie.AbstractAxis, field::RingGrids.AbstractField2D; kwargs...)
    full_field = RingGrids.interpolate(RingGrids.full_grid_type(field.grid), field.grid.nlat_half, field)
    return heatmap!(pos, full_field; kwargs...)
end

"""
$(TYPEDSIGNATURES)
Mutating variant of `heatmap` for RingGrids `Field`s. Returns both the `Axis ` as well as the
plotting object returned by `Makie.heatmap!`.
"""
function Makie.heatmap!(
        pos::Makie.GridPosition,
        field::RingGrids.AbstractFullField2D;
        title::String = default_title(field),
        aspect = 2,
        titlesize = 10,
        xticklabelsize = 10,
        yticklabelsize = 10,
        axis_kwargs = (;),
        kwargs...
    )
    ax = Axis(
        pos;
        aspect,             # 0-360˚E -90-90˚N maps have an aspect of 2:1
        title,
        titlesize,
        xticklabelsize,
        yticklabelsize,
        xticks = 0:60:360,      # label 0˚E, 60˚E, 120˚E, ...
        yticks = -60:30:60,     # label -60˚N, -30˚N, 0˚N, ...
        xtickformat = values -> ["$(round(Int, value))˚E" for value in values],
        ytickformat = values -> ["$(round(Int, value))˚N" for value in values],
        axis_kwargs...
    )
    hm = heatmap!(ax, field; kwargs...)
    return ax, hm
end

"""
$(TYPEDSIGNATURES)

Mutating variant of `heatmap` for RingGrids `Field`s that plots directly into an existing `Axis`.
Returns the plotting object created by `Makie.heatmap!`.
"""
function Makie.heatmap!(ax::Makie.AbstractAxis, field::RingGrids.AbstractFullField2D; kwargs...)
    mat = Matrix(field)                 # reshapes a full field into a matrix
    lond = RingGrids.get_lond(field)    # get lon, lat axes in degrees
    latd = RingGrids.get_latd(field)
    hm = heatmap!(ax, lond, latd, mat; kwargs...)
    return hm
end

end # module
