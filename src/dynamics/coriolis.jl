abstract type AbstractCoriolis <: AbstractModelComponent end

export Coriolis

"""Model component holding the Coriolis parameter `f` on the model grid,
a vector over the latitude rings. $(TYPEDFIELDS)"""
@kwdef struct Coriolis{VectorType} <: AbstractCoriolis
    "[DERIVED] Coriolis frequency [s^-1], scaled by radius as is vorticity = 2Ω*sin(lat)*radius"
    f::VectorType
end

Adapt.@adapt_structure Coriolis

# generator
Coriolis(SG::SpectralGrid) = Coriolis(on_architecture(SG.architecture, zeros(SG.NF, SG.nlat)))

function initialize!(coriolis::Coriolis, model::AbstractModel)
    (; radius, rotation) = model.planet
    (; sinlat) = model.geometry

    # =2Ωsin(lat) but scaled with radius as are the equations
    coriolis.f .= 2rotation * sinlat * radius
    return nothing
end

export coriolis

"""$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on the grid `Grid` of resolution `nlat_half`
on a planet of `rotation` [1/s]. Default rotation of Earth."""
function coriolis!(f::AbstractField; rotation = DEFAULT_ROTATION)
    lat = on_architecture(f, get_lat(f))     # in radians [-π/2, π/2]

    arch = architecture(f)
    (; whichring) = f.grid
    launch!(arch, RingGridWorkOrder, size(f), coriolis_kernel!, f, lat, rotation, whichring)
    return f
end

@kernel inbounds = true function coriolis_kernel!(f, lat, rotation, whichring)
    ijk = @index(Global, Cartesian)     # for 2D, 3D, ... simultaneously
    j = whichring[ijk[1]]               # latitude ring index
    f[ijk] = 2rotation * sin(lat[j])
end

"""
$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on the same grid as `field`
on a planet of kwarg `rotation` [1/s]. Default rotation of Earth."""
coriolis(field::AbstractField; kwargs...) = coriolis!(similar(field); kwargs...)
coriolis(grid::AbstractGrid, ks...; kwargs...) = coriolis!(zero(grid, ks...); kwargs...)
coriolis(::Type{Grid}, nlat_half::Integer, args...; kwargs...) where {Grid <: AbstractGrid} =
    coriolis(Grid(nlat_half), args...; kwargs...)
