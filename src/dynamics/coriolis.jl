abstract type AbstractCoriolis <: AbstractModelComponent end

export Coriolis
@kwdef struct Coriolis{VectorType} <: AbstractCoriolis
    "[DERIVED] Coriolis frequency [s^-1], scaled by radius as is vorticity = 2Ω*sin(lat)*radius"
    f::VectorType
end

Adapt.@adapt_structure Coriolis

# generator
Coriolis(SG::SpectralGrid) = Coriolis(on_architecture(SG.architecture, zeros(SG.NF, SG.nlat)))

# function barrier to unpack only model components needed
function initialize!(coriolis::Coriolis, model::PrimitiveEquation)
    model_parameters = (planet=model.planet, geometry=model.geometry)
    initialize!(coriolis, model_parameters)
end

function initialize!(coriolis::Coriolis, model)
    (; radius, rotation) = model.planet
    (; sinlat) = model.geometry

    # =2Ωsin(lat) but scaled with radius as are the equations
    coriolis.f .= 2rotation * sinlat * radius
    return nothing
end

export coriolis

"""
$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on the grid `Grid` of resolution `nlat_half`
on a planet of `rotation` [1/s]. Default rotation of Earth."""
function coriolis!(f::AbstractField; rotation = DEFAULT_ROTATION)                  
    lat = on_architecture(f, get_lat(f))     # in radians [-π/2, π/2]

    arch = architecture(f)
    (; whichring) = field.grid
    launch!(arch, RingGridWorkOrder, size(f), coriolis_kernel!, f, lat, rotation, whichring)
end

@kernel inbounds=true function coriolis_kernel!(f, lat, rotation, whichring)
    ij, k = @index(Global, NTuple)
    j = whichring[ij]
    f[ij, k] = 2rotation*sin(lat[j])
end

"""
$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on the same grid as `field`
on a planet of kwarg `rotation` [1/s]. Default rotation of Earth."""
coriolis(field::AbstractField; kwargs...) = coriolis!(similar(field); kwargs...)
coriolis(grid::AbstractGrid, ks...; kwargs...) = coriolis!(zero(grid, ks...); kwargs...)
coriolis(::Type{Grid}, nlat_half::Integer, args...; kwargs...) where {Grid<:AbstractGrid} =
    coriolis(Grid(nlat_half), args...; kwargs...)