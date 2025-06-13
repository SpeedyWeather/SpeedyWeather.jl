abstract type AbstractCoriolis <: AbstractModelComponent end

export Coriolis
@kwdef struct Coriolis{NF, VectorType} <: AbstractCoriolis
    "number of latitude rings"
    nlat::Int

    "coriolis frequency [s^-1], scaled by radius as is vorticity = 2Ω*sin(lat)*radius"
    f::VectorType = zeros(NF, nlat)
end

Coriolis(SG::SpectralGrid; kwargs...) = Coriolis{SG.NF, SG.VectorType}(nlat=SG.nlat; kwargs...)

function initialize!(coriolis::Coriolis, model::AbstractModel)
    (; rotation) = model.planet
    (; sinlat, radius) = model.geometry

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
    lat = get_lat(f)                        # in radians [-π/2, π/2]

    for (j, ring) in enumerate(eachring(f))
        fⱼ = 2rotation*sin(lat[j])
        for ij in ring
            f[ij, :] .= fⱼ                  # setindex across all ks dimensions
        end
    end
    return f
end

"""
$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on the same grid as `field`
on a planet of kwarg `rotation` [1/s]. Default rotation of Earth."""
coriolis(field::AbstractField; kwargs...) = coriolis!(similar(field); kwargs...)
coriolis(grid::AbstractGrid, ks...; kwargs...) = coriolis!(zero(grid, ks...); kwargs...)
coriolis(::Type{Grid}, nlat_half::Integer, args...; kwargs...) where {Grid<:AbstractGrid} =
    coriolis(Grid(nlat_half), args...; kwargs...)