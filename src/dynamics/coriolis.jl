abstract type AbstractCoriolis <: AbstractModelComponent end

Base.@kwdef struct Coriolis{NF} <: AbstractCoriolis
    "number of latitude rings"
    nlat::Int

    "coriolis frequency [s^-1], scaled by radius as is vorticity = 2Ω*sin(lat)*radius"
    f::Vector{NF} = zeros(NF,nlat)
end

Coriolis(SG::SpectralGrid;kwargs...) = Coriolis{SG.NF}(nlat=SG.nlat;kwargs...)

function initialize!(coriolis::Coriolis, model::ModelSetup)
    (;rotation) = model.planet
    (;sinlat, radius) = model.geometry

    # =2Ωsin(lat) but scaled with radius as are the equations
    coriolis.f .= 2rotation * sinlat * radius
end

"""
$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on the grid `Grid` of resolution `nlat_half`
on a planet of `ratation` [1/s]. Default rotation of Earth."""
function coriolis(
    ::Type{Grid},
    nlat_half::Integer;
    rotation=DEFAULT_ROTATION
) where {Grid<:AbstractGrid}
    
    f = zeros(Grid,nlat_half)
    lat = get_lat(Grid,nlat_half)

    for (j,ring) in enumerate(eachring(f))
        fⱼ = 2rotation*sin(lat[j])
        for ij in ring
            f[ij] = fⱼ
        end
    end
    return f
end

"""
$(TYPEDSIGNATURES)
Return the Coriolis parameter `f` on a grid like `grid`
on a planet of `ratation` [1/s]. Default rotation of Earth."""
function coriolis(
    grid::Grid;
    rotation=DEFAULT_ROTATION
) where {Grid<:AbstractGrid}
    return coriolis(Grid,grid.nlat_half;rotation)
end