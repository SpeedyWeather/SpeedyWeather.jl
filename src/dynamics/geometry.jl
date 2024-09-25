abstract type AbstractGeometry end
export Geometry

"""
$(TYPEDSIGNATURES)
Construct Geometry struct containing parameters and arrays describing an iso-latitude grid <:AbstractGrid
and the vertical levels. Pass on `SpectralGrid` to calculate the following fields
$(TYPEDFIELDS)
"""
@kwdef struct Geometry{NF<:AbstractFloat} <: AbstractGeometry

    "SpectralGrid that defines spectral and grid resolution"
    spectral_grid::SpectralGrid

    "grid of the dynamical core"
    Grid::Type{<:AbstractGrid} = spectral_grid.Grid

    "resolution parameter nlat_half of Grid, # of latitudes on one hemisphere (incl Equator)"
    nlat_half::Int = spectral_grid.nlat_half      


    # GRID-POINT SPACE
    "maximum number of longitudes (at/around Equator)"
    nlon_max::Int = get_nlon_max(Grid, nlat_half)

    "=nlon_max, same (used for compatibility), TODO: still needed?"
    nlon::Int = nlon_max

    "number of latitude rings"
    nlat::Int = spectral_grid.nlat
    
    "number of vertical levels"
    nlayers::Int = spectral_grid.nlayers

    "total number of grid points"
    npoints::Int = spectral_grid.npoints

    "Planet's radius [m]"
    radius::NF = spectral_grid.radius    


    # ARRAYS OF LANGITUDES/LONGITUDES
    "array of colatitudes in radians (0...π)"
    colat::Vector{Float64} = get_colat(Grid, nlat_half)

    "array of latitudes in radians (π...-π)"
    lat::Vector{NF} = get_lat(Grid, nlat_half)

    "array of latitudes in degrees (90˚...-90˚)"
    latd::Vector{Float64} = get_latd(Grid, nlat_half)

    "array of longitudes in degrees (0...360˚), empty for non-full grids"
    lond::Vector{Float64} = get_lond(Grid, nlat_half)

    "longitude (0˚...360˚) for each grid point in ring order"
    londs::Vector{NF} = get_latdlonds(Grid, nlat_half)[2]
    
    "latitude (-90˚...˚90) for each grid point in ring order"
    latds::Vector{NF} = get_latdlonds(Grid, nlat_half)[1]

    "longitude (0...2π) for each grid point in ring order"
    lons::Vector{NF} = RingGrids.get_latlons(Grid, nlat_half)[2]
    
    "latitude (-π/2...π/2) for each grid point in ring order"
    lats::Vector{NF} = RingGrids.get_latlons(Grid, nlat_half)[1]

    "sin of latitudes"
    sinlat::Vector{NF} = sind.(latd)
    
    "cos of latitudes"
    coslat::Vector{NF} = cosd.(latd)
    
    "= 1/cos(lat)"
    coslat⁻¹::Vector{NF} = 1 ./ coslat

    "= cos²(lat)"
    coslat²::Vector{NF} = coslat.^2

    "# = 1/cos²(lat)"
    coslat⁻²::Vector{NF} = 1 ./ coslat²            

    # VERTICAL SIGMA COORDINATE σ = p/p0 (fraction of surface pressure)
    "σ at half levels, σ_k+1/2"
    σ_levels_half::Vector{NF} = spectral_grid.vertical_coordinates.σ_half

    "σ at full levels, σₖ"
    σ_levels_full::Vector{NF} = 0.5*(σ_levels_half[2:end] + σ_levels_half[1:end-1])  
    
    "σ level thicknesses, σₖ₊₁ - σₖ"
    σ_levels_thick::Vector{NF} = σ_levels_half[2:end] - σ_levels_half[1:end-1]      

    "log of σ at full levels, include surface (σ=1) as last element"
    ln_σ_levels_full::Vector{NF} = log.(vcat(σ_levels_full, 1))

    "Full to half levels interpolation"
    full_to_half_interpolation::Vector{NF} = σ_interpolation_weights(σ_levels_full, σ_levels_half)
end

"""
$(TYPEDSIGNATURES)
Generator function for `Geometry` struct based on `spectral_grid`."""
function Geometry(spectral_grid::SpectralGrid)
    error_message = "nlayers=$(spectral_grid.nlayers) does not match length nlayers="*
        "$(spectral_grid.vertical_coordinates.nlayers) in spectral_grid.vertical_coordinates."
    @assert spectral_grid.nlayers == spectral_grid.vertical_coordinates.nlayers error_message
    return Geometry{spectral_grid.NF}(; spectral_grid)
end

function Base.show(io::IO, G::Geometry)
    print(io, "$(typeof(G)) for $(G.spectral_grid)")
end

"""
$(TYPEDSIGNATURES)
Interpolation weights for full to half level interpolation
on sigma coordinates. Following Fortran SPEEDY documentation eq. (1)."""
function σ_interpolation_weights(
    σ_levels_full::AbstractVector,
    σ_levels_half::AbstractVector)

    weights = zero(σ_levels_full)
    nlayers = length(weights)
    nlayers == 1 && return weights     # escape early for 1 layer to avoid out-of-bounds access

    for k in 1:nlayers-1
        weights[k] = (log(σ_levels_half[k+1]) - log(σ_levels_full[k])) /
                        (log(σ_levels_full[k+1]) - log(σ_levels_full[k]))
    end
    # was log(0.99) in Fortran SPEEDY code but doesn't make sense to me
    weights[end] =  (log(σ_levels_half[nlayers+1]) - log(σ_levels_full[nlayers])) /
                        (log(σ_levels_full[nlayers]) - log(σ_levels_full[nlayers-1]))
    
    return weights
end


"""
$(TYPEDSIGNATURES)
Given a vector in column defined at full levels, do a linear interpolation in
log(σ) to calculate its values at half-levels, skipping top (k=1/2), extrapolating to bottom (k=nlayers+1/2).
"""
function vertical_interpolate!(
    A_half::Vector,             # quantity A on half levels (excl top)
    A_full::Vector,             # quantity A on full levels
    G::Geometry,
)
    nlayers = length(A_half)
    weights = G.full_to_half_interpolation

    # full levels contain one more for surface
    # TODO this is currently confusing because the surface fluxes use full[end]
    # as surface value which is technically on half levels though!
    @boundscheck nlayers <= length(A_full) || throw(BoundsError)
    @boundscheck nlayers <= length(weights) || throw(BoundsError)

    # For A at each full level k, compute A at the half-level below, i.e. at the boundary
    # between the full levels k and k+1. Fortran SPEEDY documentation eq. (1)
    for k = 1:nlayers-1
        A_half[k] = A_full[k] + weights[k]*(A_full[k+1] - A_full[k])
    end

    # Compute the values at the surface separately
    A_half[nlayers] = A_full[nlayers] + weights[nlayers]*(A_full[nlayers] - A_full[nlayers-1])

    return nothing
end