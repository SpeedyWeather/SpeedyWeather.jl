abstract type AbstractVerticalCoordinate end

export SigmaCoordinates

"""Sigma coordinates for the vertical coordinates, defined by their half layers.
Sigma coordinates are currently hardcoded in Geometry.
$(TYPEDFIELDS)."""
struct SigmaCoordinates{NF, VectorType} <: AbstractVerticalCoordinate
    nlayers::Int
    σ_half::VectorType

    SigmaCoordinates{T, V}(nlayers::Integer, σ_half::AbstractVector) where {T, V} = sigma_okay(nlayers, σ_half) ?
    new{T, V}(nlayers, σ_half) : error("σ_half = $σ_half cannot be used for $nlayers-level SigmaCoordinates")
end

# constructor using default sigma coordinates if only nlayers provided, also collect to allow for AbstractRange
SigmaCoordinates(nlayers::Integer, σ_half::AbstractVector = default_sigma_coordinates(nlayers)) =
    SigmaCoordinates{eltype(σ_half), typeof(collect(σ_half))}(nlayers, collect(σ_half))

# constructor obtaining nlayers from σ_half
SigmaCoordinates(σ_half::AbstractVector) = SigmaCoordinates(length(σ_half)-1, σ_half)

# constructor using default nlayers if nothing provided
SigmaCoordinates() = SigmaCoordinates(DEFAULT_NLAYERS)


function Base.show(io::IO, σ::SigmaCoordinates)
    println(io, "$(σ.nlayers)-layer $(typeof(σ))")
    nchars = length(string(σ.nlayers))
    format = Printf.Format("%$(nchars)d")
    for k=1:σ.nlayers
        println(io, "├─ ", @sprintf("%1.4f", σ.σ_half[k]), "  k = ", Printf.format(format, k-1), ".5")
        σk = (σ.σ_half[k] + σ.σ_half[k+1])/2
        println(io, "│× ", @sprintf("%1.4f", σk), "  k = ", Printf.format(format, k))
    end
    print(io, "└─ ", @sprintf("%1.4f", σ.σ_half[end]), "  k = ", Printf.format(format, σ.nlayers), ".5")
end

"""
$(TYPEDSIGNATURES)
Vertical sigma coordinates defined by their nlayers+1 half levels `σ_levels_half`. Sigma coordinates are
fraction of surface pressure (p/p0) and are sorted from top (stratosphere) to bottom (surface).
The first half level is at 0 the last at 1. Default levels are equally spaced from 0 to 1 (including)."""
function default_sigma_coordinates(nlayers::Integer)

    # equi-distant in σ
    σ_half = collect(range(0, 1, nlayers+1))

    # Frierson 2006 spacing, higher resolution in surface boundary layer and in stratosphere
    # z = collect(range(1, 0, nlayers+1))
    # σ_half = @. exp(-5*(0.05*z + 0.95*z^3))
    # σ_half[1] = 0
    return σ_half
end

"""
$(TYPEDSIGNATURES)
Check that nlayers and σ_half match."""    
function sigma_okay(nlayers::Integer, σ_half::AbstractVector)
    @assert σ_half[1] >= 0 "First manually specified σ_half has to be >0"
    @assert σ_half[end] == 1 "Last manually specified σ_half has to be 1."
    @assert nlayers == (length(σ_half) - 1) "nlayers has to be length of σ_half - 1"
    @assert isincreasing(σ_half) "Vertical sigma coordinates are not increasing."
    return true
end
