abstract type VerticalCoordinates end

Base.@kwdef struct NoVerticalCoordinates <: VerticalCoordinates
    nlev::Int = 1
end

export SigmaCoordinates
Base.@kwdef struct SigmaCoordinates <: VerticalCoordinates
    nlev::Int = 8
    σ_half::Vector{Float64} = default_sigma_coordinates(nlev)

    SigmaCoordinates(nlev::Integer, σ_half::AbstractVector) = sigma_okay(nlev, σ_half) ?
    new(nlev, σ_half) : error("σ_half = $σ_half cannot be used for $nlev-level SigmaCoordinates")
end

# obtain nlev from length of predefined σ_half levels
SigmaCoordinates(σ_half::AbstractVector) = SigmaCoordinates(nlev=length(σ_half)-1;σ_half) 
SigmaCoordinates(σ_half::AbstractRange) = SigmaCoordinates(collect(σ_half))

function Base.show(io::IO, σ::SigmaCoordinates)
    println(io, "$(σ.nlev)-level SigmaCoordinates")
    nchars = length(string(σ.nlev))
    format = Printf.Format("%$(nchars)d")
    for k=1:σ.nlev
        println(io, "├─ ", @sprintf("%1.4f", σ.σ_half[k]), "  k = ", Printf.format(format, k-1), ".5")
        σk = (σ.σ_half[k] + σ.σ_half[k+1])/2
        println(io, "│× ", @sprintf("%1.4f", σk), "  k = ", Printf.format(format, k))
    end
    print(io, "└─ ", @sprintf("%1.4f", σ.σ_half[end]), "  k = ", Printf.format(format, σ.nlev), ".5")
end

"""
$(TYPEDSIGNATURES)
Vertical sigma coordinates defined by their nlev+1 half levels `σ_levels_half`. Sigma coordinates are
fraction of surface pressure (p/p0) and are sorted from top (stratosphere) to bottom (surface).
The first half level is at 0 the last at 1. Default levels are equally spaced from 0 to 1 (including)."""
function default_sigma_coordinates(nlev::Integer)

    # equi-distant in σ
    σ_half = collect(range(0, 1, nlev+1))

    # Frierson 2006 spacing, higher resolution in surface boundary layer and in stratosphere
    # z = collect(range(1, 0, nlev+1))
    # σ_half = @. exp(-5*(0.05*z + 0.95*z^3))
    # σ_half[1] = 0
    return σ_half
end

"""
$(TYPEDSIGNATURES)
Check that nlev and σ_half match."""    
function sigma_okay(nlev::Integer, σ_half::AbstractVector)
    @assert σ_half[1] >= 0 "First manually specified σ_half has to be >0"
    @assert σ_half[end] == 1 "Last manually specified σ_half has to be 1."
    @assert nlev == (length(σ_half) - 1) "nlev has to be length of σ_half - 1"
    @assert isincreasing(σ_half) "Vertical sigma coordinates are not increasing."
    return true
end

# currently not used as the grid has to be defined first
default_vertical_coordinates(::Type{<:Barotropic}) = NoVerticalCoordinates
default_vertical_coordinates(::Type{<:ShallowWater}) = NoVerticalCoordinates
default_vertical_coordinates(::Type{<:PrimitiveEquation}) = SigmaCoordinates

default_nlev(::Type{<:Barotropic}) = 1
default_nlev(::Type{<:ShallowWater}) = 1
default_nlev(::Type{<:PrimitiveEquation}) = 8
