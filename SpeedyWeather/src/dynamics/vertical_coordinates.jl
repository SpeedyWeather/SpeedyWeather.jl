abstract type AbstractVerticalCoordinate <: AbstractModelComponent end

export SigmaCoordinates

"""Sigma coordinates for the vertical coordinates, defined by their half layers.
Sigma coordinates are currently hardcoded in Geometry.
$(TYPEDFIELDS)"""
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
SigmaCoordinates(σ_half::AbstractVector) = SigmaCoordinates(length(σ_half) - 1, σ_half)

# constructor using default nlayers if nothing provided
SigmaCoordinates() = SigmaCoordinates(DEFAULT_NLAYERS)

function Base.show(io::IO, σ::SigmaCoordinates)
    println(io, "$(σ.nlayers)-layer $(typeof(σ))")
    nchars = length(string(σ.nlayers))
    σ_half = on_architecture(CPU(), σ.σ_half)
    format = Printf.Format("%$(nchars)d")
    for k in 1:σ.nlayers
        println(io, "├─ ", @sprintf("%1.4f", σ_half[k]), "  k = ", Printf.format(format, k - 1), ".5")
        σk = (σ_half[k] + σ_half[k + 1]) / 2
        println(io, "│× ", @sprintf("%1.4f", σk), "  k = ", Printf.format(format, k))
    end
    return print(io, "└─ ", @sprintf("%1.4f", σ_half[end]), "  k = ", Printf.format(format, σ.nlayers), ".5")
end

"""$(TYPEDSIGNATURES)
Vertical sigma coordinates defined by their nlayers+1 half levels `σ_levels_half`. Sigma coordinates are
fraction of surface pressure (p/p0) and are sorted from top (stratosphere) to bottom (surface).
The first half level is at 0 the last at 1. Default levels are equally spaced from 0 to 1 (including)."""
function default_sigma_coordinates(nlayers::Integer, profile = z -> z)

    # equi-distant in σ
    σ_half = collect(range(0, 1, nlayers + 1))

    # Frierson 2006 spacing, higher resolution in surface boundary layer and in stratosphere
    # z = collect(range(1, 0, nlayers+1))
    # σ_half = @. exp(-5*(0.05*z + 0.95*z^3))
    # σ_half[1] = 0
    return σ_half
end

"""$(TYPEDSIGNATURES)
Check that nlayers and σ_half match."""
function sigma_okay(nlayers::Integer, σ_half::AbstractVector)
    @assert σ_half[1] >= 0 "First manually specified σ_half has to be >0"
    @assert σ_half[end] == 1 "Last manually specified σ_half has to be 1."
    @assert nlayers == (length(σ_half) - 1) "nlayers has to be length of σ_half - 1"
    @assert Utils.isincreasing(σ_half) "Vertical sigma coordinates are not increasing."
    return true
end

export SigmaPressureCoordinates
struct SigmaPressureCoordinates{NF, IntType, VectorType, F} <: AbstractVerticalCoordinate
    reference_pressure::NF
    A_half::VectorType
    B_half::VectorType
    A_full::VectorType
    B_full::VectorType
end

function SigmaPressureCoordinates(
    spectral_grid::SpectralGrid;
    reference_pressure = 1e5,
    σ_half = range(0, 1, length=spectral_grid.nlayers+1),
    transition = σ -> σ
)
    # sigma coordinates
    σ_half = on_architecture(spectral_grid.architecture, σ_half)
    σ_full = 0.5 * (σ_half[2:end] + σ_half[1:(end - 1)])

    # hybrid coordinates
    A_half = @. σ_half * (1 - transition(σ_half))
    B_half = @. σ_half * transition(σ_half)
    A_full = @. σ_full * (1 - transition(σ_full))
    B_full = @. σ_full * transition(σ_full)

    p_ref = convert(spectral_grid.NF, reference_pressure)
    return SigmaPressureCoordinates(p_ref, A_half, B_half, A_full, B_full)
end

# function SigmaPressureCoordinateA(σ, SP::SigmaPressureCoordinates)
#     γ = SP.transition
#     return σ * (1 - γ(σ)) * SP.reference_pressure + σ

function pressure(surface_pressure::Number, k::Integer, coordinate::SigmaPressureCoordinates)
    A = coordinate.A_full
    B = coordinate.B_full
    p_ref = coordinate.reference_pressure
    return A[k]*p_ref + B[k]*surface_pressure
end

function pressure_thickness(surface_pressure::Number, k::Integer, coordinate::SigmaPressureCoordinates)
    A = coordinate.A_half
    B = coordinate.B_half
    p_ref = coordinate.reference_pressure
    return (A[k+1] - A[k])*p_ref + (B[k+1] - B[k])*surface_pressure
end