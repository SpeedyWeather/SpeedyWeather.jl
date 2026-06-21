abstract type AbstractVerticalCoordinate <: AbstractModelComponent end

export SigmaCoordinates

"""Sigma, i.e. fraction of surface pressure, vertical coordinates defined by their half layers.
First half layer has to be 0 (top of the atmosphere), last has to be 1 (surface).
$(TYPEDFIELDS)"""
@kwdef struct SigmaCoordinates{IntType, VectorType} <: AbstractVerticalCoordinate
    nlayers::IntType
    σ_half::VectorType
    σ_full::VectorType
    σ_thickness::VectorType

    SigmaCoordinates{T, V}(nlayers, σh, σf, Δσ) where {T, V} = sigma_okay(nlayers, σh) ?
        new{T, V}(nlayers, σh, σf, Δσ) : error("σ_half = $σh cannot be used for $nlayers-level SigmaCoordinates")
end

Adapt.@adapt_structure SigmaCoordinates

function SigmaCoordinates(SG::SpectralGrid, σ_half::AbstractVector = sigma_half_spacing(SG.nlayers))
    (; nlayers) = SG
    σ_half = on_architecture(SG.architecture, convert.(SG.NF, σ_half))
    σ_full = (σ_half[2:end] + σ_half[1:(end - 1)]) / 2
    σ_thickness = σ_half[2:end] - σ_half[1:(end - 1)]
    return SigmaCoordinates{typeof(nlayers), typeof(σ_half)}(; nlayers, σ_half, σ_full, σ_thickness)
end

# other constructors for convenience
SigmaCoordinates(σ_half::AbstractVector) = SigmaCoordinates(SpectralGrid(nlayers = length(σ_half) - 1), σ_half)
SigmaCoordinates() = SigmaCoordinates(SpectralGrid())

get_nlayers(σ::SigmaCoordinates) = σ.nlayers
get_σ_half(σ::SigmaCoordinates) = σ.σ_half

function Base.show(io::IO, σ::SigmaCoordinates)
    println(io, "$(σ.nlayers)-layer $(typeof(σ))")
    nchars = length(string(σ.nlayers))
    σ_half = on_architecture(CPU(), σ.σ_half)
    σ_full = on_architecture(CPU(), σ.σ_full)
    format = Printf.Format("%$(nchars)d")
    for k in 1:σ.nlayers
        println(io, "├─ ", @sprintf("%1.4f", σ_half[k]), "  k = ", Printf.format(format, k - 1), ".5")
        println(io, "│× ", @sprintf("%1.4f", σ_full[k]), "  k = ", Printf.format(format, k))
    end
    print(io, "└─ ", @sprintf("%1.4f", σ_half[end]), "  k = ", Printf.format(format, σ.nlayers), ".5")
    return nothing
end

"""$(TYPEDSIGNATURES)
Vertical sigma coordinates defined by their half levels. Sigma coordinates are
fraction of surface pressure (p/p0) and are sorted from top (stratosphere) to bottom (surface).
The first half level is at 0 the last at 1. Default `profile` is equally spaced from 0 to 1 (including)."""
function sigma_half_spacing(nlayers::Integer, profile = z -> z)
    σ_half = profile.(collect(range(0, 1, nlayers + 1)))
    σ_half[1] = 0
    return σ_half
end

"""$(TYPEDSIGNATURES)
Check that nlayers and σ_half match."""
function sigma_okay(nlayers::Integer, σ_half::AbstractVector)
    @assert σ_half[1] >= 0 "First specified σ_half has to be >=0, $(σ_half[1]) given."
    @assert σ_half[end] == 1 "Last specified σ_half has to be 1, $(σ_half[end]) given."
    @assert nlayers == (length(σ_half) - 1) "nlayers has to be length of σ_half - 1, $nlayers vs $(length(σ_half) - 1) given."
    @assert Utils.isincreasing(σ_half) "Vertical sigma coordinates are not increasing."
    return true
end

@inline function pressure(k::Integer, surface_pressure::Number, coordinate::SigmaCoordinates)
    σ = coordinate.σ_full
    return σ[k] * surface_pressure
end

@inline function pressure_thickness(k::Integer, surface_pressure::Number, coordinate::SigmaCoordinates)
    Δσ = coordinate.σ_thickness
    return Δσ[k] * surface_pressure
end

export FriersonSigmaCoordinates

"""$(TYPEDSIGNATURES)
Frierson 2006 spacing, higher resolution in surface boundary layer and in stratosphere.
Following exp(-5 * (0.05 * (1 - σ) + 0.95 * (1 - σ)^3)), σ being equi-spaced.
Originally without the 1 - σ, but then vertical ordering is reversed."""
FriersonSigmaCoordinates(SG::SpectralGrid) =
    SigmaCoordinates(SG, sigma_half_spacing(SG.nlayers, frierson_profile))

frierson_profile(σ) = exp(-5 * (0.05 * (1 - σ) + 0.95 * (1 - σ)^3))

export SigmaPressureCoordinates
struct SigmaPressureCoordinates{NF, VectorType} <: AbstractVerticalCoordinate
    reference_pressure::NF
    A_half::VectorType
    B_half::VectorType
    A_full::VectorType
    B_full::VectorType
    A_thickness::VectorType
    B_thickness::VectorType
end

Adapt.@adapt_structure SigmaPressureCoordinates

function SigmaPressureCoordinates(
        spectral_grid::SpectralGrid;
        reference_pressure = 1.0e5,
        σ_half = sigma_half_spacing(spectral_grid.nlayers),
        transition = σ -> σ
    )
    # sigma coordinates
    σ_half = on_architecture(spectral_grid.architecture, convert.(spectral_grid.NF, σ_half))
    σ_full = (σ_half[2:end] + σ_half[1:(end - 1)]) / 2

    # hybrid coordinates
    A_half = @. σ_half * (1 - transition(σ_half))
    B_half = @. σ_half * transition(σ_half)
    A_full = @. σ_full * (1 - transition(σ_full))
    B_full = @. σ_full * transition(σ_full)

    ΔA = A_half[2:end] + A_half[1:(end - 1)]
    ΔB = B_half[2:end] + B_half[1:(end - 1)]

    p_ref = convert(spectral_grid.NF, reference_pressure)
    return SigmaPressureCoordinates(p_ref, A_half, B_half, A_full, B_full, ΔA, ΔB)
end

get_nlayers(S::SigmaPressureCoordinates) = length(S.B_full)
get_σ_half(σ::SigmaPressureCoordinates) = σ.B_half

@inline function pressure(k::Integer, surface_pressure::Number, coordinate::SigmaPressureCoordinates)
    A = coordinate.A_full
    B = coordinate.B_full
    p_ref = coordinate.reference_pressure
    return A[k] * p_ref + B[k] * surface_pressure
end

@inline function pressure_thickness(k::Integer, surface_pressure::Number, coordinate::SigmaPressureCoordinates)
    ΔA = coordinate.A_thickness
    ΔB = coordinate.B_thickness
    p_ref = coordinate.reference_pressure
    return ΔA[k] * p_ref + ΔB[k] * surface_pressure
end