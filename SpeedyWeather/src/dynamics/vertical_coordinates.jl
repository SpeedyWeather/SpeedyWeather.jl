"""Abstract supertype for vertical coordinate systems. Subtypes define how model levels
are mapped to pressure given a surface pressure. Must implement `pressure`,
`pressure_half`, `pressure_thickness`, and `sigma` (see each for their signatures)."""
abstract type AbstractVerticalCoordinates <: AbstractModelComponent end

# extend for broadcasting like pressure.(k, field, vertical_coordinate)
Base.broadcastable(VC::AbstractVerticalCoordinates) = Ref(VC)

"""$(TYPEDSIGNATURES)
Pressure [Pa] at the lower interface of full level `k` (half level k+½), given `surface_pressure` [Pa].
Equivalent to `pressure_half(k+1, surface_pressure, coordinate)`."""
@inline pressure_below(k::Integer, surface_pressure::Number, coordinate::AbstractVerticalCoordinates) =
    pressure_half(k+1, surface_pressure, coordinate)

"""$(TYPEDSIGNATURES)
Pressure [Pa] at the upper interface of full level `k` (half level k-½), given `surface_pressure` [Pa].
Equivalent to `pressure_half(k, surface_pressure, coordinate)`."""
@inline pressure_above(k::Integer, surface_pressure::Number, coordinate::AbstractVerticalCoordinates) =
    pressure_half(k, surface_pressure, coordinate)

export SigmaCoordinates

"""Terrain-following sigma vertical coordinates. Pressure at full level `k` is
`p_k = σ_k * surface_pressure`, where `σ ∈ (0, 1]` is the fraction of surface pressure.
Levels are defined by `σ_half` (half levels at layer interfaces, including `σ = 0` at
the top and `σ = 1` at the surface); full levels sit at the midpoints.
$(TYPEDFIELDS)"""
struct SigmaCoordinates{IntType, VectorType} <: AbstractVerticalCoordinates
    "[DERIVED] Number of vertical layers"
    nlayers::IntType

    "[DERIVED] σ at half levels (interfaces); k = 1 is top of atmosphere (σ = 0), k = nlayers+1 is surface (σ = 1)"
    σ_half::VectorType

    "[DERIVED] σ at full levels (layer midpoints); σ_full[k] = (σ_half[k] + σ_half[k+1]) / 2"
    σ_full::VectorType

    "[DERIVED] Layer thickness in sigma; σ_thickness[k] = σ_half[k+1] - σ_half[k], sums to 1"
    σ_thickness::VectorType

    SigmaCoordinates{T, V}(nlayers, σh, σf, Δσ) where {T, V} = sigma_okay(nlayers, σh) ?
        new{T, V}(nlayers, σh, σf, Δσ) : error("σ_half = $σh cannot be used for $nlayers-level SigmaCoordinates")
end

Adapt.@adapt_structure SigmaCoordinates

function SigmaCoordinates(SG::SpectralGrid, σ_half::AbstractVector)
    (; nlayers) = SG
    σ_half = on_architecture(SG.architecture, convert.(SG.NF, σ_half))
    σ_full = (σ_half[2:end] + σ_half[1:(end - 1)]) / 2
    σ_thickness = σ_half[2:end] - σ_half[1:(end - 1)]
    return SigmaCoordinates{typeof(nlayers), typeof(σ_half)}(nlayers, σ_half, σ_full, σ_thickness)
end

# other constructors for convenience
SigmaCoordinates(SG::SpectralGrid, profile::Function) = SigmaCoordinates(SG, sigma_half_spacing(SG.nlayers, profile))
SigmaCoordinates(SG::SpectralGrid; σ_half = sigma_half_spacing(SG.nlayers)) = SigmaCoordinates(SG, σ_half)
SigmaCoordinates(σ_half::AbstractVector) = SigmaCoordinates(SpectralGrid(nlayers = length(σ_half) - 1), σ_half)
SigmaCoordinates() = SigmaCoordinates(SpectralGrid())

get_nlayers(σ::SigmaCoordinates) = σ.nlayers
get_σ_half(σ::SigmaCoordinates) = σ.σ_half
get_σ_full(σ::SigmaCoordinates) = σ.σ_full
get_σ_thickness(σ::SigmaCoordinates) = σ.σ_thickness

function Base.show(io::IO, σ::SigmaCoordinates)
    params = "{$(typeof(σ.nlayers)), $(typeof(σ.σ_half))}"
    println(io, "$(σ.nlayers)-layer ", styled"{warning:SigmaCoordinates}", styled"{note:$params}")
    nchars = length(string(σ.nlayers))
    σ_half = on_architecture(CPU(), σ.σ_half)
    σ_full = on_architecture(CPU(), σ.σ_full)
    format = Printf.Format("%$(nchars)d")
    for k in 1:σ.nlayers
        println(io, " ├─ ", @sprintf("%1.4f", σ_half[k]), "  ", "k = ", Printf.format(format, k - 1), ".5")
        println(io, " │× ", @sprintf("%1.4f", σ_full[k]), "  ", "k = ", Printf.format(format, k))
    end
    print(io, " └─ ", @sprintf("%1.4f", σ_half[end]), "  ", "k = ", Printf.format(format, σ.nlayers), ".5")
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

"""$(TYPEDSIGNATURES)
Pressure [Pa] at full level `k` given `surface_pressure` [Pa] and sigma `coordinate`."""
@inline function pressure(k::Integer, surface_pressure::Number, coordinate::SigmaCoordinates)
    σ = coordinate.σ_full
    return σ[k] * surface_pressure
end

"""$(TYPEDSIGNATURES)
Pressure [Pa] at half level `k` given `surface_pressure` [Pa] and sigma `coordinate`.
Half levels are indexed 1 (top of atmosphere, σ = 0) to nlayers+1 (surface, σ = 1)."""
@inline function pressure_half(k::Integer, surface_pressure::Number, coordinate::SigmaCoordinates)
    σ = coordinate.σ_half
    return σ[k] * surface_pressure
end

"""$(TYPEDSIGNATURES)
Pressure thickness [Pa] of full level `k` given `surface_pressure` [Pa] and sigma `coordinate`."""
@inline function pressure_thickness(k::Integer, surface_pressure::Number, coordinate::SigmaCoordinates)
    Δσ = coordinate.σ_thickness
    return Δσ[k] * surface_pressure
end

"""$(TYPEDSIGNATURES)
Sigma coordinate (fraction of surface pressure) at full level `k` for sigma `coordinate`."""
@inline sigma(k::Integer, coordinate::SigmaCoordinates) = coordinate.σ_full[k]

export FriersonSigmaCoordinates

"""$(TYPEDSIGNATURES)
Frierson 2006 spacing, higher resolution in surface boundary layer and in stratosphere.
Following exp(-5 * (0.05 * (1 - σ) + 0.95 * (1 - σ)^3)), σ being equi-spaced.
Originally without the 1 - σ, but then vertical ordering is reversed."""
FriersonSigmaCoordinates(SG::SpectralGrid) =
    SigmaCoordinates(SG, frierson_profile)

FriersonSigmaCoordinates() = FriersonSigmaCoordinates(SpectralGrid())
frierson_profile(σ) = exp(-5 * (0.05 * (1 - σ) + 0.95 * (1 - σ)^3))

export SigmaPressureCoordinates

"""Hybrid sigma-pressure vertical coordinates. Pressure at full level `k` is
`p_k = A_k * reference_pressure + B_k * surface_pressure`, where `A_k + B_k = σ_k`
(the nominal sigma level). The A/B split is controlled by a `transition` function:
`transition(σ) = 0` gives pure pressure levels (terrain-independent, fixed in space),
`transition(σ) = 1` gives pure sigma levels (terrain-following). See
`CubicSigmaPressureCoordinates` for a ready-to-use constructor with a smooth transition.
$(TYPEDFIELDS)"""
struct SigmaPressureCoordinates{NF, VectorType} <: AbstractVerticalCoordinates
    "[DERIVED] Reference pressure for the pressure component [Pa]"
    reference_pressure::NF

    "[DERIVED] Pressure component coefficient at half levels (k = 1 top, nlayers+1 surface)"
    A_half::VectorType

    "[DERIVED] Sigma component coefficient at half levels"
    B_half::VectorType

    "[DERIVED] Pressure component coefficient at full levels"
    A_full::VectorType

    "[DERIVED] Sigma component coefficient at full levels"
    B_full::VectorType

    "[DERIVED] Pressure component of layer thickness (ΔA = A_half[k+1] - A_half[k])"
    A_thickness::VectorType

    "[DERIVED] Sigma component of layer thickness (ΔB = B_half[k+1] - B_half[k])"
    B_thickness::VectorType
end

Adapt.@adapt_structure SigmaPressureCoordinates

# main constructor
function SigmaPressureCoordinates(
        spectral_grid::SpectralGrid,
        σ_half::AbstractVector = sigma_half_spacing(spectral_grid.nlayers);
        reference_pressure = 1.0e5,
        transition = σ -> σ
    )
    # sigma coordinates
    σ_half = on_architecture(spectral_grid.architecture, convert.(spectral_grid.NF, σ_half))
    # σ_full = (σ_half[2:end] + σ_half[1:(end - 1)]) / 2

    # hybrid coordinates defined via half layers
    B_half = @. σ_half * transition(σ_half)
    B_full = (B_half[2:end] + B_half[1:(end - 1)]) / 2
    
    # do not reevalute the (possibly nonlinear) transition for full layers
    # average instead to have layer centres always at mid-pressure too
    A_half = @. σ_half * (1 - transition(σ_half))
    A_full = (A_half[2:end] + A_half[1:(end - 1)]) / 2
    # A_full = maximum.(0, σ_full - B_full)     # to avoid -0
    # A_half = maximum.(0, σ_half - B_half)

    ΔA = A_half[2:end] - A_half[1:(end - 1)]
    ΔB = B_half[2:end] - B_half[1:(end - 1)]

    p_ref = convert(spectral_grid.NF, reference_pressure)
    return SigmaPressureCoordinates(p_ref, A_half, B_half, A_full, B_full, ΔA, ΔB)
end

# constructor that first evaluates the profile function to create a sigma half vector
function SigmaPressureCoordinates(
        spectral_grid::SpectralGrid,
        profile::Function;
        kwargs...
    )   
    # also allow for function to be passed on, evaluate here
    σ_half_vector = sigma_half_spacing(spectral_grid.nlayers, profile)
    return SigmaPressureCoordinates(spectral_grid, σ_half_vector; kwargs...)
end

SigmaPressureCoordinates() = SigmaPressureCoordinates(SpectralGrid())
SigmaPressureCoordinates(SG::SpectralGrid; σ_half = sigma_half_spacing(SG.nlayers), kwargs...) =
    SigmaPressureCoordinates(SG, σ_half; kwargs...)

function Base.show(io::IO, S::SigmaPressureCoordinates)
    nlayers = get_nlayers(S)
    params = "{$(typeof(S.reference_pressure)), $(typeof(S.A_half))}"
    println(io, "$nlayers-layer ", styled"{warning:SigmaPressureCoordinates}", styled"{note:$params}")
    println(io, "    ", styled"{info:p_ref}", " = $(S.reference_pressure) Pa   ░ pressure │ sigma █")
    nchars = length(string(nlayers))
    A_half = on_architecture(CPU(), S.A_half)
    B_half = on_architecture(CPU(), S.B_half)
    A_full = on_architecture(CPU(), S.A_full)
    B_full = on_architecture(CPU(), S.B_full)
    format = Printf.Format("%$(nchars)d")
    for k in 1:nlayers
        σ_half_k = A_half[k] + B_half[k]
        println(io, " ├─ ", @sprintf("%1.4f", σ_half_k), "  ", "k = ", Printf.format(format, k - 1), ".5")
        σ_full_k = A_full[k] + B_full[k]
        sigma_frac = σ_full_k > 0 ? B_full[k] / σ_full_k : 0.0
        nsigma = round(Int, 10 * sigma_frac)
        bar = "█"^nsigma * "░"^(10 - nsigma)
        println(
            io, " │× ", @sprintf("%1.4f", σ_full_k), "  ", "k = ", Printf.format(format, k),
            "   ", bar, "  ", styled"{note:A=}", @sprintf("%1.4f", A_full[k]), "  ", styled"{note:B=}", @sprintf("%1.4f", B_full[k])
        )
    end
    σ_half_end = A_half[end] + B_half[end]
    print(io, " └─ ", @sprintf("%1.4f", σ_half_end), "  ", "k = ", Printf.format(format, nlayers), ".5")
    return nothing
end

get_nlayers(S::SigmaPressureCoordinates) = length(S.B_full)
get_σ_half(σ::SigmaPressureCoordinates) = σ.B_half
get_σ_full(σ::SigmaPressureCoordinates) = σ.B_full
get_σ_thickness(σ::SigmaPressureCoordinates) = σ.B_thickness

"""$(TYPEDSIGNATURES)
Pressure [Pa] at full level `k` given `surface_pressure` [Pa] and hybrid sigma-pressure `coordinate`.
Computed as `A[k] * p_ref + B[k] * surface_pressure`, where A and B are the pressure and sigma
coefficients of the hybrid coordinate."""
@inline function pressure(k::Integer, surface_pressure::Number, coordinate::SigmaPressureCoordinates)
    A = coordinate.A_full
    B = coordinate.B_full
    p_ref = coordinate.reference_pressure
    return A[k] * p_ref + B[k] * surface_pressure
end
    
"""$(TYPEDSIGNATURES)
Pressure [Pa] at half level `k` given `surface_pressure` [Pa] and hybrid sigma-pressure `coordinate`.
Half levels are indexed 1 (top of atmosphere) to nlayers+1 (surface).
Computed as `A_half[k] * p_ref + B_half[k] * surface_pressure`."""
@inline function pressure_half(k::Integer, surface_pressure::Number, coordinate::SigmaPressureCoordinates)
    A = coordinate.A_half
    B = coordinate.B_half
    p_ref = coordinate.reference_pressure
    return A[k] * p_ref + B[k] * surface_pressure
end

"""$(TYPEDSIGNATURES)
Pressure thickness [Pa] of full level `k` given `surface_pressure` [Pa] and hybrid sigma-pressure
`coordinate`. Computed as `ΔA[k] * p_ref + ΔB[k] * surface_pressure`."""
@inline function pressure_thickness(k::Integer, surface_pressure::Number, coordinate::SigmaPressureCoordinates)
    ΔA = coordinate.A_thickness
    ΔB = coordinate.B_thickness
    p_ref = coordinate.reference_pressure
    return ΔA[k] * p_ref + ΔB[k] * surface_pressure
end

"""$(TYPEDSIGNATURES)
Sigma coordinate (fraction of surface pressure) at full level `k` for hybrid sigma-pressure
`coordinate`. Returns `A[k] + B[k]`, which equals the nominal sigma level regardless of the
pressure-sigma transition, and is independent of surface pressure."""
@inline sigma(k::Integer, coordinate::SigmaPressureCoordinates) = coordinate.A_full[k] + coordinate.B_full[k]

export CubicSigmaPressureCoordinates

"""$(TYPEDSIGNATURES)
Transition function for `CubicSigmaPressureCoordinates`. Returns 0 (pure pressure) for
σ ≤ `pressure_only_above`, 1 (pure sigma) for σ ≥ `σ_only_below`, and a cubic smoothstep in between, giving C¹
continuity at both thresholds. The smoothstep coefficients 3 and 2 are fixed by the C¹
boundary conditions on the normalised variable t ∈ [0, 1] and do not depend on the thresholds."""
function cubic_transition(σ; pressure_only_above = 0.2, σ_only_below = 0.8)
    pressure_only_above = oftype(σ, pressure_only_above)
    σ_only_below = oftype(σ, σ_only_below)
    # the actual cubic transition
    # values 3, 2 are fixed by f(0)=0, f(1)=1, f'(0)=0, f'(1)=0 on t ∈ [0,1]
    t = min(one(σ), max(zero(σ), (σ - pressure_only_above) / (σ_only_below - pressure_only_above)))
    return t * t * (3 - 2t)   
end

"""$(TYPEDSIGNATURES)
Sigma-pressure coordinate with a cubic smoothstep transition: pure pressure levels near the
top of the atmosphere (σ ≤ `pressure_only_above`), pure sigma levels near the surface (σ ≥ `σ_only_below`), and
a cubic smoothstep transition in between. The thresholds and transition shape are a pragmatic
choice and do not follow any specific operational model's formulation."""
CubicSigmaPressureCoordinates(SG::SpectralGrid, args...; pressure_only_above = 0.2, σ_only_below = 0.8, kwargs...) =
    SigmaPressureCoordinates(SG, args...; transition = σ -> cubic_transition(σ; pressure_only_above, σ_only_below), kwargs...)

CubicSigmaPressureCoordinates(args...; kwargs...) = CubicSigmaPressureCoordinates(SpectralGrid(), args...; kwargs...)