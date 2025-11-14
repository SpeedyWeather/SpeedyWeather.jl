abstract type AbstractForcing <: AbstractModelComponent end

# function barrier for all forcings to unpack model.forcing
function forcing!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        lf::Integer,
        model::AbstractModel
)
    forcing!(diagn, progn, model.forcing, lf, model)
end

## NO FORCING
forcing!(diagn, progn, forcing::Nothing, args...) = nothing

# JET STREAM FORCING
export JetStreamForcing

"""
Forcing term for the Barotropic or ShallowWaterModel with an
idealised jet stream similar to the initial conditions from
Galewsky, 2004, but mirrored for both hemispheres.

$(TYPEDFIELDS)
"""
@parameterized @kwdef mutable struct JetStreamForcing{NF, VectorType} <: AbstractForcing
    "[OPTION] jet latitude [˚N]"
    @param latitude::NF = 45 (bounds = -90..90,)

    "[OPTION] jet width [˚], default ≈ 19.29˚"
    @param width::NF = (1/4-1/7)*180 (bounds = Positive,)

    "[OPTION] sigma level [1], vertical location of jet"
    @param sigma::NF = 0.2 (bounds = Nonnegative,)

    "[OPTION] jet speed scale [m/s]"
    @param speed::NF = 85

    "[OPTION] time scale [days]"
    time_scale::Second = Day(30)

    "[DERIVED] precomputed amplitude vector [m/s²]"
    amplitude::VectorType

    "[DERIVED] precomputed vertical tapering"
    tapering::VectorType
end

function JetStreamForcing(SG::SpectralGrid; kwargs...)
    # Create arrays on the target architecture
    amplitude = on_architecture(SG.architecture, zeros(SG.NF, SG.nlat))
    tapering = on_architecture(SG.architecture, zeros(SG.NF, SG.nlayers))

    return JetStreamForcing{SG.NF, SG.VectorType}(; amplitude, tapering, kwargs...)
end

function initialize!(forcing::JetStreamForcing,
        model::AbstractModel)
    (; latitude, width, speed, time_scale, amplitude) = forcing
    (; radius) = model.planet

    # Some constants similar to Galewsky 2004
    θ₀ = (latitude-width)/360*2π        # southern boundary of jet [radians]
    θ₁ = (latitude+width)/360*2π        # northern boundary of jet
    eₙ = exp(-4/(θ₁-θ₀)^2)              # normalisation, so that speed is at max
    A₀ = speed/eₙ/time_scale.value      # amplitude [m/s²] without lat dependency
    A₀ *= radius                        # scale by radius as are the momentum equations

    (; nlat, colat) = model.geometry
    # latitude in radians, abs for north/south symmetry
    θ = @. abs(π/2 - colat)

    # Similar to u as in Galewsky, 2004 but with north/south symmetry
    @. amplitude = A₀*exp(1/(θ-θ₀)/(θ-θ₁))
    amplitude[.~(θ₀ .< θ .< θ₁)] .= 0   # apply latitude mask

    # vertical tapering
    (; sigma, tapering) = forcing
    σ = model.geometry.σ_levels_full
    tapering .= 1 .- abs.(sigma .- σ)

    return nothing
end

# function barrier
function forcing!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        forcing::JetStreamForcing,
        lf::Integer,
        model::AbstractModel
)
    forcing!(diagn, forcing)
end

"""$(TYPEDSIGNATURES)
Set for every latitude ring the tendency to the precomputed forcing
in the momentum equations following the JetStreamForcing.
The forcing is precomputed in `initialize!(::JetStreamForcing, ::AbstractModel)`."""
function forcing!(
        diagn::DiagnosticVariables,
        forcing::JetStreamForcing)
    Fu = diagn.tendencies.u_tend_grid

    (; amplitude, tapering) = forcing
    (; whichring) = Fu.grid

    arch = architecture(Fu)
    launch!(arch, RingGridWorkOrder, size(Fu), jetstream_forcing_kernel!,
        Fu, amplitude, tapering, whichring)
end

@kernel inbounds=true function jetstream_forcing_kernel!(
        Fu,
        amplitude,
        tapering,
        whichring
)
    ij, k = @index(Global, NTuple)
    j = whichring[ij]
    Fu[ij] += tapering[k] * amplitude[j]
end

export StochasticStirring
@parameterized @kwdef struct StochasticStirring{NF, VectorType} <: AbstractForcing
    "[OPTION] Stirring strength A [1/s²]"
    @param strength::NF = 2e-11

    "[OPTION] Stirring latitude [˚N]"
    @param latitude::NF = 45 (bounds = -90..90,)

    "[OPTION] Stirring width [˚]"
    @param width::NF = 24 (bounds = Positive,)

    "[DERIVED] Latitudinal mask, confined to mid-latitude storm track by default [1]"
    lat_mask::VectorType
end

function StochasticStirring(SG::SpectralGrid; kwargs...)
    lat_mask = on_architecture(SG.architecture, zeros(SG.NF, SG.nlat))
    return StochasticStirring{SG.NF, SG.VectorType}(; lat_mask, kwargs...)
end

function initialize!(
        forcing::StochasticStirring,
        model::AbstractModel)
    model.random_process isa Nothing &&
        @warn "StochasticStirring needs a random process. model.random_process is nothing."

    # precompute the latitudinal mask
    (; latd) = model.geometry

    # Gaussian centred at forcing.latitude of width forcing.width
    @. forcing.lat_mask = exp(-(forcing.latitude-latd)^2/forcing.width^2*2)
    return nothing
end

function forcing!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        forcing::StochasticStirring,
        lf::Integer,
        model::AbstractModel
)
    forcing!(diagn, forcing, model.spectral_transform)
end

function forcing!(
        diagn::DiagnosticVariables,
        forcing::StochasticStirring,
        spectral_transform::SpectralTransform
)
    # get random values from random process
    S_grid = diagn.grid.random_pattern

    # mask everything but mid-latitudes
    RingGrids._scale_lat!(S_grid, forcing.lat_mask)

    # back to spectral space
    S_masked = diagn.dynamics.a_2D
    transform!(S_masked, S_grid, diagn.dynamics.scratch_memory, spectral_transform)

    # scale by radius^2 as is the vorticity equation, and scale to forcing strength
    S_masked .*= (diagn.scale[]^2 * forcing.strength)

    # force every layer
    (; vor_tend) = diagn.tendencies
    arch = architecture(vor_tend)
    launch!(arch, SpectralWorkOrder, size(vor_tend), stochastic_stirring_kernel!,
        vor_tend, S_masked)
end

# GPU kernel for adding 2D spectral field to each layer of 3D field
@kernel inbounds=true function stochastic_stirring_kernel!(
        vor_tend,
        S_masked
)
    I = @index(Global, NTuple)
    lm = I[1]  # spectral coefficient index
    k = I[2]   # layer index

    vor_tend[lm, k] += S_masked[lm]
end

export KolmogorovFlow

"""Kolmogorov flow forcing. Fields are
$(TYPEDFIELDS)
"""
@parameterized @kwdef mutable struct KolmogorovFlow{NF} <: AbstractForcing
    "[OPTION] Strength of forcing [1/s²]"
    @param strength::NF = 1.5e-5

    "[OPTION] Wavenumber of forcing in meridional direction (pole to pole)"
    @param wavenumber::NF = 8 (bounds = Positive,)
end

KolmogorovFlow(SG::SpectralGrid; kwargs...) = KolmogorovFlow{SG.NF}(; kwargs...)

initialize!(::KolmogorovFlow, ::AbstractModel) = nothing

function forcing!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        forcing::KolmogorovFlow,
        lf::Integer,
        model::AbstractModel
)
    # scale by radius as is the vorticity equation
    s = forcing.strength * diagn.scale[]
    k = forcing.wavenumber

    Fu = diagn.tendencies.u_tend_grid
    set!(Fu, (λ, θ, σ) -> s*sind(k*θ), model.geometry)
end

export HeldSuarez

"""Temperature relaxation from Held and Suarez, 1996 BAMS
$(TYPEDFIELDS)"""
@kwdef struct HeldSuarez{NF, VectorType, MatrixType} <: AbstractForcing
    "[OPTION] sigma coordinate below which faster surface relaxation is applied"
    σb::NF = 0.7

    "[OPTION] time scale for slow global relaxation"
    relax_time_slow::Second = Day(40)

    "[OPTION] time scale for faster tropical surface relaxation"
    relax_time_fast::Second = Day(4)

    "[OPTION] minimum equilibrium temperature [K]"
    Tmin::NF = 200

    "[OPTION] maximum equilibrium temperature [K]"
    Tmax::NF = 315

    "[OPTION] meridional temperature gradient [K]"
    ΔTy::NF = 60

    "[OPTION] vertical temperature gradient [K]"
    Δθz::NF = 10

    "[DERIVED] log of sigma level per layer"
    logσ::VectorType

    "[DERIVED] relaxation time scale per layer and latitude (inverse, 1/s)"
    temp_relax_freq::MatrixType

    "[DERIVED] Term a to calculate equilibrium temperature function of latitude"
    temp_equil_a::VectorType

    "[DERIVED] Term b to calculate equilibrium temperature function of latitude"
    temp_equil_b::VectorType
end

"""
$(TYPEDSIGNATURES)
create a HeldSuarez temperature relaxation with arrays allocated given `spectral_grid`"""
function HeldSuarez(SG::SpectralGrid; kwargs...)
    (; NF, VectorType, MatrixType, nlat, nlayers) = SG

    # allocate
    logσ = on_architecture(SG.architecture, zeros(SG.NF, nlayers))
    temp_relax_freq = on_architecture(SG.architecture, zeros(SG.NF, nlayers, nlat))
    temp_equil_a = on_architecture(SG.architecture, zeros(SG.NF, nlat))
    temp_equil_b = on_architecture(SG.architecture, zeros(SG.NF, nlat))

    return HeldSuarez{NF, VectorType, MatrixType}(;
        logσ, temp_relax_freq, temp_equil_a, temp_equil_b, kwargs...)
end

"""$(TYPEDSIGNATURES)
initialize the HeldSuarez temperature relaxation by precomputing terms for the
equilibrium temperature Teq."""
function initialize!(forcing::HeldSuarez,
        model::PrimitiveEquation)
    (; coslat, sinlat) = model.geometry
    σ = model.geometry.σ_levels_full
    (; σb, ΔTy, Δθz, relax_time_slow, relax_time_fast, Tmax) = forcing
    (; logσ, temp_relax_freq, temp_equil_a, temp_equil_b) = forcing

    (; pres_ref) = model.atmosphere

    # slow relaxation everywhere, fast in the tropics
    kₐ = 1/relax_time_slow.value
    kₛ = 1/relax_time_fast.value

    logσ .= log.(σ)               # precompute log(σ) for equilibrium temperature calculation

    # Held and Suarez equation 4
    temp_relax_freq .= kₐ .+ (kₛ - kₐ)*max.(0, (σ .- σb) ./ (1-σb)) .* (coslat') .^ 4

    # Held and Suarez equation 3, split into max(Tmin, (a - b*ln(p))*(p/p₀)^κ)
    # precompute a, b to simplify online calculation
    @. temp_equil_a = Tmax - ΔTy*sinlat^2 + Δθz*log(pres_ref)*coslat^2
    @. temp_equil_b = -Δθz*coslat^2
end

"""$(TYPEDSIGNATURES)
Apply temperature relaxation following Held and Suarez 1996, BAMS."""
function forcing!(
        diagn::DiagnosticVariables,
        progn::PrognosticVariables,
        forcing::HeldSuarez,
        lf::Integer,
        model::AbstractModel
)
    temp_grid = diagn.grid.temp_grid
    pres_grid = diagn.grid.pres_grid
    temp_tend_grid = diagn.tendencies.temp_tend_grid

    (; Tmin, logσ, temp_relax_freq, temp_equil_a, temp_equil_b) = forcing
    (; κ) = model.atmosphere
    σ = model.geometry.σ_levels_full

    (; whichring) = temp_grid.grid
    launch!(architecture(temp_tend_grid), RingGridWorkOrder,
        size(temp_tend_grid), held_suarez_kernel!,
        temp_tend_grid, temp_grid, pres_grid,
        temp_relax_freq, temp_equil_a, temp_equil_b, logσ,
        Tmin, κ, σ, whichring)
end

@kernel inbounds=true function held_suarez_kernel!(
        temp_tend_grid,
        temp_grid,
        pres_grid,
        @Const(temp_relax_freq),
        @Const(temp_equil_a),
        @Const(temp_equil_b),
        @Const(logσ),
        @Const(Tmin),
        @Const(κ),
        @Const(σ),
        @Const(whichring)
)
    ij, k = @index(Global, NTuple)
    j = whichring[ij]                   # latitude ring index
    kₜ = temp_relax_freq[k, j]           # (inverse) relaxation time scale

    # Held and Suarez 1996, equation 3 with precomputed a, b during initialization
    Teq = max(Tmin, (temp_equil_a[j] + temp_equil_b[j]*logσ[k])*σ[k]^κ)
    temp_tend_grid[ij, k] -= kₜ*(temp_grid[ij, k] - Teq)  # Held and Suarez 1996, equation 2
end
