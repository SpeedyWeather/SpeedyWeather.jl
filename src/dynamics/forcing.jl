abstract type AbstractForcing <: AbstractModelComponent end

# function barrier for all forcings to unpack model.forcing
function forcing!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    model::AbstractModel,
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
@parameterized @kwdef mutable struct JetStreamForcing{NF} <: AbstractForcing
    "Number of latitude rings"
    nlat::Int = 0

    "Number of vertical layers"
    nlayers::Int = 0

    "jet latitude [˚N]"
    @param latitude::NF = 45 (bounds=-90..90,)
    
    "jet width [˚], default ≈ 19.29˚"
    @param width::NF = (1/4-1/7)*180 (bounds=Positive,)

    "sigma level [1], vertical location of jet"
    @param sigma::NF = 0.2 (bounds=Nonnegative,)

    "jet speed scale [m/s]"
    @param speed::NF = 85

    "time scale [days]"
    time_scale::Second = Day(30)

    "precomputed amplitude vector [m/s²]"
    amplitude::Vector{NF} = zeros(NF, nlat)

    "precomputed vertical tapering"
    tapering::Vector{NF} = zeros(NF, nlayers)
end

JetStreamForcing(SG::SpectralGrid; kwargs...) = JetStreamForcing{SG.NF}(
    ; nlat=SG.nlat, nlayers=SG.nlayers, kwargs...)

function initialize!(   forcing::JetStreamForcing,
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

    for j in 1:nlat
        # latitude in radians, abs for north/south symmetry
        θ = abs(π/2 - colat[j])
        if θ₀ < θ < θ₁
            # Similar to u as in Galewsky, 2004 but with north/south symmetry
            amplitude[j] = A₀*exp(1/(θ-θ₀)/(θ-θ₁))  
        else
            amplitude[j] = 0
        end
    end

    # vertical tapering
    (; nlayers, sigma, tapering) = forcing
    (; σ_levels_full) = model.geometry

    for k in 1:nlayers
        tapering[k] = 1 - abs(sigma - σ_levels_full[k])
    end

    return nothing
end

# function barrier
function forcing!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    forcing::JetStreamForcing,
    lf::Integer,
    model::AbstractModel,
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

    @inbounds for k in eachlayer(Fu) 
        for (j, ring) in enumerate(eachring(Fu))
            F = amplitude[j]
            for ij in ring
                # += to accumulate, not overwrite previous parameterizations/terms
                Fu[ij] += tapering[k]*F
            end
        end
    end
end

export StochasticStirring
@parameterized @kwdef struct StochasticStirring{NF, VectorType} <: AbstractForcing
        
    "Number of latitude rings, used for latitudinal mask"
    nlat::Int

    "[OPTION] Stirring strength A [1/s²]"
    @param strength::NF = 1e-9

    "[OPTION] Stirring latitude [˚N]"
    @param latitude::NF = 45 (bounds=-90..90,)

    "[OPTION] Stirring width [˚]"
    @param width::NF = 24 (bounds=Positive,)
    
    # TO BE INITIALISED        
    "Latitudinal mask, confined to mid-latitude storm track by default [1]"
    lat_mask::VectorType = zeros(NF, nlat)
end

function StochasticStirring(SG::SpectralGrid; kwargs...)
    return StochasticStirring{SG.NF, SG.VectorType}(; nlat=SG.nlat, kwargs...)
end

function initialize!(
    forcing::StochasticStirring,
    model::AbstractModel)
    
    model.random_process isa Nothing &&
        @warn "StochasticStirring needs a random process. model.random_process is nothing."

    # precompute the latitudinal mask
    (; latd) = model.geometry
    for j in eachindex(forcing.lat_mask)
        # Gaussian centred at forcing.latitude of width forcing.width
        forcing.lat_mask[j] = exp(-(forcing.latitude-latd[j])^2/forcing.width^2*2)
    end
end

function forcing!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    forcing::StochasticStirring,
    lf::Integer,
    model::AbstractModel,
)
    forcing!(diagn, forcing, model.spectral_transform)
end

function forcing!(
    diagn::DiagnosticVariables,
    forcing::StochasticStirring,
    spectral_transform::SpectralTransform,
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

    @inbounds for k in eachmatrix(vor_tend)
        vor_tend[:, k] .+= S_masked
    end
end

export KolmogorovFlow

"""Kolmogorov flow forcing. Fields are
$(TYPEDFIELDS)
"""
@parameterized @kwdef mutable struct KolmogorovFlow{NF} <: AbstractForcing
    "[OPTION] Strength of forcing [1/s²]"
    @param strength::NF = 3e-12

    "[OPTION] Wavenumber of forcing in meridional direction (pole to pole)"
    @param wavenumber::NF = 8 (bounds=Positive,)
end

KolmogorovFlow(SG::SpectralGrid; kwargs...) = KolmogorovFlow{SG.NF}(; kwargs...)

function forcing!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    forcing::KolmogorovFlow,
    lf::Integer,
    model::AbstractModel,
)
    # scale by radius^2 as is the vorticity equation
    s = forcing.strength * diagn.scale[]^2
    k = forcing.wavenumber

    Fu = diagn.tendencies.u_tend_grid
    set!(Fu, (λ, θ, σ) -> s*sind(k*θ), model.geometry)
end