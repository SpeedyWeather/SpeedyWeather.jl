abstract type AbstractForcing <: AbstractModelComponent end

## NO FORCING = dummy forcing
export NoForcing
struct NoForcing <: AbstractForcing end
NoForcing(SG::SpectralGrid) = NoForcing()
initialize!(::NoForcing, ::ModelSetup) = nothing

function forcing!(  diagn::DiagnosticVariablesLayer,
                    progn::PrognosticVariablesLayer,
                    forcing::NoForcing,
                    time::DateTime,
                    model::ModelSetup)
    return nothing
end

# JET STREAM FORCING
export JetStreamForcing

"""
Forcing term for the Barotropic or ShallowWaterModel with an
idealised jet stream similar to the initial conditions from
Galewsky, 2004, but mirrored for both hemispheres.

$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct JetStreamForcing{NF} <: AbstractForcing
    "Number of latitude rings"
    nlat::Int = 0

    "jet latitude [˚N]"
    latitude::NF = 45
    
    "jet width [˚], default ≈ 19.29˚"
    width::NF = (1/4-1/7)*180

    "jet speed scale [m/s]"
    speed::NF = 85

    "time scale [days]"
    time_scale::Second = Day(30)

    "precomputed amplitude vector [m/s²]"
    amplitude::Vector{NF} = zeros(NF, nlat)
end

JetStreamForcing(SG::SpectralGrid; kwargs...) = JetStreamForcing{SG.NF}(
    ; nlat=SG.nlat, kwargs...)

function initialize!(   forcing::JetStreamForcing,
                        model::ModelSetup)

    (; latitude, width, speed, time_scale, amplitude) = forcing
    (; radius) = model.spectral_grid
    
    # Some constants similar to Galewsky 2004
    θ₀ = (latitude-width)/360*2π        # southern boundary of jet [radians]
    θ₁ = (latitude+width)/360*2π        # northern boundary of jet
    eₙ = exp(-4/(θ₁-θ₀)^2)              # normalisation, so that speed is at max
    A₀ = speed/eₙ/time_scale.value      # amplitude [m/s²] without lat dependency
    A₀ *= radius                        # scale by radius as are the momentum equations

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

    return nothing
end

# function barrier
function forcing!(  diagn::DiagnosticVariablesLayer,
                    progn::PrognosticVariablesLayer,
                    forcing::JetStreamForcing,
                    time::DateTime,
                    model::ModelSetup)
    forcing!(diagn, forcing)
end

"""
$(TYPEDSIGNATURES)

Set for every latitude ring the tendency to the precomputed forcing
in the momentum equations following the JetStreamForcing.
The forcing is precomputed in `initialize!(::JetStreamForcing, ::ModelSetup)`."""
function forcing!(  diagn::DiagnosticVariablesLayer,
                    forcing::JetStreamForcing)
    Fu = diagn.tendencies.u_tend_grid
    (; amplitude) = forcing

    @inbounds for (j, ring) in enumerate(eachring(Fu))
        F = amplitude[j]
        for ij in ring
            Fu[ij] = F
        end
    end
end