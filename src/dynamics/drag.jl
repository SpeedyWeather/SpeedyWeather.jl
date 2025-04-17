abstract type AbstractDrag <: AbstractModelComponent end

# function barrier for all drags to unpack model.drag
function drag!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    lf::Integer,
    model::AbstractModel,
)
    drag!(diagn, progn, model.drag, lf, model)
end

## NO DRAG
export NoDrag
struct NoDrag <: AbstractDrag end
NoDrag(SG::SpectralGrid) = NoDrag()
initialize!(::NoDrag, ::AbstractModel) = nothing

function drag!(     diagn::DiagnosticVariables,
                    progn::PrognosticVariables,
                    drag::NoDrag,
                    args...)
    return nothing
end

# Quadratic drag
export QuadraticDrag
@kwdef mutable struct QuadraticDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1]"
    c_D::NF = 1e-5
end

QuadraticDrag(SG::SpectralGrid; kwargs...) = QuadraticDrag{SG.NF}(; kwargs...)
initialize!(::QuadraticDrag, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Quadratic drag for the momentum equations.

    F = -c_D/H*|(u, v)|*(u, v)

with c_D the non-dimensional drag coefficient as defined in `drag::QuadraticDrag`.
c_D and layer thickness `H` are precomputed in initialize!(::QuadraticDrag, ::AbstractModel)
and scaled by the radius as are the momentum equations."""
function drag!(     
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::QuadraticDrag,
    lf::Integer,
    model::AbstractModel,
)
    u = diagn.grid.u_grid
    v = diagn.grid.v_grid

    Fu = diagn.tendencies.u_tend_grid
    Fv = diagn.tendencies.v_tend_grid

    # total drag coefficient with radius scaling
    c = drag.c_D / model.atmosphere.layer_thickness
    c *= diagn.scale[]^2

    k = diagn.nlayers   # only apply to surface layer 
    @inbounds for ij in eachgridpoint(u, v, Fu, Fv)
        speed = sqrt(u[ij, k]^2 + v[ij, k]^2)
        Fu[ij, k] -= c*speed*u[ij, k]     # -= as the tendencies already contain forcing
        Fv[ij, k] -= c*speed*v[ij, k]
    end
end

export LinearVorticityDrag
@kwdef mutable struct LinearVorticityDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1/s]"
    c::NF = 1e-7
end

LinearVorticityDrag(SG::SpectralGrid; kwargs...) = LinearVorticityDrag{SG.NF}(; kwargs...)
initialize!(::LinearVorticityDrag, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Linear drag for the vorticity equations.

    F = -cξ

with c the non-dimensional drag coefficient as defined in `drag::QuadraticDrag`."""
function drag!(     
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::LinearVorticityDrag,
    lf::Integer,
    model::AbstractModel,
)
    (; vor_tend) = diagn.tendencies
    vor = progn.vor[1]
    
    # scale by radius (but only once, the second radius is in vor)
    c = drag.c * diagn.scale[]
    vor_tend .-= c*vor

    return nothing
end

export JetDrag
@kwdef struct JetDrag{NF, SpectralVariable2D} <: SpeedyWeather.AbstractDrag

    # DIMENSIONS from SpectralGrid
    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int

    "[OPTION] Relaxation time scale τ"
    time_scale::Second = Day(6)

    "[OPTION] Jet strength [m/s]"
    u₀::NF = 20

    "[OPTION] latitude of Gaussian jet [˚N]"
    latitude::NF = 30

    "[OPTION] Width of Gaussian jet [˚]"
    width::NF = 6

    # TO BE INITIALISED
    "Relaxation back to reference vorticity"
    ζ₀::SpectralVariable2D = zeros(LowerTriangularMatrix{Complex{NF}}, trunc+2, trunc+1)
end

function JetDrag(SG::SpectralGrid; kwargs...)
    return JetDrag{SG.NF, SG.SpectralVariable2D}(; SG.trunc, kwargs...)
end

function initialize!(drag::JetDrag, model::AbstractModel)
    (; spectral_grid, geometry) = model
    (; Grid, NF, nlat_half) = spectral_grid
    u = zeros(Grid{NF}, nlat_half)

    lat = geometry.latds

    for ij in eachindex(u)
        u[ij] = drag.u₀ * exp(-(lat[ij]-drag.latitude)^2/(2*drag.width^2))
    end

    û = transform(u, model.spectral_transform)
    v̂ = zero(û)
    curl!(drag.ζ₀, û, v̂, model.spectral_transform)
    return nothing
end

function drag!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::JetDrag,
    lf::Integer,
    model::AbstractModel,
)
    vor = progn.vor[lf]
    (; vor_tend) = diagn.tendencies
    (; ζ₀) = drag
    
    # scale by radius as is vorticity
    (; radius) = model.spectral_grid
    r = radius/drag.time_scale.value

    k = diagn.nlayers   # drag only on surface layer
    for lm in eachharmonic(vor_tend)
        vor_tend[lm, k] -= r*(vor[lm, k] - ζ₀[lm])
    end
end