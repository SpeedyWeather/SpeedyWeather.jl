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
@parameterized @kwdef mutable struct QuadraticDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1]"
    @param c_D::NF = 1e-12 (bounds=Nonnegative,)    # TODO is this a good default?
end

QuadraticDrag(SG::SpectralGrid; kwargs...) = QuadraticDrag{SG.NF}(; kwargs...)

initialize!(::QuadraticDrag, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Quadratic drag for the momentum equations.

    F = -c_D/H*|(u, v)|*(u, v)

with `c_D` the non-dimensional drag coefficient as defined in `drag::QuadraticDrag`.
`c_D` and layer thickness `H` are precomputed in `initialize!(::QuadraticDrag, ::AbstractModel)`
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
    
    # GPU kernel launch
    arch = architecture(u) 
    launch!(arch, LinearWorkOrder, (length(u.data),), quadratic_drag_kernel!,
            Fu, Fv, u, v, c, k)
    synchronize(arch)
end

@kernel inbounds=true function quadratic_drag_kernel!(
    Fu, Fv, u, v, c, k 
)
    ij = @index(Global, Linear)
    
    # Calculate speed at surface layer k
    speed = sqrt(u[ij, k]^2 + v[ij, k]^2)
    
    # Apply quadratic drag, -= as the tendencies already contain forcing
    Fu[ij, k] -= c * speed * u[ij, k]
    Fv[ij, k] -= c * speed * v[ij, k]
end

export LinearVorticityDrag
@parameterized @kwdef mutable struct LinearVorticityDrag{NF} <: AbstractDrag
    "[OPTION] drag coefficient [1/s]"
    @param c::NF = 1e-7 (bounds=Nonnegative,)
end

LinearVorticityDrag(SG::SpectralGrid; kwargs...) = LinearVorticityDrag{SG.NF}(; kwargs...)

initialize!(::LinearVorticityDrag, ::AbstractModel) = nothing

"""
$(TYPEDSIGNATURES)
Linear drag for the vorticity equations of the form F = -cξ
with c drag coefficient [1/s]."""


function drag!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::LinearVorticityDrag,
    lf::Integer,
    model::AbstractModel,
)
    (; vor_tend) = diagn.tendencies
    vor = get_step(progn.vor, lf)

    # scale by radius (but only once, the second radius is in vor)
    c = drag.c * diagn.scale[]
    
    # GPU kernel launch
    arch = architecture(vor_tend)
    launch!(arch, SpectralWorkOrder, size(vor_tend), linear_vorticity_drag_kernel!,
            vor_tend, vor, c)
    synchronize(arch)
end

@kernel inbounds=true function linear_vorticity_drag_kernel!(
    vor_tend, vor, c
)
    I = @index(Global, Cartesian)
    vor_tend[I] -= c * vor[I]
end

export JetDrag
@parameterized @kwdef struct JetDrag{NF, SpectralVariable2D} <: SpeedyWeather.AbstractDrag
    "[OPTION] Relaxation time scale τ"
    time_scale::Second = Day(6)

    "[OPTION] Jet strength [m/s]"
    @param u₀::NF = 20

    "[OPTION] latitude of Gaussian jet [˚N]"
    @param latitude::NF = 30 (bounds=-90..90,)

    "[OPTION] Width of Gaussian jet [˚]"
    @param width::NF = 6 (bounds=Positive,)

    # TO BE INITIALISED
    "Relaxation back to reference vorticity"
    ζ₀::SpectralVariable2D
end

function JetDrag(SG::SpectralGrid; kwargs...)
    ζ₀ = zeros(Complex{SG.NF}, SG.spectrum)
    return JetDrag{SG.NF, SG.SpectralVariable2D}(; ζ₀, kwargs...)
end

function initialize!(drag::JetDrag, model::AbstractModel)
    (; spectral_grid, geometry) = model
    (; grid, NF) = spectral_grid
    u = zeros(NF, grid)

    lat = geometry.latds

    for ij in eachindex(u)
        u[ij] = drag.u₀ * exp(-(lat[ij]-drag.latitude)^2/(2*drag.width^2))
    end

    û = transform(u, model.spectral_transform)
    v̂ = zero(û)
    curl!(drag.ζ₀, û, v̂, model.spectral_transform)
    return nothing
end

"""
function drag!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::JetDrag,
    lf::Integer,
    model::AbstractModel,
)
    vor = get_step(progn.vor, lf)
    (; vor_tend) = diagn.tendencies
    (; ζ₀) = drag

    # scale by radius as is vorticity
    s = diagn.scale[]
    r = s/drag.time_scale.value

    k = diagn.nlayers   # drag only on surface layer
    for lm in eachharmonic(vor_tend)
        vor_tend[lm, k] -= r*(vor[lm, k] - ζ₀[lm])
    end
end
"""

function drag!(
    diagn::DiagnosticVariables,
    progn::PrognosticVariables,
    drag::JetDrag,
    lf::Integer,
    model::AbstractModel,
)
    vor = get_step(progn.vor, lf)
    (; vor_tend) = diagn.tendencies
    (; ζ₀) = drag

    # scale by radius as is vorticity
    (; radius) = model.spectral_grid
    r = radius/drag.time_scale.value

    k = diagn.nlayers   # drag only on surface layer
    
    # GPU kernel launch 
    arch = architecture(diagn.grid.u_grid)
    launch!(arch, LinearWorkOrder, (size(vor_tend, 1),), jet_drag_kernel!,
            vor_tend, vor, ζ₀, r, k)
    synchronize(arch)
end

@kernel inbounds=true function jet_drag_kernel!(
    vor_tend, vor, ζ₀, r, k  
)
    lm = @index(Global, Linear)
    vor_tend[lm, k] -= r * (vor[lm, k] - ζ₀[lm])
end